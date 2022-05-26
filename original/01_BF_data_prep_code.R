# This code imports Bri's soil type qSIP data sets and prepares the data for qSIP calculations,
# including calculating tube-level WADs for all OTUs and analyzing the shifts among tubes


#Set working directory and load libraries & scripts:
  setwd("/Users/bk/Research/Projects/SIP_Modeling/qSIP")
  library(reshape2)
  library(VennDiagram)
  source("qSIP_repo/sample.vec.R")                      #sample.vec
  source("qSIP_repo/WAD.func.R")                        #WAD.func
  source("qSIP_repo/fit.norm.func.R")                   #fit.norm.func
  source("qSIP_repo/boot.WAD.func.R")                   #boot.WAD.func
  source("qSIP_repo/diff.wad.calc.R")                   #diff.wad.calc
  source("qSIP_repo/boot.diff.wad.R")                   #boot.diff.wad
  source("qSIP_repo/MW.calc.R")                         #MW.calc
  source("qSIP_repo/MW.calc.Schildkraut.R")             #MW.calc.Schildkraut
  source("qSIP_repo/comparison.message.R")              #comparison.message
  source("qSIP_repo/ape.calc.R")                        #ape.calc
  source("qSIP_repo/boot.diff.ape.R")                   #boot.diff.ape
  source("qSIP_repo/r.calc.R")                          #r.calc
  source("qSIP_repo/boot.diff.r.R")                     #boot.diff.r
  source("qSIP_repo/boot.TUBE.func.R")                  #boot.TUBE.func
  source("qSIP_repo/f.calc.R")                          #f.calc
  source("qSIP_repo/boot.diff.f.R")                     #boot.diff.f
  source("qSIP_repo/all.taxa.calcs.R")                  #all.taxa.calcs
  source("qSIP_repo/id.reps.R")                         #id.reps
  source("qSIP_repo/select.rep.R")                      #select.rep
  source("qSIP_repo/explore.filter.taxa.R")             #explore.filter.taxa
  source("qSIP_repo/filter.taxa.R")                     #filter.taxa
  source("qSIP_repo/explore.filter.fractions.taxa.R")   #explore.filter.fractions.taxa
  source("qSIP_repo/filter.fractions.taxa.R")           #filter.fractions.taxa  
  source("qSIP_repo/WAD.by.taxon.func.R")               #WAD.by.taxon.func
  source("qSIP_repo/SE.WAD.by.taxon.plot.R")            #SE.WAD.by.taxon.plot
  source("qSIP_repo/find.unlabeled.correction.R")       #find.unlabeled.correction
  source("qSIP_repo/find.labeled.correction.R")         #find.labeled.correction
  source("qSIP_repo/td.pos.resid.R")                    #td.pos.resid
  source("qSIP_repo/td.abs.resid.R")                    #td.abs.resid
  source("qSIP_repo/bu.abs.resid.R")                    #bu.abs.resid
  source("qSIP_repo/select.best.iteration.R")           #select.best.iteration
  source("qSIP_repo/find.labeled.correction.plot.R")    #find.labeled.correction.plot
  source("qSIP_repo/get.seq.taxa.nums.R")               #get.seq.taxa.nums
  source("qSIP_repo/add.lab.WAD.corr.summary.R")        #add.lab.WAD.corr.summary
  source("qSIP_repo/apply.unlabeled.correction.R")      #apply.unlabeled.correction
  source("qSIP_repo/apply.labeled.correction.R")        #apply.labeled.correction


#Import raw data & taxomic information and format raw data for analysis:
  #Read in taxonomic information
  # modification: to make melting faster, reading in taxonomy in 1 column with '_' separating phylum, kingdom...
    taxa.id <- read.table("qSIP_data/BF_L7_taxa_ids.txt", header=T, sep="\t", stringsAsFactors=T, na.strings="")
    head(taxa.id)
    names(taxa.id)
  #The number serves as a unique identifier code for each taxon to the "taxa.id" data frame:
    
  #NOTE: treatment code names cannot contain spaces
    data.all <- read.table("qSIP_data/BF_qSIP_data.txt", header=T, sep="\t", stringsAsFactors=T, check.names=F)
    dim(data.all)
    names(data.all)
    head(data.all)
    
  #Quality screening of the data: exclude samples (fractions) without sequencing data:
    data <- data.all[apply(is.na(data.all[, 16:dim(data.all)[2]]), 1, sum) != length(16:dim(data.all)[2]), ]
  #Note that some samples (fractions) have zero 16S copies for qPCR data (these are left in for now):
    data[data$avg.neat.16S.copies == 0, 1:15]
  
  #Add a column for the sum of proportional abundance of all taxa by fraction
    data$sum.abundance <- rowSums(data[,16:ncol(data)])
    data$sum.abundance
    
  #Calculate number of copies per uL, based on relative abundance and total number of copies per uL:
    ncopies <- data$avg.neat.16S.copies*data[,16:(ncol(data)-1)]
    dim(ncopies)
    ncopies <- cbind(data[,1:15], ncopies)  # add first 15 columns of data to ncopies
    dim(ncopies)
    head(ncopies)
    names(ncopies)
    
  #Melt data into long format by tube, sample, tmt, rep, fraction, DNA conc, and density;
  #Do this for copies.ul and for relative abundance, which is just our data file. Merge these to into 1 masterfile: data.melted
    ncopies.melted <- melt(ncopies, id=c("Sample",  "Density_g_ml", "DNA_ng_uL"  , "DNA_ng_fraction", "DNA_proportion", "Week", "Soil", "Substrate", "IsotopeTreat", "soil_C_isotope", "ComboTrt", "Fraction", "Tube", "avg.neat.16S.copies", "Sequenced"), measure.vars=as.character(1:length(taxa.id$taxon)), variable.name="taxon", value.name="copies.ul")
    rel.abundance.melted <- melt(data, id=c("Sample",  "Density_g_ml", "DNA_ng_uL"  , "DNA_ng_fraction", "DNA_proportion", "Week", "Soil", "Substrate", "IsotopeTreat", "soil_C_isotope", "ComboTrt", "Fraction", "Tube", "avg.neat.16S.copies", "Sequenced"), measure.vars=as.character(1:length(taxa.id$taxon)), variable.name="taxon", value.name="rel.abundance")
    data.melted <- merge(ncopies.melted, rel.abundance.melted)
    head(data.melted)
  
  #Merge taxa data and reorder data frame by taxon and SampleID and Fraction
    data.melted <- merge(data.melted, taxa.id)
    data.melted <- data.melted[order(data.melted$taxon, data.melted$Sample, data.melted$Fraction),]
    row.names(data.melted) <- 1:dim(data.melted)[1]     #rename observations to be sequential
    head(data.melted)


  #Import data frame containing the treatment comparisons to perform (split data into threee chunks according to soil type (AN, BS, GR):
    Tcompare1 <- read.table("qSIP_data/BF_qSIP_TreatmentComparisons1.txt", header=TRUE, sep="\t", colClasses=c("numeric","factor","character","character","character","numeric","character"))
    summary(Tcompare1)
    Tcompare1

    Tcompare2 <- read.table("qSIP_data/BF_qSIP_TreatmentComparisons2.txt", header=TRUE, sep="\t", colClasses=c("numeric","factor","character","character","character","numeric","character"))
    summary(Tcompare2)
    Tcompare2

    Tcompare3 <- read.table("qSIP_data/BF_qSIP_TreatmentComparisons3.txt", header=TRUE, sep="\t", colClasses=c("numeric","factor","character","character","character","numeric","character"))
    summary(Tcompare3)
    Tcompare3


#Calculate corrected WADs for each tube and taxon:
  WAD.by.taxon <- WAD.by.taxon.func(X=data.melted, vars=c("taxon", "Density_g_ml", "copies.ul", "Tube", "ComboTrt"))
  #Look at the results:
    names(WAD.by.taxon)          #names of the two data frames in the output list
    head(WAD.by.taxon$obs.wads)  #looking at the head of the first data frame
    head(WAD.by.taxon[[1]])      #another way to look at the head of the first data frame
    WAD.by.taxon$reps.by.trt     #looking at the head of the second data frame
    WAD.by.taxon[[2]]            #another way to look at the head of the second data frame


##### Perform the tube-level WAD corrections for ALL replicates of ALL treatments, and then proceed with filtering and standard qSIP analysis:


#Calculate the shift in WAD for the specified unlabeled replicates from their global mean using the taxa common to ALL replicates of ALL treatments:
  #Identify the appropriate shift for unlabeled WADs:
  #Note: this function takes a while to run (i.e., ~_?_min for ~150 taxa common to all unlabeled and all labeled treatments):
    #Include all unlabeled treatments and all labeled treatments:
      system.time(unlab.WAD.corr.list <- find.unlabeled.correction(LIST=WAD.by.taxon, unlab.tmts=c("1_AN_NoC_16O", "1_AN_Exu_16O", "1_AN_Lit_16O", "6_AN_Exu_16O", "6_AN_Lit_16O", "1_BS_NoC_16O", "1_BS_Exu_16O", "6_BS_NoC_16O", "6_BS_Exu_16O", "1_GR_Exu_16O", "1_GR_Lit_16O", "6_GR_NoC_16O", "6_GR_Exu_16O", "6_GR_Lit_16O"), lab.tmts=c("1_AN_NoC_18O", "1_BS_NoC_18O", "1_GR_NoC_18O", "1_AN_Exu_18O", "1_BS_Exu_18O", "1_GR_Exu_18O", "1_AN_Lit_18O", "1_BS_Lit_18O", "1_GR_Lit_18O", "6_AN_NoC_18O", "6_BS_NoC_18O", "6_GR_NoC_18O", "6_AN_Exu_18O", "6_BS_Exu_18O","6_GR_Exu_18O", "6_AN_Lit_18O", "6_BS_Lit_18O", "6_GR_Lit_18O"), CI=0.90))
      #Look at the results:
        names(unlab.WAD.corr.list)                           #names of the two data frames in the output list
        unlab.WAD.corr.list$WAD.norm.fit.parms               #looking at the first object -- a data frame
        unlab.WAD.corr.list[[1]]                             #another way to look at the first object
        unlab.WAD.corr.list$corr.names                       #looking at the second object in the list -- a vector of the names of the corrected unlabeled replicates
        unlab.WAD.corr.list[[2]]                             #another way to look at the second object
        head(unlab.WAD.corr.list$WAD.table.corr)             #looking at the head of the third object -- a data frame
        head(unlab.WAD.corr.list[[3]])                       #another way to look at the head of the third object


#Calculate the shift in WAD for each of the labeled replicates (18 labeled treatments x 4 replicates = 72 total labeled replicates):
  #1_AN_Exu_18O:
  system.time(R1.1.AN.Exu.18O.td.pos.resid.list <- find.labeled.correction(LIST=unlab.WAD.corr.list, lab.replicate="R1.1_AN_Exu_18O", lab.names=c("R1.1_AN_Exu_18O", "R2.1_AN_Exu_18O", "R3.1_AN_Exu_18O", "R4.1_AN_Exu_18O"), method="td.pos.resid", unlab.SD.percentile=0.5, lab.SD.percentile=0.5, min.num.nongrowers=10))
    #Find the best iteration:
    R1.1.AN.Exu.18O.td.pos.resid.best.iteration.list <- select.best.iteration(putative.nongrower.metrics=R1.1.AN.Exu.18O.td.pos.resid.list$putative.nongrower.metrics, quantile.threshold=0.10)
    find.labeled.correction.plot(find.labeled.correction.list=R1.1.AN.Exu.18O.td.pos.resid.list, filename="BF_R1_1_AN_Exu_18O_td_pos_resid_plots", path="qSIP_output/Figures/", highlight=R1.1.AN.Exu.18O.td.pos.resid.best.iteration.list$best.iteration, highlight.col="magenta")

  system.time(R2.1.AN.Exu.18O.td.pos.resid.list <- find.labeled.correction(LIST=unlab.WAD.corr.list, lab.replicate="R2.1_AN_Exu_18O", lab.names=c("R1.1_AN_Exu_18O", "R2.1_AN_Exu_18O", "R3.1_AN_Exu_18O", "R4.1_AN_Exu_18O"), method="td.pos.resid", unlab.SD.percentile=0.5, lab.SD.percentile=0.5, min.num.nongrowers=10))
    #Find the best iteration:
    R2.1.AN.Exu.18O.td.pos.resid.best.iteration.list <- select.best.iteration(putative.nongrower.metrics=R2.1.AN.Exu.18O.td.pos.resid.list$putative.nongrower.metrics, quantile.threshold=0.10)
    find.labeled.correction.plot(find.labeled.correction.list=R2.1.AN.Exu.18O.td.pos.resid.list, filename="BF_R2_1_AN_Exu_18O_td_pos_resid_plots", path="qSIP_output/Figures/", highlight=R2.1.AN.Exu.18O.td.pos.resid.best.iteration.list$best.iteration, highlight.col="magenta")

  system.time(R3.1.AN.Exu.18O.td.pos.resid.list <- find.labeled.correction(LIST=unlab.WAD.corr.list, lab.replicate="R3.1_AN_Exu_18O", lab.names=c("R1.1_AN_Exu_18O", "R2.1_AN_Exu_18O", "R3.1_AN_Exu_18O", "R4.1_AN_Exu_18O"), method="td.pos.resid", unlab.SD.percentile=0.5, lab.SD.percentile=0.5, min.num.nongrowers=10))
    #Find the best iteration:
    R3.1.AN.Exu.18O.td.pos.resid.best.iteration.list <- select.best.iteration(putative.nongrower.metrics=R3.1.AN.Exu.18O.td.pos.resid.list$putative.nongrower.metrics, quantile.threshold=0.10)
    find.labeled.correction.plot(find.labeled.correction.list=R3.1.AN.Exu.18O.td.pos.resid.list, filename="BF_R3_1_AN_Exu_18O_td_pos_resid_plots", path="qSIP_output/Figures/", highlight=R3.1.AN.Exu.18O.td.pos.resid.best.iteration.list$best.iteration, highlight.col="magenta")

  system.time(R4.1.AN.Exu.18O.td.pos.resid.list <- find.labeled.correction(LIST=unlab.WAD.corr.list, lab.replicate="R4.1_AN_Exu_18O", lab.names=c("R1.1_AN_Exu_18O", "R2.1_AN_Exu_18O", "R3.1_AN_Exu_18O", "R4.1_AN_Exu_18O"), method="td.pos.resid", unlab.SD.percentile=0.5, lab.SD.percentile=0.5, min.num.nongrowers=10))
    #Find the best iteration:
    R4.1.AN.Exu.18O.td.pos.resid.best.iteration.list <- select.best.iteration(putative.nongrower.metrics=R4.1.AN.Exu.18O.td.pos.resid.list$putative.nongrower.metrics, quantile.threshold=0.10)
    find.labeled.correction.plot(find.labeled.correction.list=R4.1.AN.Exu.18O.td.pos.resid.list, filename="BF_R4_1_AN_Exu_18O_td_pos_resid_plots", path="qSIP_output/Figures/", highlight=R4.1.AN.Exu.18O.td.pos.resid.best.iteration.list$best.iteration, highlight.col="magenta")

  # ... NEED TO ADD CORRECTION COMMANDS FOR ALL REPLICATES OF ALL LABELED TREATMENTS HERE (FOLLOW THE FORMAT FOR TREATMENT '1_AN_Exu_18O', ABOVE; NOTE: THERE ARE 4 MISSING LABELED REPLICATES ("R4.6_BS_Exu_18O", "R4.6_BS_NoC_18O", "R4.1_GR_Lit_18O", "R4.6_GR_Exu_18O"), SO DON'T RUN THE COMMANDS FOR THEM)


#Correct the summarized and raw data for tube-level shifts in WAD identified above for ALL unlabeled and labeled replicates:
  #Apply the tube-level shift to the already-summarized taxon-level labeled replicate WADs to 'correct' them and create a new summary list (analagous to 'unlab.WAD.corr.list') that includes those corrected values:
    #First, define character vectors listing the set of labeled replicates and the names of the shifts corresponding to those replicates:
      labeled.treatments <- c("1_AN_Exu_18O", "1_AN_Lit_18O", "1_AN_NoC_18O", "1_BS_Exu_18O", "1_BS_Lit_18O", "1_BS_NoC_18O", "1_GR_Exu_18O", "1_GR_Lit_18O", "1_GR_NoC_18O", "6_AN_Exu_18O", "6_AN_Lit_18O", "6_AN_NoC_18O", "6_BS_Exu_18O", "6_BS_Lit_18O", "6_BS_NoC_18O", "6_GR_Exu_18O", "6_GR_Lit_18O", "6_GR_NoC_18O")
      all.possible.labeled.replicates <- as.vector(outer(c("R1", "R2", "R3", "R4"), labeled.treatments, paste, sep="."))
      missing.labeled.replicates <- c("R4.6_BS_Exu_18O", "R4.6_BS_NoC_18O", "R4.1_GR_Lit_18O", "R4.6_GR_Exu_18O")         #there are four missing replicates from labeled treatments in this experiment
      actual.labeled.replicates <- all.possible.labeled.replicates[!(all.possible.labeled.replicates %in% missing.labeled.replicates)]
      actual.labeled.replicates.dots <- gsub(pattern="\\_", replacement="\\.", x=actual.labeled.replicates, perl=TRUE)
      actual.labeled.replicates.dots.corrections <- paste(actual.labeled.replicates.dots, "td.pos.resid.best.iteration.list$norm.correction", sep=".")
      #Evaluate the names of the corrections to get a vector of the corrections themselves:
      actual.labeled.replicates.dots.corrections.values <- unlist(lapply(as.list(actual.labeled.replicates.dots.corrections), function(x) eval(parse(text=x))))
    #Now, apply the tube-level shifts for all of those replicates:
      WAD.corr.list <- add.lab.WAD.corr.summary(summary.list=unlab.WAD.corr.list, reps.by.trt=WAD.by.taxon$reps.by.trt, lab.reps.to.add=actual.labeled.replicates, lab.shifts=actual.labeled.replicates.dots.corrections.values, CI=0.90)
      #Look at the results:
        names(WAD.corr.list)                           #names of the two data frames in the output list
        WAD.corr.list$WAD.norm.fit.parms               #looking at the first object -- a data frame
        WAD.corr.list[[1]]                             #another way to look at the first object
        WAD.corr.list$corr.names                       #looking at the second object in the list -- a vector of the names of the corrected unlabeled replicates
        WAD.corr.list[[2]]                             #another way to look at the second object
        head(WAD.corr.list$WAD.table.corr)             #looking at the head of the third object -- a data frame
        head(WAD.corr.list[[3]])                       #another way to look at the head of the third object

  #Apply the tube-level shifts to 'correct' the unlabeled WADs in 'data' (corrects all unlabeled replicates at once):
    data.corr <- apply.unlabeled.correction(raw.data=data, correction.table=unlab.WAD.corr.list$WAD.norm.fit.parms, reps.by.trt=WAD.by.taxon$reps.by.trt, vars=c("Density_g_ml", "Tube", "ComboTrt"))
  #Apply the tube-level shift to 'correct' the labeled WADs in 'data.corr' (corrects one labeled replicate at a time):
    for (i in 1:length(actual.labeled.replicates)){
      data.corr <- apply.labeled.correction(raw.data=data, raw.data.corr=data.corr, lab.replicate=actual.labeled.replicates[i], correction.value=eval(parse(text=actual.labeled.replicates.dots.corrections[i])), reps.by.trt=WAD.by.taxon$reps.by.trt, vars=c("Density_g_ml", "Tube", "ComboTrt"))
    }

  #Re-calculate number of copies per uL, based on relative abundance and total number of copies per uL:
    ncopies.corr <- data.corr$neat.avg.16S.copies*data.corr[,16:(ncol(data.corr)-1)]
    dim(ncopies.corr)
    ncopies.corr <- cbind(data.corr[,1:15], ncopies.corr)  # add first 15 columns of data.corr to ncopies.corr
    dim(ncopies.corr)
    head(ncopies.corr)

  #Melt data.corr into long format by tube, sample, tmt, rep, fraction, DNA conc, and density;
  #Do this for copies.ul and for relative abundance, which is just our data.corr file. Merge these to into 1 masterfile: data.corr.melted
    ncopies.corr.melted <- melt(ncopies.corr, id=c("Sample",  "Density_g_ml", "DNA_ng_uL"  , "DNA_ng_fraction", "DNA_proportion", "Week", "Soil", "Substrate", "IsotopeTreat", "soil_C_isotope", "ComboTrt", "Fraction", "Tube", "avg.neat.16S.copies", "Sequenced"), measure.vars=as.character(1:length(taxa.id$taxon)), variable.name="taxon", value.name="copies.ul")
    rel.abundance.corr.melted <- melt(data.corr, id=c("Sample",  "Density_g_ml", "DNA_ng_uL"  , "DNA_ng_fraction", "DNA_proportion", "Week", "Soil", "Substrate", "IsotopeTreat", "soil_C_isotope", "ComboTrt", "Fraction", "Tube", "avg.neat.16S.copies", "Sequenced"), measure.vars=as.character(1:length(taxa.id$taxon)), variable.name="taxon", value.name="rel.abundance")
    data.corr.melted <- merge(ncopies.corr.melted, rel.abundance.corr.melted)
    head(data.corr.melted)

  #Merge taxa data and reorder data frame by taxon and SampleID and fraction:
    data.corr.melted <- merge(data.corr.melted, taxa.id)
    data.corr.melted <- data.corr.melted[order(data.corr.melted$taxon, data.corr.melted$Sample, data.corr.melted$Fraction),]
    row.names(data.corr.melted) <- 1:dim(data.corr.melted)[1]   #rename observations to be sequential
    head(data.corr.melted)


#Calculate corrected WADs for each tube and taxon:
  WAD.by.taxon.corr <- WAD.by.taxon.func(X=data.corr.melted, vars=c("taxon", "Density_g_ml", "copies.ul", "Tube", "ComboTrt"))
  #Look at the results:
    names(WAD.by.taxon.corr)          #names of the two data frames in the output list
    head(WAD.by.taxon.corr$obs.wads)  #looking at the head of the first data frame
    head(WAD.by.taxon.corr[[1]])      #another way to look at the head of the first data frame
    WAD.by.taxon.corr$reps.by.trt     #looking at the head of the second data frame
    WAD.by.taxon.corr[[2]]            #another way to look at the head of the second data frame

  #Quick and dirty look at the level of replication among treatments and the occurrences of data for taxa in tubes 
    table(data.corr.melted$ComboTrt, data.corr.melted$Tube)
    table(data.corr.melted$taxon, data.corr.melted$Tube)


#Plot standard error of WADs for each taxon and treatment:
  graphics.off()
  pdf(file="qSIP_output/Figures/BF_SE_WAD_by_taxon_plots.pdf", width=6, height=11)
  WAD.SE.by.taxon.corr <- SE.WAD.by.taxon.plot(LIST=WAD.by.taxon.corr, percentile=0.95)
  dev.off()
  #look at the head of the output data frame:
    head(WAD.SE.by.taxon.corr)


#Spot-checking that SE's match up with observations:
  WAD.by.taxon.corr$reps.by.trt
  head(WAD.by.taxon.corr$obs.wads)
  head(WAD.SE.by.taxon.corr)
  sum(as.numeric(apply(WAD.by.taxon.corr$obs.wads[,16:19], 1, function(x) sd(x)/sqrt(sum(!is.na(x))))) != WAD.SE.by.taxon.corr[,2], na.rm=TRUE)
  sum(as.numeric(apply(WAD.by.taxon.corr$obs.wads[,c(11,12,14,15)], 1, function(x) sd(x)/sqrt(sum(!is.na(x))))) != WAD.SE.by.taxon.corr[,17], na.rm=TRUE)


#Identify missing replicates (tubes):
  #Look at replicates per treatment (by soil type):
    #AN:
    data.melted.AN <- data.melted[data.melted$Soil == "AN", ]
      #Convert the factor columns to factor:
      data.melted.AN <- as.data.frame(lapply(data.melted.AN, function(x) if(is.factor(x)) factor(x) else x))
    table(data.melted.AN$ComboTrt, data.melted.AN$Tube)
    #Missing replicates:
    #  6-AN-NoC-16O     (1 of 1 replicate missing)

    #BS:
    data.melted.BS <- data.melted[data.melted$Soil == "BS", ]
      #Convert the factor columns to factor:
      data.melted.BS <- as.data.frame(lapply(data.melted.BS, function(x) if(is.factor(x)) factor(x) else x))
    table(data.melted.BS$ComboTrt, data.melted.BS$Tube)
    #Missing replicates:
    #  1-BS-Lit-16O     (1 of 1 replicate missing)
    #  6-BS-Lit-16O     (1 of 1 replicate missing)
    #  6-BS-Exu-18O     (1 of 4 replicates missing)
    #  6-BS-NoC-18O     (1 of 4 replicates missing)

    #GR:
    data.melted.GR <- data.melted[data.melted$Soil == "GR", ]
      #Convert the factor columns to factor:
      data.melted.GR <- as.data.frame(lapply(data.melted.GR, function(x) if(is.factor(x)) factor(x) else x))
    table(data.melted.GR$ComboTrt, data.melted.GR$Tube)
    #Missing replicates:
    #  1-GR-NoC-16O     (1 of 1 replicate missing)
    #  1-GR-Lit-18O     (1 of 4 replicates missing)
    #  6-GR-Exu-18O     (1 of 4 replicates missing)

  #Check the number of tubes:
    x1 <- dimnames(table(data.melted.AN$ComboTrt, data.melted.AN$Tube))[[2]]
    x2 <- dimnames(table(data.melted.BS$ComboTrt, data.melted.BS$Tube))[[2]]
    x3 <- dimnames(table(data.melted.GR$ComboTrt, data.melted.GR$Tube))[[2]]
    intersect(x1, x2)
    intersect(x1, x3)
    intersect(x2, x3)
    length(x1) + length(x2) + length(x3)
    length(unique(as.character(data.melted$Tube)))
    #If no missing tubes, expect: 3 soil types * 3 carbon tmts * 2 time points * 5 reps (16O & 18O):
    3*3*2*5
  
  rm(data.melted.AN, data.melted.BS, data.melted.GR, x1, x2, x3)  


#Filter the data so that only taxa with an appropriate level of occurrence and replication among the tubes and treatments of the experiment are retained for further analysis:
  #First back up the complete data frame:
    data.melted.all <- data.melted

  #*****DO 3 OF 4/6 REPS?  OR  2 OF 4/6 REPS?*****
  #For now, we have elected to a simple, straightforward criterion that a taxon must occur in at least 3 (of 4) reps of each labeled 18O treatment within a soil type (there are 6 18O treatments within each of the three soil types) AND the taxon must also occur on at least 3 (of 6 total) reps of the 6 unlabeled treatments for that soil type (each unlabeled treatment has only 1 rep).  Note that some reps are missing, and this criterion accounts for these missing reps.  It seems to work well so far; if it does not, then we will revisit this filtering step later.
  #Note: there are 8 missing reps; one from each of these treatments (16O treatments only have 1 rep; 18O treatments have 4 reps):
    #6-AN-NoC-16O
    #1-BS-Lit-16O
    #6-BS-Lit-16O
    #6-BS-Exu-18O
    #6-BS-NoC-18O
    #1-GR-NoC-16O
    #1-GR-Lit-18O
    #6-GR-Exu-18O

  #Filtering for soil type 'AN':
  #1-AN-NoC-18O:
    Tcompare1$trt.code.2[1]
    explore.filter.taxa(DATA=data.melted, trt.code.1=NULL, trt.code.2=Tcompare1$trt.code.2[1], trt.refs=NULL, vars=c("taxon", "copies.ul", "Tube", "ComboTrt"), min.reps=3)
    dev.off()
    data.melted.1 <- filter.taxa(DATA=data.melted, trt.code.1=NULL, trt.code.2=Tcompare1$trt.code.2[1], trt.refs=NULL, vars=c("taxon", "copies.ul", "Tube", "ComboTrt"), min.reps=3)
  #1-AN-Exu-18O:
    Tcompare1$trt.code.2[2]
    explore.filter.taxa(DATA=data.melted.1, trt.code.1=NULL, trt.code.2=Tcompare1$trt.code.2[2], trt.refs=NULL, vars=c("taxon", "copies.ul", "Tube", "ComboTrt"), min.reps=3)
    dev.off()
    data.melted.1 <- filter.taxa(DATA=data.melted.1, trt.code.1=NULL, trt.code.2=Tcompare1$trt.code.2[2], trt.refs=NULL, vars=c("taxon", "copies.ul", "Tube", "ComboTrt"), min.reps=3)
  #1-AN-Lit-18O:
    Tcompare1$trt.code.2[3]
    explore.filter.taxa(DATA=data.melted.1, trt.code.1=NULL, trt.code.2=Tcompare1$trt.code.2[3], trt.refs=NULL, vars=c("taxon", "copies.ul", "Tube", "ComboTrt"), min.reps=3)
    dev.off()
    data.melted.1 <- filter.taxa(DATA=data.melted.1, trt.code.1=NULL, trt.code.2=Tcompare1$trt.code.2[3], trt.refs=NULL, vars=c("taxon", "copies.ul", "Tube", "ComboTrt"), min.reps=3)
  #6-AN-NoC-18O:
    Tcompare1$trt.code.2[4]
    explore.filter.taxa(DATA=data.melted.1, trt.code.1=NULL, trt.code.2=Tcompare1$trt.code.2[4], trt.refs=NULL, vars=c("taxon", "copies.ul", "Tube", "ComboTrt"), min.reps=3)
    dev.off()
    data.melted.1 <- filter.taxa(DATA=data.melted.1, trt.code.1=NULL, trt.code.2=Tcompare1$trt.code.2[4], trt.refs=NULL, vars=c("taxon", "copies.ul", "Tube", "ComboTrt"), min.reps=3)
  #6-AN-Exu-18O:
    Tcompare1$trt.code.2[5]
    explore.filter.taxa(DATA=data.melted.1, trt.code.1=NULL, trt.code.2=Tcompare1$trt.code.2[5], trt.refs=NULL, vars=c("taxon", "copies.ul", "Tube", "ComboTrt"), min.reps=3)
    dev.off()
    data.melted.1 <- filter.taxa(DATA=data.melted.1, trt.code.1=NULL, trt.code.2=Tcompare1$trt.code.2[5], trt.refs=NULL, vars=c("taxon", "copies.ul", "Tube", "ComboTrt"), min.reps=3)
  #6-AN-Lit-18O:
    Tcompare1$trt.code.2[6]
    explore.filter.taxa(DATA=data.melted.1, trt.code.1=NULL, trt.code.2=Tcompare1$trt.code.2[6], trt.refs=NULL, vars=c("taxon", "copies.ul", "Tube", "ComboTrt"), min.reps=3)
    dev.off()
    data.melted.1 <- filter.taxa(DATA=data.melted.1, trt.code.1=NULL, trt.code.2=Tcompare1$trt.code.2[6], trt.refs=NULL, vars=c("taxon", "copies.ul", "Tube", "ComboTrt"), min.reps=3)
  #1-AN-NoC-16O; 1-AN-Exu-16O; 1-AN-Lit-16O; 6-AN-NoC-16O; 6-AN-Exu-16O; 6-AN-Lit-16O:
    Tcompare1$trt.code.1[1]
    Tcompare1$trt.refs[1]
    explore.filter.taxa(DATA=data.melted.1, trt.code.1=NULL, trt.code.2=Tcompare1$trt.code.1[1], trt.refs=Tcompare1$trt.refs[1], vars=c("taxon", "copies.ul", "Tube", "ComboTrt"), min.reps=3)
    dev.off()
    data.melted.1 <- filter.taxa(DATA=data.melted.1, trt.code.1=NULL, trt.code.2=Tcompare1$trt.code.1[1], trt.refs=Tcompare1$trt.refs[1], vars=c("taxon", "copies.ul", "Tube", "ComboTrt"), min.reps=3)

  #Filtering for soil type 'BS':
  #1-BS-NoC-18O:
    Tcompare2$trt.code.2[1]
    explore.filter.taxa(DATA=data.melted, trt.code.1=NULL, trt.code.2=Tcompare2$trt.code.2[1], trt.refs=NULL, vars=c("taxon", "copies.ul", "Tube", "ComboTrt"), min.reps=3)
    dev.off()
    data.melted.2 <- filter.taxa(DATA=data.melted, trt.code.1=NULL, trt.code.2=Tcompare2$trt.code.2[1], trt.refs=NULL, vars=c("taxon", "copies.ul", "Tube", "ComboTrt"), min.reps=3)
  #1-BS-Exu-18O:
    Tcompare2$trt.code.2[2]
    explore.filter.taxa(DATA=data.melted.2, trt.code.1=NULL, trt.code.2=Tcompare2$trt.code.2[2], trt.refs=NULL, vars=c("taxon", "copies.ul", "Tube", "ComboTrt"), min.reps=3)
    dev.off()
    data.melted.2 <- filter.taxa(DATA=data.melted.2, trt.code.1=NULL, trt.code.2=Tcompare2$trt.code.2[2], trt.refs=NULL, vars=c("taxon", "copies.ul", "Tube", "ComboTrt"), min.reps=3)
  #1-BS-Lit-18O:
    Tcompare2$trt.code.2[3]
    explore.filter.taxa(DATA=data.melted.2, trt.code.1=NULL, trt.code.2=Tcompare2$trt.code.2[3], trt.refs=NULL, vars=c("taxon", "copies.ul", "Tube", "ComboTrt"), min.reps=3)
    dev.off()
    data.melted.2 <- filter.taxa(DATA=data.melted.2, trt.code.1=NULL, trt.code.2=Tcompare2$trt.code.2[3], trt.refs=NULL, vars=c("taxon", "copies.ul", "Tube", "ComboTrt"), min.reps=3)
  #6-BS-NoC-18O:
    Tcompare2$trt.code.2[4]
    explore.filter.taxa(DATA=data.melted.2, trt.code.1=NULL, trt.code.2=Tcompare2$trt.code.2[4], trt.refs=NULL, vars=c("taxon", "copies.ul", "Tube", "ComboTrt"), min.reps=3)
    dev.off()
    data.melted.2 <- filter.taxa(DATA=data.melted.2, trt.code.1=NULL, trt.code.2=Tcompare2$trt.code.2[4], trt.refs=NULL, vars=c("taxon", "copies.ul", "Tube", "ComboTrt"), min.reps=3)
  #6-BS-Exu-18O:
    Tcompare2$trt.code.2[5]
    explore.filter.taxa(DATA=data.melted.2, trt.code.1=NULL, trt.code.2=Tcompare2$trt.code.2[5], trt.refs=NULL, vars=c("taxon", "copies.ul", "Tube", "ComboTrt"), min.reps=3)
    dev.off()
    data.melted.2 <- filter.taxa(DATA=data.melted.2, trt.code.1=NULL, trt.code.2=Tcompare2$trt.code.2[5], trt.refs=NULL, vars=c("taxon", "copies.ul", "Tube", "ComboTrt"), min.reps=3)
  #6-BS-Lit-18O:
    Tcompare2$trt.code.2[6]
    explore.filter.taxa(DATA=data.melted.2, trt.code.1=NULL, trt.code.2=Tcompare2$trt.code.2[6], trt.refs=NULL, vars=c("taxon", "copies.ul", "Tube", "ComboTrt"), min.reps=3)
    dev.off()
    data.melted.2 <- filter.taxa(DATA=data.melted.2, trt.code.1=NULL, trt.code.2=Tcompare2$trt.code.2[6], trt.refs=NULL, vars=c("taxon", "copies.ul", "Tube", "ComboTrt"), min.reps=3)
  #1-BS-NoC-16O; 1-BS-Exu-16O; 1-BS-Lit-16O; 6-BS-NoC-16O; 6-BS-Exu-16O; 6-BS-Lit-16O:
    Tcompare2$trt.code.1[1]
    Tcompare2$trt.refs[1]
    explore.filter.taxa(DATA=data.melted.2, trt.code.1=NULL, trt.code.2=Tcompare2$trt.code.1[1], trt.refs=Tcompare2$trt.refs[1], vars=c("taxon", "copies.ul", "Tube", "ComboTrt"), min.reps=3)
    dev.off()
    data.melted.2 <- filter.taxa(DATA=data.melted.2, trt.code.1=NULL, trt.code.2=Tcompare2$trt.code.1[1], trt.refs=Tcompare2$trt.refs[1], vars=c("taxon", "copies.ul", "Tube", "ComboTrt"), min.reps=3)

  #Filtering for soil type 'GR':
  #1-GR-NoC-18O:
    Tcompare3$trt.code.2[1]
    explore.filter.taxa(DATA=data.melted, trt.code.1=NULL, trt.code.2=Tcompare3$trt.code.2[1], trt.refs=NULL, vars=c("taxon", "copies.ul", "Tube", "ComboTrt"), min.reps=3)
    dev.off()
    data.melted.3 <- filter.taxa(DATA=data.melted, trt.code.1=NULL, trt.code.2=Tcompare3$trt.code.2[1], trt.refs=NULL, vars=c("taxon", "copies.ul", "Tube", "ComboTrt"), min.reps=3)
  #1-GR-Exu-18O:
    Tcompare3$trt.code.2[2]
    explore.filter.taxa(DATA=data.melted.3, trt.code.1=NULL, trt.code.2=Tcompare3$trt.code.2[2], trt.refs=NULL, vars=c("taxon", "copies.ul", "Tube", "ComboTrt"), min.reps=3)
    dev.off()
    data.melted.3 <- filter.taxa(DATA=data.melted.3, trt.code.1=NULL, trt.code.2=Tcompare3$trt.code.2[2], trt.refs=NULL, vars=c("taxon", "copies.ul", "Tube", "ComboTrt"), min.reps=3)
  #1-GR-Lit-18O:
    Tcompare3$trt.code.2[3]
    explore.filter.taxa(DATA=data.melted.3, trt.code.1=NULL, trt.code.2=Tcompare3$trt.code.2[3], trt.refs=NULL, vars=c("taxon", "copies.ul", "Tube", "ComboTrt"), min.reps=3)
    dev.off()
    data.melted.3 <- filter.taxa(DATA=data.melted.3, trt.code.1=NULL, trt.code.2=Tcompare3$trt.code.2[3], trt.refs=NULL, vars=c("taxon", "copies.ul", "Tube", "ComboTrt"), min.reps=3)
  #6-GR-NoC-18O:
    Tcompare3$trt.code.2[4]
    explore.filter.taxa(DATA=data.melted.3, trt.code.1=NULL, trt.code.2=Tcompare3$trt.code.2[4], trt.refs=NULL, vars=c("taxon", "copies.ul", "Tube", "ComboTrt"), min.reps=3)
    dev.off()
    data.melted.3 <- filter.taxa(DATA=data.melted.3, trt.code.1=NULL, trt.code.2=Tcompare3$trt.code.2[4], trt.refs=NULL, vars=c("taxon", "copies.ul", "Tube", "ComboTrt"), min.reps=3)
  #6-GR-Exu-18O:
    Tcompare3$trt.code.2[5]
    explore.filter.taxa(DATA=data.melted.3, trt.code.1=NULL, trt.code.2=Tcompare3$trt.code.2[5], trt.refs=NULL, vars=c("taxon", "copies.ul", "Tube", "ComboTrt"), min.reps=3)
    dev.off()
    data.melted.3 <- filter.taxa(DATA=data.melted.3, trt.code.1=NULL, trt.code.2=Tcompare3$trt.code.2[5], trt.refs=NULL, vars=c("taxon", "copies.ul", "Tube", "ComboTrt"), min.reps=3)
  #6-GR-Lit-18O:
    Tcompare3$trt.code.2[6]
    explore.filter.taxa(DATA=data.melted.3, trt.code.1=NULL, trt.code.2=Tcompare3$trt.code.2[6], trt.refs=NULL, vars=c("taxon", "copies.ul", "Tube", "ComboTrt"), min.reps=3)
    dev.off()
    data.melted.3 <- filter.taxa(DATA=data.melted.3, trt.code.1=NULL, trt.code.2=Tcompare3$trt.code.2[6], trt.refs=NULL, vars=c("taxon", "copies.ul", "Tube", "ComboTrt"), min.reps=3)
  #1-GR-NoC-16O; 1-GR-Exu-16O; 1-GR-Lit-16O; 6-GR-NoC-16O; 6-GR-Exu-16O; 6-GR-Lit-16O:
    Tcompare3$trt.code.1[1]
    Tcompare3$trt.refs[1]
    explore.filter.taxa(DATA=data.melted.3, trt.code.1=NULL, trt.code.2=Tcompare3$trt.code.1[1], trt.refs=Tcompare3$trt.refs[1], vars=c("taxon", "copies.ul", "Tube", "ComboTrt"), min.reps=3)
    dev.off()
    data.melted.3 <- filter.taxa(DATA=data.melted.3, trt.code.1=NULL, trt.code.2=Tcompare3$trt.code.1[1], trt.refs=Tcompare3$trt.refs[1], vars=c("taxon", "copies.ul", "Tube", "ComboTrt"), min.reps=3)


#Save workspace:
  save(list=ls(all.names=TRUE), file="qSIP_workspaces/BF_01/.RData", envir=.GlobalEnv)



