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


#NOT DONE BECAUSE THIS TAXA ID FILE IS FOR OTUs & DATA IS FOR L7:
# #Import raw data & taxomic information and format raw data for analysis:
#   #Read in taxonomic information
#   # modification: to make melting faster, reading in taxonomy in 1 column with '_' separating phylum, kingdom...
#     taxa.id <- read.table("qSIP_data/MH_AS_taxa_ID.txt", header=T, sep="\t", stringsAsFactors=T, na.strings="")
#     head(taxa.id)
#     names(taxa.id)
#   #The number serves as a unique identifier code for each taxon to the "taxa.id" data frame:
    
  #NOTE: treatment code names cannot contain spaces
    data <- read.table("qSIP_data/MH_AS_MCLP_2500_L7_data.txt", header=T, sep="\t", stringsAsFactors=T, check.names=F)
    dim(data)
    names(data)
    head(data)

#NOT NEEDED, TAXA NAMES ARE NOT PRECEDED BY 'X's IN THIS DATAFRAME:
# #Note that columns 11 and up have a unique identifier code for taxon to match that in the "taxa.id" data frame
# names(data)[11:dim(data)[2]]
# #Rename these 'taxonID' columns to just the numbers (i.e., get rid of the leading 'X')
# names(data)[11:dim(data)[2]] <- gsub(pattern="X(\\d+)", replacement="\\1", x=names(data)[11:dim(data)[2]], perl=TRUE)

#Add a column for the sum of proportional abundance of all taxa by fraction
data$sum.abundance <- rowSums(data[,11:ncol(data)])
data$sum.abundance


#NOT NEEDED AND ALSO THE COLUMN NUMBER IS OFF; SUGGEST DELETING THIS ALL TOGETHER:
# #Do not correct the proportional abundances of all taxa by scaling by their sum
# #(they don't sum to 1, but that is presumably because proportional abundances of rare taxa are not given (because they were removed at the bioinformatics stage); even though these rare taxa are included in the total 'qPCR.16S.copies.ul' & 'DNA.ng.ul'
# needs.correcting <- F
# if (needs.correcting) {
#     data[9:(ncol(data)-1)] <- data[9:(ncol(data)-1)]/data$sum.abundance
# }

#Calculate number of copies per uL, based on relative abundance and total number of copies per uL:
ncopies <- data$qPCR.16S.copies.ul*data[,11:(ncol(data)-1)]
ncopies <- cbind(data[,1:10], ncopies)  # add first 8 columns of data to ncopies
    dim(ncopies)
    head(ncopies)
    names(ncopies)
    
#Melt data into long format by tube, sample, tmt, rep, fraction, DNA conc, and density;
#Do this for copies.ul and for relative abundance. Merge these to into 1 masterfile: data.melted
ncopies.melted <- melt(ncopies, id=c("tube", "isotope", "sample","tmt", "combo.trt.code", "rep", "fraction", "density.g.ml", "DNA.ng.ul", "qPCR.16S.copies.ul"), measure.vars=as.character(1:403), variable.name="taxon", value.name="copies.ul")
rel.abundance.melted <- melt(data, id=c("tube", "isotope", "sample",  "tmt", "combo.trt.code", "rep", "fraction", "density.g.ml", "DNA.ng.ul", "qPCR.16S.copies.ul"), measure.vars=as.character(1:403), variable.name="taxon", value.name="rel.abundance")
data.melted <- merge(ncopies.melted, rel.abundance.melted)
    head(data.melted)
  
  
#MERGING TAXA NOT DONE BECAUSE TAXA DATA WERE NOT IMPORTED ABOVE; BUT DO REORDER THE DATAFRAME:
# #Merge taxa data and reorder data frame by taxon and SampleID and Fraction
#   data.melted <- merge(data.melted, taxa.id)
  data.melted <- data.melted[order(data.melted$taxon, data.melted$sample, data.melted$fraction),]
  row.names(data.melted) <- 1:dim(data.melted)[1]     #rename observations to be sequential
  head(data.melted)

  
  #NOTE THAT YOU COULD COMBINE ALL 16O TREATMENTS INTO 'trt.code.1' and in 'trt.code.refs', ALSO, quotes and end-of-lines added to TreatmentComparison Files:
  #Import data frame containing the treatment comparisons to perform (split data into threee chunks according to soil type (AN, BS, GR):
    Tcompare1 <- read.table("qSIP_data/MH_AS_TreatmentComparisons_A.txt", header=TRUE, sep="\t", colClasses=c("numeric","factor","character","character","character","numeric","character"))
    summary(Tcompare1)
    Tcompare1

    Tcompare2 <- read.table("qSIP_data/MH_AS_TreatmentComparisons_F.txt", header=TRUE, sep="\t", colClasses=c("numeric","factor","character","character","character","numeric","character"))
    summary(Tcompare2)
    Tcompare2

    Tcompare3 <- read.table("qSIP_data/MH_AS_TreatmentComparisons_N.txt", header=TRUE, sep="\t", colClasses=c("numeric","factor","character","character","character","numeric","character"))
    summary(Tcompare3)
    Tcompare3

    Tcompare4 <- read.table("qSIP_data/MH_AS_TreatmentComparisons_O.txt", header=TRUE, sep="\t", colClasses=c("numeric","factor","character","character","character","numeric","character"))
    summary(Tcompare4)
    Tcompare4

#Calculate corrected WADs for each tube and taxon:
  WAD.by.taxon <- WAD.by.taxon.func(X=data.melted, vars=c("taxon", "density.g.ml", "copies.ul", "tube", "combo.trt.code"))
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
system.time(unlab.WAD.corr.list <- find.unlabeled.correction(LIST=WAD.by.taxon, unlab.tmts=c("H1_A_16","H2_A_16","H1_F_16","H2_F_16", "H1_N_16","H2_N_16","H1_O_16","H2_O_16"), lab.tmts=c( "H1_A_18", "H2_A_18","H1_F_18","H2_F_18","H1_N_18","H2_N_18","H1_O_18","H2_O_18"), CI=0.90))
      #Look at the results:
        names(unlab.WAD.corr.list)                           #names of the two data frames in the output list
        unlab.WAD.corr.list$WAD.norm.fit.parms               #looking at the first object -- a data frame
        unlab.WAD.corr.list[[1]]                             #another way to look at the first object
        unlab.WAD.corr.list$corr.names                       #looking at the second object in the list -- a vector of the names of the corrected unlabeled replicates
        unlab.WAD.corr.list[[2]]                             #another way to look at the second object
        head(unlab.WAD.corr.list$WAD.table.corr)             #looking at the head of the third object -- a data frame
        head(unlab.WAD.corr.list[[3]])                       #another way to look at the head of the third object
        
        
#Calculate the shift in WAD for each of the labeled replicates (8 labeled treatments x 4 replicates = 32 total labeled replicates):
  #1_AN_Exu_18O:
  #Can also try method="td.abs.resid"
  system.time(R1.H1_A_18.td.pos.resid.list <- find.labeled.correction(LIST=unlab.WAD.corr.list, lab.replicate="R1.H1_A_18", lab.names=c("R1.H1_A_18", "R2.H1_A_18", "R3.H1_A_18", "R4.H1_A_18"), method="td.abs.resid", unlab.SD.percentile=0.70, lab.SD.percentile=0.70, min.num.nongrowers=10))
    #Find the best iteration:
    R1.H1_A_18.td.pos.resid.best.iteration.list <- select.best.iteration(putative.nongrower.metrics=R1.H1_A_18.td.pos.resid.list$putative.nongrower.metrics, quantile.threshold=0.10)
    find.labeled.correction.plot(find.labeled.correction.list=R1.H1_A_18.td.pos.resid.list, filename="MH_AS_R1_H1_A_18_td_pos_resid_plots", path="qSIP_output/Figures/", highlight=R1.H1_A_18.td.pos.resid.best.iteration.list$best.iteration, highlight.col="magenta")

        


#Save workspace:
  save(list=ls(all.names=TRUE), file="qSIP_workspaces/MH_AS_01/.RData", envir=.GlobalEnv)



