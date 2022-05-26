# This code imports the BH qSIP data sets and performs the basic calculations on the data


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


#Import raw data & format raw data for analysis:
  #Read in raw uclust sequencing and taxonomy output file obtained from Greg Caporaso & format it properly (separate import steps for the taxonomy/names and the data itself):
  #NOTE: the raw data are in a genus-level BIOM table (sample x taxon abundance table); the file contains all genera assigned by both taxonomy assignment methods, so there are many taxa that have abundances of 0 because they were only assigned in the RDP-assignment
  #NOTE: according to Greg, there were 790 genera (for the uclust-analysis) and 979 (for the RDP-analysis); the uclust assignments are represented in this data file
    #First, read in the raw data with column headings:
      Gdata1 <- read.table("qSIP_data/BH_table_mc4049_sorted_L6_sorted_and_filled.txt", header=T, sep="\t", stringsAsFactors=F)
      dim(Gdata1)
    #Transpose the data so that each taxon is a column:
      Gdata1 <- t(Gdata1)
      dim(Gdata1)
    #Extract and format the row names; these are the sample IDs:
      SampleID <- gsub(pattern="X(.+)", replacement="\\1", x=row.names(Gdata1)[2:dim(Gdata1)[1]], perl=TRUE)
    #Store the taxon names as a vector:
      column.names <- Gdata1[1,]
    #Second, read in the raw data without headers so that columns are formatted as numeric:
      Gdata2 <- read.table("qSIP_data/BH_table_mc4049_sorted_L6_sorted_and_filled.txt", skip=1, header=F, sep="\t", stringsAsFactors=F)
      dim(Gdata2)
    #Transpose the data so that each taxon is a column (skip the first column, since it contains the taxon names):
      Gdata2 <- t(Gdata2[,2:dim(Gdata2)[2]])
      dim(Gdata2)
      Gdata <- Gdata2
    #Append the column of sample IDs to the correctly formatted data and assign the taxon names as the names for all other columns:
      Gdata <- data.frame(SampleID, Gdata)
      names(Gdata)[2:dim(Gdata)[2]] <- column.names
      dim(Gdata)
    
  #Read in raw data (obtained from Bruce) to get the qPCR and tube-level data (skip first few lines, including header); taxonomic relative abundances are incomplete in this file:
  #NOTE: treatment code names cannot contain spaces
    data <- read.table("qSIP_data/BH_Dual_SIP_IRMS_Data.txt", skip=1, header=F, sep="\t", stringsAsFactors=T)
  #Name columns 1 through 5; columns 6 and up will be dropped below (they have incomplete taxonomic abundances):
    names(data)[1:5] <- c("SampleID", "trt.code", "DNA.ng.ul", "density.g.ml", "copy.number")
  
  #Remove a duplicate sample in Gdata to avoid introducing bias in calculating WAD:
    #According to Becky's email on Jan 29, 2015: "Looks like 5.20 (119.5.20) was just sequenced twice, but is the same sample, so either can be used."
    #Figure out which one is represented in 'data':
      Gdata[Gdata$SampleID %in% "5.20",]
      Gdata[Gdata$SampleID %in% "119.5.20",]
      data[data$SampleID == "5.2",]
      Gdata.5.20.abunds <- as.numeric(Gdata[Gdata$SampleID %in% "5.20", 2:dim(Gdata)[2]])
      Gdata.119.5.20.abunds <- as.numeric(Gdata[Gdata$SampleID %in% "119.5.20", 2:dim(Gdata)[2]])
      data.5.2.abunds <- as.numeric(data[data$SampleID == "5.2", 6:dim(data)[2]])
      #See which one has relative abundances matching those of sample "5.2" in data (exclude genera with abundances of zero):
      sum(round(data.5.2.abunds[data.5.2.abunds != 0], 5) %in% round(Gdata.5.20.abunds[Gdata.5.20.abunds != 0], 5))
      sum(round(data.5.2.abunds[data.5.2.abunds != 0], 5) %in% round(Gdata.119.5.20.abunds[Gdata.119.5.20.abunds != 0], 5))
    #Sample "119.5.20" (in Gdata) is the one that is represented in 'data', so keep that one:
      dim(Gdata)
      Gdata <- Gdata[Gdata$SampleID != "5.20",]
      Gdata$SampleID <- factor(Gdata$SampleID)
      row.names(Gdata) <- 1:dim(Gdata)[1]
      dim(Gdata)

  #Drop columns 6 and up (they have incomplete taxonomic abundances):
    data <- data[,1:5]
  
  #Fix SampleID codes in 'Gdata' to match those in 'data':
    #Pull out the correct SampleID numbers from the obscure codings:
      SampleID1 <- Gdata$SampleID
      SampleID2 <- SampleID3 <- character(length(SampleID1))
      for (i in 1:length(SampleID1)){
        if (grepl(pattern="\\d+\\.\\d+\\..+", x=SampleID1[i], perl=TRUE)){
          temp <- as.numeric(gsub(pattern="(\\d+)\\.\\d+\\..+", replacement="\\1", x=SampleID1[i], perl=TRUE))
          if (temp > 36){
            SampleID2[i] <- gsub(pattern="\\d+\\.(\\d+\\..+)", replacement="\\1", x=SampleID1[i], perl=TRUE)
          }
          if (temp <= 36){
            SampleID2[i] <- gsub(pattern="(\\d+\\.)\\d+\\.(.+)", replacement="\\1\\2", x=SampleID1[i], perl=TRUE)
          }
        }
        else {
          SampleID2[i] <- as.character(SampleID1)[i]
        }
      }
      for (i in 1:length(SampleID2)){
        if (grepl(pattern="\\d+\\.2\\..+", x=SampleID2[i], perl=TRUE)){
          SampleID3[i] <- gsub(pattern="(\\d+\\.)2\\.(.+)", replacement="\\1\\2", x=SampleID2[i], perl=TRUE)
        }
        else {
          SampleID3[i] <- as.character(SampleID2)[i]
        }
      }
      for (i in 1:length(SampleID3)){
        temp1 <- gsub(pattern="(\\d+\\.)\\d+", replacement="\\1", x=SampleID3[i], perl=TRUE)
        temp2 <- as.numeric(gsub(pattern="\\d+\\.(\\d+)", replacement="\\1", x=SampleID3[i], perl=TRUE))
        if (!is.na(temp2) & temp2 >= 7 & temp2 < 10){
          SampleID3[i] <- paste(temp1, "0", temp2, sep="")
        }
      }
      SampleID.new <- data.frame(SampleID1=SampleID1, SampleID2=SampleID2, SampleID3=as.numeric(as.character(SampleID3)), ID3=as.numeric(as.character(SampleID3)))
      SampleID.new

    Gdata$SampleID <- as.numeric(Gdata$SampleID)
    Gdata$SampleID <- SampleID.new$ID3
  
  #Identify the remaining samples in Gdata that were not represented in data:
    Gdata$SampleID[!is.element(Gdata$SampleID, data$SampleID)]

  #And limit Gdata to just those samples that were shared with data:
    Gdata <- Gdata[Gdata$SampleID %in% Gdata$SampleID[is.element(Gdata$SampleID, data$SampleID)],]
    Gdata$SampleID <- factor(Gdata$SampleID)
    row.names(Gdata) <- 1:dim(Gdata)[1]
    dim(Gdata)
    dim(data)
  
  #Combine Gdata and data (first sort by SampleID so that both dataframes are in the same order):
    Gdata <- Gdata[order(Gdata$SampleID),]
    data <- data[order(data$SampleID),]
    Gdata$SampleID == data$SampleID
    sum(Gdata$SampleID == data$SampleID)
    data <- data.frame(data, Gdata[,2:dim(Gdata)[2]])
    names(data)[6:dim(data)[2]] <- names(Gdata)[2:dim(Gdata)[2]]
  
  #Save taxon names to create a dataframe of numeric taxonomic codes and taxonomic classifications:
    taxa.names <- names(data)[6:ncol(data)]

  #Create a dataframe of numeric taxonomic codes and taxonomic classifications:
    length(taxa.names)
    taxa.id <- data.frame(taxon=seq(1, length(taxa.names)), code=taxa.names)
    taxa.id$kingdom <- factor(gsub(pattern="[k__]*(.*)\\;[p__]*(.*)\\;[c__]*(.*)\\;[o__]*(.*)\\;[f__]*(.*)\\;[g__]*(.*)", replacement="\\1", x=taxa.names, perl=TRUE))
    taxa.id$phylum <- factor(gsub(pattern="[k__]*(.*)\\;[p__]*(.*)\\;[c__]*(.*)\\;[o__]*(.*)\\;[f__]*(.*)\\;[g__]*(.*)", replacement="\\2", x=taxa.names, perl=TRUE))
    taxa.id$class <- factor(gsub(pattern="[k__]*(.*)\\;[p__]*(.*)\\;[c__]*(.*)\\;[o__]*(.*)\\;[f__]*(.*)\\;[g__]*(.*)", replacement="\\3", x=taxa.names, perl=TRUE))
    taxa.id$order <- factor(gsub(pattern="[k__]*(.*)\\;[p__]*(.*)\\;[c__]*(.*)\\;[o__]*(.*)\\;[f__]*(.*)\\;[g__]*(.*)", replacement="\\4", x=taxa.names, perl=TRUE))
    taxa.id$family <- factor(gsub(pattern="[k__]*(.*)\\;[p__]*(.*)\\;[c__]*(.*)\\;[o__]*(.*)\\;[f__]*(.*)\\;[g__]*(.*)", replacement="\\5", x=taxa.names, perl=TRUE))
    taxa.id$genus <- factor(gsub(pattern="[k__]*(.*)\\;[p__]*(.*)\\;[c__]*(.*)\\;[o__]*(.*)\\;[f__]*(.*)\\;[g__]*(.*)", replacement="\\6", x=taxa.names, perl=TRUE))
    #Replace blanks with NAs:
    taxa.id[taxa.id == ""] <- NA
      #Re-convert the factor columns to factor (to eliminate the level for the blanks):
      taxa.id$kingdom <- factor(taxa.id$kingdom)
      taxa.id$phylum <- factor(taxa.id$phylum)
      taxa.id$class <- factor(taxa.id$class)
      taxa.id$order <- factor(taxa.id$order)
      taxa.id$family <- factor(taxa.id$family)
      taxa.id$genus <- factor(taxa.id$genus)

  #Rename taxa columns in 'data'; columns 6 and up are given the unique numeric identifier code for taxon from the "taxa.id" data frame:
    names(data)[6:ncol(data)] <- taxa.id$taxon

  #Add a column for the sum of proportional abundance of all taxa by fraction:
    data$sum.abundance <- rowSums(data[,6:ncol(data)])

  #Do not correct the proportional abundances of all taxa by scaling by their sum 
  #(they already sum to 1 because proportional abundances of singeltons and unassigned OTUs are not included in this dataset (those were removed at the bioinformatics stage) even though the 16S sequences representing these excluded rare 'taxa' are included in the total 'copy.number'
    needs.correcting <- F
    if (needs.correcting) {
      data[6:(ncol(data)-1)] <- data[6:(ncol(data)-1)]/data$sum.abundance
    }

  #Calculate number of copies per uL, based on relative abundance and total number of copies per uL:
    ncopies <- data$copy.number*data[,6:(ncol(data)-1)]
    ncopies <- cbind(data[,1:5], ncopies)  # add first 5 columns of data to ncopies

  #Add tube ID to ncopies:
    dotpos <- regexpr("\\.[^\\.]*$", ncopies$SampleID)	# find ".", (number prior to "." is tube number)
    tube <- substr(ncopies$SampleID, 1, dotpos-1)
    ncopies <- cbind(tube, ncopies)

  #Melt data into long format by tube, trt.code, DNA conc, density and copy number;
  #Do this for ncopies and for relative abundance. Merge these to into 1 masterfile: data.melted
    ncopies.melted <- melt(ncopies, variable.name="taxon", id=c("tube","SampleID", "trt.code", "DNA.ng.ul", "density.g.ml", "copy.number"), value.name="copies")
    rel.abundance.melted <- melt(data, variable.name="taxon", id=c("SampleID", "trt.code", "DNA.ng.ul", "density.g.ml", "copy.number"), value.name="rel.abundance")
    data.melted <- merge(ncopies.melted, rel.abundance.melted)

  #Merge taxa data and reorder data frame by taxon and SampleID:
    data.melted <- merge(data.melted, taxa.id)
    data.melted <- data.melted[order(data.melted$taxon, data.melted$SampleID),]


#Import soil extraction data & calculate mass of soil in each tube:
  Sdat <- read.table("qSIP_data/BH_SoilExtractionData.txt", header=TRUE, sep="\t")
  Sdat$tube <- factor(Sdat$tube)
  Sdat$g.soil <- Sdat$ug.DNA.added.to.tube / Sdat$ug.DNA.g.soil
  summary(Sdat)


#Import data frame containing the treatment comparisons to perform:
  Tcompare <- read.table("qSIP_data/BH_TreatmentComparisons.txt", header=TRUE, sep="\t", colClasses=c("numeric","factor","character","character","character","numeric","character"))
  summary(Tcompare)
  Tcompare


#Filtering: keep only those genera that occurred in all tubes of the experiment:
  #First back up the complete data frame:
    data.melted.all <- data.melted
    
  #Subset the melted data to grab only those taxon-fractions with copies present:
    data.melted.occurrences <- data.melted[data.melted$copies > 0,]
      #Convert the factor columns to factor:
        data.melted.occurrences$taxon <- factor(data.melted.occurrences$taxon)
        data.melted.occurrences$SampleID <- factor(data.melted.occurrences$SampleID)
        data.melted.occurrences$trt.code <- factor(data.melted.occurrences$trt.code)
        data.melted.occurrences$tube <- factor(data.melted.occurrences$tube)
        data.melted.occurrences$kingdom <- factor(data.melted.occurrences$kingdom)
        data.melted.occurrences$phylum <- factor(data.melted.occurrences$phylum)
        data.melted.occurrences$class <- factor(data.melted.occurrences$class)
        data.melted.occurrences$order <- factor(data.melted.occurrences$order)
        data.melted.occurrences$family <- factor(data.melted.occurrences$family)
        data.melted.occurrences$genus <- factor(data.melted.occurrences$genus)
    dim(data.melted)
    dim(data.melted.occurrences)

  #Create a function to pass to tapply for calculating the number of unique tubes having copies for each taxon (genus):
    length.unique <- function(X){
      length(unique(X))
    }

  #Calculate the number of unique tubes with copies present for each taxon (genus):
    tubes.per.genus <- tapply(data.melted.occurrences$tube, data.melted.occurrences$taxon, length.unique)
    tubes.per.genus <- sort(tubes.per.genus)

  #Keep only those genera that occur in all 36 tubes of the experiment:
    dev.off()
    dev.new(width=3.5, height=7)
    par(mfrow=c(2,1))
    hist(tubes.per.genus, main="")
    abline(v=35.5, col="red")
    plot(x=1:length(tubes.per.genus), y=as.numeric(tubes.per.genus), xlab="Number of genera", ylab="Tubes per genus", type="l")
    abline(h=35.5, col="red")
    par(mfrow=c(1,1))
    # tubes.per.genus[as.numeric(tubes.per.genus) == 36]
    length(tubes.per.genus[as.numeric(tubes.per.genus) >= 36])

    dev.off()

  #Subset the data.melted dataframe so that it only contains taxon-fractions for genera occurring in all 36 tubes:
    data.melted <- data.melted[data.melted$taxon %in% names(tubes.per.genus[as.numeric(tubes.per.genus) == 36]),]
    row.names(data.melted) <- 1:dim(data.melted)[1]
      #Re-convert the factor columns to factor:
        data.melted$taxon <- factor(data.melted$taxon)
        data.melted$SampleID <- factor(data.melted$SampleID)
        data.melted$trt.code <- factor(data.melted$trt.code)
        data.melted$tube <- factor(data.melted$tube)
        data.melted$kingdom <- factor(data.melted$kingdom)
        data.melted$phylum <- factor(data.melted$phylum)
        data.melted$class <- factor(data.melted$class)
        data.melted$order <- factor(data.melted$order)
        data.melted$family <- factor(data.melted$family)
        data.melted$genus <- factor(data.melted$genus)
    dim(data.melted.all)
    dim(data.melted)
    

#Some examples of how to use the functions:
  #Create a subset of the data containing only data for taxon 104 in tube 2 of treatment '1_C_16O'
    T104T1T2 <- data.melted[data.melted$taxon==104 & data.melted$trt.code=="1_C_16O" & data.melted$tube=="2",]
  #Calculate the weighted average density for taxon 104 in tube 2:
    T104T1T2.wad.out <- WAD.func(y=T104T1T2$copies, x=T104T1T2$density.g.ml)
    T104T1T2.wad.out
  #Create a subset of the data containing only data for taxon 104 in all tubes of treatment '1_C_16O'
    T104T1 <- data.melted[data.melted$taxon==104 & data.melted$trt.code=="1_C_16O",]
  #Calculate the bootstrapped (and observed) weighted average density estimates for taxon 104 in treatment '1_C_16O'
    T104T1.boot.out <- boot.WAD.func(X=T104T1, vars=c("density.g.ml", "copies", "tube"), CI=0.90, draws=1000)
    T104T1.boot.out
  #Create a subset of the data containing only data for taxon 104 in all tubes of treatment '1_C_18O'
    T104T2 <- data.melted[data.melted$taxon==104 & data.melted$trt.code=="1_C_18O",]
  #Calculate the difference in weighted average density (with bootstrapped estimates of uncertainty) between the '1_C_18O' and '1_C_16O' treatments for taxon 104
    T104T1v2.out <- boot.diff.wad(T104T1, T104T2, vars=c("density.g.ml", "copies", "tube", "trt.code"), CI=0.90, draws=1000, tailed.test=1)
    T104T1v2.out
  #Create a subset of the data containing only data for taxon 104 in all tubes of the 'reference' treatments (i.e., '1_C_16O' & '1_P_12Cplus16O')
    T104ref <- data.melted[data.melted$taxon==104 & (data.melted$trt.code=="1_C_16O" | data.melted$trt.code=="1_P_12Cplus16O"),]
  #Calculate the molecular weight, GC content, and average number of carbon and oxygen atoms per DNA nucleotide for taxon 104:
    T104ref.MW.out <- MW.calc(X=T104ref, vars=c("density.g.ml", "copies", "tube"))
    T104ref.MW.out
  #Calculate the excess atom fraction of 18O (with bootstrapped estimates of uncertainty) in the '1_C_18O' vs. the '1_C_16O' treatment for taxon 104
    T104T1v2.ape.out <- boot.diff.ape(X.light=T104T1, X.heavy=T104T2, X.reference=T104ref, iso.compare="18O", vars=c("density.g.ml", "copies", "tube", "trt.code"), CI=0.90, draws=1000)
    T104T1v2.ape.out
  #Calculate the intrinsic rate of increase (r, with bootstrapped estimates of uncertainty) for taxon 104 as determined by 18O incorporation
    T104T1v2.r.out <- boot.diff.r(X.light=T104T1, X.heavy=T104T2, X.reference=T104ref, M.soil=Sdat, iso.compare="18O", days=7, vars=c("density.g.ml", "copies", "tube", "trt.code", "g.soil"), growth.model="exponential", prop.O.from.water=0.33, v.frac=50, CI=0.90, draws=1000)
    T104T1v2.r.out
  #Calculate the bootstrapped (and observed) weighted average density estimates for taxon 104 in treatment '1_C_16O' along with the corresponding values for soil mass and total 16S copies from the tubes chosen for each bootstrap sample
    T104T1.tube.out <- boot.TUBE.func(X=T104T1, M.soil=Sdat, vars=c("density.g.ml", "copies", "tube", "g.soil"), v.frac=50, CI=0.90, draws=1000)
    T104T1.tube.out
  #Calculate the flux of carbon into biomass (with bootstrapped estimates of uncertainty) for taxon 104 as determined by 18O incorporation
    T104T1v2.f.out <- boot.diff.f(X.light=T104T1, X.heavy=T104T2, X.reference=T104ref, M.soil=Sdat, iso.compare="18O", days=7, vars=c("density.g.ml", "copies", "tube", "trt.code", "g.soil"), growth.model="exponential", prop.O.from.water=0.33, v.frac=50, copies.cell=6, pgC.cell=0.1, CI=0.90, draws=1000)
    T104T1v2.f.out
  #Note that the 'message' in the output of most of the functions warns if a taxon was not detected in a tube(s), however after the filtering step applied above, this is never the case for this particular analysis:


#Run all wad.diff, ape, r, & flux calculations for all taxa and all comparisons:
  #NOTES: the reference treatments for the comparisons are restricted to those in the same week of the experiment
  #       using set.seed ensures that subsequent function calls utilize the same random number seed, so results are replicated exactly)
  #       using system.time calculates how long it took the function to run
  #       set prop.O.from.water=0.33 according to McHugh et al. unpublished data from E. coli pure cultures incubated with 18O-H2O
  #       rename the bootstrapped output files so that they are not overwritten if 'all.taxa.calcs' is re-run
  set.seed(100)
  system.time(all.comparisons <- all.taxa.calcs(X.all=data.melted, comparisons=Tcompare, M.soil=Sdat, vars=c("taxon", "density.g.ml", "copies", "tube", "trt.code", "g.soil"), growth.model="exponential", prop.O.from.water=0.33, v.frac=50, copies.cell=6, pgC.cell=0.1, CI=0.90, draws=1000, tailed.test=1))
  bootstrapped.filenames <- paste("qSIP_output/", c("bootstrapped_wad1.txt", "bootstrapped_wad2.txt", "bootstrapped_wad_diff.txt", "bootstrapped_ape.txt", "bootstrapped_r.txt", "bootstrapped_C_fluxes.txt", "bootstrapped_N1.txt", "bootstrapped_N2.txt", "bootstrapped_N.txt"), sep="")
  new.bootstrapped.filenames <- paste("qSIP_output/", "BH_", c("bootstrapped_wad1.txt", "bootstrapped_wad2.txt", "bootstrapped_wad_diff.txt", "bootstrapped_ape.txt", "bootstrapped_r.txt", "bootstrapped_C_fluxes.txt", "bootstrapped_N1.txt", "bootstrapped_N2.txt", "bootstrapped_N.txt"), sep="")
  file.rename(from=bootstrapped.filenames, to=new.bootstrapped.filenames)
  summary(all.comparisons)
  dim(all.comparisons)


#Write the results (all.comparisons) to a text file:
  dir.create(path=paste(getwd(), "/qSIP_output", sep=""), showWarnings=FALSE)
  write.table(all.comparisons, "qSIP_output/BH_all_comparisons.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)


#Write the taxa.id dataframe to a text file:
  write.table(taxa.id, "qSIP_output/BH_taxa_ID.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
