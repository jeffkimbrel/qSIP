# This code imports the TM qSIP data sets, formats them, and performs quality control measures on the data


graphics.off()	#close all graphics windows


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
  source("qSIP_repo/boot.pop.R")                        #boot.pop
  source("qSIP_repo/boot.TUBE.pop.R")                   #boot.TUBE.pop
  source("qSIP_repo/comparison.message.pop.R")          #comparison.message.pop
  source("qSIP_repo/r.calc.pop.R")                      #r.calc.pop
  source("qSIP_repo/f.calc.pop.R")                      #f.calc.pop
  source("qSIP_repo/boot.r.pop.R")                      #boot.r.pop
  source("qSIP_repo/boot.f.pop.R")                      #boot.f.pop
  source("qSIP_repo/all.taxa.calcs.pop.R")              #all.taxa.calcs.pop
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


#Read in lab qSIP data (obtained from Theresa) to get the qPCR and tube-level information:
  #NOTE: treatment code names cannot contain spaces
  labdata <- read.table("qSIP_data/TM_qSIP_data.txt", header=T, sep="\t", stringsAsFactors=T)
    #Create tube IDs and add fraction numbers to 'data'
      #Assign tube numbers according to replicate ID and treatment ID:
        tube.key <- data.frame(tube=1:15, rep=rep(c(1,2,3,4,5), each=3), tmt=rep(c("Time0", "16O", "18O"), times=5))
      tube <- numeric(dim(labdata)[1])
      for (i in 1:dim(labdata)[1]){
        tube[i] <- tube.key$tube[tube.key$rep == labdata$rep[i] & tube.key$tmt == labdata$tmt[i]]
      }
      fraction <- gsub(pattern="^([^\\-])$", replacement=NA, x=labdata$sample, perl=TRUE)      #samples without a '-' in thier 'sample' name get an 'NA' for 'fraction'
      fraction <- gsub(pattern="(.*)(\\d)\\-(.+)", replacement="\\3", x=fraction, perl=TRUE)   #all other samples get the last digits in their 'sample' name for 'fraction'
      fraction <- as.numeric(as.character(fraction))
      labdata <- data.frame(tube, fraction, labdata)
      summary(labdata)


#Import soil extraction data & calculate mass of soil represented by the DNA in each tube:
  Sdat <- read.table("qSIP_data/TM_SoilExtractionData.txt", header=TRUE, sep="\t")
  #Add tube IDs, rep, and tmt to 'Sdat':
    tmt <- gsub(pattern="^([^\\-])$", replacement="Time0", x=Sdat$sample, perl=TRUE)           #samples without a '-' in thier 'sample' name get are 'Time0' treatments
    tmt <- gsub(pattern="(.*)\\-(\\d+)", replacement="\\1", x=tmt, perl=TRUE)                  #all other samples get the prefix in 'sample' name as their treatment code
    tmt <- factor(as.character(tmt))
    rep <- gsub(pattern="(.*)\\-(\\d+)", replacement="\\2", x=Sdat$sample, perl=TRUE)          #the second number in the 'sample' name is the rep; samples without a '-' in thier 'sample' name get that full string for 'rep'
    rep <- as.numeric(as.character(rep))
    Sdat <- data.frame(rep, tmt, Sdat)
    tube <- numeric(dim(Sdat)[1])
    for (i in 1:dim(Sdat)[1]){
      tube[i] <- tube.key$tube[tube.key$rep == Sdat$rep[i] & tube.key$tmt == Sdat$tmt[i]]
    }
    Sdat <- data.frame(tube, Sdat)
  #Calculate mass of soil represented by the DNA in each tube:
    Sdat$g.soil <- (Sdat$mean.qubit.ng.ul * Sdat$ul.used.to.make.800.ng) * (1/Sdat$total.ng.DNA) * (Sdat$g.soil.extracted)
    Sdat$g.soil - Sdat$g.soil.tube
    Sdat
    summary(Sdat)


#Read in the mapping file (obtained from Michaela Hayer) to translate the sample codes in the sequence data into meaningful treatment, tube, and fraction information:
  map <- read.table("qSIP_data/TM_mapping_file.txt", header=T, sep="\t", comment.char="", stringsAsFactors=F)
  #Edit the mapping file so that variable names and codes match those in the lab qSIP data:
    names(map) <- c("SampleID", "sample.rep.fraction", "tmt", "description")
    map$tmt[map$tmt == "C"] <- "Time0"
    map$sample.rep.fraction <- factor(map$sample.rep.fraction)
    map$tmt <- factor(map$tmt)
    map$description <- factor(map$description)
    #Add tube IDs, rep, and fraction numbers to 'map'
      rep <- gsub(pattern="^(\\d+)\\.(\\d+)\\.(\\d+)$", replacement="\\2", x=map$sample.rep.fraction, perl=TRUE)   #the second number in the 'sample' name is the rep; samples without a '.' in thier 'sample' name get that full string for 'rep'
      rep <- as.numeric(as.character(rep))
      map <- data.frame(map, rep)
      tube <- numeric(dim(map)[1])
      for (i in 1:dim(map)[1]){
        tube[i] <- tube.key$tube[tube.key$rep == map$rep[i] & tube.key$tmt == map$tmt[i]]
      }
      fraction <- gsub(pattern="^([^\\.]$)", replacement=NA, x=map$sample.rep.fraction, perl=TRUE)                 #samples without a '.' in thier 'sample' name get an 'NA' for 'fraction'
      fraction <- gsub(pattern="(.*)(\\d)\\.(.+)", replacement="\\3", x=fraction, perl=TRUE)                       #all other samples get the last digits in their 'sample' name for 'fraction'
      fraction <- as.numeric(as.character(fraction))
      map <- data.frame(map, tube, fraction)
      summary(map)
  

#Import raw sequencing data & save taxon names to create a dataframe of taxonomic codes and data:
  #NOTE: there are 5 versions of the sequencing data:
  # 1) 'Gdata'  genus-level (L6) BIOM table (sample x taxon abundance table) with relative abundances calculated after filtering out taxa using the 0.005% threshold
  # 2) 'Sdata'  species-level (L7) BIOM table (sample x taxon abundance table) with relative abundances calculated after filtering out taxa using the 0.005% threshold
  # 3) 'O5data'  OTU-level BIOM table (sample x taxon abundance table) with raw sequence counts; this dataset was filtered to remove OTUs below the 0.005% threshold
  # 4) 'O12data'  OTU-level BIOM table (sample x taxon abundance table) with raw sequence counts; this dataset was filtered to remove singleton and doubleton OTUs; note that 15 samples (#121-132,147,149,178) from a separate experiment are included in this data file becasue they were sequenced in the same run; those samples are not in the mapping file; also because this data was not filtered at the 0.005% threshold, it contains 5 samples (#24,106,108,110,115) that do not appear in any of the 0.005%-filtered datasets because those samples had very few sequence reads and were therefore filtered out when the 0.005% threshold was applied
  # 5) 'O1data'  OTU-level BIOM table (sample x taxon abundance table) with raw sequence counts; this dataset was filtered to remove only singleton OTUs; note that 15 samples (#121-132,147,149,178) from a separate experiment are included in this data file becasue they were sequenced in the same run; those samples are not in the mapping file; also because this data was not filtered at the 0.005% threshold, it contains 5 samples (#24,106,108,110,115) that do not appear in any of the 0.005%-filtered datasets because those samples had very few sequence reads and were therefore filtered out when the 0.005% threshold was applied
  
  #Genus-level relative abundance data after 0.005% filtering:
  #Read in raw sequencing and taxonomy output file obtained from Michaela Hayer & format it properly (separate import steps for the taxonomy/names and the data itself):
    #First, read in the raw data with column headings:
      Gdata1 <- read.table("qSIP_data/TM_uc_97_non_chimeras_otu_table_no_low_sample_filtered_.005%filtered_L6.txt", header=T, sep="\t", stringsAsFactors=F)
      dim(Gdata1)
    #Transpose the data so that each taxon is a column:
      Gdata1 <- t(Gdata1)
      dim(Gdata1)
    #Extract and format the row names; these are the sample IDs:
      SampleID <- gsub(pattern="X(.+)", replacement="\\1", x=row.names(Gdata1)[2:dim(Gdata1)[1]], perl=TRUE)
    #Store the taxon names as a vector:
      taxa.names <- Gdata1[1,]
    #Second, read in the raw data without headers so that columns are formatted as numeric:
      Gdata2 <- read.table("qSIP_data/TM_uc_97_non_chimeras_otu_table_no_low_sample_filtered_.005%filtered_L6.txt", skip=1, header=F, sep="\t", stringsAsFactors=F)
      dim(Gdata2)
    #Transpose the data so that each taxon is a column (skip the first column, since it contains the taxon names):
      Gdata2 <- t(Gdata2[,2:dim(Gdata2)[2]])
      dim(Gdata2)
      Gdata <- Gdata2
    #Append the column of sample IDs to the correctly formatted data and assign the taxon names as the names for all other columns:
      Gdata <- data.frame(SampleID, Gdata)
      names(Gdata)[2:dim(Gdata)[2]] <- taxa.names
      Gdata$SampleID <- as.numeric(as.character(Gdata$SampleID))
      row.names(Gdata) <- 1:dim(Gdata)[1]
      dim(Gdata)

  #Create a dataframe of numeric taxonomic codes and taxonomic classifications:
    length(taxa.names)
    taxa.id <- data.frame(taxon=seq(1, length(taxa.names)), code=taxa.names)
    taxa.id$kingdom <- factor(gsub(pattern="[k__]*(.*)\\;[p__]*(.*)\\;[c__]*(.*)\\;[o__]*(.*)\\;[f__]*(.*)\\;[g__]*(.*)", replacement="\\1", x=taxa.names, perl=TRUE))
    taxa.id$phylum <- factor(gsub(pattern="[k__]*(.*)\\;[p__]*(.*)\\;[c__]*(.*)\\;[o__]*(.*)\\;[f__]*(.*)\\;[g__]*(.*)", replacement="\\2", x=taxa.names, perl=TRUE))
    taxa.id$class <- factor(gsub(pattern="[k__]*(.*)\\;[p__]*(.*)\\;[c__]*(.*)\\;[o__]*(.*)\\;[f__]*(.*)\\;[g__]*(.*)", replacement="\\3", x=taxa.names, perl=TRUE))
    taxa.id$order <- factor(gsub(pattern="[k__]*(.*)\\;[p__]*(.*)\\;[c__]*(.*)\\;[o__]*(.*)\\;[f__]*(.*)\\;[g__]*(.*)", replacement="\\4", x=taxa.names, perl=TRUE))
    taxa.id$family <- factor(gsub(pattern="[k__]*(.*)\\;[p__]*(.*)\\;[c__]*(.*)\\;[o__]*(.*)\\;[f__]*(.*)\\;[g__]*(.*)", replacement="\\5", x=taxa.names, perl=TRUE))
    taxa.id$genus <- factor(gsub(pattern="[k__]*(.*)\\;[p__]*(.*)\\;[c__]*(.*)\\;[o__]*(.*)\\;[f__]*(.*)\\;[g__]*(.*)", replacement="\\6", x=taxa.names, perl=TRUE))
    #Replace blank cells in 'taxa.id' with NAs:
    taxa.id[taxa.id == ""] <- NA
      #Re-convert the factor columns to factor (to eliminate the level for the blanks):
      taxa.id$kingdom <- factor(taxa.id$kingdom)
      taxa.id$phylum <- factor(taxa.id$phylum)
      taxa.id$class <- factor(taxa.id$class)
      taxa.id$order <- factor(taxa.id$order)
      taxa.id$family <- factor(taxa.id$family)
      taxa.id$genus <- factor(taxa.id$genus)

  #Rename taxa columns in 'data'; columns 2 and up are given the unique numeric identifier code for taxon from the "taxa.id" data frame:
    names(Gdata)[2:ncol(Gdata)] <- taxa.id$taxon


#Combine the lab qSIP data and the sequencing data into one master data file by using the mapping file:
    data1 <- merge(map, Gdata, all.x=TRUE)
    row.names(data1) <- 1:dim(data1)[1]
    data.all <- merge(labdata, data1, all.x=TRUE)
    row.names(data.all) <- 1:dim(data.all)[1]


##########{___________Quality screening of the data: identify & exclude density outliers_____________________________

  #Identify & exclude outlier density fractions that should be excluded because water mixed in with the CsCl at the light (top) end:
    par(mfrow=c(3,2))
    for (i in 1:length(levels(factor(as.character(data.all$rep))))){
      plot(x=data.all$density.g.ml[data.all$tmt == "16O" & data.all$rep == i], y=data.all$DNA.ng.ul[data.all$tmt == "16O" & data.all$rep == i], pch=21, bg="blue", col="black", xlab="density (g/ml)", ylab="[DNA] (ng/ul)", main=paste("rep ", i, sep=""))
      points(x=data.all$density.g.ml[data.all$tmt == "18O" & data.all$rep == i], y=data.all$DNA.ng.ul[data.all$tmt == "18O" & data.all$rep == i], pch=21, bg="red", col="black", cex=0.6)
      abline(v=1.65) 
    }
    par(mfrow=c(1,1))
    #View on a shorter axis as DNA concentration:
      par(mfrow=c(3,2))
      for (i in 1:length(levels(factor(as.character(data.all$rep))))){
        plot(x=data.all$density.g.ml[data.all$tmt == "16O" & data.all$rep == i], y=data.all$DNA.ng.ul[data.all$tmt == "16O" & data.all$rep == i], pch=21, bg="blue", col="black", xlim=c(1.575, max(data.all$density.g.ml, na.rm=TRUE)), xlab="density (g/ml)", ylab="[DNA] (ng/ul)", main=paste("rep ", i, sep=""))
        points(x=data.all$density.g.ml[data.all$tmt == "18O" & data.all$rep == i], y=data.all$DNA.ng.ul[data.all$tmt == "18O" & data.all$rep == i], pch=21, bg="red", col="black", cex=0.6)
        abline(v=1.65) 
      }
      par(mfrow=c(1,1))
    #View on a shorter axis as number of 16S copies:
      par(mfrow=c(3,2))
      for (i in 1:length(levels(factor(as.character(data.all$rep))))){
        plot(x=data.all$density.g.ml[data.all$tmt == "16O" & data.all$rep == i], y=data.all$qPCR.16S.copies.ul[data.all$tmt == "16O" & data.all$rep == i], pch=21, bg="blue", col="black", xlim=c(1.575, max(data.all$density.g.ml, na.rm=TRUE)), xlab="density (g/ml)", ylab="copies/ul", main=paste("rep ", i, sep=""))
        points(x=data.all$density.g.ml[data.all$tmt == "18O" & data.all$rep == i], y=data.all$qPCR.16S.copies.ul[data.all$tmt == "18O" & data.all$rep == i], pch=21, bg="red", col="black", cex=0.6)
        abline(v=1.65) 
      }
      par(mfrow=c(1,1))

      graphics.off()

    #Exclude fractions that are less than 1.65 g/ml:
      data <- data.all[data.all$density.g.ml >= 1.65 | is.na(data.all$density.g.ml),]

##########}__________________________________________________________________________________________________________


##########{___________Calculate 16S copies per uL & 'melt' data into long format; merge in taxa data_________________
  
  #Add a column for the sum of proportional abundance of all taxa by fraction
    data$sum.abundance <- rowSums(data[,17:ncol(data)])

  #Do not correct the proportional abundances of all taxa by scaling by their sum because all 'true' taxa are included and they already sum to 1
  #OTUs were filtered to include only those ≥0.005% of total sequences.
  #The excluded OTUs can be considered erroneous and therefore it is proper to exclude them before calculating relative abundances, as was done in this case.
    data$sum.abundance

  #Calculate number of copies per uL, based on relative abundance and total number of copies per uL
    ncopies <- data$qPCR.16S.copies.ul*data[,17:(ncol(data)-1)]
    ncopies <- cbind(data[,1:16], ncopies)  # add first 16 columns of data to ncopies

  #Melt data into long format by rep, tmt, density , DNA.ng.ul, etc.;
  #Do this for ncopies and for relative abundance, then merge these to into 1 masterfile: data.melted
    ncopies.melted <- melt(data=ncopies, id.vars=c("tube", "fraction", "rep", "tmt", "sample", "nD", "density.g.ml", "DNA.ng.ul", "DNA.percent.tube", "qPCR.16S.copies.ul", "qPCR.DNA.ng.ul", "ng.DNA.g.soil.calc", "copies.g.soil.calc", "SampleID"), measure.vars=as.character(1:364), variable.name="taxon", value.name="copies")
    rel.abundance.melted <- melt(data=data, id.vars=c("tube", "fraction", "rep", "tmt", "sample", "nD", "density.g.ml", "DNA.ng.ul", "DNA.percent.tube", "qPCR.16S.copies.ul", "qPCR.DNA.ng.ul", "ng.DNA.g.soil.calc", "copies.g.soil.calc", "SampleID"), measure.vars=as.character(1:364), variable.name="taxon", value.name="rel.abundance")
    data.melted <- merge(ncopies.melted, rel.abundance.melted)

  #Merge taxa data (not including the 'code column') and reorder data frame by taxon, tube (rep & treatment), and fraction:
    data.melted <- merge(data.melted, taxa.id[,c(1,3:8)])
    data.melted <- data.melted[order(data.melted$taxon, data.melted$tube, data.melted$fraction),]
    rownames(data.melted) <- 1:dim(data.melted)[1]

##########}__________________________________________________________________________________________________________


##########{___________Calculate tube-level abundance of 16S genes for all replicates_________________________________
  
#Create a tube-level 16S dataset using tube-level data from the Time0 samples and by creating analagous tube-level data from the fraction-based dataset (assume tube = sum of fractions; see the check on this assumption below):
  #Create an empty data frame:
    ncopies.tube <- data.frame(matrix(NA, nrow=length(levels(factor(ncopies$rep))) * length(levels(ncopies$tmt)), ncol=5+dim(ncopies)[2]-16))
    names(ncopies.tube) <- c("tube", "rep", "tmt", "DNA.ng.ul", "qPCR.16S.copies.ul", names(ncopies)[17:dim(ncopies)[2]])
    # ncopies.tube[,] <- NA

  #Calculate number of copies per uL:
    for (k in unique(ncopies$tube)){
      ncopies.tube$tube[k] <- k
      if (unique(ncopies$tmt[ncopies$tube == k]) == "Time0"){
        ncopies.tube$rep[k] <- unique(ncopies$rep[ncopies$tube == k & is.na(ncopies$fraction)])
        ncopies.tube$tmt[k] <- as.character(unique(ncopies$tmt[ncopies$tube == k & is.na(ncopies$fraction)]))
        ncopies.tube$qPCR.16S.copies.ul[k] <- unique(ncopies$qPCR.16S.copies.ul[ncopies$tube == k & is.na(ncopies$fraction)])
        ncopies.tube$DNA.ng.ul[k] <- unique(ncopies$DNA.ng.ul[ncopies$tube == k & is.na(ncopies$fraction)])
        ncopies.tube[k,6:dim(ncopies.tube)[2]] <- apply(ncopies[ncopies$tube == k & is.na(ncopies$fraction), 17:dim(ncopies)[2]], 2, sum, na.rm=TRUE)
      }
      else {
        ncopies.tube$rep[k] <- unique(ncopies$rep[ncopies$tube == k])
        ncopies.tube$tmt[k] <- as.character(unique(ncopies$tmt[ncopies$tube == k]))
        ncopies.tube[k,6:dim(ncopies.tube)[2]] <- apply(ncopies[ncopies$tube == k, 17:dim(ncopies)[2]], 2, sum, na.rm=TRUE)
        ncopies.tube$qPCR.16S.copies.ul[k] <- sum(ncopies.tube[k,6:dim(ncopies.tube)[2]])
        ncopies.tube$DNA.ng.ul[k] <- Sdat$mean.qubit.ng.ul[Sdat$tube == k]
      }
    }
    ncopies.tube$tmt <- factor(ncopies.tube$tmt)
    ncopies.tube[,1:10]
    
  #Check that the sum of copies/uL across all taxa is the same as the tube-level estimate of total copies/uL:
    data.frame(ncopies.tube$qPCR.16S.copies.ul, apply(ncopies.tube[,6:dim(ncopies.tube)[2]], 1, sum))

  #Melt tube-level 16S data into long format by tube, rep, tmt:
    ncopies.tube.melted <- melt(data=ncopies.tube, id.vars=c("tube", "rep", "tmt"), measure.vars=as.character(1:364), variable.name="taxon", value.name="copies")

  #Merge with taxa data (not including the 'code column') and reorder data frame by taxon and tube (rep & treatment):
    ncopies.tube.melted <- merge(ncopies.tube.melted, taxa.id[,c(1,3:8)])
    ncopies.tube.melted <- ncopies.tube.melted[order(ncopies.tube.melted$taxon, ncopies.tube.melted$tube),]
    rownames(ncopies.tube.melted) <- 1:dim(ncopies.tube.melted)[1]

##########}__________________________________________________________________________________________________________


##########{___________Quality screening of the data: exclude replicates likely to have inadequate 16S data___________
      
  #Remove observations (fractions) where 16S copies were not measured:
    #First look at all the data and note that tubes #6 & #14 did not have any unmeasured fractions:
    tapply(data.melted$copies, INDEX=list(data.melted$tube, data.melted$tmt), function(x) sum(is.na(x)))
    #Collect those fractions where 16S copies were not measured to assess the proportion of total DNA represented in those fractions
    data.melted.no.16S <- data.melted[is.na(data.melted$copies),]
    data.melted <- data.melted[!is.na(data.melted$copies),]
    #Proportion of DNA (for each tube and taxon) that was in fractions lacking 16S copy data:
    #(so long as these proportions are low, there shouldn't be a problem estimating total 16S copies per tube based on the sum of fractions with 16S copies)
    #(Note again that tubes #6 & #14 did not have any unmeasured fractions, therefore we need to create rows for these tubes below)
    table(data.melted.no.16S$tube, data.melted.no.16S$tmt)
    prop.DNA.missing <- tapply(data.melted.no.16S$DNA.percent.tube, INDEX=list(data.melted.no.16S$tube, data.melted.no.16S$taxon), sum)
    prop.DNA.missing <- rbind(prop.DNA.missing[1:5,], rep(0, dim(prop.DNA.missing)[2]), prop.DNA.missing[6:12,], rep(0, dim(prop.DNA.missing)[2]), prop.DNA.missing[13,])
    dimnames(prop.DNA.missing)[[1]][c(6,14:15)] <- c("6", "14", "15")
    #Look at the proportion of DNA (for each tube and taxon) that was in fractions lacking 16S copy data
    #remember, there is no taxon-level variation (because DNA is quantified in bulk); therefore look at the values per tube across all taxa:
    apply(prop.DNA.missing, 1, max)
    tube.key
    #For all tubes except 5 & 9 & 15 (16O-rep2 & 18O-rep3 & 18O-rep5), there is ≤1% of total DNA in fractions that lack 16S data, indicating that for those tubes, the sum of 16S copies across all fractions should adequately represent tube-level 16S abundances
    #Tube 5 (16O-rep2) has 17% of total DNA in missing 16S fractions
    #Tube 9 (18O-rep3) has 10% of total DNA in missing 16S fractions
    #Tube 15 (18O-rep5) has 24% of total DNA in missing 16S fractions
    #(Note that tubes 1,4,7,10,13 are the Time0 treatments and lack 16S copy data only because individual fractions were not sequenced)

##########}__________________________________________________________________________________________________________


#BASED ON THE ABOVE ANALYSIS, DROP TUBES 5, 9, AND 15 (16O-REP2, 18O-REP3, & 18O-REP5) FROM FURTHER ANALYSIS:
  #qSIP data:
    data.melted <- data.melted[!is.element(data.melted$tube, c(5,9,15)),]
    #(reorder data frame by taxon, tube (rep & treatment), and fraction)
    data.melted <- data.melted[order(data.melted$taxon, data.melted$tube, data.melted$fraction),]
    rownames(data.melted) <- 1:dim(data.melted)[1]
    #remove any empty factor levels:
    data.melted$taxon <- factor(data.melted$taxon)
    data.melted$tmt <- factor(data.melted$tmt)
    data.melted$sample <- factor(data.melted$sample)
    data.melted$kingdom <- factor(data.melted$kingdom)
    data.melted$phylum <- factor(data.melted$phylum)
    data.melted$class <- factor(data.melted$class)
    data.melted$order <- factor(data.melted$order)
    data.melted$family <- factor(data.melted$family)
    data.melted$genus <- factor(data.melted$genus)

  #tube-level abundance data:
    ncopies.tube.melted <- ncopies.tube.melted[!is.element(ncopies.tube.melted$tube, c(5,9,15)),]
    #(reorder data frame by taxon and tube (rep & treatment))
    ncopies.tube.melted <- ncopies.tube.melted[order(ncopies.tube.melted$taxon, ncopies.tube.melted$tube),]
    rownames(ncopies.tube.melted) <- 1:dim(ncopies.tube.melted)[1]
    #remove any empty factor levels:
    ncopies.tube.melted$taxon <- factor(ncopies.tube.melted$taxon)
    ncopies.tube.melted$tmt <- factor(ncopies.tube.melted$tmt)
    ncopies.tube.melted$kingdom <- factor(ncopies.tube.melted$kingdom)
    ncopies.tube.melted$phylum <- factor(ncopies.tube.melted$phylum)
    ncopies.tube.melted$class <- factor(ncopies.tube.melted$class)
    ncopies.tube.melted$order <- factor(ncopies.tube.melted$order)
    ncopies.tube.melted$family <- factor(ncopies.tube.melted$family)
    ncopies.tube.melted$genus <- factor(ncopies.tube.melted$genus)


  #Note that it could be defensible to do this screening-out step above differently.
    #Here are some possibilities:
      #Exclude: 0 (leaves 5 18O, 5 16O, 5 Time0)
      #         2 (leaves 4 18O, 4 16O, 5 Time0) (uses a 15% cutoff)
      #         3 (leaves 3 180, 4 16O, 5 Time0) (uses a 10% cutoff (effectively a 1% cutoff in this case))
      #If paired bootstrapping is implemented (and it probably should be on this data set because total copies per tube seem to vary mostly among reps, not independently among tubes or isotope treatments), then we should probably go with the option of not excluding any tubes at this step based on inadequate 16S data.
      #Then, we would filter rare taxa such that a taxon has to occur in at least 3 tubes of each treatment and after meeting that criterion, then it must also occur in all 3 treatments within a rep (to allow paired calculations to be performed); if not, it is filtered out

  #Depending on which of these options is chosen, the filtering of rare taxa may need to be done differently below.
    #Here are some options for the case where paired bootstrapping is NOT implemented:
      #If 0 tubes are excluded in the step above, then impose one of these critera: a taxon must be in ≥13 of 15 total tubes  OR  must be in ≥3 reps of each treatment;  (these have the net result of 13of15 OR 9of15)
      #If 2 tubes are excluded in the step above, then impose one of these critera: a taxon must be in ≥11 of 13 remaining tubes  OR  must be in ≥3 of the 4 remaining 18O reps AND in ≥3 of the 4 remaining 16O reps AND in ≥3or4 of the 5 Time0 reps;  (these have the net result of 11of13 OR 9or10of13)
      #If 3 tubes are excluded in the step above, then impose one of these critera: a taxon must be in all 3 remaining 18O reps; in ≥3 of the 4 remaining 16O reps; in ≥3or4 of the 5 Time0 reps;  (these have the net result of 9or10of15)
    #Below, I have elected to use the simplest, most straightforward criterion (one that can also be implemented regardless of the screening option chosen above) that a taxon must occur in at least 3 reps (post-screening) of all three treatments (18O, 16O, Time0).  It seems to work well so far; if it does not, then I will revisit this filtering step later.


##########{___________Filter out rare taxa that do not occur in a sufficient number of replicates____________________

#Filter the data so that only taxa with an appropriate level of occurrence and replication among the tubes and treatments of the experiment are retained for further analysis:
  #{
  #First back up the complete data frames:
    data.melted.all <- data.melted
    ncopies.tube.melted.all <- ncopies.tube.melted

  #Keep only those taxa that occurred in nearly all tubes of the experiment:
    #Because two tubes of the 18O treatment and one tube of the 16O treatment were already removed from analysis, do the filtering like this (3 STEPS):
    #To remain in the analysis, a taxon must occur in all three remaining 18O tubes.
    #That taxon must also occur in at least 3 of the 4 remaining 16O tubes and in 3 of the 5 Time0 tubes.
  
    #First, we'll determine which taxa meet those criteria, using the 'ncopies.tube.melted' data frame, 
    #because it contains tube-level taxon abundances across all treatments (Time0 included)
    #then, we'll also keep only those same taxa in the 'data.melted' data frame


    #STEP 1:
      #Subset the data into only those taxon-tubes with copies present in the 18O treatment:
      ncopies.tube.melted.occurrences.18O <- ncopies.tube.melted[ncopies.tube.melted$copies > 0 & ncopies.tube.melted$tmt == "18O",]
          #Convert the factor columns to factor:
            ncopies.tube.melted.occurrences.18O$taxon <- factor(ncopies.tube.melted.occurrences.18O$taxon)
            ncopies.tube.melted.occurrences.18O$kingdom <- factor(ncopies.tube.melted.occurrences.18O$kingdom)
            ncopies.tube.melted.occurrences.18O$phylum <- factor(ncopies.tube.melted.occurrences.18O$phylum)
            ncopies.tube.melted.occurrences.18O$class <- factor(ncopies.tube.melted.occurrences.18O$class)
            ncopies.tube.melted.occurrences.18O$order <- factor(ncopies.tube.melted.occurrences.18O$order)
            ncopies.tube.melted.occurrences.18O$family <- factor(ncopies.tube.melted.occurrences.18O$family)
            ncopies.tube.melted.occurrences.18O$genus <- factor(ncopies.tube.melted.occurrences.18O$genus)
        dim(ncopies.tube.melted)
        dim(ncopies.tube.melted.occurrences.18O)

      #Create a function to pass to tapply for calculating the number of unique tubes having copies for each taxon:
        length.unique <- function(X){
          length(unique(X))
        }

      #Calculate the number of unique tubes in the 18O treatment with copies present for each taxon:
        tubes.per.taxon.18O <- tapply(ncopies.tube.melted.occurrences.18O$tube, ncopies.tube.melted.occurrences.18O$taxon, length.unique)
        tubes.per.taxon.18O <- sort(tubes.per.taxon.18O)
  
      #Keep only those taxa that occur in all three of the 3 remaining tubes of the 18O treatment:
        dev.off()
        dev.new(width=3.5, height=7)
        par(mfrow=c(2,1))
        hist(tubes.per.taxon.18O, xlab="# of tubes per taxon", ylab="# of taxa", main="")
        abline(v=2.75, col="red")
        plot(x=1:length(tubes.per.taxon.18O), y=as.numeric(tubes.per.taxon.18O), xlab="taxon (sorted by # tubes per taxon)", ylab="# of tubes per taxon", type="l")
        abline(h=2.75, col="red")
        par(mfrow=c(1,1))
        length(tubes.per.taxon.18O)   #total number of taxa occuring in the 18O treatment (after removing 2 'bad' tubes above)
        # tubes.per.taxon.18O[as.numeric(tubes.per.taxon.18O) >= 3]
        length(tubes.per.taxon.18O[as.numeric(tubes.per.taxon.18O) >= 3])   #number of taxa that occurred in all 3 of the 3 remaining tubes of the 18O treatment

        dev.off()

      #Now, subset the ncopies.tube.melted dataframe FOR REAL so that it only contains taxon-tubes for taxa occurring in all 3 of the remaining 18O reps:
        ncopies.tube.melted <- ncopies.tube.melted[ncopies.tube.melted$taxon %in% names(tubes.per.taxon.18O[as.numeric(tubes.per.taxon.18O) >= 3]),]
        row.names(ncopies.tube.melted) <- 1:dim(ncopies.tube.melted)[1]
          #Re-convert the factor columns to factor:
            ncopies.tube.melted$taxon <- factor(ncopies.tube.melted$taxon)
            ncopies.tube.melted$kingdom <- factor(ncopies.tube.melted$kingdom)
            ncopies.tube.melted$phylum <- factor(ncopies.tube.melted$phylum)
            ncopies.tube.melted$class <- factor(ncopies.tube.melted$class)
            ncopies.tube.melted$order <- factor(ncopies.tube.melted$order)
            ncopies.tube.melted$family <- factor(ncopies.tube.melted$family)
            ncopies.tube.melted$genus <- factor(ncopies.tube.melted$genus)
        dim(ncopies.tube.melted.all)
        dim(ncopies.tube.melted)


    #STEP 2:
      #Next, subset the step-1 filtered data above into only those taxon-tubes with copies present in the 16O treatment:
      ncopies.tube.melted.occurrences.16O <- ncopies.tube.melted[ncopies.tube.melted$copies > 0 & ncopies.tube.melted$tmt == "16O",]
          #Convert the factor columns to factor:
            ncopies.tube.melted.occurrences.16O$taxon <- factor(ncopies.tube.melted.occurrences.16O$taxon)
            ncopies.tube.melted.occurrences.16O$kingdom <- factor(ncopies.tube.melted.occurrences.16O$kingdom)
            ncopies.tube.melted.occurrences.16O$phylum <- factor(ncopies.tube.melted.occurrences.16O$phylum)
            ncopies.tube.melted.occurrences.16O$class <- factor(ncopies.tube.melted.occurrences.16O$class)
            ncopies.tube.melted.occurrences.16O$order <- factor(ncopies.tube.melted.occurrences.16O$order)
            ncopies.tube.melted.occurrences.16O$family <- factor(ncopies.tube.melted.occurrences.16O$family)
            ncopies.tube.melted.occurrences.16O$genus <- factor(ncopies.tube.melted.occurrences.16O$genus)
        dim(ncopies.tube.melted.all)
        dim(ncopies.tube.melted)
        dim(ncopies.tube.melted.occurrences.16O)

      #Calculate the number of unique tubes in the 16O treatment with copies present for each taxon:
        tubes.per.taxon.16O <- tapply(ncopies.tube.melted.occurrences.16O$tube, ncopies.tube.melted.occurrences.16O$taxon, length.unique)
        tubes.per.taxon.16O <- sort(tubes.per.taxon.16O)

      #Keep only those taxa that occur in at least 3 of the 4 remaining 16O tubes of the experiment:
        dev.off()
        dev.new(width=3.5, height=7)
        par(mfrow=c(2,1))
        hist(tubes.per.taxon.16O, xlab="# of tubes per taxon", ylab="# of taxa", main="")
        abline(v=2.75, col="red")
        plot(x=1:length(tubes.per.taxon.16O), y=as.numeric(tubes.per.taxon.16O), xlab="taxon (sorted by # tubes per taxon)", ylab="# of tubes per taxon", type="l")
        abline(h=2.75, col="red")
        par(mfrow=c(1,1))
        length(tubes.per.taxon.16O)   #total number of taxa occuring in the 16O treatment
        # tubes.per.taxon.16O[as.numeric(tubes.per.taxon.16O) >= 3]
        length(tubes.per.taxon.16O[as.numeric(tubes.per.taxon.16O) >= 3])   #number of taxa that occurred in at least 3 of the 4 remaining tubes from the 16O treatment

        dev.off()

      #Now, subset the ncopies.tube.melted dataframe FOR REAL AGAIN so that it only contains taxon-tubes for taxa occurring in at least 3 of the 4 remaining 16O tubes of the experiment:
      #(This is in addition to the STEP 1 filtering step above, which kept only taxa that occurred in all 3 of the remaining 18O reps)
        ncopies.tube.melted <- ncopies.tube.melted[ncopies.tube.melted$taxon %in% names(tubes.per.taxon.16O[as.numeric(tubes.per.taxon.16O) >= 3]),]
        row.names(ncopies.tube.melted) <- 1:dim(ncopies.tube.melted)[1]
          #Re-convert the factor columns to factor:
            ncopies.tube.melted$taxon <- factor(ncopies.tube.melted$taxon)
            ncopies.tube.melted$kingdom <- factor(ncopies.tube.melted$kingdom)
            ncopies.tube.melted$phylum <- factor(ncopies.tube.melted$phylum)
            ncopies.tube.melted$class <- factor(ncopies.tube.melted$class)
            ncopies.tube.melted$order <- factor(ncopies.tube.melted$order)
            ncopies.tube.melted$family <- factor(ncopies.tube.melted$family)
            ncopies.tube.melted$genus <- factor(ncopies.tube.melted$genus)
        dim(ncopies.tube.melted.all)
        dim(ncopies.tube.melted)
        length(levels(ncopies.tube.melted.all$taxon))
        length(levels(ncopies.tube.melted$taxon))
        #Taxa excluded according to the filtering steps above:
        length(levels(ncopies.tube.melted.all$taxon)[!is.element(levels(ncopies.tube.melted.all$taxon), levels(ncopies.tube.melted$taxon))])        
        levels(ncopies.tube.melted.all$taxon)[!is.element(levels(ncopies.tube.melted.all$taxon), levels(ncopies.tube.melted$taxon))]


    #STEP 3:
      #Next, subset the step-1 & step-2 filtered data above into only those taxon-tubes with copies present in the Time0 treatment:
      ncopies.tube.melted.occurrences.Time0 <- ncopies.tube.melted[ncopies.tube.melted$copies > 0 & ncopies.tube.melted$tmt == "Time0",]
          #Convert the factor columns to factor:
            ncopies.tube.melted.occurrences.Time0$taxon <- factor(ncopies.tube.melted.occurrences.Time0$taxon)
            ncopies.tube.melted.occurrences.Time0$kingdom <- factor(ncopies.tube.melted.occurrences.Time0$kingdom)
            ncopies.tube.melted.occurrences.Time0$phylum <- factor(ncopies.tube.melted.occurrences.Time0$phylum)
            ncopies.tube.melted.occurrences.Time0$class <- factor(ncopies.tube.melted.occurrences.Time0$class)
            ncopies.tube.melted.occurrences.Time0$order <- factor(ncopies.tube.melted.occurrences.Time0$order)
            ncopies.tube.melted.occurrences.Time0$family <- factor(ncopies.tube.melted.occurrences.Time0$family)
            ncopies.tube.melted.occurrences.Time0$genus <- factor(ncopies.tube.melted.occurrences.Time0$genus)
        dim(ncopies.tube.melted.all)
        dim(ncopies.tube.melted)
        dim(ncopies.tube.melted.occurrences.Time0)

      #Calculate the number of unique tubes in the Time0 treatment with copies present for each taxon:
        tubes.per.taxon.Time0 <- tapply(ncopies.tube.melted.occurrences.Time0$tube, ncopies.tube.melted.occurrences.Time0$taxon, length.unique)
        tubes.per.taxon.Time0 <- sort(tubes.per.taxon.Time0)

      #Keep only those taxa that occur in at least 3 of the Time0 tubes of the experiment:
        dev.off()
        dev.new(width=3.5, height=7)
        par(mfrow=c(2,1))
        hist(tubes.per.taxon.Time0, xlab="# of tubes per taxon", ylab="# of taxa", main="")
        abline(v=2.75, col="red")
        plot(x=1:length(tubes.per.taxon.Time0), y=as.numeric(tubes.per.taxon.Time0), xlab="taxon (sorted by # tubes per taxon)", ylab="# of tubes per taxon", type="l")
        abline(h=2.75, col="red")
        par(mfrow=c(1,1))
        length(tubes.per.taxon.Time0)   #total number of taxa occuring in the Time0 treatment
        # tubes.per.taxon.Time0[as.numeric(tubes.per.taxon.Time0) >= 3]
        length(tubes.per.taxon.Time0[as.numeric(tubes.per.taxon.Time0) >= 3])   #number of taxa that occurred in at least 3 of the 5 tubes from the Time0 treatment

        dev.off()

      #Now, subset the ncopies.tube.melted dataframe FOR REAL AGAIN so that it only contains taxon-tubes for taxa occurring in at least 3 of the 5 Time0 tubes of the experiment:
      #(This is in addition to the STEP 1 filtering step above, which kept only taxa that occurred in all 3 of the remaining 18O reps, and the STEP 2 filtering step above, which kept only taxa that occurred in 3 of the remaining 4 16O reps)
        ncopies.tube.melted <- ncopies.tube.melted[ncopies.tube.melted$taxon %in% names(tubes.per.taxon.Time0[as.numeric(tubes.per.taxon.Time0) >= 3]),]
        row.names(ncopies.tube.melted) <- 1:dim(ncopies.tube.melted)[1]
          #Re-convert the factor columns to factor:
            ncopies.tube.melted$taxon <- factor(ncopies.tube.melted$taxon)
            ncopies.tube.melted$kingdom <- factor(ncopies.tube.melted$kingdom)
            ncopies.tube.melted$phylum <- factor(ncopies.tube.melted$phylum)
            ncopies.tube.melted$class <- factor(ncopies.tube.melted$class)
            ncopies.tube.melted$order <- factor(ncopies.tube.melted$order)
            ncopies.tube.melted$family <- factor(ncopies.tube.melted$family)
            ncopies.tube.melted$genus <- factor(ncopies.tube.melted$genus)
        dim(ncopies.tube.melted.all)
        dim(ncopies.tube.melted)
        length(levels(ncopies.tube.melted.all$taxon))
        length(levels(ncopies.tube.melted$taxon))
        #Taxa excluded according to the filtering steps above:
        length(levels(ncopies.tube.melted.all$taxon)[!is.element(levels(ncopies.tube.melted.all$taxon), levels(ncopies.tube.melted$taxon))])        
        levels(ncopies.tube.melted.all$taxon)[!is.element(levels(ncopies.tube.melted.all$taxon), levels(ncopies.tube.melted$taxon))]


    #FINAL STEP:
      #Filter the 'data.melted' qSIP data frame in the same way as was done above for the 'ncopies.tube.melted' dataframe so that they both contain the same set of taxa:
        data.melted  <- data.melted[data.melted$taxon %in% levels(ncopies.tube.melted$taxon),]
        row.names(data.melted) <- 1:dim(data.melted)[1]
          #Re-convert the factor columns to factor:
            data.melted$taxon <- factor(data.melted$taxon)
            data.melted$kingdom <- factor(data.melted$kingdom)
            data.melted$phylum <- factor(data.melted$phylum)
            data.melted$class <- factor(data.melted$class)
            data.melted$order <- factor(data.melted$order)
            data.melted$family <- factor(data.melted$family)
            data.melted$genus <- factor(data.melted$genus)
        dim(data.melted.all)
        dim(data.melted)
        length(levels(data.melted.all$taxon))
        length(levels(data.melted$taxon))
        #Taxa excluded according to the filtering steps above:
        length(levels(data.melted.all$taxon)[!is.element(levels(data.melted.all$taxon), levels(data.melted$taxon))])        
        levels(data.melted.all$taxon)[!is.element(levels(data.melted.all$taxon), levels(data.melted$taxon))]
    

    #Calculate the proportion of sequence reads excluded (after QIIME filtering & after exclusion of outlier density fractions) when the threshold for inclusion of a taxon is set such that it must occur according to the criteria above (in at least 3 reps of each treatment):
      #First, calculate the sum of all relative abundances of all taxa prior to filtering:
        #Do this for the Time0 and 16O and 18O data (but remember to exclude tubes 5 & 9 & 15 (16O-rep2 & 18O-rep3 & 18O-rep5)):
          taxa.rel.abundance <- apply(data[!is.element(data$tube, c(5,9,15)), 17:(dim(data)[2]-1)], 2, sum, na.rm=TRUE)
          sum(as.numeric(taxa.rel.abundance)) == sum(data[!is.element(data$tube, c(5,9,15)), 17:(dim(data)[2]-1)], na.rm=TRUE)
          tot.rel.abundance <- sum(as.numeric(taxa.rel.abundance))
          tot.rel.abundance
      #Next, calculate the sum of all relative abundances of all taxa that remain after filtering:
        #Do this for the Time0 and 16O and 18O data (but remember to exclude tubes 5 & 9 & 15 (16O-rep2 & 18O-rep3 & 18O-rep5)):
          data.rel.abundance.filtered <- data[!is.element(data$tube, c(5,9,15)), names(data)[names(data) %in% levels(ncopies.tube.melted$taxon)]]
          sum(data.rel.abundance.filtered, na.rm=TRUE)
      #Excluded reads are Y proportion of total reads:
        (tot.rel.abundance - sum(data.rel.abundance.filtered, na.rm=TRUE)) / tot.rel.abundance
      #Included reads are Y proportion of total reads:
        sum(data.rel.abundance.filtered, na.rm=TRUE) / tot.rel.abundance
        #~99% of reads are kept using this filtering approach (after the filtering that was done in QIIME & after exclusion of outlier density fractions)
  #}

##########}__________________________________________________________________________________________________________


##########{___________Some examples of how to use the qSIP functions_________________________________________________

  #Create a subset of the data containing only data for taxon 104 in tube 2 of treatment '16O'
    T104T1T2 <- data.melted[data.melted$taxon==104 & data.melted$tmt=="16O" & data.melted$tube==2,]
  #Calculate the weighted average density for taxon 104 in tube 2:
    T104T1T2.wad.out <- WAD.func(y=T104T1T2$copies, x=T104T1T2$density.g.ml)
    T104T1T2.wad.out
  #Create a subset of the data containing only data for taxon 104 in all tubes of treatment '16O'
    T104T1 <- data.melted[data.melted$taxon==104 & data.melted$tmt=="16O",]
  #Calculate the bootstrapped (and observed) weighted average density estimates for taxon 104 in treatment '16O'
    T104T1.boot.out <- boot.WAD.func(X=T104T1, vars=c("density.g.ml", "copies", "tube"), CI=0.90, draws=1000)
    T104T1.boot.out
  #Create a subset of the data containing only data for taxon 104 in all tubes of treatment '18O'
    T104T2 <- data.melted[data.melted$taxon==104 & data.melted$tmt=="18O",]
  #Calculate the difference in weighted average density (with bootstrapped estimates of uncertainty) between the '18O' and '16O' treatments for taxon 104
    T104T1v2.out <- boot.diff.wad(T104T1, T104T2, vars=c("density.g.ml", "copies", "tube", "tmt"), CI=0.90, draws=1000, tailed.test=1)
    T104T1v2.out
  #Create a subset of the data containing only data for taxon 104 in all tubes of the 'reference' treatments (i.e., '16O')
    T104ref <- data.melted[data.melted$taxon==104 & (data.melted$tmt=="16O"),]
  #Calculate the molecular weight, GC content, and average number of carbon and oxygen atoms per DNA nucleotide for taxon 104:
    T104ref.MW.out <- MW.calc(X=T104ref, vars=c("density.g.ml", "copies", "tube"))
    T104ref.MW.out
  #Calculate the excess atom fraction of 18O (with bootstrapped estimates of uncertainty) in the '18O' vs. the '16O' treatment for taxon 104
    T104T1v2.ape.out <- boot.diff.ape(X.light=T104T1, X.heavy=T104T2, X.reference=T104ref, iso.compare="18O", vars=c("density.g.ml", "copies", "tube", "tmt"), CI=0.90, draws=1000)
    T104T1v2.ape.out
  #Calculate the intrinsic rate of increase (r, with bootstrapped estimates of uncertainty) for taxon 104 as determined by 18O incorporation; set prop.O.from.water=0.33 according to McHugh et al. unpublished data from E. coli pure cultures incubated with 18O-H2O
    T104T1v2.r.out <- boot.diff.r(X.light=T104T1, X.heavy=T104T2, X.reference=T104ref, M.soil=Sdat, iso.compare="18O", days=10, vars=c("density.g.ml", "copies", "tube", "tmt", "g.soil"), growth.model="exponential", prop.O.from.water=0.33, v.frac=50, CI=0.90, draws=1000)
    T104T1v2.r.out
  #Calculate the bootstrapped (and observed) weighted average density estimates for taxon 104 in treatment '16O' along with the corresponding values for soil mass and total 16S copies from the tubes chosen for each bootstrap sample
    T104T1.tube.out <- boot.TUBE.func(X=T104T1, M.soil=Sdat, vars=c("density.g.ml", "copies", "tube", "g.soil"), v.frac=50, CI=0.90, draws=1000)
    T104T1.tube.out
  #Calculate the flux of carbon into biomass (with bootstrapped estimates of uncertainty) for taxon 104 as determined by 18O incorporation; set prop.O.from.water=0.33 according to McHugh et al. unpublished data from E. coli pure cultures incubated with 18O-H2O
    T104T1v2.f.out <- boot.diff.f(X.light=T104T1, X.heavy=T104T2, X.reference=T104ref, M.soil=Sdat, iso.compare="18O", days=10, vars=c("density.g.ml", "copies", "tube", "tmt", "g.soil"), growth.model="exponential", prop.O.from.water=0.33, v.frac=50, copies.cell=6, pgC.cell=0.1, CI=0.90, draws=1000)
    T104T1v2.f.out
  #Note that the 'message' in the output of most of the functions warns if a taxon was not detected in a tube(s) (as was the case for taxon 341):
    T341T1 <- data.melted[data.melted$taxon==341 & data.melted$tmt=="16O",]
    T341T1.boot.out <- boot.WAD.func(X=T341T1, vars=c("density.g.ml", "copies", "tube"), CI=0.90, draws=1000)
    T341T1.boot.out
    T341T2 <- data.melted[data.melted$taxon==341 & data.melted$tmt=="18O",]
    T341T2.boot.out <- boot.WAD.func(X=T341T2, vars=c("density.g.ml", "copies", "tube"), CI=0.90, draws=1000)
    T341T2.boot.out
    T341T1v2.out <- boot.diff.wad(T341T1, T341T2, vars=c("density.g.ml", "copies", "tube", "tmt"), CI=0.90, draws=1000, tailed.test=1)
    T341T1v2.out

##########}__________________________________________________________________________________________________________


##########{___________Some examples of how to use the population-level functions_____________________________________

  #Create a subset of the data containing only data for taxon 74 in all tubes of treatment 'Time0'
    T74T0.tube <- ncopies.tube.melted[ncopies.tube.melted$taxon==74 & ncopies.tube.melted$tmt=="Time0",]
  #Calculate the bootstrapped (and observed) copies/uL for taxon 74 in treatment 'Time0'
    T74T0.tube.copies.out <- boot.pop(X=T74T0.tube, vars=c("copies", "tube"), CI=0.90, draws=1000)
    T74T0.tube.copies.out
  #Create a subset of the data containing only data for taxon 74 in all tubes of treatment '16O'
    T74T1.tube <- ncopies.tube.melted[ncopies.tube.melted$taxon==74 & ncopies.tube.melted$tmt=="16O",]
  #Calculate the bootstrapped (and observed) copies/uL for taxon 74 in treatment '16O'
    T74T1.tube.copies.out <- boot.pop(X=T74T1.tube, vars=c("copies", "tube"), CI=0.90, draws=1000)
    T74T1.tube.copies.out
  #Create a subset of the data containing only data for taxon 74 in all tubes of treatment '18O'
    T74T2.tube <- ncopies.tube.melted[ncopies.tube.melted$taxon==74 & ncopies.tube.melted$tmt=="18O",]
  #Calculate the bootstrapped (and observed) copies/uL for taxon 74 in treatment '18O'
    T74T2.tube.copies.out <- boot.pop(X=T74T2.tube, vars=c("copies", "tube"), CI=0.90, draws=1000)
    T74T2.tube.copies.out
  #Calculate the bootstrapped (and observed) copies/uL for taxon 74 in treatments '16O' and '18O' (both treatments at day 10 combined)
    T74T12.tube.copies.out <- boot.pop(X=rbind(T74T1.tube, T74T2.tube), vars=c("copies", "tube"), CI=0.90, draws=1000)
    T74T12.tube.copies.out
  #Calculate the bootstrapped (and observed) copies/uL for taxon 74 in treatment 'Time0' along with the corresponding values for soil mass and total 16S copies from the tubes chosen for each bootstrap sample
    T74T0.tube.copies.out <- boot.TUBE.pop(X=T74T0.tube, M.soil=Sdat, vars=c("copies", "tube", "g.soil"), vol=100, CI=0.90, draws=1000)
    T74T0.tube.copies.out
  #Calculate the bootstrapped (and observed) copies/uL for taxon 74 in treatment '16O' along with the corresponding values for soil mass and total 16S copies from the tubes chosen for each bootstrap sample
    T74T1.tube.copies.out <- boot.TUBE.pop(X=T74T1.tube, M.soil=Sdat, vars=c("copies", "tube", "g.soil"), vol=50, CI=0.90, draws=1000)
    T74T1.tube.copies.out
  #Calculate the bootstrapped (and observed) copies/uL for taxon 74 in treatment '18O' along with the corresponding values for soil mass and total 16S copies from the tubes chosen for each bootstrap sample
    T74T2.tube.copies.out <- boot.TUBE.pop(X=T74T2.tube, M.soil=Sdat, vars=c("copies", "tube", "g.soil"), vol=50, CI=0.90, draws=1000)
    T74T2.tube.copies.out
  #Calculate the bootstrapped (and observed) copies/uL for taxon 74 in treatments '16O' & '18O' (both treatments at day 10 combined) along with the corresponding values for soil mass and total 16S copies from the tubes chosen for each bootstrap sample
    T74T12.tube.copies.out <- boot.TUBE.pop(X=rbind(T74T1.tube, T74T2.tube), M.soil=Sdat, vars=c("copies", "tube", "g.soil"), vol=50, CI=0.90, draws=1000)
    T74T12.tube.copies.out
  #Calculate the rate of population increase (r, with bootstrapped estimates of uncertainty) for taxon 74 as determined by the change in the number of 16S copies over time (using a linear growth model for 'Time0' and '16O' treatments only)
    T74T0v1.r.out.lin <- boot.r.pop(T0=T74T0.tube, Tt=T74T1.tube, M.soil=Sdat, days=10, vars=c("copies", "tube", "tmt", "g.soil"), growth.model="linear", vol=c(100, 50), CI=0.90, draws=1000)
    T74T0v1.r.out.lin
  #Calculate the rate of population increase (r, with bootstrapped estimates of uncertainty) for taxon 74 as determined by the change in the number of 16S copies over time (using an exponential growth model for 'Time0' and '16O' treatments only)
    T74T0v1.r.out <- boot.r.pop(T0=T74T0.tube, Tt=T74T1.tube, M.soil=Sdat, days=10, vars=c("copies", "tube", "tmt", "g.soil"), growth.model="exponential", vol=c(100, 50), CI=0.90, draws=1000)
    T74T0v1.r.out
  #Calculate the rate of population increase (r, with bootstrapped estimates of uncertainty) for taxon 74 as determined by the change in the number of 16S copies over time (using an exponential growth model for 'Time0' and '18O' treatments only)
    T74T0v2.r.out <- boot.r.pop(T0=T74T0.tube, Tt=T74T2.tube, M.soil=Sdat, days=10, vars=c("copies", "tube", "tmt", "g.soil"), growth.model="exponential", vol=c(100, 50), CI=0.90, draws=1000)
    T74T0v2.r.out
  #Calculate the rate of population increase (r, with bootstrapped estimates of uncertainty) for taxon 74 as determined by the change in the number of 16S copies over time (using an exponential growth model for 'Time0' to '16O' & '18O' treatments combined at day 10)
    T74T0v12.r.out <- boot.r.pop(T0=T74T0.tube, Tt=rbind(T74T1.tube, T74T2.tube), M.soil=Sdat, days=10, vars=c("copies", "tube", "tmt", "g.soil"), growth.model="exponential", vol=c(100, 50), CI=0.90, draws=1000)
    T74T0v12.r.out
  #Calculate the flux of carbon into biomass (f, with bootstrapped estimates of uncertainty) for taxon 74 as determined by the change in the number of 16S copies over time (using an exponential growth model for 'Time0' and '16O' treatments only)
    T74T0v1.f.out <- boot.f.pop(T0=T74T0.tube, Tt=T74T1.tube, M.soil=Sdat, days=10, vars=c("copies", "tube", "tmt", "g.soil"), growth.model="exponential", vol=c(100, 50), copies.cell=6, pgC.cell=0.1, CI=0.90, draws=1000)
    T74T0v1.f.out
  #Calculate the flux of carbon into biomass (f, with bootstrapped estimates of uncertainty) for taxon 74 as determined by the change in the number of 16S copies over time (using an exponential growth model for 'Time0' and '18O' treatments only)
    T74T0v2.f.out <- boot.f.pop(T0=T74T0.tube, Tt=T74T2.tube, M.soil=Sdat, days=10, vars=c("copies", "tube", "tmt", "g.soil"), growth.model="exponential", vol=c(100, 50), copies.cell=6, pgC.cell=0.1, CI=0.90, draws=1000)
    T74T0v2.f.out
  #Calculate the flux of carbon into biomass (f, with bootstrapped estimates of uncertainty) for taxon 74 as determined by the change in the number of 16S copies over time (using an exponential growth model for 'Time0' to '16O' & '18O' treatments combined at day 10)
    T74T0v12.f.out <- boot.f.pop(T0=T74T0.tube, Tt=rbind(T74T1.tube, T74T2.tube), M.soil=Sdat, days=10, vars=c("copies", "tube", "tmt", "g.soil"), growth.model="exponential", vol=c(100, 50), copies.cell=6, pgC.cell=0.1, CI=0.90, draws=1000)
    T74T0v12.f.out

##########}__________________________________________________________________________________________________________




#Save workspace:
  save(list=ls(all.names=TRUE), file="qSIP_workspaces/TM_01/.RData", envir=.GlobalEnv)



