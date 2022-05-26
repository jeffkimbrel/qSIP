# This code imports the BH OTU qSIP data sets and performs the basic calculations on the data


#Set working directory and load libraries & scripts:
  #Set working directory:
    setwd("/Users/bk/Research/Projects/SIP_Modeling/qSIP")
  #Set path for the qSIP functions:
    func.path <- "/Users/bk/Research/Projects/SIP_Modeling/qSIP/qSIP_repo"
  #Set path for the qSIP input data:
    data.path <- "/Users/bk/Research/Projects/SIP_Modeling/qSIP/qSIP_data"
  library(reshape2)
  library(VennDiagram)
  source(paste(func.path, "/sample.vec.R", sep=""))                     #sample.vec
  source(paste(func.path, "/WAD.func.R", sep=""))                       #WAD.func
  source(paste(func.path, "/fit.norm.func.R", sep=""))                  #fit.norm.func
  source(paste(func.path, "/boot.WAD.func.R", sep=""))                  #boot.WAD.func
  source(paste(func.path, "/diff.wad.calc.R", sep=""))                  #diff.wad.calc
  source(paste(func.path, "/boot.diff.wad.R", sep=""))                  #boot.diff.wad
  source(paste(func.path, "/MW.calc.R", sep=""))                        #MW.calc
  source(paste(func.path, "/MW.calc.Schildkraut.R", sep=""))            #MW.calc.Schildkraut
  source(paste(func.path, "/comparison.message.R", sep=""))             #comparison.message
  source(paste(func.path, "/ape.calc.R", sep=""))                       #ape.calc
  source(paste(func.path, "/boot.diff.ape.R", sep=""))                  #boot.diff.ape
  source(paste(func.path, "/r.calc.R", sep=""))                         #r.calc
  source(paste(func.path, "/boot.diff.r.R", sep=""))                    #boot.diff.r
  source(paste(func.path, "/boot.TUBE.func.R", sep=""))                 #boot.TUBE.func
  source(paste(func.path, "/f.calc.R", sep=""))                         #f.calc
  source(paste(func.path, "/boot.diff.f.R", sep=""))                    #boot.diff.f
  source(paste(func.path, "/all.taxa.calcs.R", sep=""))                 #all.taxa.calcs
  source(paste(func.path, "/id.reps.R", sep=""))                        #id.reps
  source(paste(func.path, "/select.rep.R", sep=""))                     #select.rep
  source(paste(func.path, "/explore.filter.taxa.R", sep=""))            #explore.filter.taxa
  source(paste(func.path, "/filter.taxa.R", sep=""))                    #filter.taxa
  source(paste(func.path, "/explore.filter.fractions.taxa.R", sep=""))  #explore.filter.fractions.taxa
  source(paste(func.path, "/filter.fractions.taxa.R", sep=""))          #filter.fractions.taxa
  source(paste(func.path, "/WAD.by.taxon.func.R", sep=""))              #WAD.by.taxon.func
  source(paste(func.path, "/SE.WAD.by.taxon.plot.R", sep=""))           #SE.WAD.by.taxon.plot
  source(paste(func.path, "/find.unlabeled.correction.R", sep=""))      #find.unlabeled.correction
  source(paste(func.path, "/find.labeled.correction.R", sep=""))        #find.labeled.correction
  source(paste(func.path, "/td.pos.resid.R", sep=""))                   #td.pos.resid
  source(paste(func.path, "/td.abs.resid.R", sep=""))                   #td.abs.resid
  source(paste(func.path, "/bu.abs.resid.R", sep=""))                   #bu.abs.resid
  source(paste(func.path, "/select.best.iteration.R", sep=""))          #select.best.iteration
  source(paste(func.path, "/find.labeled.correction.plot.R", sep=""))   #find.labeled.correction.plot
  source(paste(func.path, "/get.seq.taxa.nums.R", sep=""))              #get.seq.taxa.nums
  source(paste(func.path, "/add.lab.WAD.corr.summary.R", sep=""))       #add.lab.WAD.corr.summary
  source(paste(func.path, "/apply.unlabeled.correction.R", sep=""))     #apply.unlabeled.correction
  source(paste(func.path, "/apply.labeled.correction.R", sep=""))       #apply.labeled.correction


#Import raw data & taxomic information and format raw data for analysis:
  #Read in taxonomic information (currently represented as long column names in the raw data file)
    taxa.id <- read.table(paste(data.path, "/BH_OTU_taxa_id.txt", sep=""), header=T, sep="\t", stringsAsFactors=T, na.strings="")
  #The OTU number serves as a unique identifier code for each taxon to the "taxa.id" data frame:
    
  #Read in raw data (skip header)
  #NOTE: this data is still in units of copies, not relative abundance
  #NOTE: treatment code names cannot contain spaces
    data <- read.table(paste(data.path, "/BH_OTU_qSIP_data.txt", sep=""), header=F, sep="\t", stringsAsFactors=T)
  #Name columns; columns 7 and up are given the unique OTU identifier code for taxon to match that in the "taxa.id" data frame
  # (OTUs in the original data frame are in the same order as those in the imported taxa.id data frame)
    names(data)[1:6] <- c("SampleID", "trt.code", "DNA.ng.ul", "density.g.ml", "copy.number", "g.soil.extracted")
    names(data)[7:ncol(data)] <- taxa.id$taxon
    
  #Add a column for the sum of sequence abundance of all taxa by fraction
    data$sum.abundance <- rowSums(data[,7:ncol(data)])

  #Calculate relative abundance of each OTU in each fraction:
    data.rel <- (1/data$sum.abundance)*data[,7:(ncol(data)-1)]
    data.rel <- cbind(data[,1:6], data.rel)  # add first 6 columns of data to data.rel

  #Add a column for the sum of proportional abundance of all taxa by fraction
    data.rel$sum.abundance <- rowSums(data.rel[,7:ncol(data.rel)])

  #Do not correct the proportional abundances of all taxa by scaling by their sum because all OTUs are included and they already sum to 1
  #NOTE: becasue singletons & doubletons were filtered out prior to calculating relative abundance, the relative abundances for the remaining OTUs may be overestimated
    data.rel$sum.abundance[data.rel$sum.abundance != 1]

  #Calculate number of copies per uL, based on relative abundance and total number of copies per uL
    ncopies <- data.rel$copy.number*data.rel[,7:(ncol(data.rel)-1)]
    ncopies <- cbind(data.rel[,1:6], ncopies)  # add first 6 columns of data to ncopies

  #Add tube ID to ncopies:
    #First, pull out the correct SampleID numbers from the obscure codings:
      SampleID1 <- ncopies$SampleID
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
    #Add the tube ID to ncopies:
      dotpos <- regexpr("\\.[^\\.]*$", SampleID.new$ID3)	# find ".", (number prior to "." is tube number)
      tube <- substr(SampleID.new$ID3, 1, dotpos-1)
      ncopies <- cbind(tube, ncopies)
      ncopies$SampleID <- factor(ncopies$SampleID)

  #Remove a duplicate sample to avoid introducing bias in calculating WAD:
    #according to Becky on Jan 29, 2015: "Looks like 5.20 (119.5.20) was just sequenced twice, but is the same sample, so either can be used"
    frac.5.20 <- t(ncopies[ncopies$SampleID %in% c("5.2", "119.5.20"),8:dim(ncopies)[2]])
    plot(log10(frac.5.20[,1]) ~ log10(frac.5.20[,2]))
    abline(a=0, b=1)
    #Arbitrarily decided to keep the second case of sample 5.20 (119.5.20):
    dev.off()
    dim(ncopies)
    ncopies <- ncopies[ncopies$SampleID != "5.2",]
    ncopies$SampleID <- factor(ncopies$SampleID)
    row.names(ncopies) <- 1:dim(ncopies)[1]
    dim(ncopies)
    dim(data.rel)
    data.rel <- data.rel[data.rel$SampleID != "5.2",]
    data.rel$SampleID <- factor(data.rel$SampleID)
    row.names(data.rel) <- 1:dim(data.rel)[1]
    dim(data.rel)

  #Melt data into long format by tube, trt.code, DNA conc, density and copy number;
  #Do this for ncopies and for relative abundance. Merge these to into 1 masterfile: data.melted
    ncopies.melted <- melt(ncopies, id=c("tube","SampleID", "trt.code", "DNA.ng.ul", "density.g.ml", "copy.number"), measure.vars=names(ncopies)[8:dim(ncopies)[2]], variable.name="taxon", value.name="copies")
    rel.abundance.melted <- melt(data.rel, id=c("SampleID", "trt.code", "DNA.ng.ul", "density.g.ml", "copy.number"), measure.vars=names(data.rel)[7:(dim(data.rel)[2]-1)], variable.name="taxon", value.name="rel.abundance")
    data.melted <- merge(ncopies.melted, rel.abundance.melted)

  #Merge taxa data and reorder data frame by taxon and SampleID
    data.melted <- merge(data.melted, taxa.id)
    data.melted <- data.melted[order(data.melted$taxon, data.melted$SampleID),]


#Import soil extraction data & calculate mass of soil representative of each tube:
  #(Use this calculation instead of using the mass of soil extracted ('g.soil.extracted') in the qSIP data)
  Sdat <- read.table(paste(data.path, "/BH_SoilExtractionData.txt", sep=""), header=TRUE, sep="\t")
  Sdat$tube <- factor(Sdat$tube)
  Sdat$g.soil <- Sdat$ug.DNA.added.to.tube / Sdat$ug.DNA.g.soil
  summary(Sdat)


#Import data frame containing the treatment comparisons to perform:
  Tcompare <- read.table(paste(data.path, "/BH_TreatmentComparisons.txt", sep=""), header=TRUE, sep="\t", colClasses=c("numeric","factor","character","character","character","numeric","character"))
  summary(Tcompare)
  Tcompare


#Filter the data so that only OTUs with an appropriate level of occurrence and replication among the tubes and treatments of the experiment are retained for further analysis:
  #First back up the complete data frame:
    data.melted.all <- data.melted

  #Examining the 0.005% filtering criterion:
    #Calculate the number of OTUs that meet the threshold of making up at least 0.005% of total sequence reads (after the singletons and doubletons were removed):
      OTU.tot.reads <- apply(data[,7:(dim(data)[2]-1)], 2, sum)
      sum(as.numeric(OTU.tot.reads)) == sum(data[,7:(dim(data)[2]-1)])
      tot.reads <- sum(as.numeric(OTU.tot.reads))
      OTU.prop.tot.reads <- OTU.tot.reads / tot.reads
      # OTU.tot.reads[OTU.prop.tot.reads < 0.00005]
      length(OTU.tot.reads[OTU.prop.tot.reads < 0.00005])
      length(OTU.tot.reads[OTU.prop.tot.reads >= 0.00005])    
    
    #For comparison, the total number of OTUs (after the singletons and doubletons were removed) was:
      length(levels(data.melted.all$taxon))
    
    #Calculate the number (and proportion) of sequence reads excluded (after the exclusion of singletons & doubletons) when 
      #the threshold for inclusion of an OTU is set at at least 0.005% of total reads across all reps and treatments: 
      sum(OTU.tot.reads[OTU.prop.tot.reads < 0.00005])
      sum(OTU.tot.reads[OTU.prop.tot.reads < 0.00005]) / tot.reads
      #it is only ~6% of reads excluded using the 0.00005 threshold level
  
  
  #An alternative, and perhaps preferred filtering method is to keep only those OTUs that occurred in all (or nearly all) tubes of the experiment:
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
          data.melted.occurrences$species <- factor(data.melted.occurrences$species)
      dim(data.melted)
      dim(data.melted.occurrences)
  
    #Create a function to pass to tapply for calculating the number of unique tubes having copies for each taxon (OTU):
      length.unique <- function(X){
        length(unique(X))
      }
  
    #Calculate the number of unique tubes with copies present for each taxon (OTU):
      tubes.per.OTU <- tapply(data.melted.occurrences$tube, data.melted.occurrences$taxon, length.unique)
      tubes.per.OTU <- sort(tubes.per.OTU)
  
    #Keep only those OTUs that occur in at least 34 of the 36 tubes of the experiment:
      dev.off()
      dev.new(width=3.5, height=7)
      par(mfrow=c(2,1))
      hist(tubes.per.OTU, main="")
      abline(v=33.5, col="red")
      plot(x=1:length(tubes.per.OTU), y=as.numeric(tubes.per.OTU), type="l")
      abline(h=33.5, col="red")
      par(mfrow=c(1,1))
      # tubes.per.OTU[as.numeric(tubes.per.OTU) >= 34]
      length(tubes.per.OTU[as.numeric(tubes.per.OTU) >= 34])
  
      dev.off()
  
    #Subset the data.melted dataframe so that it only contains taxon-fractions for OTUs occurring in at least 34 of the 36 tubes:
      data.melted <- data.melted[data.melted$taxon %in% names(tubes.per.OTU[as.numeric(tubes.per.OTU) >= 34]),]
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
          data.melted$species <- factor(data.melted$species)
      dim(data.melted.all)
      dim(data.melted)
      
    #Calculate the number (and proportion) of sequence reads excluded (after the exclusion of singletons & doubletons) when 
      #the threshold for inclusion of an OTU is set such that it must occur in at least 34 of the 36 tubes of the experiment: 
      data.reads.filtered <- data[, names(data)[names(data) %in% levels(data.melted$taxon)]]
      sum(data.reads.filtered)
      tot.reads
      #Excluded reads:
        tot.reads - sum(data.reads.filtered)
        (tot.reads - sum(data.reads.filtered)) / tot.reads
        #it is only ~7% of reads excluded using this filtering approach
      #Included reads:
        sum(data.reads.filtered)
        sum(data.reads.filtered) / tot.reads
        #~93% of reads are kept using this filtering approach
  

#Some examples of how to use the functions:
  #Create a subset of the data containing only data for taxon 170339 in tube 2 of treatment '1_C_16O'
    T170339T1T2 <- data.melted[data.melted$taxon==170339 & data.melted$trt.code=="1_C_16O" & data.melted$tube=="2",]
  #Calculate the weighted average density for taxon 170339 in tube 2:
    T170339T1T2.wad.out <- WAD.func(y=T170339T1T2$copies, x=T170339T1T2$density.g.ml)
    T170339T1T2.wad.out
  #Create a subset of the data containing only data for taxon 170339 in all tubes of treatment '1_C_16O'
    T170339T1 <- data.melted[data.melted$taxon==170339 & data.melted$trt.code=="1_C_16O",]
  #Calculate the bootstrapped (and observed) weighted average density estimates for taxon 170339 in treatment '1_C_16O'
    T170339T1.boot.out <- boot.WAD.func(X=T170339T1, vars=c("density.g.ml", "copies", "tube"), CI=0.90, draws=1000)
    T170339T1.boot.out
  #Create a subset of the data containing only data for taxon 170339 in all tubes of treatment '1_C_18O'
    T170339T2 <- data.melted[data.melted$taxon==170339 & data.melted$trt.code=="1_C_18O",]
  #Calculate the difference in weighted average density (with bootstrapped estimates of uncertainty) between the '1_C_18O' and '1_C_16O' treatments for taxon 170339
    T170339T1v2.out <- boot.diff.wad(T170339T1, T170339T2, vars=c("density.g.ml", "copies", "tube", "trt.code"), CI=0.90, draws=1000, tailed.test=1)
    T170339T1v2.out
  #Create a subset of the data containing only data for taxon 170339 in all tubes of the 'reference' treatments (i.e., '1_C_16O' & '1_P_12Cplus16O')
    T170339ref <- data.melted[data.melted$taxon==170339 & (data.melted$trt.code=="1_C_16O" | data.melted$trt.code=="1_P_12Cplus16O"),]
  #Calculate the molecular weight, GC content, and average number of carbon and oxygen atoms per DNA nucleotide for taxon 170339:
    T170339ref.MW.out <- MW.calc(X=T170339ref, vars=c("density.g.ml", "copies", "tube"))
    T170339ref.MW.out
  #Calculate the excess atom fraction of 18O (with bootstrapped estimates of uncertainty) in the '1_C_18O' vs. the '1_C_16O' treatment for taxon 170339
    T170339T1v2.ape.out <- boot.diff.ape(X.light=T170339T1, X.heavy=T170339T2, X.reference=T170339ref, iso.compare="18O", vars=c("density.g.ml", "copies", "tube", "trt.code"), CI=0.90, draws=1000)
    T170339T1v2.ape.out
  #Calculate the bootstrapped (and observed) weighted average density estimates for taxon 170339 in treatment '1_C_16O' along with the corresponding values for soil mass and total 16S copies from the tubes chosen for each bootstrap sample
    T170339T1.tube.out <- boot.TUBE.func(X=T170339T1, M.soil=Sdat, vars=c("density.g.ml", "copies", "tube", "g.soil"), v.frac=50, CI=0.90, draws=1000)
    T170339T1.tube.out
  #Note that the 'message' in the output of most of the functions warns if a taxon was not detected in a tube(s) (as was the case for taxon 4335437):
    T4335437T1 <- data.melted[data.melted$taxon==4335437 & data.melted$trt.code=="1_C_16O",]
    T4335437T1.boot.out <- boot.WAD.func(X=T4335437T1, vars=c("density.g.ml", "copies", "tube"), CI=0.90, draws=1000)
    T4335437T1.boot.out
    T4335437T2 <- data.melted[data.melted$taxon==4335437 & data.melted$trt.code=="1_C_18O",]
    T4335437T1v2.out <- boot.diff.wad(T4335437T1, T4335437T2, vars=c("density.g.ml", "copies", "tube", "trt.code"), CI=0.90, draws=1000, tailed.test=1)
    T4335437T1v2.out


#First create a small test data set for the all.taxa.calcs function that only has the first 20 OTUs:
  data.melted.1.20 <- data.melted[data.melted$taxon %in% levels(data.melted$taxon)[1:20],]
    #Re-convert the factor columns to factor:
      data.melted.1.20$taxon <- factor(data.melted.1.20$taxon)
      data.melted.1.20$SampleID <- factor(data.melted.1.20$SampleID)
      data.melted.1.20$trt.code <- factor(data.melted.1.20$trt.code)
      data.melted.1.20$tube <- factor(data.melted.1.20$tube)
      data.melted.1.20$kingdom <- factor(data.melted.1.20$kingdom)
      data.melted.1.20$phylum <- factor(data.melted.1.20$phylum)
      data.melted.1.20$class <- factor(data.melted.1.20$class)
      data.melted.1.20$order <- factor(data.melted.1.20$order)
      data.melted.1.20$family <- factor(data.melted.1.20$family)
      data.melted.1.20$genus <- factor(data.melted.1.20$genus)
      data.melted.1.20$species <- factor(data.melted.1.20$species)


#To run the all.taxa.calcs function on the small test data set created above, replace 'data.melted' with 'data.melted.1.20' below:


#Run all wad.diff, ape, r, & flux calculations for all taxa and all comparisons:
  #NOTES: r and f quantities are not returned when M.soil is unspecified
  #       the reference treatments for the comparisons are restricted to those in the same week of the experiment
  #       using set.seed ensures that subsequent function calls utilize the same random number seed, so results are replicated exactly)
  #       using system.time calculates how long it took the function to run
  #       rename the bootstrapped output files so that they are not overwritten if 'all.taxa.calcs' is re-run
  set.seed(100)
  system.time(all.comparisons <- all.taxa.calcs(X.all=data.melted, comparisons=Tcompare, vars=c("taxon", "density.g.ml", "copies", "tube", "trt.code", "g.soil"), growth.model="exponential", v.frac=50, CI=0.90, draws=1000, tailed.test=1))
  bootstrapped.filenames <- paste("qSIP_output/", c("bootstrapped_wad1.txt", "bootstrapped_wad2.txt", "bootstrapped_wad_diff.txt", "bootstrapped_ape.txt"), sep="")
  new.bootstrapped.filenames <- paste("qSIP_output/", "BH_OTU_", c("bootstrapped_wad1.txt", "bootstrapped_wad2.txt", "bootstrapped_wad_diff.txt", "bootstrapped_ape.txt"), sep="")
  file.rename(from=bootstrapped.filenames, to=new.bootstrapped.filenames)
  summary(all.comparisons)
  dim(all.comparisons)


#Write the results (all.comparisons) to a text file:
  dir.create(path=paste(getwd(), "/qSIP_output", sep=""), showWarnings=FALSE)
  write.table(all.comparisons, "qSIP_output/BH_OTU_all_comparisons.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)


#Write the taxa.id dataframe to a text file:
  write.table(taxa.id, "qSIP_output/BH_OTU_taxa_ID.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
