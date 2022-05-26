# This code imports the MH ANT qSIP data sets and performs the basic calculations on the data


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
    taxa.id <- read.table("qSIP_data/MH_ANT_taxa_id.txt", header=T, sep="\t", stringsAsFactors=T, na.strings="")
  #The number serves as a unique identifier code for each taxon to the "taxa.id" data frame:
    
  #NOTE: treatment code names cannot contain spaces
    data <- read.table("qSIP_data/MH_ANT_qSIP_data.txt", header=T, sep="\t", stringsAsFactors=T)
    names(data)
  #Note that columns 10 and up have a unique identifier code for taxon to match that in the "taxa.id" data frame
    names(data)[10:dim(data)[2]]
  #Rename these 'taxonID' columns to just the numbers (i.e., get rid of the leading 'X')
    names(data)[10:dim(data)[2]] <- gsub(pattern="X(\\d+)", replacement="\\1", x=names(data)[10:dim(data)[2]], perl=TRUE)

  #Identify & exclude outlier density fractions that should be excluded because water mixed in with the CsCl at the light (top) end:
    #Treatment "C":
      par(mfrow=c(2,2))
      for (i in 1:length(levels(factor(as.character(data$rep))))){
        plot(x=data$density.g.ml[data$combo.trt.code == "C_16O" & data$rep == i], y=data$DNA.ng.ul[data$combo.trt.code == "C_16O" & data$rep == i], pch=21, bg="blue", col="black", xlab="density (g/ml)", ylab="[DNA] (ng/ul)", main=paste("C_", "rep ", i, sep=""))
        points(x=data$density.g.ml[data$combo.trt.code == "C_18O" & data$rep == i], y=data$DNA.ng.ul[data$combo.trt.code == "C_18O" & data$rep == i], pch=21, bg="red", col="black", cex=0.6)
        abline(v=1.65) 
      }
      par(mfrow=c(1,1))

    #Treatment "Carb":
      par(mfrow=c(2,2))
      for (i in 1:length(levels(factor(as.character(data$rep))))){
        plot(x=data$density.g.ml[data$combo.trt.code == "Carb_16O" & data$rep == i], y=data$DNA.ng.ul[data$combo.trt.code == "Carb_16O" & data$rep == i], pch=21, bg="blue", col="black", xlab="density (g/ml)", ylab="[DNA] (ng/ul)", main=paste("Carb_", "rep ", i, sep=""))
        points(x=data$density.g.ml[data$combo.trt.code == "Carb_18O" & data$rep == i], y=data$DNA.ng.ul[data$combo.trt.code == "Carb_18O" & data$rep == i], pch=21, bg="red", col="black", cex=0.6)
        abline(v=1.65) 
      }
      par(mfrow=c(1,1))

    #Treatment "N":
      par(mfrow=c(2,2))
      for (i in 1:length(levels(factor(as.character(data$rep))))){
        plot(x=data$density.g.ml[data$combo.trt.code == "N_16O" & data$rep == i], y=data$DNA.ng.ul[data$combo.trt.code == "N_16O" & data$rep == i], pch=21, bg="blue", col="black", xlab="density (g/ml)", ylab="[DNA] (ng/ul)", main=paste("N_", "rep ", i, sep=""))
        points(x=data$density.g.ml[data$combo.trt.code == "N_18O" & data$rep == i], y=data$DNA.ng.ul[data$combo.trt.code == "N_18O" & data$rep == i], pch=21, bg="red", col="black", cex=0.6)
        abline(v=1.65) 
      }
      par(mfrow=c(1,1))

    #Treatment "pH":
      par(mfrow=c(2,2))
      for (i in 1:length(levels(factor(as.character(data$rep))))){
        plot(x=data$density.g.ml[data$combo.trt.code == "pH_16O" & data$rep == i], y=data$DNA.ng.ul[data$combo.trt.code == "pH_16O" & data$rep == i], pch=21, bg="blue", col="black", xlab="density (g/ml)", ylab="[DNA] (ng/ul)", main=paste("pH_", "rep ", i, sep=""))
        points(x=data$density.g.ml[data$combo.trt.code == "pH_18O" & data$rep == i], y=data$DNA.ng.ul[data$combo.trt.code == "pH_18O" & data$rep == i], pch=21, bg="red", col="black", cex=0.6)
        abline(v=1.65) 
      }
      par(mfrow=c(1,1))

    #Treatment "S":
      par(mfrow=c(2,2))
      for (i in 1:length(levels(factor(as.character(data$rep))))){
        plot(x=data$density.g.ml[data$combo.trt.code == "S_16O" & data$rep == i], y=data$DNA.ng.ul[data$combo.trt.code == "S_16O" & data$rep == i], pch=21, bg="blue", col="black", xlab="density (g/ml)", ylab="[DNA] (ng/ul)", main=paste("S_", "rep ", i, sep=""))
        points(x=data$density.g.ml[data$combo.trt.code == "S_18O" & data$rep == i], y=data$DNA.ng.ul[data$combo.trt.code == "S_18O" & data$rep == i], pch=21, bg="red", col="black", cex=0.6)
        abline(v=1.65) 
      }
      par(mfrow=c(1,1))

      graphics.off()

    #Exclude fractions that are less than 1.65 g/ml:
      data <- data[data$density.g.ml >= 1.65,]

  #Add a column for the sum of proportional abundance of all taxa by fraction
    data$sum.abundance <- rowSums(data[,10:ncol(data)])
    data$sum.abundance

  #Do not correct the proportional abundances of all taxa by scaling by their sum 
  #(they don't sum to 1, but that is presumably because proportional abundances of rare taxa are not given (because they were removed at the bioinformatics stage); even though these rare taxa are included in the total 'DNA.ng.ul'
    needs.correcting <- F
    if (needs.correcting) {
      data[10:(ncol(data)-1)] <- data[10:(ncol(data)-1)]/data$sum.abundance
    }

  #Calculate mass of DNA per uL, based on relative abundance and total mass of DNA per uL
  #NOTE: THIS CALCULATION MAY BE PROBLEMATIC; IT ASSUMES THAT THE NUMBER OF 16S COPIES IS PERFECTLY PROPORTIONAL TO THE MASS OF DNA
  #ALTHOUGH THESE TWO QUANTITIES ARE CORRELATED, THERE IS GREATER THAN ORDER-OF MAGNITUDE VARIATION IN DNA MASS FOR A GIVEN NUMBER OF 16S COPIES
    mass.ul <- data$DNA.ng.ul*data[,10:(ncol(data)-1)]
    mass.ul <- cbind(data[,1:9], mass.ul)  # add first 9 columns of data to mass.ul

  #Melt data into long format by tube, sample, tmt, isotope, combo.trt.code, rep, fraction, DNA conc, and density;
  #Do this for mass.ul and for relative abundance. Merge these to into 1 masterfile: data.melted
    mass.ul.melted <- melt(mass.ul, id=c("tube", "sample", "tmt", "isotope", "combo.trt.code", "rep", "fraction", "density.g.ml", "DNA.ng.ul"), measure.vars=as.character(1:310), variable.name="taxon", value.name="mass.ul")
    rel.abundance.melted <- melt(data, id=c("tube", "sample", "tmt", "isotope", "combo.trt.code", "rep", "fraction", "density.g.ml", "DNA.ng.ul"), measure.vars=as.character(1:310), variable.name="taxon", value.name="rel.abundance")
    data.melted <- merge(mass.ul.melted, rel.abundance.melted)

  #Merge taxa data and reorder data frame by taxon and SampleID
    data.melted <- merge(data.melted, taxa.id)
    data.melted <- data.melted[order(data.melted$taxon, data.melted$sample),]


#Import data frame containing the treatment comparisons to perform:
#NOTE: THE DURATION OF THE INCUBATION 'days' IS MADE-UP TO ENABLE RUNNING 'all.taxa.calcs' BELOW:
  Tcompare <- read.table("qSIP_data/MH_ANT_TreatmentComparisons.txt", header=TRUE, sep="\t", colClasses=c("numeric","factor","character","character","character","numeric","character"))
  summary(Tcompare)
  Tcompare


#Run all wad.diff, ape, r, & flux calculations for all taxa and all comparisons:
  #NOTES: only wad and ape results are valid; r, f (C fluxes) results are not calculated here (soil data is not specified); note that incubation duration values in Tcompare are made-up to enable running the 'all.taxa.calcs' script (but the made-up durations do not affect any of the results that are returned)
  #       using set.seed ensures that subsequent function calls utilize the same random number seed, so results are replicated exactly)
  #       using system.time calculates how long it took the function to run
  #       set prop.O.from.water=0.33 according to McHugh et al. unpublished data from E. coli pure cultures incubated with 18O-H2O
  #       rename the bootstrapped output files so that they are not overwritten if 'all.taxa.calcs' is re-run
  set.seed(100)
  system.time(all.comparisons <- all.taxa.calcs(X.all=data.melted, comparisons=Tcompare, vars=c("taxon", "density.g.ml", "mass.ul", "tube", "combo.trt.code", "g.soil"), growth.model="exponential", prop.O.from.water=0.33, v.frac=50, copies.cell=6, pgC.cell=0.1, CI=0.90, draws=1000, tailed.test=1))
  bootstrapped.filenames <- paste("qSIP_output/", c("bootstrapped_wad1.txt", "bootstrapped_wad2.txt", "bootstrapped_wad_diff.txt", "bootstrapped_ape.txt"), sep="")
  new.bootstrapped.filenames <- paste("qSIP_output/", "MH_ANT_", c("bootstrapped_wad1.txt", "bootstrapped_wad2.txt", "bootstrapped_wad_diff.txt", "bootstrapped_ape.txt"), sep="")
  file.rename(from=bootstrapped.filenames, to=new.bootstrapped.filenames)
  summary(all.comparisons)
  dim(all.comparisons)


#Write the results (all.comparisons) to a text file:
  dir.create(path=paste(getwd(), "/qSIP_output", sep=""), showWarnings=FALSE)
  write.table(all.comparisons, "qSIP_output/MH_ANT_all_comparisons.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)


#Write the taxa.id dataframe to a text file:
  write.table(taxa.id, "qSIP_output/MH_ANT_taxa_ID.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
