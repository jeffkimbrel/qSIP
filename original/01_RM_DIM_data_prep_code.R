# This code imports the Dimensions week 1 qSIP data sets and prepares the data for qSIP calculations,
# including calculating tube-level WADs for all OTUs and analyzing the shifts among tubes
# There are 20 comparisons, so we will subset the comparisons per ecosystem (5 comparisons each)
# This code only performs grassland (GL) calculations


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
    taxa.id <- read.table("qSIP_data/RM_DIM_L7_taxa_ids.txt", header=T, sep="\t", stringsAsFactors=T, na.strings="")
    head(taxa.id)
    names(taxa.id)
  #The number serves as a unique identifier code for each taxon to the "taxa.id" data frame:
    
  #NOTE: treatment code names cannot contain spaces
    data <- read.table("qSIP_data/RM_DIM_final_qSIP_data_file.txt", header=T, sep="\t", stringsAsFactors=T, check.names=F)
    dim(data)
    names(data)
    head(data)
    
  #Add a column for the sum of proportional abundance of all taxa by fraction
    data$sum.abundance <- rowSums(data[,7:ncol(data)])
    data$sum.abundance

  #Calculate number of copies per uL, based on relative abundance and total number of copies per uL:
    ncopies <- data$neat.avg.16S.copies*data[,7:(ncol(data)-1)]
    dim(ncopies)
    ncopies <- cbind(data[,1:6], ncopies)  # add first 6 columns of data to ncopies
    dim(ncopies)
    head(ncopies)
    
  #Melt data into long format by tube, sample, tmt, rep, fraction, DNA conc, and density;
  #Do this for copies.ul and for relative abundance, which is just our data file. Merge these to into 1 masterfile: data.melted
    ncopies.melted <- melt(ncopies, id=c("Sample", "fraction", "iso.treat.eco", "neat.avg.16S.copies", "Density.g.ml", "DNA.ng.ul"), measure.vars=as.character(1:length(taxa.id$taxon)), variable.name="taxon", value.name="copies.ul")
    rel.abundance.melted <- melt(data, id=c("Sample", "fraction", "iso.treat.eco", "neat.avg.16S.copies", "Density.g.ml", "DNA.ng.ul"), measure.vars=as.character(1:length(taxa.id$taxon)), variable.name="taxon", value.name="rel.abundance")
    data.melted <- merge(ncopies.melted, rel.abundance.melted)
    head(data.melted)

  #Merge taxa data and reorder data frame by taxon and SampleID AND FRACTION
    data.melted <- merge(data.melted, taxa.id)
    data.melted <- data.melted[order(data.melted$taxon, data.melted$Sample, data.melted$fraction),]
    row.names(data.melted) <- 1:dim(data.melted)[1]     ####BJK - RENAME OBSERVATIONS TO BE SEQUENTIAL
    head(data.melted)


  #Import data frame containing the treatment comparisons to perform:
  #NOTE: THE DURATION OF THE INCUBATION 'days' IS MADE-UP TO ENABLE RUNNING 'all.taxa.calcs' BELOW:
    Tcompare.GL <- read.table("qSIP_data/RM_DIM_TreatmentComparisons_GL.txt", header=TRUE, sep="\t", colClasses=c("numeric","factor","character","character","character","numeric","character"))
    summary(Tcompare.GL)
    Tcompare.GL


  #Read in the list of all the sample-fractions prepped for sequencing:
    seqd.data <- read.table("qSIP_data/RM_DIM_Sequenced_Fractions.txt", header=T, sep="\t", stringsAsFactors=F, check.names=F)


#Calculate the number of missing fractions per replicate based on sequencing failures:
  #Create a vector of sample-fractions that were submitted for sequencing:
    seqd <- seqd.data$sample_fraction

  #Create a vector of sample-fractions that passed sequencing (i.e., those that ended up in the qSIP data set):
    passed <- paste(data$Sample, data$fraction, sep="_")

  #Create a vector of sample-fractions that failed sequencing:
    failed <- setdiff(seqd, passed)

  #Numbers of sample-fractions that were sequenced, passed, and failed:
    length(seqd)
    length(passed)
    length(failed)

  #Create vectors of tubes (replicates) that were sequenced, passed, and failed:
    seqd.tube <- gsub(pattern="(.*)\\_(.*)", replacement="\\1", x=seqd, perl=TRUE)
    passed.tube <- gsub(pattern="(.*)\\_(.*)", replacement="\\1", x=passed, perl=TRUE)
    failed.tube <- gsub(pattern="(.*)\\_(.*)", replacement="\\1", x=failed, perl=TRUE)

  #Put all of this information together into a data.frame that lists the number of missing fractions per replicate based on sequencing failures:
    fracs.per.rep <- data.frame(matrix(NA, nrow=length(levels(data$Sample)), ncol=5))
    names(fracs.per.rep) <- c("tube", "iso.treat.eco", "num.seqd.fracs", "num.passed.fracs", "num.missing.fracs")
    for(i in 1: length(levels(data$Sample))){
      fracs.per.rep$tube[i] <- levels(data$Sample)[i]
      fracs.per.rep$iso.treat.eco[i] <- unique(as.character(data$iso.treat.eco[data$Sample == levels(data$Sample)[i]]))
      fracs.per.rep$num.seqd.fracs[i] <- sum(seqd.tube == levels(data$Sample)[i])
      fracs.per.rep$num.passed.fracs[i] <- sum(passed.tube == levels(data$Sample)[i])
      fracs.per.rep$num.missing.fracs[i] <- sum(failed.tube == levels(data$Sample)[i])
    }
    fracs.per.rep$tube <- factor(fracs.per.rep$tube)
    fracs.per.rep$iso.treat.eco <- factor(fracs.per.rep$iso.treat.eco)
    fracs.per.rep


#Write the fracs.per.rep dataframe to a text file:
  write.table(fracs.per.rep, "qSIP_output/RM_DIM_fractions_per_replicate.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)


#Plot tube-level density curves (all taxa lumped together):
  #####{
  #Isolate the data for each replicate of the "16O.GL" treatment:
    id.reps(DATA=ncopies, focal.tmt="16O.GL", vars=c("Sample", "iso.treat.eco"))
    unlab.1 <- select.rep(DATA=ncopies, focal.tmt="16O.GL", replicate.index=1, vars=c("Density.g.ml", "neat.avg.16S.copies", "Sample", "iso.treat.eco"))
    unlab.2 <- select.rep(DATA=ncopies, focal.tmt="16O.GL", replicate.index=2, vars=c("Density.g.ml", "neat.avg.16S.copies", "Sample", "iso.treat.eco"))
    unlab.3 <- select.rep(DATA=ncopies, focal.tmt="16O.GL", replicate.index=3, vars=c("Density.g.ml", "neat.avg.16S.copies", "Sample", "iso.treat.eco"))
  #Isolate the data for each replicate of the "12C.GL" treatment:
    id.reps(DATA=ncopies, focal.tmt="12C.GL", vars=c("Sample", "iso.treat.eco"))
    unlab.4 <- select.rep(DATA=ncopies, focal.tmt="12C.GL", replicate.index=1, vars=c("Density.g.ml", "neat.avg.16S.copies", "Sample", "iso.treat.eco"))
    unlab.5 <- select.rep(DATA=ncopies, focal.tmt="12C.GL", replicate.index=2, vars=c("Density.g.ml", "neat.avg.16S.copies", "Sample", "iso.treat.eco"))
    unlab.6 <- select.rep(DATA=ncopies, focal.tmt="12C.GL", replicate.index=3, vars=c("Density.g.ml", "neat.avg.16S.copies", "Sample", "iso.treat.eco"))
  #Isolate the data for each replicate of the "12C_N.GL" treatment:
    id.reps(DATA=ncopies, focal.tmt="12C_N.GL", vars=c("Sample", "iso.treat.eco"))
    unlab.7 <- select.rep(DATA=ncopies, focal.tmt="12C_N.GL", replicate.index=1, vars=c("Density.g.ml", "neat.avg.16S.copies", "Sample", "iso.treat.eco"))
    unlab.8 <- select.rep(DATA=ncopies, focal.tmt="12C_N.GL", replicate.index=2, vars=c("Density.g.ml", "neat.avg.16S.copies", "Sample", "iso.treat.eco"))
    unlab.9 <- select.rep(DATA=ncopies, focal.tmt="12C_N.GL", replicate.index=3, vars=c("Density.g.ml", "neat.avg.16S.copies", "Sample", "iso.treat.eco"))
  #Isolate the data for each replicate of the "18O.GL" treatment:
    id.reps(DATA=ncopies, focal.tmt="18O.GL", vars=c("Sample", "iso.treat.eco"))
    lab.1 <- select.rep(DATA=ncopies, focal.tmt="18O.GL", replicate.index=1, vars=c("Density.g.ml", "neat.avg.16S.copies", "Sample", "iso.treat.eco"))
    lab.2 <- select.rep(DATA=ncopies, focal.tmt="18O.GL", replicate.index=2, vars=c("Density.g.ml", "neat.avg.16S.copies", "Sample", "iso.treat.eco"))
    lab.3 <- select.rep(DATA=ncopies, focal.tmt="18O.GL", replicate.index=3, vars=c("Density.g.ml", "neat.avg.16S.copies", "Sample", "iso.treat.eco"))

  #Plot density curves for all unlabeled replicates:
    par(mai=c(1.02, 1.02, 0.22, 0.22))
    XLIM <- c(min(c(unlab.1$X, unlab.2$X, unlab.3$X, unlab.4$X, unlab.5$X, unlab.6$X, unlab.7$X, unlab.8$X, unlab.9$X, lab.1$X, lab.2$X, lab.3$X)), max(c(unlab.1$X, unlab.2$X, unlab.3$X, unlab.4$X, unlab.5$X, unlab.6$X, unlab.7$X, unlab.8$X, unlab.9$X, lab.1$X, lab.2$X, lab.3$X)))
    YLIM <- c(min(c(unlab.1$Y, unlab.2$Y, unlab.3$Y, unlab.4$Y, unlab.5$Y, unlab.6$Y, unlab.7$Y, unlab.8$Y, unlab.9$Y, lab.1$Y, lab.2$Y, lab.3$Y)), max(c(unlab.1$Y, unlab.2$Y, unlab.3$Y, unlab.4$Y, unlab.5$Y, unlab.6$Y, unlab.7$Y, unlab.8$Y, unlab.9$Y, lab.1$Y, lab.2$Y, lab.3$Y)))
    plot(unlab.1$Y~unlab.1$X, type="n", pch=21, col="black", bg="green", xlim=XLIM, ylim=YLIM, xlab="Density g/ml", ylab="16S copies", cex=1.25, cex.axis=1.5, cex.lab=1.8)
    points(unlab.1$Y~unlab.1$X, pch=21, cex=1.25, col="black", bg="green")
    points(unlab.2$Y~unlab.2$X, pch=21, cex=1.25, col="black", bg="red")
    points(unlab.3$Y~unlab.3$X, pch=21, cex=1.25, col="black", bg="yellow")
    points(unlab.4$Y~unlab.4$X, pch=21, cex=1.25, col="black", bg="purple")
    points(unlab.5$Y~unlab.5$X, pch=21, cex=1.25, col="black", bg="blue")
    points(unlab.6$Y~unlab.6$X, pch=21, cex=1.25, col="black", bg="orange")
    points(unlab.7$Y~unlab.7$X, pch=21, cex=1.25, col="black", bg="cyan")
    points(unlab.8$Y~unlab.8$X, pch=21, cex=1.25, col="black", bg="magenta")
    points(unlab.9$Y~unlab.9$X, pch=21, cex=1.25, col="black", bg="wheat")
    abline(v=WAD.func(unlab.1$Y, unlab.1$X), cex=1.25, col="green", lwd=2)
    abline(v=WAD.func(unlab.2$Y, unlab.2$X), cex=1.25, col="red", lwd=2)
    abline(v=WAD.func(unlab.3$Y, unlab.3$X), cex=1.25, col="yellow", lwd=2)
    abline(v=WAD.func(unlab.4$Y, unlab.4$X), cex=1.25, col="purple", lwd=2)
    abline(v=WAD.func(unlab.5$Y, unlab.5$X), cex=1.25, col="blue", lwd=2)
    abline(v=WAD.func(unlab.6$Y, unlab.6$X), cex=1.25, col="orange", lwd=2)
    abline(v=WAD.func(unlab.7$Y, unlab.7$X), cex=1.25, col="cyan", lwd=2)
    abline(v=WAD.func(unlab.8$Y, unlab.8$X), cex=1.25, col="magenta", lwd=2)
    abline(v=WAD.func(unlab.9$Y, unlab.9$X), cex=1.25, col="wheat", lwd=2)

  #Plot a density curve for one unlabeled replicate:
    par(mai=c(1.02, 1.02, 0.22, 0.22))
    XLIM <- c(min(c(unlab.1$X, unlab.2$X, unlab.3$X, unlab.4$X, unlab.5$X, unlab.6$X, unlab.7$X, unlab.8$X, unlab.9$X, lab.1$X, lab.2$X, lab.3$X)), max(c(unlab.1$X, unlab.2$X, unlab.3$X, unlab.4$X, unlab.5$X, unlab.6$X, unlab.7$X, unlab.8$X, unlab.9$X, lab.1$X, lab.2$X, lab.3$X)))
    YLIM <- c(min(c(unlab.1$Y, unlab.2$Y, unlab.3$Y, unlab.4$Y, unlab.5$Y, unlab.6$Y, unlab.7$Y, unlab.8$Y, unlab.9$Y, lab.1$Y, lab.2$Y, lab.3$Y)), max(c(unlab.1$Y, unlab.2$Y, unlab.3$Y, unlab.4$Y, unlab.5$Y, unlab.6$Y, unlab.7$Y, unlab.8$Y, unlab.9$Y, lab.1$Y, lab.2$Y, lab.3$Y)))
    plot(unlab.1$Y~unlab.1$X, type="n", pch=21, col="black", bg="green", xlim=XLIM, ylim=YLIM, xlab="Density g/ml", ylab="16S copies", cex=1.25, cex.axis=1.5, cex.lab=1.8)
    points(unlab.8$Y~unlab.8$X, pch=21, cex=1.25, col="black", bg="magenta")
    abline(v=WAD.func(unlab.8$Y, unlab.8$X), cex=1.25, col="magenta", lwd=2)

  #Plot a density curve for one labeled replicate:
    par(mai=c(1.02, 1.02, 0.22, 0.22))
    XLIM <- c(min(c(unlab.1$X, unlab.2$X, unlab.3$X, unlab.4$X, unlab.5$X, unlab.6$X, unlab.7$X, unlab.8$X, unlab.9$X, lab.1$X, lab.2$X, lab.3$X)), max(c(unlab.1$X, unlab.2$X, unlab.3$X, unlab.4$X, unlab.5$X, unlab.6$X, unlab.7$X, unlab.8$X, unlab.9$X, lab.1$X, lab.2$X, lab.3$X)))
    YLIM <- c(min(c(unlab.1$Y, unlab.2$Y, unlab.3$Y, unlab.4$Y, unlab.5$Y, unlab.6$Y, unlab.7$Y, unlab.8$Y, unlab.9$Y, lab.1$Y, lab.2$Y, lab.3$Y)), max(c(unlab.1$Y, unlab.2$Y, unlab.3$Y, unlab.4$Y, unlab.5$Y, unlab.6$Y, unlab.7$Y, unlab.8$Y, unlab.9$Y, lab.1$Y, lab.2$Y, lab.3$Y)))
    plot(unlab.1$Y~unlab.1$X, type="n", pch=21, col="black", bg="green", xlim=XLIM, ylim=YLIM, xlab="Density g/ml", ylab="16S copies", cex=1.25, cex.axis=1.5, cex.lab=1.8)
    points(lab.3$Y~lab.3$X, pch=21, cex=1.25, col="black", bg="brown")
    abline(v=WAD.func(lab.2$Y, lab.2$X), cex=1.25, col="brown", lwd=2)

  #Plot a density curve for one labeled and one unlabeled replicate:
    par(mai=c(1.02, 1.02, 0.22, 0.22))
    XLIM <- c(min(c(unlab.1$X, unlab.2$X, unlab.3$X, unlab.4$X, unlab.5$X, unlab.6$X, unlab.7$X, unlab.8$X, unlab.9$X, lab.1$X, lab.2$X, lab.3$X)), max(c(unlab.1$X, unlab.2$X, unlab.3$X, unlab.4$X, unlab.5$X, unlab.6$X, unlab.7$X, unlab.8$X, unlab.9$X, lab.1$X, lab.2$X, lab.3$X)))
    YLIM <- c(min(c(unlab.1$Y, unlab.2$Y, unlab.3$Y, unlab.4$Y, unlab.5$Y, unlab.6$Y, unlab.7$Y, unlab.8$Y, unlab.9$Y, lab.1$Y, lab.2$Y, lab.3$Y)), max(c(unlab.1$Y, unlab.2$Y, unlab.3$Y, unlab.4$Y, unlab.5$Y, unlab.6$Y, unlab.7$Y, unlab.8$Y, unlab.9$Y, lab.1$Y, lab.2$Y, lab.3$Y)))
    plot(unlab.1$Y~unlab.1$X, type="n", pch=21, col="black", bg="green", xlim=XLIM, ylim=YLIM, xlab="Density g/ml", ylab="16S copies", cex=1.25, cex.axis=1.5, cex.lab=1.8)
    points(unlab.8$Y~unlab.8$X, pch=21, cex=1.25, col="black", bg="magenta")
    points(lab.3$Y~lab.3$X, pch=21, cex=1.25, col="black", bg="brown")
    abline(v=WAD.func(unlab.8$Y, unlab.8$X), cex=1.25, col="magenta", lwd=2)
    abline(v=WAD.func(lab.2$Y, lab.2$X), cex=1.25, col="brown", lwd=2)

  #Plot density curves for two unlabeled replicates:
    par(mai=c(1.02, 1.02, 0.22, 0.22))
    XLIM <- c(min(c(unlab.1$X, unlab.2$X, unlab.3$X, unlab.4$X, unlab.5$X, unlab.6$X, unlab.7$X, unlab.8$X, unlab.9$X, lab.1$X, lab.2$X, lab.3$X)), max(c(unlab.1$X, unlab.2$X, unlab.3$X, unlab.4$X, unlab.5$X, unlab.6$X, unlab.7$X, unlab.8$X, unlab.9$X, lab.1$X, lab.2$X, lab.3$X)))
    YLIM <- c(min(c(unlab.1$Y, unlab.2$Y, unlab.3$Y, unlab.4$Y, unlab.5$Y, unlab.6$Y, unlab.7$Y, unlab.8$Y, unlab.9$Y, lab.1$Y, lab.2$Y, lab.3$Y)), max(c(unlab.1$Y, unlab.2$Y, unlab.3$Y, unlab.4$Y, unlab.5$Y, unlab.6$Y, unlab.7$Y, unlab.8$Y, unlab.9$Y, lab.1$Y, lab.2$Y, lab.3$Y)))
    plot(unlab.1$Y~unlab.1$X, type="n", pch=21, col="black", bg="green", xlim=XLIM, ylim=YLIM, xlab="Density g/ml", ylab="16S copies", cex=1.25, cex.axis=1.5, cex.lab=1.8)
    points(unlab.4$Y~unlab.4$X, pch=21, cex=1.25, col="black", bg="purple")
    points(unlab.8$Y~unlab.8$X, pch=21, cex=1.25, col="black", bg="magenta")
    abline(v=WAD.func(unlab.4$Y, unlab.4$X), cex=1.25, col="purple", lwd=2)
    abline(v=WAD.func(unlab.8$Y, unlab.8$X), cex=1.25, col="magenta", lwd=2)

  #Plot density curves for all labeled replicates:
    par(mai=c(1.02, 1.02, 0.22, 0.22))
    XLIM <- c(min(c(unlab.1$X, unlab.2$X, unlab.3$X, unlab.4$X, unlab.5$X, unlab.6$X, unlab.7$X, unlab.8$X, unlab.9$X, lab.1$X, lab.2$X, lab.3$X)), max(c(unlab.1$X, unlab.2$X, unlab.3$X, unlab.4$X, unlab.5$X, unlab.6$X, unlab.7$X, unlab.8$X, unlab.9$X, lab.1$X, lab.2$X, lab.3$X)))
    YLIM <- c(min(c(unlab.1$Y, unlab.2$Y, unlab.3$Y, unlab.4$Y, unlab.5$Y, unlab.6$Y, unlab.7$Y, unlab.8$Y, unlab.9$Y, lab.1$Y, lab.2$Y, lab.3$Y)), max(c(unlab.1$Y, unlab.2$Y, unlab.3$Y, unlab.4$Y, unlab.5$Y, unlab.6$Y, unlab.7$Y, unlab.8$Y, unlab.9$Y, lab.1$Y, lab.2$Y, lab.3$Y)))
    plot(unlab.1$Y~unlab.1$X, type="n", pch=21, col="black", bg="green", xlim=XLIM, ylim=YLIM, xlab="Density g/ml", ylab="16S copies", cex=1.25, cex.axis=1.5, cex.lab=1.8)
    points(lab.1$Y~lab.1$X, pch=21, cex=1.25, col="black", bg="gray")
    points(lab.2$Y~lab.2$X, pch=21, cex=1.25, col="black", bg="brown")
    points(lab.3$Y~lab.3$X, pch=21, cex=1.25, col="black", bg="darkgoldenrod2")
    abline(v=WAD.func(lab.1$Y, lab.1$X), cex=1.25, col="gray", lwd=2)
    abline(v=WAD.func(lab.2$Y, lab.2$X), cex=1.25, col="brown", lwd=2)
    abline(v=WAD.func(lab.3$Y, lab.3$X), cex=1.25, col="darkgoldenrod2", lwd=2)

  #Plot density curves for two labeled replicates:
    par(mai=c(1.02, 1.02, 0.22, 0.22))
    XLIM <- c(min(c(unlab.1$X, unlab.2$X, unlab.3$X, unlab.4$X, unlab.5$X, unlab.6$X, unlab.7$X, unlab.8$X, unlab.9$X, lab.1$X, lab.2$X, lab.3$X)), max(c(unlab.1$X, unlab.2$X, unlab.3$X, unlab.4$X, unlab.5$X, unlab.6$X, unlab.7$X, unlab.8$X, unlab.9$X, lab.1$X, lab.2$X, lab.3$X)))
    YLIM <- c(min(c(unlab.1$Y, unlab.2$Y, unlab.3$Y, unlab.4$Y, unlab.5$Y, unlab.6$Y, unlab.7$Y, unlab.8$Y, unlab.9$Y, lab.1$Y, lab.2$Y, lab.3$Y)), max(c(unlab.1$Y, unlab.2$Y, unlab.3$Y, unlab.4$Y, unlab.5$Y, unlab.6$Y, unlab.7$Y, unlab.8$Y, unlab.9$Y, lab.1$Y, lab.2$Y, lab.3$Y)))
    plot(unlab.1$Y~unlab.1$X, type="n", pch=21, col="black", bg="green", xlim=XLIM, ylim=YLIM, xlab="Density g/ml", ylab="16S copies", cex=1.25, cex.axis=1.5, cex.lab=1.8)
    points(lab.2$Y~lab.2$X, pch=21, cex=1.25, col="black", bg="gray")
    points(lab.3$Y~lab.3$X, pch=21, cex=1.25, col="black", bg="brown")
    abline(v=WAD.func(lab.1$Y, lab.1$X), cex=1.25, col="gray", lwd=2)
    abline(v=WAD.func(lab.2$Y, lab.2$X), cex=1.25, col="brown", lwd=2)

  #Plot a density curve for one labeled and one unlabaled replicate with labeled tube less than unlabeled tube!:
    par(mai=c(1.02, 1.02, 0.22, 0.22))
    XLIM <- c(min(c(unlab.1$X, unlab.2$X, unlab.3$X, unlab.4$X, unlab.5$X, unlab.6$X, unlab.7$X, unlab.8$X, unlab.9$X, lab.1$X, lab.2$X, lab.3$X)), max(c(unlab.1$X, unlab.2$X, unlab.3$X, unlab.4$X, unlab.5$X, unlab.6$X, unlab.7$X, unlab.8$X, unlab.9$X, lab.1$X, lab.2$X, lab.3$X)))
    YLIM <- c(min(c(unlab.1$Y, unlab.2$Y, unlab.3$Y, unlab.4$Y, unlab.5$Y, unlab.6$Y, unlab.7$Y, unlab.8$Y, unlab.9$Y, lab.1$Y, lab.2$Y, lab.3$Y)), max(c(unlab.1$Y, unlab.2$Y, unlab.3$Y, unlab.4$Y, unlab.5$Y, unlab.6$Y, unlab.7$Y, unlab.8$Y, unlab.9$Y, lab.1$Y, lab.2$Y, lab.3$Y)))
    plot(unlab.1$Y~unlab.1$X, type="n", pch=21, col="black", bg="green", xlim=XLIM, ylim=YLIM, xlab="Density g/ml", ylab="16S copies", cex=1.25, cex.axis=1.5, cex.lab=1.8)
    points(unlab.4$Y~unlab.4$X, pch=21, cex=1.25, col="black", bg="purple")
    points(lab.2$Y~lab.2$X, pch=21, cex=1.25, col="black", bg="gray")
    abline(v=WAD.func(lab.1$Y, lab.1$X), cex=1.25, col="gray", lwd=2)
    abline(v=WAD.func(unlab.4$Y, unlab.4$X), cex=1.25, col="purple", lwd=2)
  #####}

graphics.off()

#Calculate WADs for each tube and taxon (Grassland only):
  data.melted.GL <- data.melted[data.melted$iso.treat.eco %in% c("12C_18O_N.GL", "12C_18O.GL", "12C_N.GL", "12C.GL", "13C_N.GL", "13C.GL", "16O.GL", "18O.GL"),]
  data.melted.GL$taxon <- factor(data.melted.GL$taxon)
  data.melted.GL$Sample <- factor(data.melted.GL$Sample)
  data.melted.GL$iso.treat.eco <- factor(data.melted.GL$iso.treat.eco)
  row.names(data.melted.GL) <- 1:dim(data.melted.GL)[1]
  WAD.by.taxon.GL <- WAD.by.taxon.func(X=data.melted.GL, vars=c("taxon", "Density.g.ml", "copies.ul", "Sample", "iso.treat.eco"))
  #Look at the results:
    names(WAD.by.taxon.GL)          #names of the two data frames in the output list
    head(WAD.by.taxon.GL$obs.wads)  #looking at the head of the first data frame
    head(WAD.by.taxon.GL[[1]])      #another way to look at the head of the first data frame
    WAD.by.taxon.GL$reps.by.trt     #looking at the second data frame
    WAD.by.taxon.GL[[2]]            #another way to look at the second data frame
  
  #Quick and dirty look at the level of replication among treatments and the occurrences of data for taxa in tubes 
    table(data.melted.GL$iso.treat.eco, data.melted.GL$Sample)
    table(data.melted.GL$taxon, data.melted.GL$Sample)


#Unsorted line graph of WADs -- taxon numbers are the same across tubes -- (all unlabeled & 18O replicates: 4 trts * 3 reps/trt = 12 tubes):
#(note: because taxa are unsorted, excluded taxa (i.e., ones that did not occur in all replicates) are skipped when drawing lines between the taxa that are present in all replicates)
#####{
  graphics.off()
  dev.new(width=6, height=7.74)
  par(mai=c(0.44,0.57,0.22,1.60), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
    #Select only the data for taxa common to all replicates of the specified treatments:
    obs.wads.by.taxon <- WAD.by.taxon.GL$obs.wads
    reps.by.trt <- WAD.by.taxon.GL$reps.by.trt
    treatments=c("16O.GL", "12C.GL", "12C_N.GL", "18O.GL")
    reps.by.trt <- reps.by.trt[reps.by.trt[,1] %in% treatments,]
      reps.by.trt[,1] <- factor(reps.by.trt[,1])
      reps.by.trt[,2] <- factor(reps.by.trt[,2])
      reps.by.trt[,3] <- factor(reps.by.trt[,3])
      reps.by.trt[,4] <- factor(reps.by.trt[,4])
    curr.data <- cbind(taxon=obs.wads.by.taxon[,1], obs.wads.by.taxon[,names(obs.wads.by.taxon) %in% as.character(unlist(reps.by.trt[,2:4]))])
      #Remove any taxon that does not occur in all 12 replicates:
      curr.data <- curr.data[apply(is.na(curr.data[2:dim(curr.data)[2]]), 1, sum) == 0, ]
      curr.data$taxon <- factor(curr.data$taxon)
    #Set the x-limits for the graph:
    x.min <- min(curr.data[,2:dim(curr.data)[2]], na.rm=TRUE)
    x.max <- max(curr.data[,2:dim(curr.data)[2]], na.rm=TRUE)
    AT.y <- seq(0, max(as.numeric(as.character(obs.wads.by.taxon$taxon))), 50)
  plot(x=curr.data[,2], y=as.numeric(as.character(curr.data$taxon)), type="n", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min, x.max))
  par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
  axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
  mtext(expression(paste("Observed WAD (g ml"^-1, ")", sep="")), side=1, line=1.35, cex=0.75)
  par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
  axis(side=2, at=AT.y, labels=AT.y, tck=-0.015, las=1, cex.axis=0.6)
  mtext("Taxon", side=2, line=2.35, cex=0.75)
  #16O.GL
    # points(x=curr.data[,2], y=as.numeric(as.character(curr.data$taxon)), pch=21, cex=0.7, col="black", bg="green")           #R1
    # points(x=curr.data[,9], y=as.numeric(as.character(curr.data$taxon)), pch=21, cex=0.7, col="black", bg="red")             #R2
    # points(x=curr.data[,10], y=as.numeric(as.character(curr.data$taxon)), pch=21, cex=0.7, col="black", bg="yellow")         #R3
    lines(x=curr.data[,2], y=as.numeric(as.character(curr.data$taxon)), col="green")
    lines(x=curr.data[,9], y=as.numeric(as.character(curr.data$taxon)), col="red")
    lines(x=curr.data[,10], y=as.numeric(as.character(curr.data$taxon)), col="yellow")
  #12C.GL
    # points(x=curr.data[,11], y=as.numeric(as.character(curr.data$taxon)), pch=21, cex=0.7, col="black", bg="purple")         #R1
    # points(x=curr.data[,12], y=as.numeric(as.character(curr.data$taxon)), pch=21, cex=0.7, col="black", bg="blue")           #R2
    # points(x=curr.data[,13], y=as.numeric(as.character(curr.data$taxon)), pch=21, cex=0.7, col="black", bg="orange")         #R3
    lines(x=curr.data[,11], y=as.numeric(as.character(curr.data$taxon)), col="purple")
    lines(x=curr.data[,12], y=as.numeric(as.character(curr.data$taxon)), col="blue")
    lines(x=curr.data[,13], y=as.numeric(as.character(curr.data$taxon)), col="orange")
  #12C_N.GL
    # points(x=curr.data[,3], y=as.numeric(as.character(curr.data$taxon)), pch=21, cex=0.7, col="black", bg="cyan")            #R1
    # points(x=curr.data[,4], y=as.numeric(as.character(curr.data$taxon)), pch=21, cex=0.7, col="black", bg="magenta")         #R2
    # points(x=curr.data[,5], y=as.numeric(as.character(curr.data$taxon)), pch=21, cex=0.7, col="black", bg="wheat")           #R3
    lines(x=curr.data[,3], y=as.numeric(as.character(curr.data$taxon)), col="cyan")
    lines(x=curr.data[,4], y=as.numeric(as.character(curr.data$taxon)), col="magenta")
    lines(x=curr.data[,5], y=as.numeric(as.character(curr.data$taxon)), col="wheat")
  #18O.GL
    # points(x=curr.data[,6], y=as.numeric(as.character(curr.data$taxon)), pch=24, cex=0.7, col="black", bg="gray")            #R1
    # points(x=curr.data[,7], y=as.numeric(as.character(curr.data$taxon)), pch=24, cex=0.7, col="black", bg="brown")           #R2
    # points(x=curr.data[,8], y=as.numeric(as.character(curr.data$taxon)), pch=24, cex=0.7, col="black", bg="darkgoldenrod2")  #R3
    lines(x=curr.data[,6], y=as.numeric(as.character(curr.data$taxon)), col="gray")
    lines(x=curr.data[,7], y=as.numeric(as.character(curr.data$taxon)), col="brown")
    lines(x=curr.data[,8], y=as.numeric(as.character(curr.data$taxon)), col="darkgoldenrod2")
    par(xpd=NA)
    ink <- c("green", "red", "yellow", "purple", "blue", "orange", "cyan", "magenta", "wheat", "gray", "brown", "darkgoldenrod2")
    legend(x=x.max*1.0055, y=max(as.numeric(as.character(curr.data$taxon)))+(0.04*(max(as.numeric(as.character(curr.data$taxon)))-1)), legend=c("R1.16O.GL", "R2.16O.GL", "R3.16O.GL", "R1.12C.GL", "R2.12C.GL", "R3.12C.GL", "R1.12C_N.GL", "R2.12C_N.GL", "R3.12C_N.GL", "R1.18O.GL", "R2.18O.GL", "R3.18O.GL"), col=ink, pt.bg=ink, pch=c(rep(21,9), 24, 24, 24), pt.cex=0.7, lwd=1, bty="o", cex=0.6, xjust=0, yjust=1)
    par(xpd=FALSE)
    rm(obs.wads.by.taxon, reps.by.trt, treatments, curr.data, ink)
#####}

#Unsorted line graph of WADs -- taxon numbers are the same across tubes -- (only 16O & 18O replicates: 2 trts * 3 reps/trt = 6 tubes):
#(note: because taxa are unsorted, excluded taxa (i.e., ones that did not occur in all replicates) are skipped when drawing lines between the taxa that are present in all replicates)
#(note: the three 16O reps are all black in this plot)
#####{
  graphics.off()
  dev.new(width=6, height=7.74)
  par(mai=c(0.44,0.57,0.22,1.60), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
    #Select only the data for taxa common to all replicates of the specified treatments:
    obs.wads.by.taxon <- WAD.by.taxon.GL$obs.wads
    reps.by.trt <- WAD.by.taxon.GL$reps.by.trt
    treatments=c("16O.GL", "12C.GL", "12C_N.GL", "18O.GL")
    reps.by.trt <- reps.by.trt[reps.by.trt[,1] %in% treatments,]
      reps.by.trt[,1] <- factor(reps.by.trt[,1])
      reps.by.trt[,2] <- factor(reps.by.trt[,2])
      reps.by.trt[,3] <- factor(reps.by.trt[,3])
      reps.by.trt[,4] <- factor(reps.by.trt[,4])
    curr.data <- cbind(taxon=obs.wads.by.taxon[,1], obs.wads.by.taxon[,names(obs.wads.by.taxon) %in% as.character(unlist(reps.by.trt[,2:4]))])
      #Remove any taxon that does not occur in all 12 replicates:
      curr.data <- curr.data[apply(is.na(curr.data[2:dim(curr.data)[2]]), 1, sum) == 0, ]
      curr.data$taxon <- factor(curr.data$taxon)
    #Set the x-limits for the graph:
    x.min <- min(curr.data[,2:dim(curr.data)[2]], na.rm=TRUE)
    x.max <- max(curr.data[,2:dim(curr.data)[2]], na.rm=TRUE)
    AT.y <- seq(0, max(as.numeric(as.character(obs.wads.by.taxon$taxon))), 50)
  plot(x=curr.data[,2], y=as.numeric(as.character(curr.data$taxon)), type="n", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min, x.max))
  par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
  axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
  mtext(expression(paste("Observed WAD (g ml"^-1, ")", sep="")), side=1, line=1.35, cex=0.75)
  par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
  axis(side=2, at=AT.y, labels=AT.y, tck=-0.015, las=1, cex.axis=0.6)
  mtext("Taxon", side=2, line=2.35, cex=0.75)
  #16O.GL
    # points(x=curr.data[,2], y=as.numeric(as.character(curr.data$taxon)), pch=21, cex=0.7, col="black", bg="black")         #R1
    # points(x=curr.data[,9], y=as.numeric(as.character(curr.data$taxon)), pch=21, cex=0.7, col="black", bg="black")         #R2
    # points(x=curr.data[,10], y=as.numeric(as.character(curr.data$taxon)), pch=21, cex=0.7, col="black", bg="black")        #R3
    lines(x=curr.data[,2], y=as.numeric(as.character(curr.data$taxon)), col="black")
    lines(x=curr.data[,9], y=as.numeric(as.character(curr.data$taxon)), col="black")
    lines(x=curr.data[,10], y=as.numeric(as.character(curr.data$taxon)), col="black")
  #18O.GL
    # points(x=curr.data[,6], y=as.numeric(as.character(curr.data$taxon)), pch=24, cex=0.7, col="black", bg="gray")            #R1
    # points(x=curr.data[,7], y=as.numeric(as.character(curr.data$taxon)), pch=24, cex=0.7, col="black", bg="brown")           #R2
    # points(x=curr.data[,8], y=as.numeric(as.character(curr.data$taxon)), pch=24, cex=0.7, col="black", bg="darkgoldenrod2")  #R3
    lines(x=curr.data[,6], y=as.numeric(as.character(curr.data$taxon)), col="gray")
    lines(x=curr.data[,7], y=as.numeric(as.character(curr.data$taxon)), col="brown")
    lines(x=curr.data[,8], y=as.numeric(as.character(curr.data$taxon)), col="darkgoldenrod2")
    par(xpd=NA)
    ink <- c("black", "black", "black", "gray", "brown", "darkgoldenrod2")
    legend(x=x.max*1.0055, y=max(as.numeric(as.character(curr.data$taxon)))+(0.04*(max(as.numeric(as.character(curr.data$taxon)))-1)), legend=c("R1.16O.GL", "R2.16O.GL", "R3.16O.GL", "R1.18O.GL", "R2.18O.GL", "R3.18O.GL"), col=ink, pt.bg=ink, pch=c(rep(21,3), 24, 24, 24), pt.cex=0.7, lwd=1, bty="o", cex=0.6, xjust=0, yjust=1)
    par(xpd=FALSE)
    rm(obs.wads.by.taxon, reps.by.trt, treatments, curr.data, ink)
#####}


graphics.off()


#Calculate the shift in WAD for the specified unlabeled replicates from their global mean:

  #Identifying the appropriate shift in unlabeled WADs invloves first filtering the data to include only those taxa occuring in all replicates of all specified treatments, 
  #then fitting a normal distribution to the mean WADs of the unlabeled replicates ranked by mean WAD, 
  #then fitting normal distributions to the individual replicates, and then calculating the shift
  #Note: this function takes a few minutes to run (i.e., ~4min for ~250 taxa common to 3 unlabeled and 1 labeled treatment):
    #Include all four unlabeled treatments and the "18O" treatment from the Grassland ecosystem:
      system.time(unlab.WAD.corr.18O.GL.list <- find.unlabeled.correction(LIST=WAD.by.taxon.GL, unlab.tmts=c("16O.GL", "12C.GL", "12C_N.GL"), lab.tmts=c("18O.GL"), CI=0.90))
      #Look at the results:
        names(unlab.WAD.corr.18O.GL.list)                           #names of the two data frames in the output list
        unlab.WAD.corr.18O.GL.list$WAD.norm.fit.parms               #looking at the first object -- a data frame
        unlab.WAD.corr.18O.GL.list[[1]]                             #another way to look at the first object
        unlab.WAD.corr.18O.GL.list$corr.names                       #looking at the second object in the list -- a vector of the names of the corrected unlabeled replicates
        unlab.WAD.corr.18O.GL.list[[2]]                             #another way to look at the second object
        head(unlab.WAD.corr.18O.GL.list$WAD.table.corr)             #looking at the head of the third object -- a data frame
        head(unlab.WAD.corr.18O.GL.list[[3]])                       #another way to look at the head of the third object


#Resorted graph of WADs -- taxon numbers are NOT the same across tubes -- (all replicates: 4 trts * 3 reps/trt = 12 tubes):
#####{
  graphics.off()
  # dev.new(width=6, height=7.74)
  pdf(file="qSIP_output/Figures/RM_DIM_GL_WAD_by_taxon_resorted_all_reps_raw.pdf", width=6, height=7.74)
  par(mai=c(0.44,0.57,0.22,1.60), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
    #Select the data to plot:
    obs.wads.by.taxon <- WAD.by.taxon.GL$obs.wads
    curr.data.orders <- unlab.WAD.corr.18O.GL.list$WAD.table.corr
    #Set the x-limits for the graph:
    x.min <- min(curr.data.orders[,c(19:21,27:35)], na.rm=TRUE)
    x.max <- max(curr.data.orders[,c(19:21,27:35)], na.rm=TRUE)
    AT.y <- seq(0, max(as.numeric(as.character(obs.wads.by.taxon$taxon))), 50)
    plot(x=curr.data.orders$R1.16O.GL, y=curr.data.orders$R1.16O.GL.order, type="n", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min, x.max))
  par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
  axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
  mtext(expression(paste("Observed WAD (g ml"^-1, ")", sep="")), side=1, line=1.35, cex=0.75)
  par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
  axis(side=2, at=AT.y, labels=AT.y, tck=-0.015, las=1, cex.axis=0.6)
  mtext("Taxon (each replicate ranked independently)", side=2, line=2.35, cex=0.75)
  #16O.GL
    points(x=curr.data.orders$R1.16O.GL, y=curr.data.orders$R1.16O.GL.order, pch=21, cex=0.7, col="black", bg="green")
    points(x=curr.data.orders$R2.16O.GL, y=curr.data.orders$R2.16O.GL.order, pch=21, cex=0.7, col="black", bg="red")
    points(x=curr.data.orders$R3.16O.GL, y=curr.data.orders$R3.16O.GL.order, pch=21, cex=0.7, col="black", bg="yellow")
  #12C.GL
    points(x=curr.data.orders$R1.12C.GL, y=curr.data.orders$R1.12C.GL.order, pch=21, cex=0.7, col="black", bg="purple")
    points(x=curr.data.orders$R2.12C.GL, y=curr.data.orders$R2.12C.GL.order, pch=21, cex=0.7, col="black", bg="blue")
    points(x=curr.data.orders$R3.12C.GL, y=curr.data.orders$R3.12C.GL.order, pch=21, cex=0.7, col="black", bg="orange")
  #12C_N.GL
    points(x=curr.data.orders$R1.12C_N.GL, y=curr.data.orders$R1.12C_N.GL.order, pch=21, cex=0.7, col="black", bg="cyan")
    points(x=curr.data.orders$R2.12C_N.GL, y=curr.data.orders$R2.12C_N.GL.order, pch=21, cex=0.7, col="black", bg="magenta")
    points(x=curr.data.orders$R3.12C_N.GL, y=curr.data.orders$R3.12C_N.GL.order, pch=21, cex=0.7, col="black", bg="wheat")
  #18.O.GL
    points(x=curr.data.orders$R1.18O.GL, y=curr.data.orders$R1.18O.GL.order, pch=21, cex=0.7, col="black", bg="gray")
    points(x=curr.data.orders$R2.18O.GL, y=curr.data.orders$R2.18O.GL.order, pch=21, cex=0.7, col="black", bg="brown")
    points(x=curr.data.orders$R3.18O.GL, y=curr.data.orders$R3.18O.GL.order, pch=21, cex=0.7, col="black", bg="darkgoldenrod2")
    par(xpd=NA)
    ink <- c("green", "red", "yellow", "purple", "blue", "orange", "cyan", "magenta", "wheat", "gray", "brown", "darkgoldenrod2")
    legend(x=x.max*1.0055, y=max(as.numeric(as.character(curr.data.orders$mean.unlab.WAD.order)))+(0.04*(max(as.numeric(as.character(curr.data.orders$mean.unlab.WAD.order)))-1)), legend=c("R1.16O.GL", "R2.16O.GL", "R3.16O.GL", "R1.12C.GL", "R2.12C.GL", "R3.12C.GL", "R1.12C_N.GL", "R2.12C_N.GL", "R3.12C_N.GL", "R1.18O.GL", "R2.18O.GL", "R3.18O.GL"), col=ink, pt.bg=ink, pch=rep(21,12), pt.cex=0.7, lwd=1, bty="o", cex=0.6, xjust=0, yjust=1)
    par(xpd=FALSE)
    rm(obs.wads.by.taxon, curr.data.orders, ink)
#####}
dev.off()   #(run graph code through this line only if creating a pdf version of the plot)

#Sorted graph of WADs -- taxon numbers are the same across tubes -- (only 16O replicates: 3 trts * 3 reps/trt = 9 tubes):
#White line = mean of unlabeled data
#Green line = normal cdf fit to the mean of all unlabeled data
#####{
  graphics.off()
  # dev.new(width=6, height=7.74)
  pdf(file="qSIP_output/Figures/RM_DIM_GL_WAD_by_taxon_unlabeled_reps_raw.pdf", width=6, height=7.74)
  par(mai=c(0.44,0.57,0.22,1.60), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
    #Select the data to plot:
    obs.wads.by.taxon <- WAD.by.taxon.GL$obs.wads
    unlab.WAD.norm.fit.parms <- unlab.WAD.corr.18O.GL.list$WAD.norm.fit.parms
    curr.data.orders <- unlab.WAD.corr.18O.GL.list$WAD.table.corr
    unlab.reps <- c("R1.16O.GL", "R2.16O.GL", "R3.16O.GL", "R1.12C.GL", "R2.12C.GL", "R3.12C.GL", "R1.12C_N.GL", "R2.12C_N.GL", "R3.12C_N.GL")
    unlab.WADs <- curr.data.orders[, names(curr.data.orders) %in% unlab.reps]
    #Set the x-limits for the graph:
    x.min <- min(curr.data.orders[,c(19:21,27:35)], na.rm=TRUE)
    x.max <- max(curr.data.orders[,c(19:21,27:35)], na.rm=TRUE)
    AT.y <- seq(0, max(as.numeric(as.character(obs.wads.by.taxon$taxon))), 50)
  plot(x=curr.data.orders$R1.16O.GL, y=curr.data.orders$R1.16O.GL.order, type="n", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min, x.max))
  par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
  axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
  mtext(expression(paste("Observed WAD (g ml"^-1, ")", sep="")), side=1, line=1.35, cex=0.75)
  par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
  axis(side=2, at=AT.y, labels=AT.y, tck=-0.015, las=1, cex.axis=0.6)
  mtext("Taxon", side=2, line=2.4, cex=0.75)
  #16O.GL
    points(x=curr.data.orders$R1.16O.GL, y=curr.data.orders$mean.unlab.WAD.order, pch=21, cex=0.7, col="black", bg="green")
    points(x=curr.data.orders$R2.16O.GL, y=curr.data.orders$mean.unlab.WAD.order, pch=21, cex=0.7, col="black", bg="red")
    points(x=curr.data.orders$R3.16O.GL, y=curr.data.orders$mean.unlab.WAD.order, pch=21, cex=0.7, col="black", bg="yellow")
    # new.order <- order(curr.data.orders$mean.unlab.WAD.order)
    # lines(x=curr.data.orders$R1.16O.GL[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=2, col="green")
    # lines(x=curr.data.orders$R2.16O.GL[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=2, col="red")
    # lines(x=curr.data.orders$R3.16O.GL[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=2, col="black")
    # lines(x=curr.data.orders$R3.16O.GL[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=1.5, col="yellow")
  #12C.GL
    points(x=curr.data.orders$R1.12C.GL, y=curr.data.orders$mean.unlab.WAD.order, pch=21, cex=0.7, col="black", bg="purple")
    points(x=curr.data.orders$R2.12C.GL, y=curr.data.orders$mean.unlab.WAD.order, pch=21, cex=0.7, col="black", bg="blue")
    points(x=curr.data.orders$R3.12C.GL, y=curr.data.orders$mean.unlab.WAD.order, pch=21, cex=0.7, col="black", bg="orange")
    # new.order <- order(curr.data.orders$mean.unlab.WAD.order)
    # lines(x=curr.data.orders$R1.12C.GL[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=2, col="purple")
    # lines(x=curr.data.orders$R2.12C.GL[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=2, col="blue")
    # lines(x=curr.data.orders$R3.12C.GL[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=2, col="orange")
  #12C_N.GL
    points(x=curr.data.orders$R1.12C_N.GL, y=curr.data.orders$mean.unlab.WAD.order, pch=21, cex=0.7, col="black", bg="cyan")
    points(x=curr.data.orders$R2.12C_N.GL, y=curr.data.orders$mean.unlab.WAD.order, pch=21, cex=0.7, col="black", bg="magenta")
    points(x=curr.data.orders$R3.12C_N.GL, y=curr.data.orders$mean.unlab.WAD.order, pch=21, cex=0.7, col="black", bg="wheat")
    # new.order <- order(curr.data.orders$mean.unlab.WAD.order)
    # lines(x=curr.data.orders$R1.12C_N.GL[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=2, col="cyan")
    # lines(x=curr.data.orders$R2.12C_N.GL[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=2, col="magenta")
    # lines(x=curr.data.orders$R3.12C_N.GL[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=2, col="wheat")
  #Unlabeled average:
    # points(x=apply(unlab.WADs, 1, mean, na.rm=TRUE), y=curr.data.orders$mean.unlab.WAD.order, pch=21, cex=0.35, col="black", bg="black")
    new.order <- order(curr.data.orders$mean.unlab.WAD.order)
    lines(x=apply(unlab.WADs, 1, mean, na.rm=TRUE)[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=3, col="black")
    lines(x=apply(unlab.WADs, 1, mean, na.rm=TRUE)[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=2, col="white")
  #Normal cdf fit to the mean of all unlabeled data:
    x <- seq(min(unlab.WADs), max(unlab.WADs), 0.000001)
    y <- pnorm(q=x, mean=unlab.WAD.norm.fit.parms$mean[unlab.WAD.norm.fit.parms$trt == "unlabeled"], sd=unlab.WAD.norm.fit.parms$stdev[unlab.WAD.norm.fit.parms$trt == "unlabeled"])*max(curr.data.orders$mean.unlab.WAD.order)
    lines(x, y, lwd=3, col="black")
    lines(x, y, lwd=2, col="lawngreen")
    par(xpd=NA)
    ink <- c("green", "red", "yellow", "purple", "blue", "orange", "cyan", "magenta", "wheat", "white", "lawngreen")
    legend(x=x.max*1.0040, y=max(as.numeric(as.character(curr.data.orders$mean.unlab.WAD.order)))+(0.04*(max(as.numeric(as.character(curr.data.orders$mean.unlab.WAD.order)))-1)), legend=c("R1.16O.GL", "R2.16O.GL", "R3.16O.GL", "R1.12C.GL", "R2.12C.GL", "R3.12C.GL", "R1.12C_N.GL", "R2.12C_N.GL", "R3.12C_N.GL", "unlab empirical mean", "unlabeled normal fit"), col=c(ink[1:9], "black", "black"), pt.bg=ink, pch=rep(21,12), pt.cex=0.7, lwd=c(rep(1, 9), 3, 3), bty="o", cex=0.6, xjust=0, yjust=1)
    legend(x=x.max*1.0040, y=max(as.numeric(as.character(curr.data.orders$mean.unlab.WAD.order)))+(0.04*(max(as.numeric(as.character(curr.data.orders$mean.unlab.WAD.order)))-1)), legend=c("R1.16O.GL", "R2.16O.GL", "R3.16O.GL", "R1.12C.GL", "R2.12C.GL", "R3.12C.GL", "R1.12C_N.GL", "R2.12C_N.GL", "R3.12C_N.GL", "unlab empirical mean", "unlabeled normal fit"), col=ink, pt.bg=ink, pch=rep(21,12), pt.cex=0.7, lwd=c(rep(1, 9), 2, 2), bty="o", cex=0.6, xjust=0, yjust=1)
    par(xpd=FALSE)
    rm(obs.wads.by.taxon, unlab.WAD.norm.fit.parms, curr.data.orders, unlab.reps, unlab.WADs, ink)
#####}
dev.off()   #(run graph code through this line only if creating a pdf version of the plot)

#Sorted graph of WADs -- taxon numbers are the same across tubes -- (only 18O replicates: 1 trt * 3 reps/trt = 3 tubes):
#White line = mean of unlabeled data
#Green line = normal cdf fit to the mean of all unlabeled data
#####{
  graphics.off()
  # dev.new(width=6, height=7.74)
  pdf(file="qSIP_output/Figures/RM_DIM_GL_WAD_by_taxon_labeled_reps_raw.pdf", width=6, height=7.74)
  par(mai=c(0.44,0.57,0.22,1.60), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
    #Select the data to plot:
    obs.wads.by.taxon <- WAD.by.taxon.GL$obs.wads
    unlab.WAD.norm.fit.parms <- unlab.WAD.corr.18O.GL.list$WAD.norm.fit.parms
    curr.data.orders <- unlab.WAD.corr.18O.GL.list$WAD.table.corr
    unlab.reps <- c("R1.16O.GL", "R2.16O.GL", "R3.16O.GL", "R1.12C.GL", "R2.12C.GL", "R3.12C.GL", "R1.12C_N.GL", "R2.12C_N.GL", "R3.12C_N.GL")
    unlab.WADs <- curr.data.orders[, names(curr.data.orders) %in% unlab.reps]
    #Set the x-limits for the graph:
    x.min <- min(curr.data.orders[,c(19:21,27:35)], na.rm=TRUE)
    x.max <- max(curr.data.orders[,c(19:21,27:35)], na.rm=TRUE)
    AT.y <- seq(0, max(as.numeric(as.character(obs.wads.by.taxon$taxon))), 50)
  plot(x=curr.data.orders$R1.16O.GL, y=curr.data.orders$R1.16O.GL.order, type="n", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min, x.max))
  par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
  axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
  mtext(expression(paste("Observed WAD (g ml"^-1, ")", sep="")), side=1, line=1.35, cex=0.75)
  par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
  axis(side=2, at=AT.y, labels=AT.y, tck=-0.015, las=1, cex.axis=0.6)
  mtext("Taxon", side=2, line=2.4, cex=0.75)
  #18.O.GL
    points(x=curr.data.orders$R1.18O.GL, y=curr.data.orders$mean.unlab.WAD.order, pch=21, cex=0.7, col="black", bg="gray")
    points(x=curr.data.orders$R2.18O.GL, y=curr.data.orders$mean.unlab.WAD.order, pch=21, cex=0.7, col="black", bg="brown")
    points(x=curr.data.orders$R3.18O.GL, y=curr.data.orders$mean.unlab.WAD.order, pch=21, cex=0.7, col="black", bg="darkgoldenrod2")
    # new.order <- order(curr.data.orders$mean.unlab.WAD.order)
    # lines(x=curr.data.orders$R1.18O.GL[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=2, col="gray")
    # lines(x=curr.data.orders$R2.18O.GL[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=2, col="brown")
    # lines(x=curr.data.orders$R3.18O.GL[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=2, col="darkgoldenrod2")
  #Unlabeled average:
    # points(x=apply(unlab.WADs, 1, mean, na.rm=TRUE), y=curr.data.orders$mean.unlab.WAD.order, pch=21, cex=0.35, col="black", bg="black")
    new.order <- order(curr.data.orders$mean.unlab.WAD.order)
    lines(x=apply(unlab.WADs, 1, mean, na.rm=TRUE)[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=3, col="black")
    lines(x=apply(unlab.WADs, 1, mean, na.rm=TRUE)[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=2, col="white")
  #Normal cdf fit to the mean of all unlabeled data:
    x <- seq(min(unlab.WADs), max(unlab.WADs), 0.000001)
    y <- pnorm(q=x, mean=unlab.WAD.norm.fit.parms$mean[unlab.WAD.norm.fit.parms$trt == "unlabeled"], sd=unlab.WAD.norm.fit.parms$stdev[unlab.WAD.norm.fit.parms$trt == "unlabeled"])*max(curr.data.orders$mean.unlab.WAD.order)
    lines(x, y, lwd=3, col="black")
    lines(x, y, lwd=2, col="lawngreen")
    par(xpd=NA)
    ink <- c("gray", "brown", "darkgoldenrod2", "white", "lawngreen")
    legend(x=x.max*1.0040, y=max(as.numeric(as.character(curr.data.orders$mean.unlab.WAD.order)))+(0.04*(max(as.numeric(as.character(curr.data.orders$mean.unlab.WAD.order)))-1)), legend=c("R1.18O.GL", "R2.18O.GL", "R3.18O.GL", "unlab empirical mean", "unlabeled normal fit"), col=c(ink[1:3], "black", "black"), pt.bg=ink, pch=rep(21,5), pt.cex=0.7, lwd=c(rep(1, 3), 3, 3), bty="o", cex=0.6, xjust=0, yjust=1)
    legend(x=x.max*1.0040, y=max(as.numeric(as.character(curr.data.orders$mean.unlab.WAD.order)))+(0.04*(max(as.numeric(as.character(curr.data.orders$mean.unlab.WAD.order)))-1)), legend=c("R1.18O.GL", "R2.18O.GL", "R3.18O.GL", "unlab empirical mean", "unlabeled normal fit"), col=ink, pt.bg=ink, pch=rep(21,5), pt.cex=0.7, lwd=c(rep(1, 3), 2, 2), bty="o", cex=0.6, xjust=0, yjust=1)
    par(xpd=FALSE)
    rm(obs.wads.by.taxon, unlab.WAD.norm.fit.parms, curr.data.orders, unlab.reps, unlab.WADs, ink)
#####}
dev.off()   #(run graph code through this line only if creating a pdf version of the plot)

#Sorted graph of WADs -- taxon numbers are the same across tubes -- (16O & 18O replicates: 4 trts * 3 reps/trt = 12 tubes):
#White line = mean of unlabeled data
#Green line = normal cdf fit to the mean of all unlabeled data
#####{
  graphics.off()
  # dev.new(width=6, height=7.74)
  pdf(file="qSIP_output/Figures/RM_DIM_GL_WAD_by_taxon_all_reps_raw.pdf", width=6, height=7.74)
  par(mai=c(0.44,0.57,0.22,1.60), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
    #Select the data to plot:
    obs.wads.by.taxon <- WAD.by.taxon.GL$obs.wads
    unlab.WAD.norm.fit.parms <- unlab.WAD.corr.18O.GL.list$WAD.norm.fit.parms
    curr.data.orders <- unlab.WAD.corr.18O.GL.list$WAD.table.corr
    unlab.reps <- c("R1.16O.GL", "R2.16O.GL", "R3.16O.GL", "R1.12C.GL", "R2.12C.GL", "R3.12C.GL", "R1.12C_N.GL", "R2.12C_N.GL", "R3.12C_N.GL")
    unlab.WADs <- curr.data.orders[, names(curr.data.orders) %in% unlab.reps]
    #Set the x-limits for the graph:
    x.min <- min(curr.data.orders[,c(19:21,27:35)], na.rm=TRUE)
    x.max <- max(curr.data.orders[,c(19:21,27:35)], na.rm=TRUE)
    AT.y <- seq(0, max(as.numeric(as.character(obs.wads.by.taxon$taxon))), 50)
  plot(x=curr.data.orders$R1.16O.GL, y=curr.data.orders$R1.16O.GL.order, type="n", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min, x.max))
  par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
  axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
  mtext(expression(paste("Observed WAD (g ml"^-1, ")", sep="")), side=1, line=1.35, cex=0.75)
  par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
  axis(side=2, at=AT.y, labels=AT.y, tck=-0.015, las=1, cex.axis=0.6)
  mtext("Taxon", side=2, line=2.4, cex=0.75)
  #16O.GL
    points(x=curr.data.orders$R1.16O.GL, y=curr.data.orders$mean.unlab.WAD.order, pch=21, cex=0.7, col="black", bg="green")
    points(x=curr.data.orders$R2.16O.GL, y=curr.data.orders$mean.unlab.WAD.order, pch=21, cex=0.7, col="black", bg="red")
    points(x=curr.data.orders$R3.16O.GL, y=curr.data.orders$mean.unlab.WAD.order, pch=21, cex=0.7, col="black", bg="yellow")
    # new.order <- order(curr.data.orders$mean.unlab.WAD.order)
    # lines(x=curr.data.orders$R1.16O.GL[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=2, col="green")
    # lines(x=curr.data.orders$R2.16O.GL[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=2, col="red")
    # lines(x=curr.data.orders$R3.16O.GL[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=2, col="black")
    # lines(x=curr.data.orders$R3.16O.GL[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=1.5, col="yellow")
  #12C.GL
    points(x=curr.data.orders$R1.12C.GL, y=curr.data.orders$mean.unlab.WAD.order, pch=21, cex=0.7, col="black", bg="purple")
    points(x=curr.data.orders$R2.12C.GL, y=curr.data.orders$mean.unlab.WAD.order, pch=21, cex=0.7, col="black", bg="blue")
    points(x=curr.data.orders$R3.12C.GL, y=curr.data.orders$mean.unlab.WAD.order, pch=21, cex=0.7, col="black", bg="orange")
    # new.order <- order(curr.data.orders$mean.unlab.WAD.order)
    # lines(x=curr.data.orders$R1.12C.GL[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=2, col="purple")
    # lines(x=curr.data.orders$R2.12C.GL[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=2, col="blue")
    # lines(x=curr.data.orders$R3.12C.GL[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=2, col="orange")
  #12C_N.GL
    points(x=curr.data.orders$R1.12C_N.GL, y=curr.data.orders$mean.unlab.WAD.order, pch=21, cex=0.7, col="black", bg="cyan")
    points(x=curr.data.orders$R2.12C_N.GL, y=curr.data.orders$mean.unlab.WAD.order, pch=21, cex=0.7, col="black", bg="magenta")
    points(x=curr.data.orders$R3.12C_N.GL, y=curr.data.orders$mean.unlab.WAD.order, pch=21, cex=0.7, col="black", bg="wheat")
    # new.order <- order(curr.data.orders$mean.unlab.WAD.order)
    # lines(x=curr.data.orders$R1.12C_N.GL[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=2, col="cyan")
    # lines(x=curr.data.orders$R2.12C_N.GL[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=2, col="magenta")
    # lines(x=curr.data.orders$R3.12C_N.GL[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=2, col="wheat")
  #18.O.GL
    points(x=curr.data.orders$R1.18O.GL, y=curr.data.orders$mean.unlab.WAD.order, pch=21, cex=0.7, col="black", bg="gray")
    points(x=curr.data.orders$R2.18O.GL, y=curr.data.orders$mean.unlab.WAD.order, pch=21, cex=0.7, col="black", bg="brown")
    points(x=curr.data.orders$R3.18O.GL, y=curr.data.orders$mean.unlab.WAD.order, pch=21, cex=0.7, col="black", bg="darkgoldenrod2")
    # new.order <- order(curr.data.orders$mean.unlab.WAD.order)
    # lines(x=curr.data.orders$R1.18O.GL[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=2, col="gray")
    # lines(x=curr.data.orders$R2.18O.GL[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=2, col="brown")
    # lines(x=curr.data.orders$R3.18O.GL[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=2, col="darkgoldenrod2")
  #Unlabeled average:
    # points(x=apply(unlab.WADs, 1, mean, na.rm=TRUE), y=curr.data.orders$mean.unlab.WAD.order, pch=21, cex=0.35, col="black", bg="black")
    new.order <- order(curr.data.orders$mean.unlab.WAD.order)
    lines(x=apply(unlab.WADs, 1, mean, na.rm=TRUE)[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=3, col="black")
    lines(x=apply(unlab.WADs, 1, mean, na.rm=TRUE)[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=2, col="white")
  #Normal cdf fit to the mean of all unlabeled data:
    x <- seq(min(unlab.WADs), max(unlab.WADs), 0.000001)
    y <- pnorm(q=x, mean=unlab.WAD.norm.fit.parms$mean[unlab.WAD.norm.fit.parms$trt == "unlabeled"], sd=unlab.WAD.norm.fit.parms$stdev[unlab.WAD.norm.fit.parms$trt == "unlabeled"])*max(curr.data.orders$mean.unlab.WAD.order)
    lines(x, y, lwd=3, col="black")
    lines(x, y, lwd=2, col="lawngreen")
    par(xpd=NA)
    ink <- c("green", "red", "yellow", "purple", "blue", "orange", "cyan", "magenta", "wheat", "gray", "brown", "darkgoldenrod2", "white", "lawngreen")
    legend(x=x.max*1.0040, y=max(as.numeric(as.character(curr.data.orders$mean.unlab.WAD.order)))+(0.04*(max(as.numeric(as.character(curr.data.orders$mean.unlab.WAD.order)))-1)), legend=c("R1.16O.GL", "R2.16O.GL", "R3.16O.GL", "R1.12C.GL", "R2.12C.GL", "R3.12C.GL", "R1.12C_N.GL", "R2.12C_N.GL", "R3.12C_N.GL", "R1.18O.GL", "R2.18O.GL", "R3.18O.GL", "unlab empirical mean", "unlabeled normal fit"), col=c(ink[1:12], "black", "black"), pt.bg=ink, pch=rep(21,12), pt.cex=0.7, lwd=c(rep(1, 12), 3, 3), bty="o", cex=0.6, xjust=0, yjust=1)
    legend(x=x.max*1.0040, y=max(as.numeric(as.character(curr.data.orders$mean.unlab.WAD.order)))+(0.04*(max(as.numeric(as.character(curr.data.orders$mean.unlab.WAD.order)))-1)), legend=c("R1.16O.GL", "R2.16O.GL", "R3.16O.GL", "R1.12C.GL", "R2.12C.GL", "R3.12C.GL", "R1.12C_N.GL", "R2.12C_N.GL", "R3.12C_N.GL", "R1.18O.GL", "R2.18O.GL", "R3.18O.GL", "unlab empirical mean", "unlabeled normal fit"), col=ink, pt.bg=ink, pch=rep(21,12), pt.cex=0.7, lwd=c(rep(1, 12), 2, 2), bty="o", cex=0.6, xjust=0, yjust=1)
    par(xpd=FALSE)
    rm(obs.wads.by.taxon, unlab.WAD.norm.fit.parms, curr.data.orders, unlab.reps, unlab.WADs, ink)
#####}
dev.off()   #(run graph code through this line only if creating a pdf version of the plot)


# Compare replicate #1 of the 18O Grassland treatment to the mean of all (9) unlabeled replicates:

  #Iterate through the remaining taxa to create a data frame of metrics associated with different subsets of putative nongrowing taxa in the labeled treatment for use in assessing a 'match' with unlabeled taxa.
  #And simultaneously create three more data frames: one that lists the taxa included and one the taxa excluded for each iteration, the function also returns the raw data filtered according to the specifiec SD thresholds.
  #Note: the 'fit.norm.func' function calls in these functions slow them process down substantially (i.e., for ~60 taxa: from <1min without them to ~10min with them)
    #Start with all data after imposing SD thresholds; iteratively remove the data point with the largest positive residual:
      system.time(R1.18O.GL.td.pos.resid.list <- find.labeled.correction(LIST=unlab.WAD.corr.18O.GL.list, lab.replicate="R1.18O.GL", lab.names=c("R1.18O.GL", "R2.18O.GL", "R3.18O.GL"), method="td.pos.resid", unlab.SD.percentile=0.5, lab.SD.percentile=0.5, min.num.nongrowers=10))
    #Start with all data after imposing SD thresholds; iteratively remove the data point with the largest absolute value residual:
      system.time(R1.18O.GL.td.abs.resid.list <- find.labeled.correction(LIST=unlab.WAD.corr.18O.GL.list, lab.replicate="R1.18O.GL", lab.names=c("R1.18O.GL", "R2.18O.GL", "R3.18O.GL"), method="td.abs.resid", unlab.SD.percentile=0.5, lab.SD.percentile=0.5, min.num.nongrowers=10))
    #Start with only the set of taxa having the lowest delta-WAD after imposing SD thresholds; iteratively add taxa with the next-smallest absolute value residual (this is effectively the same as 'td.pos.resid'):
      system.time(R1.18O.GL.bu.abs.resid.list <- find.labeled.correction(LIST=unlab.WAD.corr.18O.GL.list, lab.replicate="R1.18O.GL", lab.names=c("R1.18O.GL", "R2.18O.GL", "R3.18O.GL"), method="bu.abs.resid", unlab.SD.percentile=0.5, lab.SD.percentile=0.5, min.num.nongrowers=10))

  #Write the find.labeled.correction output tables to text files:
    #td.pos.resid results:
      write.table(R1.18O.GL.td.pos.resid.list$putative.nongrower.metrics, "qSIP_output/RM_DIM_R1_18O_GL_td_pos_resid_putative_nongrower_metrics.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R1.18O.GL.td.pos.resid.list$taxa.in, "qSIP_output/RM_DIM_R1_18O_GL_td_pos_resid_taxa_in.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R1.18O.GL.td.pos.resid.list$taxa.out, "qSIP_output/RM_DIM_R1_18O_GL_td_pos_resid_taxa_out.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R1.18O.GL.td.pos.resid.list$data.low.SD, "qSIP_output/RM_DIM_R1_18O_GL_td_pos_resid_data_low_SD.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
    #td.abs.resid results:
      write.table(R1.18O.GL.td.abs.resid.list$putative.nongrower.metrics, "qSIP_output/RM_DIM_R1_18O_GL_td_abs_resid_putative_nongrower_metrics.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R1.18O.GL.td.abs.resid.list$taxa.in, "qSIP_output/RM_DIM_R1_18O_GL_td_abs_resid_taxa_in.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R1.18O.GL.td.abs.resid.list$taxa.out, "qSIP_output/RM_DIM_R1_18O_GL_td_abs_resid_taxa_out.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R1.18O.GL.td.abs.resid.list$data.low.SD, "qSIP_output/RM_DIM_R1_18O_GL_td_abs_resid_data_low_SD.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
    #bu.abs.resid results:
      write.table(R1.18O.GL.bu.abs.resid.list$putative.nongrower.metrics, "qSIP_output/RM_DIM_R1_18O_GL_bu_abs_resid_putative_nongrower_metrics.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R1.18O.GL.bu.abs.resid.list$taxa.in, "qSIP_output/RM_DIM_R1_18O_GL_bu_abs_resid_taxa_in.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R1.18O.GL.bu.abs.resid.list$taxa.out, "qSIP_output/RM_DIM_R1_18O_GL_bu_abs_resid_taxa_out.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R1.18O.GL.bu.abs.resid.list$data.low.SD, "qSIP_output/RM_DIM_R1_18O_GL_bu_abs_resid_data_low_SD.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)

  #Some examples of plotting the results of the search for the 'ideal' set of 'nongrowers' for R1.18O.GL:
  #(NOTE: plots are written and immediately opened as pdfs (uncomment-out the 'system' commands to achieve this); as currently named, these files will be overwritten later)
    #Plotting all iterations (no best iteration highlighted):
      graphics.off()
      find.labeled.correction.plot(find.labeled.correction.list=R1.18O.GL.td.pos.resid.list, filename="RM_DIM_GL_R1_18O_GL_td_pos_resid_plots", path="qSIP_output/Figures/")
      # system('open qSIP_output/Figures/RM_DIM_GL_R1_18O_GL_td_pos_resid_plots.pdf')
    #Finding and plotting the best set of iterations (highlighting those iterations that appeared most frequently among the fit assessment metrics using a quantile threshold of 0.20 to identify the 'best' iterations for each metric):
      R1.18O.GL.td.pos.resid.best.iteration.list <- select.best.iteration(putative.nongrower.metrics=R1.18O.GL.td.pos.resid.list$putative.nongrower.metrics, quantile.threshold=0.20)
      R1.18O.GL.td.pos.resid.best.iteration.list
      graphics.off()
      find.labeled.correction.plot(find.labeled.correction.list=R1.18O.GL.td.pos.resid.list, filename="RM_DIM_GL_R1_18O_GL_td_pos_resid_plots", path="qSIP_output/Figures/", highlight=R1.18O.GL.td.pos.resid.best.iteration.list$freq.iterations, highlight.col=c("gold", "red", "magenta"))
      # system('open qSIP_output/Figures/RM_DIM_GL_R1_18O_GL_td_pos_resid_plots.pdf')
    #Finding and plotting the single best iteration using a quantile threshold of 0.10 to identify the 'best' iterations for each metric:
      R1.18O.GL.td.pos.resid.best.iteration.list <- select.best.iteration(putative.nongrower.metrics=R1.18O.GL.td.pos.resid.list$putative.nongrower.metrics, quantile.threshold=0.10)
      R1.18O.GL.td.pos.resid.best.iteration.list
      graphics.off()
      find.labeled.correction.plot(find.labeled.correction.list=R1.18O.GL.td.pos.resid.list, filename="RM_DIM_GL_R1_18O_GL_td_pos_resid_plots", path="qSIP_output/Figures/", highlight=R1.18O.GL.td.pos.resid.best.iteration.list$best.iteration, highlight.col="magenta")
      # system('open qSIP_output/Figures/RM_DIM_GL_R1_18O_GL_td_pos_resid_plots.pdf')

  graphics.off()

  #Find the best iteration (using all three methods):
    R1.18O.GL.td.pos.resid.best.iteration.list <- select.best.iteration(putative.nongrower.metrics=R1.18O.GL.td.pos.resid.list$putative.nongrower.metrics, quantile.threshold=0.10)
    R1.18O.GL.td.pos.resid.best.iteration.list
    find.labeled.correction.plot(find.labeled.correction.list=R1.18O.GL.td.pos.resid.list, filename="RM_DIM_GL_R1_18O_GL_td_pos_resid_plots", path="qSIP_output/Figures/", highlight=R1.18O.GL.td.pos.resid.best.iteration.list$best.iteration, highlight.col="magenta")

    R1.18O.GL.td.abs.resid.best.iteration.list <- select.best.iteration(putative.nongrower.metrics=R1.18O.GL.td.abs.resid.list$putative.nongrower.metrics, quantile.threshold=0.10)
    R1.18O.GL.td.abs.resid.best.iteration.list
    find.labeled.correction.plot(find.labeled.correction.list=R1.18O.GL.td.abs.resid.list, filename="RM_DIM_GL_R1_18O_GL_td_abs_resid_plots", path="qSIP_output/Figures/", highlight=R1.18O.GL.td.abs.resid.best.iteration.list$best.iteration, highlight.col="magenta")

    R1.18O.GL.bu.abs.resid.best.iteration.list <- select.best.iteration(putative.nongrower.metrics=R1.18O.GL.bu.abs.resid.list$putative.nongrower.metrics, quantile.threshold=0.10)
    R1.18O.GL.bu.abs.resid.best.iteration.list
    find.labeled.correction.plot(find.labeled.correction.list=R1.18O.GL.bu.abs.resid.list, filename="RM_DIM_GL_R1_18O_GL_bu_abs_resid_plots", path="qSIP_output/Figures/", highlight=R1.18O.GL.bu.abs.resid.best.iteration.list$best.iteration, highlight.col="magenta")

  #Get the sequence of taxa added at each iteration:
    get.seq.taxa.nums(R1.18O.GL.td.pos.resid.list$taxa.in)
    get.seq.taxa.nums(R1.18O.GL.td.abs.resid.list$taxa.in)
    get.seq.taxa.nums(R1.18O.GL.bu.abs.resid.list$taxa.in)

  #Sanity check: sequence of taxa added or removed is the same (as it should be) within each algorithm:
    sum(get.seq.taxa.nums(R1.18O.GL.td.pos.resid.list$taxa.in) != get.seq.taxa.nums(R1.18O.GL.td.pos.resid.list$taxa.out))
    sum(get.seq.taxa.nums(R1.18O.GL.td.abs.resid.list$taxa.in) != get.seq.taxa.nums(R1.18O.GL.td.abs.resid.list$taxa.out))

  #Compare the sequence of taxa added or removed between the 'td.pos.resid' and 'bu.abs.resid' algorithms:
    sum(get.seq.taxa.nums(R1.18O.GL.td.pos.resid.list$taxa.in)[2:length(get.seq.taxa.nums(R1.18O.GL.td.pos.resid.list$taxa.in))] != get.seq.taxa.nums(R1.18O.GL.bu.abs.resid.list$taxa.in)[length(get.seq.taxa.nums(R1.18O.GL.bu.abs.resid.list$taxa.in)):2])
    plot(x=get.seq.taxa.nums(R1.18O.GL.td.pos.resid.list$taxa.in)[2:length(get.seq.taxa.nums(R1.18O.GL.td.pos.resid.list$taxa.in))], y=get.seq.taxa.nums(R1.18O.GL.bu.abs.resid.list$taxa.in)[length(get.seq.taxa.nums(R1.18O.GL.bu.abs.resid.list$taxa.in)):2], xlab="Top-down sequence", ylab="Bottom-up sequence")
    abline(a=0, b=1, col="red")
    #They are the same

  #Compare the sequence of taxa removed between the 'td.pos.resid' and 'td.abs.resid' algorithms:
    sum(get.seq.taxa.nums(R1.18O.GL.td.pos.resid.list$taxa.in)[1:length(get.seq.taxa.nums(R1.18O.GL.td.pos.resid.list$taxa.in))] != get.seq.taxa.nums(R1.18O.GL.td.abs.resid.list$taxa.in)[1:length(get.seq.taxa.nums(R1.18O.GL.td.abs.resid.list$taxa.in))])
    plot(x=get.seq.taxa.nums(R1.18O.GL.td.pos.resid.list$taxa.in)[1:length(get.seq.taxa.nums(R1.18O.GL.td.pos.resid.list$taxa.in))], y=get.seq.taxa.nums(R1.18O.GL.td.abs.resid.list$taxa.in)[1:length(get.seq.taxa.nums(R1.18O.GL.td.abs.resid.list$taxa.in))], xlab="Top-down sequence", ylab="Bottom-up sequence")
    abline(a=0, b=1, col="red")
    #They are not the same


graphics.off()


# Compare replicate #2 of the 18O Grassland treatment to the mean of all (9) unlabeled replicates:

  #Iterate through the remaining taxa to create a data frame of metrics associated with different subsets of putative nongrowing taxa in the labeled treatment for use in assessing a 'match' with unlabeled taxa.
  #And simultaneously create three more data frames: one that lists the taxa included and one the taxa excluded for each iteration, the function also returns the raw data filtered according to the specifiec SD thresholds.
  #Note: the 'fit.norm.func' function calls in these functions slow them process down substantially (i.e., for ~60 taxa: from <1min without them to ~10min with them)
    #Start with all data after imposing SD thresholds; iteratively remove the data point with the largest positive residual:
      system.time(R2.18O.GL.td.pos.resid.list <- find.labeled.correction(LIST=unlab.WAD.corr.18O.GL.list, lab.replicate="R2.18O.GL", lab.names=c("R1.18O.GL", "R2.18O.GL", "R3.18O.GL"), method="td.pos.resid", unlab.SD.percentile=0.5, lab.SD.percentile=0.5, min.num.nongrowers=10))
    #Start with all data after imposing SD thresholds; iteratively remove the data point with the largest absolute value residual:
      system.time(R2.18O.GL.td.abs.resid.list <- find.labeled.correction(LIST=unlab.WAD.corr.18O.GL.list, lab.replicate="R2.18O.GL", lab.names=c("R1.18O.GL", "R2.18O.GL", "R3.18O.GL"), method="td.abs.resid", unlab.SD.percentile=0.5, lab.SD.percentile=0.5, min.num.nongrowers=10))
    #Start with only the set of taxa having the lowest delta-WAD after imposing SD thresholds; iteratively add taxa with the next-smallest absolute value residual (this is effectively the same as 'td.pos.resid'):
      system.time(R2.18O.GL.bu.abs.resid.list <- find.labeled.correction(LIST=unlab.WAD.corr.18O.GL.list, lab.replicate="R2.18O.GL", lab.names=c("R1.18O.GL", "R2.18O.GL", "R3.18O.GL"), method="bu.abs.resid", unlab.SD.percentile=0.5, lab.SD.percentile=0.5, min.num.nongrowers=10))

  #Write the find.labeled.correction output tables to text files:
    #td.pos.resid results:
      write.table(R2.18O.GL.td.pos.resid.list$putative.nongrower.metrics, "qSIP_output/RM_DIM_R2_18O_GL_td_pos_resid_putative_nongrower_metrics.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R2.18O.GL.td.pos.resid.list$taxa.in, "qSIP_output/RM_DIM_R2_18O_GL_td_pos_resid_taxa_in.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R2.18O.GL.td.pos.resid.list$taxa.out, "qSIP_output/RM_DIM_R2_18O_GL_td_pos_resid_taxa_out.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R2.18O.GL.td.pos.resid.list$data.low.SD, "qSIP_output/RM_DIM_R2_18O_GL_td_pos_resid_data_low_SD.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
    #td.abs.resid results:
      write.table(R2.18O.GL.td.abs.resid.list$putative.nongrower.metrics, "qSIP_output/RM_DIM_R2_18O_GL_td_abs_resid_putative_nongrower_metrics.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R2.18O.GL.td.abs.resid.list$taxa.in, "qSIP_output/RM_DIM_R2_18O_GL_td_abs_resid_taxa_in.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R2.18O.GL.td.abs.resid.list$taxa.out, "qSIP_output/RM_DIM_R2_18O_GL_td_abs_resid_taxa_out.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R2.18O.GL.td.abs.resid.list$data.low.SD, "qSIP_output/RM_DIM_R2_18O_GL_td_abs_resid_data_low_SD.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
    #bu.abs.resid results:
      write.table(R2.18O.GL.bu.abs.resid.list$putative.nongrower.metrics, "qSIP_output/RM_DIM_R2_18O_GL_bu_abs_resid_putative_nongrower_metrics.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R2.18O.GL.bu.abs.resid.list$taxa.in, "qSIP_output/RM_DIM_R2_18O_GL_bu_abs_resid_taxa_in.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R2.18O.GL.bu.abs.resid.list$taxa.out, "qSIP_output/RM_DIM_R2_18O_GL_bu_abs_resid_taxa_out.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R2.18O.GL.bu.abs.resid.list$data.low.SD, "qSIP_output/RM_DIM_R2_18O_GL_bu_abs_resid_data_low_SD.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)

  #Find the best iteration:
    R2.18O.GL.td.pos.resid.best.iteration.list <- select.best.iteration(putative.nongrower.metrics=R2.18O.GL.td.pos.resid.list$putative.nongrower.metrics, quantile.threshold=0.10)
    R2.18O.GL.td.pos.resid.best.iteration.list
    find.labeled.correction.plot(find.labeled.correction.list=R2.18O.GL.td.pos.resid.list, filename="RM_DIM_GL_R2_18O_GL_td_pos_resid_plots", path="qSIP_output/Figures/", highlight=R2.18O.GL.td.pos.resid.best.iteration.list$best.iteration, highlight.col="magenta")

    R2.18O.GL.td.abs.resid.best.iteration.list <- select.best.iteration(putative.nongrower.metrics=R2.18O.GL.td.abs.resid.list$putative.nongrower.metrics, quantile.threshold=0.10)
    R2.18O.GL.td.abs.resid.best.iteration.list
    find.labeled.correction.plot(find.labeled.correction.list=R2.18O.GL.td.abs.resid.list, filename="RM_DIM_GL_R2_18O_GL_td_abs_resid_plots", path="qSIP_output/Figures/", highlight=R2.18O.GL.td.abs.resid.best.iteration.list$best.iteration, highlight.col="magenta")

    R2.18O.GL.bu.abs.resid.best.iteration.list <- select.best.iteration(putative.nongrower.metrics=R2.18O.GL.bu.abs.resid.list$putative.nongrower.metrics, quantile.threshold=0.10)
    R2.18O.GL.bu.abs.resid.best.iteration.list
    find.labeled.correction.plot(find.labeled.correction.list=R2.18O.GL.bu.abs.resid.list, filename="RM_DIM_GL_R2_18O_GL_bu_abs_resid_plots", path="qSIP_output/Figures/", highlight=R2.18O.GL.bu.abs.resid.best.iteration.list$best.iteration, highlight.col="magenta")

  #Get the sequence of taxa added at each iteration:
    get.seq.taxa.nums(R2.18O.GL.td.pos.resid.list$taxa.in)
    get.seq.taxa.nums(R2.18O.GL.td.abs.resid.list$taxa.in)
    get.seq.taxa.nums(R2.18O.GL.bu.abs.resid.list$taxa.in)


# Compare replicate #3 of the 18O Grassland treatment to the mean of all (9) unlabeled replicates:

  #Iterate through the remaining taxa to create a data frame of metrics associated with different subsets of putative nongrowing taxa in the labeled treatment for use in assessing a 'match' with unlabeled taxa.
  #And simultaneously create three more data frames: one that lists the taxa included and one the taxa excluded for each iteration, the function also returns the raw data filtered according to the specifiec SD thresholds.
  #Note: the 'fit.norm.func' function calls in these functions slow them process down substantially (i.e., for ~60 taxa: from <1min without them to ~10min with them)
    #Start with all data after imposing SD thresholds; iteratively remove the data point with the largest positive residual:
      system.time(R3.18O.GL.td.pos.resid.list <- find.labeled.correction(LIST=unlab.WAD.corr.18O.GL.list, lab.replicate="R3.18O.GL", lab.names=c("R1.18O.GL", "R2.18O.GL", "R3.18O.GL"), method="td.pos.resid", unlab.SD.percentile=0.5, lab.SD.percentile=0.5, min.num.nongrowers=10))
    #Start with all data after imposing SD thresholds; iteratively remove the data point with the largest absolute value residual:
      system.time(R3.18O.GL.td.abs.resid.list <- find.labeled.correction(LIST=unlab.WAD.corr.18O.GL.list, lab.replicate="R3.18O.GL", lab.names=c("R1.18O.GL", "R2.18O.GL", "R3.18O.GL"), method="td.abs.resid", unlab.SD.percentile=0.5, lab.SD.percentile=0.5, min.num.nongrowers=10))
    #Start with only the set of taxa having the lowest delta-WAD after imposing SD thresholds; iteratively add taxa with the next-smallest absolute value residual (this is effectively the same as 'td.pos.resid'):
      system.time(R3.18O.GL.bu.abs.resid.list <- find.labeled.correction(LIST=unlab.WAD.corr.18O.GL.list, lab.replicate="R3.18O.GL", lab.names=c("R1.18O.GL", "R2.18O.GL", "R3.18O.GL"), method="bu.abs.resid", unlab.SD.percentile=0.5, lab.SD.percentile=0.5, min.num.nongrowers=10))

  #Write the find.labeled.correction output tables to text files:
    #td.pos.resid results:
      write.table(R3.18O.GL.td.pos.resid.list$putative.nongrower.metrics, "qSIP_output/RM_DIM_R3_18O_GL_td_pos_resid_putative_nongrower_metrics.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R3.18O.GL.td.pos.resid.list$taxa.in, "qSIP_output/RM_DIM_R3_18O_GL_td_pos_resid_taxa_in.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R3.18O.GL.td.pos.resid.list$taxa.out, "qSIP_output/RM_DIM_R3_18O_GL_td_pos_resid_taxa_out.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R3.18O.GL.td.pos.resid.list$data.low.SD, "qSIP_output/RM_DIM_R3_18O_GL_td_pos_resid_data_low_SD.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
    #td.abs.resid results:
      write.table(R3.18O.GL.td.abs.resid.list$putative.nongrower.metrics, "qSIP_output/RM_DIM_R3_18O_GL_td_abs_resid_putative_nongrower_metrics.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R3.18O.GL.td.abs.resid.list$taxa.in, "qSIP_output/RM_DIM_R3_18O_GL_td_abs_resid_taxa_in.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R3.18O.GL.td.abs.resid.list$taxa.out, "qSIP_output/RM_DIM_R3_18O_GL_td_abs_resid_taxa_out.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R3.18O.GL.td.abs.resid.list$data.low.SD, "qSIP_output/RM_DIM_R3_18O_GL_td_abs_resid_data_low_SD.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
    #bu.abs.resid results:
      write.table(R3.18O.GL.bu.abs.resid.list$putative.nongrower.metrics, "qSIP_output/RM_DIM_R3_18O_GL_bu_abs_resid_putative_nongrower_metrics.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R3.18O.GL.bu.abs.resid.list$taxa.in, "qSIP_output/RM_DIM_R3_18O_GL_bu_abs_resid_taxa_in.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R3.18O.GL.bu.abs.resid.list$taxa.out, "qSIP_output/RM_DIM_R3_18O_GL_bu_abs_resid_taxa_out.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R3.18O.GL.bu.abs.resid.list$data.low.SD, "qSIP_output/RM_DIM_R3_18O_GL_bu_abs_resid_data_low_SD.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)

  #Find the best iteration:
    R3.18O.GL.td.pos.resid.best.iteration.list <- select.best.iteration(putative.nongrower.metrics=R3.18O.GL.td.pos.resid.list$putative.nongrower.metrics, quantile.threshold=0.10)
    R3.18O.GL.td.pos.resid.best.iteration.list
    find.labeled.correction.plot(find.labeled.correction.list=R3.18O.GL.td.pos.resid.list, filename="RM_DIM_GL_R3_18O_GL_td_pos_resid_plots", path="qSIP_output/Figures/", highlight=R3.18O.GL.td.pos.resid.best.iteration.list$best.iteration, highlight.col="magenta")

    R3.18O.GL.td.abs.resid.best.iteration.list <- select.best.iteration(putative.nongrower.metrics=R3.18O.GL.td.abs.resid.list$putative.nongrower.metrics, quantile.threshold=0.10)
    R3.18O.GL.td.abs.resid.best.iteration.list
    find.labeled.correction.plot(find.labeled.correction.list=R3.18O.GL.td.abs.resid.list, filename="RM_DIM_GL_R3_18O_GL_td_abs_resid_plots", path="qSIP_output/Figures/", highlight=R3.18O.GL.td.abs.resid.best.iteration.list$best.iteration, highlight.col="magenta")

    R3.18O.GL.bu.abs.resid.best.iteration.list <- select.best.iteration(putative.nongrower.metrics=R3.18O.GL.bu.abs.resid.list$putative.nongrower.metrics, quantile.threshold=0.10)
    R3.18O.GL.bu.abs.resid.best.iteration.list
    find.labeled.correction.plot(find.labeled.correction.list=R3.18O.GL.bu.abs.resid.list, filename="RM_DIM_GL_R3_18O_GL_bu_abs_resid_plots", path="qSIP_output/Figures/", highlight=R3.18O.GL.bu.abs.resid.best.iteration.list$best.iteration, highlight.col="magenta")

  #Get the sequence of taxa added at each iteration:
    get.seq.taxa.nums(R3.18O.GL.td.pos.resid.list$taxa.in)
    get.seq.taxa.nums(R3.18O.GL.td.abs.resid.list$taxa.in)
    get.seq.taxa.nums(R3.18O.GL.bu.abs.resid.list$taxa.in)


# Compare the mean of the three 18O replicates of the Grassland ecosystem to the mean of all (9) unlabeled replicates:
# (NOTE: using the mean of three replicates doesn't work as well as correcting each replicate one at a time; it is shown here only for illustrative purposes)

  #Iterate through the remaining taxa to create a data frame of metrics associated with different subsets of putative nongrowing taxa in the labeled treatment for use in assessing a 'match' with unlabeled taxa.
  #And simultaneously create three more data frames: one that lists the taxa included and one the taxa excluded for each iteration, the function also returns the raw data filtered according to the specifiec SD thresholds.
  #Note: the 'fit.norm.func' function calls in these functions slow them process down substantially (i.e., for ~60 taxa: from <1min without them to ~10min with them)
    #Start with all data after imposing SD thresholds; iteratively remove the data point with the largest positive residual:
      system.time(mean.18O.GL.td.pos.resid.list <- find.labeled.correction(LIST=unlab.WAD.corr.18O.GL.list, lab.replicate="lab.mean", lab.names=c("R1.18O.GL", "R2.18O.GL", "R3.18O.GL"), method="td.pos.resid", unlab.SD.percentile=0.5, lab.SD.percentile=0.5, min.num.nongrowers=10))

  #Write the find.labeled.correction output tables to text files:
    #td.pos.resid results:
      write.table(mean.18O.GL.td.pos.resid.list$putative.nongrower.metrics, "qSIP_output/RM_DIM_mean_18O_GL_td_pos_resid_putative_nongrower_metrics.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(mean.18O.GL.td.pos.resid.list$taxa.in, "qSIP_output/RM_DIM_mean_18O_GL_td_pos_resid_taxa_in.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(mean.18O.GL.td.pos.resid.list$taxa.out, "qSIP_output/RM_DIM_mean_18O_GL_td_pos_resid_taxa_out.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(mean.18O.GL.td.pos.resid.list$data.low.SD, "qSIP_output/RM_DIM_mean_18O_GL_td_pos_resid_data_low_SD.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)

  #Find the best iteration:
    mean.18O.GL.td.pos.resid.best.iteration.list <- select.best.iteration(putative.nongrower.metrics=mean.18O.GL.td.pos.resid.list$putative.nongrower.metrics, quantile.threshold=0.10)
    mean.18O.GL.td.pos.resid.best.iteration.list
    find.labeled.correction.plot(find.labeled.correction.list=mean.18O.GL.td.pos.resid.list, filename="RM_DIM_GL_mean_18O_GL_td_pos_resid_plots", path="qSIP_output/Figures/", highlight=mean.18O.GL.td.pos.resid.best.iteration.list$best.iteration, highlight.col="magenta")

  #Get the sequence of taxa added at each iteration:
    get.seq.taxa.nums(mean.18O.GL.td.pos.resid.list$taxa.in)


#Compare the iteration number and the set of 'nongrowers' identified for each of the three 18O.GL replicates (only focus on 'td.pos.resid' results:
  #Best iteration identified for each replicate:
    R1.18O.GL.td.pos.resid.best.iteration.list$best.iteration
    R2.18O.GL.td.pos.resid.best.iteration.list$best.iteration
    R3.18O.GL.td.pos.resid.best.iteration.list$best.iteration
  #Define the set of 'nongrowers' for each 18O.GL replicate:
    R1.18O.GL.td.pos.nongrowers <- as.numeric(as.character(R1.18O.GL.td.pos.resid.list$taxa.in[, R1.18O.GL.td.pos.resid.best.iteration.list$best.iteration]))
    R1.18O.GL.td.pos.nongrowers <- R1.18O.GL.td.pos.nongrowers[!is.na(R1.18O.GL.td.pos.nongrowers)]
    R2.18O.GL.td.pos.nongrowers <- as.numeric(as.character(R2.18O.GL.td.pos.resid.list$taxa.in[, R2.18O.GL.td.pos.resid.best.iteration.list$best.iteration]))
    R2.18O.GL.td.pos.nongrowers <- R2.18O.GL.td.pos.nongrowers[!is.na(R2.18O.GL.td.pos.nongrowers)]
    R3.18O.GL.td.pos.nongrowers <- as.numeric(as.character(R3.18O.GL.td.pos.resid.list$taxa.in[, R3.18O.GL.td.pos.resid.best.iteration.list$best.iteration]))
    R3.18O.GL.td.pos.nongrowers <- R3.18O.GL.td.pos.nongrowers[!is.na(R3.18O.GL.td.pos.nongrowers)]
  #Look at the set of 'nongrowers' for each 18O.GL replicate:
    R1.18O.GL.td.pos.nongrowers
    R2.18O.GL.td.pos.nongrowers
    R3.18O.GL.td.pos.nongrowers
  #The set of 'nongrowing' taxa common to all three 18O.GL replicates:
    intersect(intersect(R1.18O.GL.td.pos.nongrowers, R2.18O.GL.td.pos.nongrowers), R3.18O.GL.td.pos.nongrowers)
  #View a Venn diagram of the 'nongrowers' for each 18O.GL replicate:
    #NOTE: the 'VennDiagram' package must be loaded first
    graphics.off()
    n12 <- intersect(R1.18O.GL.td.pos.nongrowers, R2.18O.GL.td.pos.nongrowers)
    n23 <- intersect(R2.18O.GL.td.pos.nongrowers, R3.18O.GL.td.pos.nongrowers)
    n13 <- intersect(R1.18O.GL.td.pos.nongrowers, R3.18O.GL.td.pos.nongrowers)
    n123 <- intersect(intersect(R1.18O.GL.td.pos.nongrowers, R2.18O.GL.td.pos.nongrowers), R3.18O.GL.td.pos.nongrowers)
    draw.triple.venn( area1=length(R1.18O.GL.td.pos.nongrowers), 
                      area2=length(R2.18O.GL.td.pos.nongrowers), 
                      area3=length(R3.18O.GL.td.pos.nongrowers), 
                      n12=length(n12), 
                      n23=length(n23), 
                      n13=length(n13), 
                      n123=length(n123),
                      category=c("R1.18O.GL", "R2.18O.GL", "R3.18O.GL"))


graphics.off()


#Correct the summarized and raw data for tube-level shifts in WAD identified above for the unlabeled and 18O replicates:

  #Apply the tube-level shift to the already-summarized taxon-level labeled replicate WADs to 'correct' them and create a new summary list (analagous to 'unlab.WAD.corr.18O.GL.list') that includes those corrected values:
    WAD.corr.18O.GL.list <- add.lab.WAD.corr.summary(summary.list=unlab.WAD.corr.18O.GL.list, reps.by.trt=WAD.by.taxon.GL$reps.by.trt, lab.reps.to.add=c("R1.18O.GL", "R2.18O.GL", "R3.18O.GL"), lab.shifts=c(R1.18O.GL.td.pos.resid.best.iteration.list$norm.correction, R2.18O.GL.td.pos.resid.best.iteration.list$norm.correction, R3.18O.GL.td.pos.resid.best.iteration.list$norm.correction), CI=0.90)
    #Look at the results:
      names(WAD.corr.18O.GL.list)                           #names of the two data frames in the output list
      WAD.corr.18O.GL.list$WAD.norm.fit.parms               #looking at the first object -- a data frame
      WAD.corr.18O.GL.list[[1]]                             #another way to look at the first object
      WAD.corr.18O.GL.list$corr.names                       #looking at the second object in the list -- a vector of the names of the corrected unlabeled replicates
      WAD.corr.18O.GL.list[[2]]                             #another way to look at the second object
      head(WAD.corr.18O.GL.list$WAD.table.corr)             #looking at the head of the third object -- a data frame
      head(WAD.corr.18O.GL.list[[3]])                       #another way to look at the head of the third object

  #Apply the tube-level shifts to 'correct' the unlabeled WADs in 'data' (corrects all unlabeled replicates at once):
    data.corr <- apply.unlabeled.correction(raw.data=data, correction.table=unlab.WAD.corr.18O.GL.list$WAD.norm.fit.parms, reps.by.trt=WAD.by.taxon.GL$reps.by.trt, vars=c("Density.g.ml", "Sample", "iso.treat.eco"))
  #Apply the tube-level shift to 'correct' the labeled WADs in 'data.corr' (corrects one labeled replicate at a time):
    data.corr <- apply.labeled.correction(raw.data=data, raw.data.corr=data.corr, lab.replicate="R1.18O.GL", correction.value=R1.18O.GL.td.pos.resid.best.iteration.list$norm.correction, reps.by.trt=WAD.by.taxon.GL$reps.by.trt, vars=c("Density.g.ml", "Sample", "iso.treat.eco"))
    data.corr <- apply.labeled.correction(raw.data=data, raw.data.corr=data.corr, lab.replicate="R2.18O.GL", correction.value=R2.18O.GL.td.pos.resid.best.iteration.list$norm.correction, reps.by.trt=WAD.by.taxon.GL$reps.by.trt, vars=c("Density.g.ml", "Sample", "iso.treat.eco"))
    data.corr <- apply.labeled.correction(raw.data=data, raw.data.corr=data.corr, lab.replicate="R3.18O.GL", correction.value=R3.18O.GL.td.pos.resid.best.iteration.list$norm.correction, reps.by.trt=WAD.by.taxon.GL$reps.by.trt, vars=c("Density.g.ml", "Sample", "iso.treat.eco"))

  #Re-calculate number of copies per uL, based on relative abundance and total number of copies per uL:
    ncopies.corr <- data.corr$neat.avg.16S.copies*data.corr[,7:(ncol(data.corr)-1)]
    dim(ncopies.corr)
    ncopies.corr <- cbind(data.corr[,1:6], ncopies.corr)  # add first 6 columns of data.corr to ncopies.corr
    dim(ncopies.corr)
    head(ncopies.corr)

  #Melt data.corr into long format by tube, sample, tmt, rep, fraction, DNA conc, and density;
  #Do this for copies.ul and for relative abundance, which is just our data.corr file. Merge these to into 1 masterfile: data.corr.melted
    ncopies.corr.melted <- melt(ncopies.corr, id=c("Sample", "fraction", "iso.treat.eco", "neat.avg.16S.copies", "Density.g.ml", "DNA.ng.ul"), measure.vars=as.character(1:length(taxa.id$taxon)), variable.name="taxon", value.name="copies.ul")
    rel.abundance.corr.melted <- melt(data.corr, id=c("Sample", "fraction", "iso.treat.eco", "neat.avg.16S.copies", "Density.g.ml", "DNA.ng.ul"), measure.vars=as.character(1:length(taxa.id$taxon)), variable.name="taxon", value.name="rel.abundance")
    data.corr.melted <- merge(ncopies.corr.melted, rel.abundance.corr.melted)
    head(data.corr.melted)

  #Merge taxa data and reorder data frame by taxon and SampleID AND FRACTION
    data.corr.melted <- merge(data.corr.melted, taxa.id)
    data.corr.melted <- data.corr.melted[order(data.corr.melted$taxon, data.corr.melted$Sample, data.corr.melted$fraction),]
    row.names(data.corr.melted) <- 1:dim(data.corr.melted)[1]   #rename observations to be sequential
    head(data.corr.melted)


#Plot tube-level density curves (all taxa lumped together):
  #####{
  #Isolate the data for each replicate of the "16O.GL" treatment:
    id.reps(DATA=ncopies.corr, focal.tmt="16O.GL", vars=c("Sample", "iso.treat.eco"))
    unlab.corr.1 <- select.rep(DATA=ncopies.corr, focal.tmt="16O.GL", replicate.index=1, vars=c("Density.g.ml", "neat.avg.16S.copies", "Sample", "iso.treat.eco"))
    unlab.corr.2 <- select.rep(DATA=ncopies.corr, focal.tmt="16O.GL", replicate.index=2, vars=c("Density.g.ml", "neat.avg.16S.copies", "Sample", "iso.treat.eco"))
    unlab.corr.3 <- select.rep(DATA=ncopies.corr, focal.tmt="16O.GL", replicate.index=3, vars=c("Density.g.ml", "neat.avg.16S.copies", "Sample", "iso.treat.eco"))
  #Isolate the data for each replicate of the "12C.GL" treatment:
    id.reps(DATA=ncopies, focal.tmt="12C.GL", vars=c("Sample", "iso.treat.eco"))
    unlab.corr.4 <- select.rep(DATA=ncopies.corr, focal.tmt="12C.GL", replicate.index=1, vars=c("Density.g.ml", "neat.avg.16S.copies", "Sample", "iso.treat.eco"))
    unlab.corr.5 <- select.rep(DATA=ncopies.corr, focal.tmt="12C.GL", replicate.index=2, vars=c("Density.g.ml", "neat.avg.16S.copies", "Sample", "iso.treat.eco"))
    unlab.corr.6 <- select.rep(DATA=ncopies.corr, focal.tmt="12C.GL", replicate.index=3, vars=c("Density.g.ml", "neat.avg.16S.copies", "Sample", "iso.treat.eco"))
  #Isolate the data for each replicate of the "12C_N.GL" treatment:
    id.reps(DATA=ncopies, focal.tmt="12C_N.GL", vars=c("Sample", "iso.treat.eco"))
    unlab.corr.7 <- select.rep(DATA=ncopies.corr, focal.tmt="12C_N.GL", replicate.index=1, vars=c("Density.g.ml", "neat.avg.16S.copies", "Sample", "iso.treat.eco"))
    unlab.corr.8 <- select.rep(DATA=ncopies.corr, focal.tmt="12C_N.GL", replicate.index=2, vars=c("Density.g.ml", "neat.avg.16S.copies", "Sample", "iso.treat.eco"))
    unlab.corr.9 <- select.rep(DATA=ncopies.corr, focal.tmt="12C_N.GL", replicate.index=3, vars=c("Density.g.ml", "neat.avg.16S.copies", "Sample", "iso.treat.eco"))
  #Isolate the data for each replicate of the "18O.GL" treatment:
    id.reps(DATA=ncopies, focal.tmt="18O.GL", vars=c("Sample", "iso.treat.eco"))
    lab.corr.1 <- select.rep(DATA=ncopies.corr, focal.tmt="18O.GL", replicate.index=1, vars=c("Density.g.ml", "neat.avg.16S.copies", "Sample", "iso.treat.eco"))
    lab.corr.2 <- select.rep(DATA=ncopies.corr, focal.tmt="18O.GL", replicate.index=2, vars=c("Density.g.ml", "neat.avg.16S.copies", "Sample", "iso.treat.eco"))
    lab.corr.3 <- select.rep(DATA=ncopies.corr, focal.tmt="18O.GL", replicate.index=3, vars=c("Density.g.ml", "neat.avg.16S.copies", "Sample", "iso.treat.eco"))

  #Plot density curves for all unlabeled replicates:
    par(mai=c(1.02, 1.02, 0.22, 0.22))
    XLIM <- c(min(c(unlab.corr.1$X, unlab.corr.2$X, unlab.corr.3$X, unlab.corr.4$X, unlab.corr.5$X, unlab.corr.6$X, unlab.corr.7$X, unlab.corr.8$X, unlab.corr.9$X, lab.corr.1$X, lab.corr.2$X, lab.corr.3$X)), max(c(unlab.corr.1$X, unlab.corr.2$X, unlab.corr.3$X, unlab.corr.4$X, unlab.corr.5$X, unlab.corr.6$X, unlab.corr.7$X, unlab.corr.8$X, unlab.corr.9$X, lab.corr.1$X, lab.corr.2$X, lab.corr.3$X)))
    YLIM <- c(min(c(unlab.corr.1$Y, unlab.corr.2$Y, unlab.corr.3$Y, unlab.corr.4$Y, unlab.corr.5$Y, unlab.corr.6$Y, unlab.corr.7$Y, unlab.corr.8$Y, unlab.corr.9$Y, lab.corr.1$Y, lab.corr.2$Y, lab.corr.3$Y)), max(c(unlab.corr.1$Y, unlab.corr.2$Y, unlab.corr.3$Y, unlab.corr.4$Y, unlab.corr.5$Y, unlab.corr.6$Y, unlab.corr.7$Y, unlab.corr.8$Y, unlab.corr.9$Y, lab.corr.1$Y, lab.corr.2$Y, lab.corr.3$Y)))
    plot(unlab.corr.1$Y~unlab.corr.1$X, type="n", pch=21, col="black", bg="green", xlim=XLIM, ylim=YLIM, xlab="Density g/ml", ylab="16S copies", cex=1.25, cex.axis=1.5, cex.lab=1.8)
    points(unlab.corr.1$Y~unlab.corr.1$X, pch=21, cex=1.25, col="black", bg="green")
    points(unlab.corr.2$Y~unlab.corr.2$X, pch=21, cex=1.25, col="black", bg="red")
    points(unlab.corr.3$Y~unlab.corr.3$X, pch=21, cex=1.25, col="black", bg="yellow")
    points(unlab.corr.4$Y~unlab.corr.4$X, pch=21, cex=1.25, col="black", bg="purple")
    points(unlab.corr.5$Y~unlab.corr.5$X, pch=21, cex=1.25, col="black", bg="blue")
    points(unlab.corr.6$Y~unlab.corr.6$X, pch=21, cex=1.25, col="black", bg="orange")
    points(unlab.corr.7$Y~unlab.corr.7$X, pch=21, cex=1.25, col="black", bg="cyan")
    points(unlab.corr.8$Y~unlab.corr.8$X, pch=21, cex=1.25, col="black", bg="magenta")
    points(unlab.corr.9$Y~unlab.corr.9$X, pch=21, cex=1.25, col="black", bg="wheat")
    abline(v=WAD.func(unlab.corr.1$Y, unlab.corr.1$X), cex=1.25, col="green", lwd=2)
    abline(v=WAD.func(unlab.corr.2$Y, unlab.corr.2$X), cex=1.25, col="red", lwd=2)
    abline(v=WAD.func(unlab.corr.3$Y, unlab.corr.3$X), cex=1.25, col="yellow", lwd=2)
    abline(v=WAD.func(unlab.corr.4$Y, unlab.corr.4$X), cex=1.25, col="purple", lwd=2)
    abline(v=WAD.func(unlab.corr.5$Y, unlab.corr.5$X), cex=1.25, col="blue", lwd=2)
    abline(v=WAD.func(unlab.corr.6$Y, unlab.corr.6$X), cex=1.25, col="orange", lwd=2)
    abline(v=WAD.func(unlab.corr.7$Y, unlab.corr.7$X), cex=1.25, col="cyan", lwd=2)
    abline(v=WAD.func(unlab.corr.8$Y, unlab.corr.8$X), cex=1.25, col="magenta", lwd=2)
    abline(v=WAD.func(unlab.corr.9$Y, unlab.corr.9$X), cex=1.25, col="wheat", lwd=2)

  #Plot a density curve for one unlabeled replicate:
    par(mai=c(1.02, 1.02, 0.22, 0.22))
    XLIM <- c(min(c(unlab.corr.1$X, unlab.corr.2$X, unlab.corr.3$X, unlab.corr.4$X, unlab.corr.5$X, unlab.corr.6$X, unlab.corr.7$X, unlab.corr.8$X, unlab.corr.9$X, lab.corr.1$X, lab.corr.2$X, lab.corr.3$X)), max(c(unlab.corr.1$X, unlab.corr.2$X, unlab.corr.3$X, unlab.corr.4$X, unlab.corr.5$X, unlab.corr.6$X, unlab.corr.7$X, unlab.corr.8$X, unlab.corr.9$X, lab.corr.1$X, lab.corr.2$X, lab.corr.3$X)))
    YLIM <- c(min(c(unlab.corr.1$Y, unlab.corr.2$Y, unlab.corr.3$Y, unlab.corr.4$Y, unlab.corr.5$Y, unlab.corr.6$Y, unlab.corr.7$Y, unlab.corr.8$Y, unlab.corr.9$Y, lab.corr.1$Y, lab.corr.2$Y, lab.corr.3$Y)), max(c(unlab.corr.1$Y, unlab.corr.2$Y, unlab.corr.3$Y, unlab.corr.4$Y, unlab.corr.5$Y, unlab.corr.6$Y, unlab.corr.7$Y, unlab.corr.8$Y, unlab.corr.9$Y, lab.corr.1$Y, lab.corr.2$Y, lab.corr.3$Y)))
    plot(unlab.corr.1$Y~unlab.corr.1$X, type="n", pch=21, col="black", bg="green", xlim=XLIM, ylim=YLIM, xlab="Density g/ml", ylab="16S copies", cex=1.25, cex.axis=1.5, cex.lab=1.8)
    points(unlab.corr.8$Y~unlab.corr.8$X, pch=21, cex=1.25, col="black", bg="magenta")
    abline(v=WAD.func(unlab.corr.8$Y, unlab.corr.8$X), cex=1.25, col="magenta", lwd=2)

  #Plot a density curve for one labeled replicate:
    par(mai=c(1.02, 1.02, 0.22, 0.22))
    XLIM <- c(min(c(unlab.corr.1$X, unlab.corr.2$X, unlab.corr.3$X, unlab.corr.4$X, unlab.corr.5$X, unlab.corr.6$X, unlab.corr.7$X, unlab.corr.8$X, unlab.corr.9$X, lab.corr.1$X, lab.corr.2$X, lab.corr.3$X)), max(c(unlab.corr.1$X, unlab.corr.2$X, unlab.corr.3$X, unlab.corr.4$X, unlab.corr.5$X, unlab.corr.6$X, unlab.corr.7$X, unlab.corr.8$X, unlab.corr.9$X, lab.corr.1$X, lab.corr.2$X, lab.corr.3$X)))
    YLIM <- c(min(c(unlab.corr.1$Y, unlab.corr.2$Y, unlab.corr.3$Y, unlab.corr.4$Y, unlab.corr.5$Y, unlab.corr.6$Y, unlab.corr.7$Y, unlab.corr.8$Y, unlab.corr.9$Y, lab.corr.1$Y, lab.corr.2$Y, lab.corr.3$Y)), max(c(unlab.corr.1$Y, unlab.corr.2$Y, unlab.corr.3$Y, unlab.corr.4$Y, unlab.corr.5$Y, unlab.corr.6$Y, unlab.corr.7$Y, unlab.corr.8$Y, unlab.corr.9$Y, lab.corr.1$Y, lab.corr.2$Y, lab.corr.3$Y)))
    plot(unlab.corr.1$Y~unlab.corr.1$X, type="n", pch=21, col="black", bg="green", xlim=XLIM, ylim=YLIM, xlab="Density g/ml", ylab="16S copies", cex=1.25, cex.axis=1.5, cex.lab=1.8)
    points(lab.corr.3$Y~lab.corr.3$X, pch=21, cex=1.25, col="black", bg="brown")
    abline(v=WAD.func(lab.corr.2$Y, lab.corr.2$X), cex=1.25, col="brown", lwd=2)

  #Plot a density curve for one labeled and one unlabeled replicate:
    par(mai=c(1.02, 1.02, 0.22, 0.22))
    XLIM <- c(min(c(unlab.corr.1$X, unlab.corr.2$X, unlab.corr.3$X, unlab.corr.4$X, unlab.corr.5$X, unlab.corr.6$X, unlab.corr.7$X, unlab.corr.8$X, unlab.corr.9$X, lab.corr.1$X, lab.corr.2$X, lab.corr.3$X)), max(c(unlab.corr.1$X, unlab.corr.2$X, unlab.corr.3$X, unlab.corr.4$X, unlab.corr.5$X, unlab.corr.6$X, unlab.corr.7$X, unlab.corr.8$X, unlab.corr.9$X, lab.corr.1$X, lab.corr.2$X, lab.corr.3$X)))
    YLIM <- c(min(c(unlab.corr.1$Y, unlab.corr.2$Y, unlab.corr.3$Y, unlab.corr.4$Y, unlab.corr.5$Y, unlab.corr.6$Y, unlab.corr.7$Y, unlab.corr.8$Y, unlab.corr.9$Y, lab.corr.1$Y, lab.corr.2$Y, lab.corr.3$Y)), max(c(unlab.corr.1$Y, unlab.corr.2$Y, unlab.corr.3$Y, unlab.corr.4$Y, unlab.corr.5$Y, unlab.corr.6$Y, unlab.corr.7$Y, unlab.corr.8$Y, unlab.corr.9$Y, lab.corr.1$Y, lab.corr.2$Y, lab.corr.3$Y)))
    plot(unlab.corr.1$Y~unlab.corr.1$X, type="n", pch=21, col="black", bg="green", xlim=XLIM, ylim=YLIM, xlab="Density g/ml", ylab="16S copies", cex=1.25, cex.axis=1.5, cex.lab=1.8)
    points(unlab.corr.8$Y~unlab.corr.8$X, pch=21, cex=1.25, col="black", bg="magenta")
    points(lab.corr.3$Y~lab.corr.3$X, pch=21, cex=1.25, col="black", bg="brown")
    abline(v=WAD.func(unlab.corr.8$Y, unlab.corr.8$X), cex=1.25, col="magenta", lwd=2)
    abline(v=WAD.func(lab.corr.2$Y, lab.corr.2$X), cex=1.25, col="brown", lwd=2)

  #Plot density curves for two unlabeled replicates:
    par(mai=c(1.02, 1.02, 0.22, 0.22))
    XLIM <- c(min(c(unlab.corr.1$X, unlab.corr.2$X, unlab.corr.3$X, unlab.corr.4$X, unlab.corr.5$X, unlab.corr.6$X, unlab.corr.7$X, unlab.corr.8$X, unlab.corr.9$X, lab.corr.1$X, lab.corr.2$X, lab.corr.3$X)), max(c(unlab.corr.1$X, unlab.corr.2$X, unlab.corr.3$X, unlab.corr.4$X, unlab.corr.5$X, unlab.corr.6$X, unlab.corr.7$X, unlab.corr.8$X, unlab.corr.9$X, lab.corr.1$X, lab.corr.2$X, lab.corr.3$X)))
    YLIM <- c(min(c(unlab.corr.1$Y, unlab.corr.2$Y, unlab.corr.3$Y, unlab.corr.4$Y, unlab.corr.5$Y, unlab.corr.6$Y, unlab.corr.7$Y, unlab.corr.8$Y, unlab.corr.9$Y, lab.corr.1$Y, lab.corr.2$Y, lab.corr.3$Y)), max(c(unlab.corr.1$Y, unlab.corr.2$Y, unlab.corr.3$Y, unlab.corr.4$Y, unlab.corr.5$Y, unlab.corr.6$Y, unlab.corr.7$Y, unlab.corr.8$Y, unlab.corr.9$Y, lab.corr.1$Y, lab.corr.2$Y, lab.corr.3$Y)))
    plot(unlab.corr.1$Y~unlab.corr.1$X, type="n", pch=21, col="black", bg="green", xlim=XLIM, ylim=YLIM, xlab="Density g/ml", ylab="16S copies", cex=1.25, cex.axis=1.5, cex.lab=1.8)
    points(unlab.corr.4$Y~unlab.corr.4$X, pch=21, cex=1.25, col="black", bg="purple")
    points(unlab.corr.8$Y~unlab.corr.8$X, pch=21, cex=1.25, col="black", bg="magenta")
    abline(v=WAD.func(unlab.corr.4$Y, unlab.corr.4$X), cex=1.25, col="purple", lwd=2)
    abline(v=WAD.func(unlab.corr.8$Y, unlab.corr.8$X), cex=1.25, col="magenta", lwd=2)

  #Plot density curves for all labeled replicates:
    par(mai=c(1.02, 1.02, 0.22, 0.22))
    XLIM <- c(min(c(unlab.corr.1$X, unlab.corr.2$X, unlab.corr.3$X, unlab.corr.4$X, unlab.corr.5$X, unlab.corr.6$X, unlab.corr.7$X, unlab.corr.8$X, unlab.corr.9$X, lab.corr.1$X, lab.corr.2$X, lab.corr.3$X)), max(c(unlab.corr.1$X, unlab.corr.2$X, unlab.corr.3$X, unlab.corr.4$X, unlab.corr.5$X, unlab.corr.6$X, unlab.corr.7$X, unlab.corr.8$X, unlab.corr.9$X, lab.corr.1$X, lab.corr.2$X, lab.corr.3$X)))
    YLIM <- c(min(c(unlab.corr.1$Y, unlab.corr.2$Y, unlab.corr.3$Y, unlab.corr.4$Y, unlab.corr.5$Y, unlab.corr.6$Y, unlab.corr.7$Y, unlab.corr.8$Y, unlab.corr.9$Y, lab.corr.1$Y, lab.corr.2$Y, lab.corr.3$Y)), max(c(unlab.corr.1$Y, unlab.corr.2$Y, unlab.corr.3$Y, unlab.corr.4$Y, unlab.corr.5$Y, unlab.corr.6$Y, unlab.corr.7$Y, unlab.corr.8$Y, unlab.corr.9$Y, lab.corr.1$Y, lab.corr.2$Y, lab.corr.3$Y)))
    plot(unlab.corr.1$Y~unlab.corr.1$X, type="n", pch=21, col="black", bg="green", xlim=XLIM, ylim=YLIM, xlab="Density g/ml", ylab="16S copies", cex=1.25, cex.axis=1.5, cex.lab=1.8)
    points(lab.corr.1$Y~lab.corr.1$X, pch=21, cex=1.25, col="black", bg="gray")
    points(lab.corr.2$Y~lab.corr.2$X, pch=21, cex=1.25, col="black", bg="brown")
    points(lab.corr.3$Y~lab.corr.3$X, pch=21, cex=1.25, col="black", bg="darkgoldenrod2")
    abline(v=WAD.func(lab.corr.1$Y, lab.corr.1$X), cex=1.25, col="gray", lwd=2)
    abline(v=WAD.func(lab.corr.2$Y, lab.corr.2$X), cex=1.25, col="brown", lwd=2)
    abline(v=WAD.func(lab.corr.3$Y, lab.corr.3$X), cex=1.25, col="darkgoldenrod2", lwd=2)

  #Plot density curves for two labeled replicates:
    par(mai=c(1.02, 1.02, 0.22, 0.22))
    XLIM <- c(min(c(unlab.corr.1$X, unlab.corr.2$X, unlab.corr.3$X, unlab.corr.4$X, unlab.corr.5$X, unlab.corr.6$X, unlab.corr.7$X, unlab.corr.8$X, unlab.corr.9$X, lab.corr.1$X, lab.corr.2$X, lab.corr.3$X)), max(c(unlab.corr.1$X, unlab.corr.2$X, unlab.corr.3$X, unlab.corr.4$X, unlab.corr.5$X, unlab.corr.6$X, unlab.corr.7$X, unlab.corr.8$X, unlab.corr.9$X, lab.corr.1$X, lab.corr.2$X, lab.corr.3$X)))
    YLIM <- c(min(c(unlab.corr.1$Y, unlab.corr.2$Y, unlab.corr.3$Y, unlab.corr.4$Y, unlab.corr.5$Y, unlab.corr.6$Y, unlab.corr.7$Y, unlab.corr.8$Y, unlab.corr.9$Y, lab.corr.1$Y, lab.corr.2$Y, lab.corr.3$Y)), max(c(unlab.corr.1$Y, unlab.corr.2$Y, unlab.corr.3$Y, unlab.corr.4$Y, unlab.corr.5$Y, unlab.corr.6$Y, unlab.corr.7$Y, unlab.corr.8$Y, unlab.corr.9$Y, lab.corr.1$Y, lab.corr.2$Y, lab.corr.3$Y)))
    plot(unlab.corr.1$Y~unlab.corr.1$X, type="n", pch=21, col="black", bg="green", xlim=XLIM, ylim=YLIM, xlab="Density g/ml", ylab="16S copies", cex=1.25, cex.axis=1.5, cex.lab=1.8)
    points(lab.corr.2$Y~lab.corr.2$X, pch=21, cex=1.25, col="black", bg="gray")
    points(lab.corr.3$Y~lab.corr.3$X, pch=21, cex=1.25, col="black", bg="brown")
    abline(v=WAD.func(lab.corr.1$Y, lab.corr.1$X), cex=1.25, col="gray", lwd=2)
    abline(v=WAD.func(lab.corr.2$Y, lab.corr.2$X), cex=1.25, col="brown", lwd=2)

  #Plot a density curve for one labeled and one unlabaled replicate with labeled tube less than unlabeled tube!:
    par(mai=c(1.02, 1.02, 0.22, 0.22))
    XLIM <- c(min(c(unlab.corr.1$X, unlab.corr.2$X, unlab.corr.3$X, unlab.corr.4$X, unlab.corr.5$X, unlab.corr.6$X, unlab.corr.7$X, unlab.corr.8$X, unlab.corr.9$X, lab.corr.1$X, lab.corr.2$X, lab.corr.3$X)), max(c(unlab.corr.1$X, unlab.corr.2$X, unlab.corr.3$X, unlab.corr.4$X, unlab.corr.5$X, unlab.corr.6$X, unlab.corr.7$X, unlab.corr.8$X, unlab.corr.9$X, lab.corr.1$X, lab.corr.2$X, lab.corr.3$X)))
    YLIM <- c(min(c(unlab.corr.1$Y, unlab.corr.2$Y, unlab.corr.3$Y, unlab.corr.4$Y, unlab.corr.5$Y, unlab.corr.6$Y, unlab.corr.7$Y, unlab.corr.8$Y, unlab.corr.9$Y, lab.corr.1$Y, lab.corr.2$Y, lab.corr.3$Y)), max(c(unlab.corr.1$Y, unlab.corr.2$Y, unlab.corr.3$Y, unlab.corr.4$Y, unlab.corr.5$Y, unlab.corr.6$Y, unlab.corr.7$Y, unlab.corr.8$Y, unlab.corr.9$Y, lab.corr.1$Y, lab.corr.2$Y, lab.corr.3$Y)))
    plot(unlab.corr.1$Y~unlab.corr.1$X, type="n", pch=21, col="black", bg="green", xlim=XLIM, ylim=YLIM, xlab="Density g/ml", ylab="16S copies", cex=1.25, cex.axis=1.5, cex.lab=1.8)
    points(unlab.corr.4$Y~unlab.corr.4$X, pch=21, cex=1.25, col="black", bg="purple")
    points(lab.corr.2$Y~lab.corr.2$X, pch=21, cex=1.25, col="black", bg="gray")
    abline(v=WAD.func(lab.corr.1$Y, lab.corr.1$X), cex=1.25, col="gray", lwd=2)
    abline(v=WAD.func(unlab.corr.4$Y, unlab.corr.4$X), cex=1.25, col="purple", lwd=2)
  #####}

graphics.off()

#Calculate corrected WADs for each tube and taxon (Grassland only; this actually only calculates the corrected WADS for the unlabeled and 18O treatments):
  data.corr.melted.GL <- data.corr.melted[data.corr.melted$iso.treat.eco %in% c("12C_18O_N.GL", "12C_18O.GL", "12C_N.GL", "12C.GL", "13C_N.GL", "13C.GL", "16O.GL", "18O.GL"),]
  data.corr.melted.GL$taxon <- factor(data.corr.melted.GL$taxon)
  data.corr.melted.GL$Sample <- factor(data.corr.melted.GL$Sample)
  data.corr.melted.GL$iso.treat.eco <- factor(data.corr.melted.GL$iso.treat.eco)
  row.names(data.corr.melted.GL) <- 1:dim(data.corr.melted.GL)[1]
  WAD.by.taxon.GL.corr <- WAD.by.taxon.func(X=data.corr.melted.GL, vars=c("taxon", "Density.g.ml", "copies.ul", "Sample", "iso.treat.eco"))
  #Look at the results:
    names(WAD.by.taxon.GL.corr)          #names of the two data frames in the output list
    head(WAD.by.taxon.GL.corr$obs.wads)  #looking at the head of the first data frame
    head(WAD.by.taxon.GL.corr[[1]])      #another way to look at the head of the first data frame
    WAD.by.taxon.GL.corr$reps.by.trt     #looking at the head of the second data frame
    WAD.by.taxon.GL.corr[[2]]            #another way to look at the head of the second data frame

  #Quick and dirty look at the level of replication among treatments and the occurrences of data for taxa in tubes 
    table(data.corr.melted.GL$iso.treat.eco, data.corr.melted.GL$Sample)
    table(data.corr.melted.GL$taxon, data.corr.melted.GL$Sample)


#Unsorted line graph of corrected WADs -- taxon numbers are the same across tubes -- (all unlabeled & 18O replicates: 4 trts * 3 reps/trt = 12 tubes):
#(note: because taxa are unsorted, excluded taxa (i.e., ones that did not occur in all replicates) are skipped when drawing lines between the taxa that are present in all replicates)
#####{
  graphics.off()
  dev.new(width=6, height=7.74)
  par(mai=c(0.44,0.57,0.22,1.60), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
    #Select only the data for taxa common to all replicates of the specified treatments:
    obs.wads.by.taxon <- WAD.by.taxon.GL.corr$obs.wads
    reps.by.trt <- WAD.by.taxon.GL.corr$reps.by.trt
    treatments=c("16O.GL", "12C.GL", "12C_N.GL", "18O.GL")
    reps.by.trt <- reps.by.trt[reps.by.trt[,1] %in% treatments,]
      reps.by.trt[,1] <- factor(reps.by.trt[,1])
      reps.by.trt[,2] <- factor(reps.by.trt[,2])
      reps.by.trt[,3] <- factor(reps.by.trt[,3])
      reps.by.trt[,4] <- factor(reps.by.trt[,4])
    curr.data <- cbind(taxon=obs.wads.by.taxon[,1], obs.wads.by.taxon[,names(obs.wads.by.taxon) %in% as.character(unlist(reps.by.trt[,2:4]))])
      #Remove any taxon that does not occur in all 12 replicates:
      curr.data <- curr.data[apply(is.na(curr.data[2:dim(curr.data)[2]]), 1, sum) == 0, ]
      curr.data$taxon <- factor(curr.data$taxon)
    #Set the x-limits for the graph:
    x.min <- min(curr.data[,2:dim(curr.data)[2]], na.rm=TRUE)
    x.max <- max(curr.data[,2:dim(curr.data)[2]], na.rm=TRUE)
    AT.y <- seq(0, max(as.numeric(as.character(obs.wads.by.taxon$taxon))), 50)
  plot(x=curr.data[,2], y=as.numeric(as.character(curr.data$taxon)), type="n", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min, x.max))
  par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
  axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
  mtext(expression(paste("Observed WAD (g ml"^-1, ")", sep="")), side=1, line=1.35, cex=0.75)
  par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
  axis(side=2, at=AT.y, labels=AT.y, tck=-0.015, las=1, cex.axis=0.6)
  mtext("Taxon", side=2, line=2.35, cex=0.75)
  #16O.GL
    # points(x=curr.data[,2], y=as.numeric(as.character(curr.data$taxon)), pch=21, cex=0.7, col="black", bg="green")           #R1
    # points(x=curr.data[,9], y=as.numeric(as.character(curr.data$taxon)), pch=21, cex=0.7, col="black", bg="red")             #R2
    # points(x=curr.data[,10], y=as.numeric(as.character(curr.data$taxon)), pch=21, cex=0.7, col="black", bg="yellow")         #R3
    lines(x=curr.data[,2], y=as.numeric(as.character(curr.data$taxon)), col="green")
    lines(x=curr.data[,9], y=as.numeric(as.character(curr.data$taxon)), col="red")
    lines(x=curr.data[,10], y=as.numeric(as.character(curr.data$taxon)), col="yellow")
  #12C.GL
    # points(x=curr.data[,11], y=as.numeric(as.character(curr.data$taxon)), pch=21, cex=0.7, col="black", bg="purple")         #R1
    # points(x=curr.data[,12], y=as.numeric(as.character(curr.data$taxon)), pch=21, cex=0.7, col="black", bg="blue")           #R2
    # points(x=curr.data[,13], y=as.numeric(as.character(curr.data$taxon)), pch=21, cex=0.7, col="black", bg="orange")         #R3
    lines(x=curr.data[,11], y=as.numeric(as.character(curr.data$taxon)), col="purple")
    lines(x=curr.data[,12], y=as.numeric(as.character(curr.data$taxon)), col="blue")
    lines(x=curr.data[,13], y=as.numeric(as.character(curr.data$taxon)), col="orange")
  #12C_N.GL
    # points(x=curr.data[,3], y=as.numeric(as.character(curr.data$taxon)), pch=21, cex=0.7, col="black", bg="cyan")            #R1
    # points(x=curr.data[,4], y=as.numeric(as.character(curr.data$taxon)), pch=21, cex=0.7, col="black", bg="magenta")         #R2
    # points(x=curr.data[,5], y=as.numeric(as.character(curr.data$taxon)), pch=21, cex=0.7, col="black", bg="wheat")           #R3
    lines(x=curr.data[,3], y=as.numeric(as.character(curr.data$taxon)), col="cyan")
    lines(x=curr.data[,4], y=as.numeric(as.character(curr.data$taxon)), col="magenta")
    lines(x=curr.data[,5], y=as.numeric(as.character(curr.data$taxon)), col="wheat")
  #18O.GL
    # points(x=curr.data[,6], y=as.numeric(as.character(curr.data$taxon)), pch=24, cex=0.7, col="black", bg="gray")            #R1
    # points(x=curr.data[,7], y=as.numeric(as.character(curr.data$taxon)), pch=24, cex=0.7, col="black", bg="brown")           #R2
    # points(x=curr.data[,8], y=as.numeric(as.character(curr.data$taxon)), pch=24, cex=0.7, col="black", bg="darkgoldenrod2")  #R3
    lines(x=curr.data[,6], y=as.numeric(as.character(curr.data$taxon)), col="gray")
    lines(x=curr.data[,7], y=as.numeric(as.character(curr.data$taxon)), col="brown")
    lines(x=curr.data[,8], y=as.numeric(as.character(curr.data$taxon)), col="darkgoldenrod2")
    par(xpd=NA)
    ink <- c("green", "red", "yellow", "purple", "blue", "orange", "cyan", "magenta", "wheat", "gray", "brown", "darkgoldenrod2")
    legend(x=x.max*1.0055, y=max(as.numeric(as.character(curr.data$taxon)))+(0.04*(max(as.numeric(as.character(curr.data$taxon)))-1)), legend=c("R1.16O.GL", "R2.16O.GL", "R3.16O.GL", "R1.12C.GL", "R2.12C.GL", "R3.12C.GL", "R1.12C_N.GL", "R2.12C_N.GL", "R3.12C_N.GL", "R1.18O.GL", "R2.18O.GL", "R3.18O.GL"), col=ink, pt.bg=ink, pch=c(rep(21,9), 24, 24, 24), pt.cex=0.7, lwd=1, bty="o", cex=0.6, xjust=0, yjust=1)
    par(xpd=FALSE)
    rm(obs.wads.by.taxon, reps.by.trt, treatments, curr.data, ink)
#####}

#Unsorted line graph of corrected WADs -- taxon numbers are the same across tubes -- (only 16O & 18O replicates: 2 trts * 3 reps/trt = 6 tubes):
#(note: because taxa are unsorted, excluded taxa (i.e., ones that did not occur in all replicates) are skipped when drawing lines between the taxa that are present in all replicates)
#(note: the three 16O reps are all black in this plot)
#####{
  graphics.off()
  dev.new(width=6, height=7.74)
  par(mai=c(0.44,0.57,0.22,1.60), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
    #Select only the data for taxa common to all replicates of the specified treatments:
    obs.wads.by.taxon <- WAD.by.taxon.GL.corr$obs.wads
    reps.by.trt <- WAD.by.taxon.GL.corr$reps.by.trt
    treatments=c("16O.GL", "12C.GL", "12C_N.GL", "18O.GL")
    reps.by.trt <- reps.by.trt[reps.by.trt[,1] %in% treatments,]
      reps.by.trt[,1] <- factor(reps.by.trt[,1])
      reps.by.trt[,2] <- factor(reps.by.trt[,2])
      reps.by.trt[,3] <- factor(reps.by.trt[,3])
      reps.by.trt[,4] <- factor(reps.by.trt[,4])
    curr.data <- cbind(taxon=obs.wads.by.taxon[,1], obs.wads.by.taxon[,names(obs.wads.by.taxon) %in% as.character(unlist(reps.by.trt[,2:4]))])
      #Remove any taxon that does not occur in all 12 replicates:
      curr.data <- curr.data[apply(is.na(curr.data[2:dim(curr.data)[2]]), 1, sum) == 0, ]
      curr.data$taxon <- factor(curr.data$taxon)
    #Set the x-limits for the graph:
    x.min <- min(curr.data[,2:dim(curr.data)[2]], na.rm=TRUE)
    x.max <- max(curr.data[,2:dim(curr.data)[2]], na.rm=TRUE)
    AT.y <- seq(0, max(as.numeric(as.character(obs.wads.by.taxon$taxon))), 50)
  plot(x=curr.data[,2], y=as.numeric(as.character(curr.data$taxon)), type="n", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min, x.max))
  par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
  axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
  mtext(expression(paste("Observed WAD (g ml"^-1, ")", sep="")), side=1, line=1.35, cex=0.75)
  par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
  axis(side=2, at=AT.y, labels=AT.y, tck=-0.015, las=1, cex.axis=0.6)
  mtext("Taxon", side=2, line=2.35, cex=0.75)
  #16O.GL
    # points(x=curr.data[,2], y=as.numeric(as.character(curr.data$taxon)), pch=21, cex=0.7, col="black", bg="black")         #R1
    # points(x=curr.data[,9], y=as.numeric(as.character(curr.data$taxon)), pch=21, cex=0.7, col="black", bg="black")         #R2
    # points(x=curr.data[,10], y=as.numeric(as.character(curr.data$taxon)), pch=21, cex=0.7, col="black", bg="black")        #R3
    lines(x=curr.data[,2], y=as.numeric(as.character(curr.data$taxon)), col="black")
    lines(x=curr.data[,9], y=as.numeric(as.character(curr.data$taxon)), col="black")
    lines(x=curr.data[,10], y=as.numeric(as.character(curr.data$taxon)), col="black")
  #18O.GL
    # points(x=curr.data[,6], y=as.numeric(as.character(curr.data$taxon)), pch=24, cex=0.7, col="black", bg="gray")            #R1
    # points(x=curr.data[,7], y=as.numeric(as.character(curr.data$taxon)), pch=24, cex=0.7, col="black", bg="brown")           #R2
    # points(x=curr.data[,8], y=as.numeric(as.character(curr.data$taxon)), pch=24, cex=0.7, col="black", bg="darkgoldenrod2")  #R3
    lines(x=curr.data[,6], y=as.numeric(as.character(curr.data$taxon)), col="gray")
    lines(x=curr.data[,7], y=as.numeric(as.character(curr.data$taxon)), col="brown")
    lines(x=curr.data[,8], y=as.numeric(as.character(curr.data$taxon)), col="darkgoldenrod2")
    par(xpd=NA)
    ink <- c("black", "black", "black", "gray", "brown", "darkgoldenrod2")
    legend(x=x.max*1.0055, y=max(as.numeric(as.character(curr.data$taxon)))+(0.04*(max(as.numeric(as.character(curr.data$taxon)))-1)), legend=c("R1.16O.GL", "R2.16O.GL", "R3.16O.GL", "R1.18O.GL", "R2.18O.GL", "R3.18O.GL"), col=ink, pt.bg=ink, pch=c(rep(21,3), 24, 24, 24), pt.cex=0.7, lwd=1, bty="o", cex=0.6, xjust=0, yjust=1)
    par(xpd=FALSE)
    rm(obs.wads.by.taxon, reps.by.trt, treatments, curr.data, ink)
#####}

#Resorted graph of corrected WADs -- taxon numbers are NOT the same across tubes -- (all replicates: 4 trts * 3 reps/trt = 12 tubes):
#####{
  graphics.off()
  # dev.new(width=6, height=7.74)
  pdf(file="qSIP_output/Figures/RM_DIM_GL_WAD_by_taxon_resorted_all_reps_corr.pdf", width=6, height=7.74)
  par(mai=c(0.44,0.57,0.22,1.60), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
    #Select the data to plot:
    obs.wads.by.taxon <- WAD.by.taxon.GL.corr$obs.wads
    curr.data.orders <- WAD.corr.18O.GL.list$WAD.table.corr
    #Set the x-limits for the graph:
    x.min <- min(curr.data.orders[,c(19:21,27:35)], na.rm=TRUE)
    x.max <- max(curr.data.orders[,c(19:21,27:35)], na.rm=TRUE)
    AT.y <- seq(0, max(as.numeric(as.character(obs.wads.by.taxon$taxon))), 50)
    plot(x=curr.data.orders$R1.16O.GL, y=curr.data.orders$R1.16O.GL.order, type="n", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min, x.max))
  par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
  axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
  mtext(expression(paste("Observed WAD (g ml"^-1, ")", sep="")), side=1, line=1.35, cex=0.75)
  par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
  axis(side=2, at=AT.y, labels=AT.y, tck=-0.015, las=1, cex.axis=0.6)
  mtext("Taxon (each replicate ranked independently)", side=2, line=2.35, cex=0.75)
  #16O.GL
    points(x=curr.data.orders$R1.16O.GL.corr, y=curr.data.orders$R1.16O.GL.order, pch=21, cex=0.7, col="black", bg="green")
    points(x=curr.data.orders$R2.16O.GL.corr, y=curr.data.orders$R2.16O.GL.order, pch=21, cex=0.7, col="black", bg="red")
    points(x=curr.data.orders$R3.16O.GL.corr, y=curr.data.orders$R3.16O.GL.order, pch=21, cex=0.7, col="black", bg="yellow")
  #12C.GL
    points(x=curr.data.orders$R1.12C.GL.corr, y=curr.data.orders$R1.12C.GL.order, pch=21, cex=0.7, col="black", bg="purple")
    points(x=curr.data.orders$R2.12C.GL.corr, y=curr.data.orders$R2.12C.GL.order, pch=21, cex=0.7, col="black", bg="blue")
    points(x=curr.data.orders$R3.12C.GL.corr, y=curr.data.orders$R3.12C.GL.order, pch=21, cex=0.7, col="black", bg="orange")
  #12C_N.GL
    points(x=curr.data.orders$R1.12C_N.GL.corr, y=curr.data.orders$R1.12C_N.GL.order, pch=21, cex=0.7, col="black", bg="cyan")
    points(x=curr.data.orders$R2.12C_N.GL.corr, y=curr.data.orders$R2.12C_N.GL.order, pch=21, cex=0.7, col="black", bg="magenta")
    points(x=curr.data.orders$R3.12C_N.GL.corr, y=curr.data.orders$R3.12C_N.GL.order, pch=21, cex=0.7, col="black", bg="wheat")
  #18.O.GL
    points(x=curr.data.orders$R1.18O.GL.corr, y=curr.data.orders$R1.18O.GL.order, pch=21, cex=0.7, col="black", bg="gray")
    points(x=curr.data.orders$R2.18O.GL.corr, y=curr.data.orders$R2.18O.GL.order, pch=21, cex=0.7, col="black", bg="brown")
    points(x=curr.data.orders$R3.18O.GL.corr, y=curr.data.orders$R3.18O.GL.order, pch=21, cex=0.7, col="black", bg="darkgoldenrod2")
    par(xpd=NA)
    ink <- c("green", "red", "yellow", "purple", "blue", "orange", "cyan", "magenta", "wheat", "gray", "brown", "darkgoldenrod2")
    legend(x=x.max*1.0055, y=max(as.numeric(as.character(curr.data.orders$mean.unlab.WAD.order)))+(0.04*(max(as.numeric(as.character(curr.data.orders$mean.unlab.WAD.order)))-1)), legend=c("R1.16O.GL", "R2.16O.GL", "R3.16O.GL", "R1.12C.GL", "R2.12C.GL", "R3.12C.GL", "R1.12C_N.GL", "R2.12C_N.GL", "R3.12C_N.GL", "R1.18O.GL", "R2.18O.GL", "R3.18O.GL"), col=ink, pt.bg=ink, pch=rep(21,12), pt.cex=0.7, lwd=1, bty="o", cex=0.6, xjust=0, yjust=1)
    par(xpd=FALSE)
    rm(obs.wads.by.taxon, curr.data.orders, ink)
#####}
dev.off()   #(run graph code through this line only if creating a pdf version of the plot)

#Sorted graph of corrected WADs -- taxon numbers are the same across tubes -- (only 16O replicates: 3 trts * 3 reps/trt = 9 tubes):
#White line = mean of (uncorrected) unlabeled data
#Green line = normal cdf fit to the mean of all (uncorrected) unlabeled data
#####{
  graphics.off()
  # dev.new(width=6, height=7.74)
  pdf(file="qSIP_output/Figures/RM_DIM_GL_WAD_by_taxon_unlabeled_reps_corr.pdf", width=6, height=7.74)
  par(mai=c(0.44,0.57,0.22,1.60), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
    #Select the data to plot:
    obs.wads.by.taxon <- WAD.by.taxon.GL.corr$obs.wads
    WAD.norm.fit.parms <- WAD.corr.18O.GL.list$WAD.norm.fit.parms
    curr.data.orders <- WAD.corr.18O.GL.list$WAD.table.corr
    unlab.reps <- c("R1.16O.GL", "R2.16O.GL", "R3.16O.GL", "R1.12C.GL", "R2.12C.GL", "R3.12C.GL", "R1.12C_N.GL", "R2.12C_N.GL", "R3.12C_N.GL")
    unlab.WADs <- curr.data.orders[, names(curr.data.orders) %in% unlab.reps]
    #Set the x-limits for the graph:
    x.min <- min(curr.data.orders[,c(19:21,27:35)], na.rm=TRUE)
    x.max <- max(curr.data.orders[,c(19:21,27:35)], na.rm=TRUE)
    AT.y <- seq(0, max(as.numeric(as.character(obs.wads.by.taxon$taxon))), 50)
  plot(x=curr.data.orders$R1.16O.GL, y=curr.data.orders$R1.16O.GL.order, type="n", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min, x.max))
  par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
  axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
  mtext(expression(paste("Observed WAD (g ml"^-1, ")", sep="")), side=1, line=1.35, cex=0.75)
  par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
  axis(side=2, at=AT.y, labels=AT.y, tck=-0.015, las=1, cex.axis=0.6)
  mtext("Taxon", side=2, line=2.4, cex=0.75)
  #16O.GL
    points(x=curr.data.orders$R1.16O.GL.corr, y=curr.data.orders$mean.unlab.WAD.order, pch=21, cex=0.7, col="black", bg="green")
    points(x=curr.data.orders$R2.16O.GL.corr, y=curr.data.orders$mean.unlab.WAD.order, pch=21, cex=0.7, col="black", bg="red")
    points(x=curr.data.orders$R3.16O.GL.corr, y=curr.data.orders$mean.unlab.WAD.order, pch=21, cex=0.7, col="black", bg="yellow")
    # new.order <- order(curr.data.orders$mean.unlab.WAD.order)
    # lines(x=curr.data.orders$R1.16O.GL.corr[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=2, col="green")
    # lines(x=curr.data.orders$R2.16O.GL.corr[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=2, col="red")
    # lines(x=curr.data.orders$R3.16O.GL.corr[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=2, col="black")
    # lines(x=curr.data.orders$R3.16O.GL.corr[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=1.5, col="yellow")
  #12C.GL
    points(x=curr.data.orders$R1.12C.GL.corr, y=curr.data.orders$mean.unlab.WAD.order, pch=21, cex=0.7, col="black", bg="purple")
    points(x=curr.data.orders$R2.12C.GL.corr, y=curr.data.orders$mean.unlab.WAD.order, pch=21, cex=0.7, col="black", bg="blue")
    points(x=curr.data.orders$R3.12C.GL.corr, y=curr.data.orders$mean.unlab.WAD.order, pch=21, cex=0.7, col="black", bg="orange")
    # new.order <- order(curr.data.orders$mean.unlab.WAD.order)
    # lines(x=curr.data.orders$R1.12C.GL.corr[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=2, col="purple")
    # lines(x=curr.data.orders$R2.12C.GL.corr[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=2, col="blue")
    # lines(x=curr.data.orders$R3.12C.GL.corr[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=2, col="orange")
  #12C_N.GL
    points(x=curr.data.orders$R1.12C_N.GL.corr, y=curr.data.orders$mean.unlab.WAD.order, pch=21, cex=0.7, col="black", bg="cyan")
    points(x=curr.data.orders$R2.12C_N.GL.corr, y=curr.data.orders$mean.unlab.WAD.order, pch=21, cex=0.7, col="black", bg="magenta")
    points(x=curr.data.orders$R3.12C_N.GL.corr, y=curr.data.orders$mean.unlab.WAD.order, pch=21, cex=0.7, col="black", bg="wheat")
    # new.order <- order(curr.data.orders$mean.unlab.WAD.order)
    # lines(x=curr.data.orders$R1.12C_N.GL.corr[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=2, col="cyan")
    # lines(x=curr.data.orders$R2.12C_N.GL.corr[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=2, col="magenta")
    # lines(x=curr.data.orders$R3.12C_N.GL.corr[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=2, col="wheat")
  #Unlabeled average:
    # points(x=apply(unlab.WADs, 1, mean, na.rm=TRUE), y=curr.data.orders$mean.unlab.WAD.order, pch=21, cex=0.35, col="black", bg="black")
    new.order <- order(curr.data.orders$mean.unlab.WAD.order)
    lines(x=apply(unlab.WADs, 1, mean, na.rm=TRUE)[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=3, col="black")
    lines(x=apply(unlab.WADs, 1, mean, na.rm=TRUE)[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=2, col="white")
  #Normal cdf fit to the mean of all unlabeled data:
    x <- seq(min(unlab.WADs), max(unlab.WADs), 0.000001)
    y <- pnorm(q=x, mean=WAD.norm.fit.parms$mean[WAD.norm.fit.parms$trt == "unlabeled"], sd=WAD.norm.fit.parms$stdev[WAD.norm.fit.parms$trt == "unlabeled"])*max(curr.data.orders$mean.unlab.WAD.order)
    lines(x, y, lwd=3, col="black")
    lines(x, y, lwd=2, col="lawngreen")
    par(xpd=NA)
    ink <- c("green", "red", "yellow", "purple", "blue", "orange", "cyan", "magenta", "wheat", "white", "lawngreen")
    legend(x=x.max*1.0040, y=max(as.numeric(as.character(curr.data.orders$mean.unlab.WAD.order)))+(0.04*(max(as.numeric(as.character(curr.data.orders$mean.unlab.WAD.order)))-1)), legend=c("R1.16O.GL", "R2.16O.GL", "R3.16O.GL", "R1.12C.GL", "R2.12C.GL", "R3.12C.GL", "R1.12C_N.GL", "R2.12C_N.GL", "R3.12C_N.GL", "unlab empirical mean", "unlabeled normal fit"), col=c(ink[1:9], "black", "black"), pt.bg=ink, pch=rep(21,12), pt.cex=0.7, lwd=c(rep(1, 9), 3, 3), bty="o", cex=0.6, xjust=0, yjust=1)
    legend(x=x.max*1.0040, y=max(as.numeric(as.character(curr.data.orders$mean.unlab.WAD.order)))+(0.04*(max(as.numeric(as.character(curr.data.orders$mean.unlab.WAD.order)))-1)), legend=c("R1.16O.GL", "R2.16O.GL", "R3.16O.GL", "R1.12C.GL", "R2.12C.GL", "R3.12C.GL", "R1.12C_N.GL", "R2.12C_N.GL", "R3.12C_N.GL", "unlab empirical mean", "unlabeled normal fit"), col=ink, pt.bg=ink, pch=rep(21,12), pt.cex=0.7, lwd=c(rep(1, 9), 2, 2), bty="o", cex=0.6, xjust=0, yjust=1)
    par(xpd=FALSE)
    rm(obs.wads.by.taxon, WAD.norm.fit.parms, curr.data.orders, unlab.reps, unlab.WADs, ink)
#####}
dev.off()   #(run graph code through this line only if creating a pdf version of the plot)

#Sorted graph of corrected WADs -- taxon numbers are the same across tubes -- (only 18O replicates: 1 trt * 3 reps/trt = 3 tubes):
#White line = mean of (uncorrected) unlabeled data
#Green line = normal cdf fit to the mean of all (uncorrected) unlabeled data
#####{
  graphics.off()
  # dev.new(width=6, height=7.74)
  pdf(file="qSIP_output/Figures/RM_DIM_GL_WAD_by_taxon_labeled_reps_corr.pdf", width=6, height=7.74)
  par(mai=c(0.44,0.57,0.22,1.60), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
    #Select the data to plot:
    obs.wads.by.taxon <- WAD.by.taxon.GL.corr$obs.wads
    WAD.norm.fit.parms <- WAD.corr.18O.GL.list$WAD.norm.fit.parms
    curr.data.orders <- WAD.corr.18O.GL.list$WAD.table.corr
    unlab.reps <- c("R1.16O.GL", "R2.16O.GL", "R3.16O.GL", "R1.12C.GL", "R2.12C.GL", "R3.12C.GL", "R1.12C_N.GL", "R2.12C_N.GL", "R3.12C_N.GL")
    unlab.WADs <- curr.data.orders[, names(curr.data.orders) %in% unlab.reps]
    #Set the x-limits for the graph:
    x.min <- min(curr.data.orders[,c(19:21,27:35)], na.rm=TRUE)
    x.max <- max(curr.data.orders[,c(19:21,27:35)], na.rm=TRUE)
    AT.y <- seq(0, max(as.numeric(as.character(obs.wads.by.taxon$taxon))), 50)
  plot(x=curr.data.orders$R1.16O.GL, y=curr.data.orders$R1.16O.GL.order, type="n", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min, x.max))
  par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
  axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
  mtext(expression(paste("Observed WAD (g ml"^-1, ")", sep="")), side=1, line=1.35, cex=0.75)
  par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
  axis(side=2, at=AT.y, labels=AT.y, tck=-0.015, las=1, cex.axis=0.6)
  mtext("Taxon", side=2, line=2.4, cex=0.75)
  #18.O.GL
    points(x=curr.data.orders$R1.18O.GL.corr, y=curr.data.orders$mean.unlab.WAD.order, pch=21, cex=0.7, col="black", bg="gray")
    points(x=curr.data.orders$R2.18O.GL.corr, y=curr.data.orders$mean.unlab.WAD.order, pch=21, cex=0.7, col="black", bg="brown")
    points(x=curr.data.orders$R3.18O.GL.corr, y=curr.data.orders$mean.unlab.WAD.order, pch=21, cex=0.7, col="black", bg="darkgoldenrod2")
    # new.order <- order(curr.data.orders$mean.unlab.WAD.order)
    # lines(x=curr.data.orders$R1.18O.GL.corr[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=2, col="gray")
    # lines(x=curr.data.orders$R2.18O.GL.corr[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=2, col="brown")
    # lines(x=curr.data.orders$R3.18O.GL.corr[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=2, col="darkgoldenrod2")
  #Unlabeled average:
    # points(x=apply(unlab.WADs, 1, mean, na.rm=TRUE), y=curr.data.orders$mean.unlab.WAD.order, pch=21, cex=0.35, col="black", bg="black")
    new.order <- order(curr.data.orders$mean.unlab.WAD.order)
    lines(x=apply(unlab.WADs, 1, mean, na.rm=TRUE)[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=3, col="black")
    lines(x=apply(unlab.WADs, 1, mean, na.rm=TRUE)[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=2, col="white")
  #Normal cdf fit to the mean of all unlabeled data:
    x <- seq(min(unlab.WADs), max(unlab.WADs), 0.000001)
    y <- pnorm(q=x, mean=WAD.norm.fit.parms$mean[WAD.norm.fit.parms$trt == "unlabeled"], sd=WAD.norm.fit.parms$stdev[WAD.norm.fit.parms$trt == "unlabeled"])*max(curr.data.orders$mean.unlab.WAD.order)
    lines(x, y, lwd=3, col="black")
    lines(x, y, lwd=2, col="lawngreen")
    par(xpd=NA)
    ink <- c("gray", "brown", "darkgoldenrod2", "white", "lawngreen")
    legend(x=x.max*1.0040, y=max(as.numeric(as.character(curr.data.orders$mean.unlab.WAD.order)))+(0.04*(max(as.numeric(as.character(curr.data.orders$mean.unlab.WAD.order)))-1)), legend=c("R1.18O.GL", "R2.18O.GL", "R3.18O.GL", "unlab empirical mean", "unlabeled normal fit"), col=c(ink[1:3], "black", "black"), pt.bg=ink, pch=rep(21,5), pt.cex=0.7, lwd=c(rep(1, 3), 3, 3), bty="o", cex=0.6, xjust=0, yjust=1)
    legend(x=x.max*1.0040, y=max(as.numeric(as.character(curr.data.orders$mean.unlab.WAD.order)))+(0.04*(max(as.numeric(as.character(curr.data.orders$mean.unlab.WAD.order)))-1)), legend=c("R1.18O.GL", "R2.18O.GL", "R3.18O.GL", "unlab empirical mean", "unlabeled normal fit"), col=ink, pt.bg=ink, pch=rep(21,5), pt.cex=0.7, lwd=c(rep(1, 3), 2, 2), bty="o", cex=0.6, xjust=0, yjust=1)
    par(xpd=FALSE)
    rm(obs.wads.by.taxon, WAD.norm.fit.parms, curr.data.orders, unlab.reps, unlab.WADs, ink)
#####}
dev.off()   #(run graph code through this line only if creating a pdf version of the plot)

#Sorted graph of corrected WADs -- taxon numbers are the same across tubes -- (16O & 18O replicates: 4 trts * 3 reps/trt = 12 tubes):
#White line = mean of (uncorrected) unlabeled data
#Green line = normal cdf fit to the mean of all (uncorrected) unlabeled data
#####{
  graphics.off()
  # dev.new(width=6, height=7.74)
  pdf(file="qSIP_output/Figures/RM_DIM_GL_WAD_by_taxon_all_reps_corr.pdf", width=6, height=7.74)
  par(mai=c(0.44,0.57,0.22,1.60), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
    #Select the data to plot:
    obs.wads.by.taxon <- WAD.by.taxon.GL.corr$obs.wads
    WAD.norm.fit.parms <- WAD.corr.18O.GL.list$WAD.norm.fit.parms
    curr.data.orders <- WAD.corr.18O.GL.list$WAD.table.corr
    unlab.reps <- c("R1.16O.GL", "R2.16O.GL", "R3.16O.GL", "R1.12C.GL", "R2.12C.GL", "R3.12C.GL", "R1.12C_N.GL", "R2.12C_N.GL", "R3.12C_N.GL")
    unlab.WADs <- curr.data.orders[, names(curr.data.orders) %in% unlab.reps]
    #Set the x-limits for the graph:
    x.min <- min(curr.data.orders[,c(19:21,27:35)], na.rm=TRUE)
    x.max <- max(curr.data.orders[,c(19:21,27:35)], na.rm=TRUE)
    AT.y <- seq(0, max(as.numeric(as.character(obs.wads.by.taxon$taxon))), 50)
  plot(x=curr.data.orders$R1.16O.GL, y=curr.data.orders$R1.16O.GL.order, type="n", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min, x.max))
  par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
  axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
  mtext(expression(paste("Observed WAD (g ml"^-1, ")", sep="")), side=1, line=1.35, cex=0.75)
  par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
  axis(side=2, at=AT.y, labels=AT.y, tck=-0.015, las=1, cex.axis=0.6)
  mtext("Taxon", side=2, line=2.4, cex=0.75)
  #16O.GL
    points(x=curr.data.orders$R1.16O.GL.corr, y=curr.data.orders$mean.unlab.WAD.order, pch=21, cex=0.7, col="black", bg="green")
    points(x=curr.data.orders$R2.16O.GL.corr, y=curr.data.orders$mean.unlab.WAD.order, pch=21, cex=0.7, col="black", bg="red")
    points(x=curr.data.orders$R3.16O.GL.corr, y=curr.data.orders$mean.unlab.WAD.order, pch=21, cex=0.7, col="black", bg="yellow")
    # new.order <- order(curr.data.orders$mean.unlab.WAD.order)
    # lines(x=curr.data.orders$R1.16O.GL.corr[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=2, col="green")
    # lines(x=curr.data.orders$R2.16O.GL.corr[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=2, col="red")
    # lines(x=curr.data.orders$R3.16O.GL.corr[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=2, col="black")
    # lines(x=curr.data.orders$R3.16O.GL.corr[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=1.5, col="yellow")
  #12C.GL
    points(x=curr.data.orders$R1.12C.GL.corr, y=curr.data.orders$mean.unlab.WAD.order, pch=21, cex=0.7, col="black", bg="purple")
    points(x=curr.data.orders$R2.12C.GL.corr, y=curr.data.orders$mean.unlab.WAD.order, pch=21, cex=0.7, col="black", bg="blue")
    points(x=curr.data.orders$R3.12C.GL.corr, y=curr.data.orders$mean.unlab.WAD.order, pch=21, cex=0.7, col="black", bg="orange")
    # new.order <- order(curr.data.orders$mean.unlab.WAD.order)
    # lines(x=curr.data.orders$R1.12C.GL.corr[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=2, col="purple")
    # lines(x=curr.data.orders$R2.12C.GL.corr[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=2, col="blue")
    # lines(x=curr.data.orders$R3.12C.GL.corr[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=2, col="orange")
  #12C_N.GL
    points(x=curr.data.orders$R1.12C_N.GL.corr, y=curr.data.orders$mean.unlab.WAD.order, pch=21, cex=0.7, col="black", bg="cyan")
    points(x=curr.data.orders$R2.12C_N.GL.corr, y=curr.data.orders$mean.unlab.WAD.order, pch=21, cex=0.7, col="black", bg="magenta")
    points(x=curr.data.orders$R3.12C_N.GL.corr, y=curr.data.orders$mean.unlab.WAD.order, pch=21, cex=0.7, col="black", bg="wheat")
    # new.order <- order(curr.data.orders$mean.unlab.WAD.order)
    # lines(x=curr.data.orders$R1.12C_N.GL.corr[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=2, col="cyan")
    # lines(x=curr.data.orders$R2.12C_N.GL.corr[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=2, col="magenta")
    # lines(x=curr.data.orders$R3.12C_N.GL.corr[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=2, col="wheat")
  #18.O.GL
    points(x=curr.data.orders$R1.18O.GL.corr, y=curr.data.orders$mean.unlab.WAD.order, pch=21, cex=0.7, col="black", bg="gray")
    points(x=curr.data.orders$R2.18O.GL.corr, y=curr.data.orders$mean.unlab.WAD.order, pch=21, cex=0.7, col="black", bg="brown")
    points(x=curr.data.orders$R3.18O.GL.corr, y=curr.data.orders$mean.unlab.WAD.order, pch=21, cex=0.7, col="black", bg="darkgoldenrod2")
    # new.order <- order(curr.data.orders$mean.unlab.WAD.order)
    # lines(x=curr.data.orders$R1.18O.GL.corr[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=2, col="gray")
    # lines(x=curr.data.orders$R2.18O.GL.corr[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=2, col="brown")
    # lines(x=curr.data.orders$R3.18O.GL.corr[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=2, col="darkgoldenrod2")
  #Unlabeled average:
    # points(x=apply(unlab.WADs, 1, mean, na.rm=TRUE), y=curr.data.orders$mean.unlab.WAD.order, pch=21, cex=0.35, col="black", bg="black")
    new.order <- order(curr.data.orders$mean.unlab.WAD.order)
    lines(x=apply(unlab.WADs, 1, mean, na.rm=TRUE)[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=3, col="black")
    lines(x=apply(unlab.WADs, 1, mean, na.rm=TRUE)[new.order], y=curr.data.orders$mean.unlab.WAD.order[new.order], lwd=2, col="white")
  #Normal cdf fit to the mean of all unlabeled data:
    x <- seq(min(unlab.WADs), max(unlab.WADs), 0.000001)
    y <- pnorm(q=x, mean=WAD.norm.fit.parms$mean[WAD.norm.fit.parms$trt == "unlabeled"], sd=WAD.norm.fit.parms$stdev[WAD.norm.fit.parms$trt == "unlabeled"])*max(curr.data.orders$mean.unlab.WAD.order)
    lines(x, y, lwd=3, col="black")
    lines(x, y, lwd=2, col="lawngreen")
    par(xpd=NA)
    ink <- c("green", "red", "yellow", "purple", "blue", "orange", "cyan", "magenta", "wheat", "gray", "brown", "darkgoldenrod2", "white", "lawngreen")
    legend(x=x.max*1.0040, y=max(as.numeric(as.character(curr.data.orders$mean.unlab.WAD.order)))+(0.04*(max(as.numeric(as.character(curr.data.orders$mean.unlab.WAD.order)))-1)), legend=c("R1.16O.GL", "R2.16O.GL", "R3.16O.GL", "R1.12C.GL", "R2.12C.GL", "R3.12C.GL", "R1.12C_N.GL", "R2.12C_N.GL", "R3.12C_N.GL", "R1.18O.GL", "R2.18O.GL", "R3.18O.GL", "unlab empirical mean", "unlabeled normal fit"), col=c(ink[1:12], "black", "black"), pt.bg=ink, pch=rep(21,12), pt.cex=0.7, lwd=c(rep(1, 12), 3, 3), bty="o", cex=0.6, xjust=0, yjust=1)
    legend(x=x.max*1.0040, y=max(as.numeric(as.character(curr.data.orders$mean.unlab.WAD.order)))+(0.04*(max(as.numeric(as.character(curr.data.orders$mean.unlab.WAD.order)))-1)), legend=c("R1.16O.GL", "R2.16O.GL", "R3.16O.GL", "R1.12C.GL", "R2.12C.GL", "R3.12C.GL", "R1.12C_N.GL", "R2.12C_N.GL", "R3.12C_N.GL", "R1.18O.GL", "R2.18O.GL", "R3.18O.GL", "unlab empirical mean", "unlabeled normal fit"), col=ink, pt.bg=ink, pch=rep(21,12), pt.cex=0.7, lwd=c(rep(1, 12), 2, 2), bty="o", cex=0.6, xjust=0, yjust=1)
    par(xpd=FALSE)
    rm(obs.wads.by.taxon, WAD.norm.fit.parms, curr.data.orders, unlab.reps, unlab.WADs, ink)
#####}
dev.off()   #(run graph code through this line only if creating a pdf version of the plot)






##### Now perform the tube-level WAD corrections for ALL GRASSLAND replicates of ALL treatments (not just unlabeled and 18O), and then proceed with standard qSIP analysis:


#Calculate the shift in WAD for the specified unlabeled replicates from their global mean using the taxa common to ALL replicates of ALL treatments:
  #Identify the appropriate shift for unlabeled WADs:
  #Note: this function takes a few minutes to run (i.e., ~4min for ~250 taxa common to 3 unlabeled and 5 labeled treatment):
    #Include all four unlabeled treatments and all five labeled treatments from the Grassland ecosystem:
      system.time(unlab.WAD.corr.GL.list <- find.unlabeled.correction(LIST=WAD.by.taxon.GL, unlab.tmts=c("16O.GL", "12C.GL", "12C_N.GL"), lab.tmts=c("12C_18O_N.GL", "12C_18O.GL", "13C_N.GL", "13C.GL", "18O.GL"), CI=0.90))
      #Look at the results:
        names(unlab.WAD.corr.GL.list)                           #names of the two data frames in the output list
        unlab.WAD.corr.GL.list$WAD.norm.fit.parms               #looking at the first object -- a data frame
        unlab.WAD.corr.GL.list[[1]]                             #another way to look at the first object
        unlab.WAD.corr.GL.list$corr.names                       #looking at the second object in the list -- a vector of the names of the corrected unlabeled replicates
        unlab.WAD.corr.GL.list[[2]]                             #another way to look at the second object
        head(unlab.WAD.corr.GL.list$WAD.table.corr)             #looking at the head of the third object -- a data frame
        head(unlab.WAD.corr.GL.list[[3]])                       #another way to look at the head of the third object

      #Note that the unlabeled shifts are slightly different when all treatments are included compared to when only the 18O treatment is included:
        unlab.WAD.corr.18O.GL.list$WAD.norm.fit.parms
        unlab.WAD.corr.GL.list$WAD.norm.fit.parms


#Calculate the shift in WAD for each of the labeled replicates (5 labeled treatments x 3 replicates = 15 total labeled replicates in the Grassland ecosystem):
  #18O.GL:
    system.time(R1.18O.GL.td.pos.resid.list <- find.labeled.correction(LIST=unlab.WAD.corr.GL.list, lab.replicate="R1.18O.GL", lab.names=c("R1.18O.GL", "R2.18O.GL", "R3.18O.GL"), method="td.pos.resid", unlab.SD.percentile=0.5, lab.SD.percentile=0.5, min.num.nongrowers=10))
      #Write the find.labeled.correction output tables to text files:
      write.table(R1.18O.GL.td.pos.resid.list$putative.nongrower.metrics, "qSIP_output/RM_DIM_R1_18O_GL_td_pos_resid_putative_nongrower_metrics.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R1.18O.GL.td.pos.resid.list$taxa.in, "qSIP_output/RM_DIM_R1_18O_GL_td_pos_resid_taxa_in.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R1.18O.GL.td.pos.resid.list$taxa.out, "qSIP_output/RM_DIM_R1_18O_GL_td_pos_resid_taxa_out.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R1.18O.GL.td.pos.resid.list$data.low.SD, "qSIP_output/RM_DIM_R1_18O_GL_td_pos_resid_data_low_SD.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      #Find the best iteration:
      R1.18O.GL.td.pos.resid.best.iteration.list <- select.best.iteration(putative.nongrower.metrics=R1.18O.GL.td.pos.resid.list$putative.nongrower.metrics, quantile.threshold=0.10)
      find.labeled.correction.plot(find.labeled.correction.list=R1.18O.GL.td.pos.resid.list, filename="RM_DIM_GL_R1_18O_GL_td_pos_resid_plots", path="qSIP_output/Figures/", highlight=R1.18O.GL.td.pos.resid.best.iteration.list$best.iteration, highlight.col="magenta")
    system.time(R2.18O.GL.td.pos.resid.list <- find.labeled.correction(LIST=unlab.WAD.corr.GL.list, lab.replicate="R2.18O.GL", lab.names=c("R1.18O.GL", "R2.18O.GL", "R3.18O.GL"), method="td.pos.resid", unlab.SD.percentile=0.5, lab.SD.percentile=0.5, min.num.nongrowers=10))
      #Write the find.labeled.correction output tables to text files:
      write.table(R2.18O.GL.td.pos.resid.list$putative.nongrower.metrics, "qSIP_output/RM_DIM_R2_18O_GL_td_pos_resid_putative_nongrower_metrics.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R2.18O.GL.td.pos.resid.list$taxa.in, "qSIP_output/RM_DIM_R2_18O_GL_td_pos_resid_taxa_in.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R2.18O.GL.td.pos.resid.list$taxa.out, "qSIP_output/RM_DIM_R2_18O_GL_td_pos_resid_taxa_out.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R2.18O.GL.td.pos.resid.list$data.low.SD, "qSIP_output/RM_DIM_R2_18O_GL_td_pos_resid_data_low_SD.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      #Find the best iteration:
      R2.18O.GL.td.pos.resid.best.iteration.list <- select.best.iteration(putative.nongrower.metrics=R2.18O.GL.td.pos.resid.list$putative.nongrower.metrics, quantile.threshold=0.10)
      find.labeled.correction.plot(find.labeled.correction.list=R2.18O.GL.td.pos.resid.list, filename="RM_DIM_GL_R2_18O_GL_td_pos_resid_plots", path="qSIP_output/Figures/", highlight=R2.18O.GL.td.pos.resid.best.iteration.list$best.iteration, highlight.col="magenta")
    system.time(R3.18O.GL.td.pos.resid.list <- find.labeled.correction(LIST=unlab.WAD.corr.GL.list, lab.replicate="R3.18O.GL", lab.names=c("R1.18O.GL", "R2.18O.GL", "R3.18O.GL"), method="td.pos.resid", unlab.SD.percentile=0.5, lab.SD.percentile=0.5, min.num.nongrowers=10))
      #Write the find.labeled.correction output tables to text files:
      write.table(R3.18O.GL.td.pos.resid.list$putative.nongrower.metrics, "qSIP_output/RM_DIM_R3_18O_GL_td_pos_resid_putative_nongrower_metrics.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R3.18O.GL.td.pos.resid.list$taxa.in, "qSIP_output/RM_DIM_R3_18O_GL_td_pos_resid_taxa_in.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R3.18O.GL.td.pos.resid.list$taxa.out, "qSIP_output/RM_DIM_R3_18O_GL_td_pos_resid_taxa_out.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R3.18O.GL.td.pos.resid.list$data.low.SD, "qSIP_output/RM_DIM_R3_18O_GL_td_pos_resid_data_low_SD.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      #Find the best iteration:
      R3.18O.GL.td.pos.resid.best.iteration.list <- select.best.iteration(putative.nongrower.metrics=R3.18O.GL.td.pos.resid.list$putative.nongrower.metrics, quantile.threshold=0.10)
      find.labeled.correction.plot(find.labeled.correction.list=R3.18O.GL.td.pos.resid.list, filename="RM_DIM_GL_R3_18O_GL_td_pos_resid_plots", path="qSIP_output/Figures/", highlight=R3.18O.GL.td.pos.resid.best.iteration.list$best.iteration, highlight.col="magenta")
  #13C.GL
    system.time(R1.13C.GL.td.pos.resid.list <- find.labeled.correction(LIST=unlab.WAD.corr.GL.list, lab.replicate="R1.13C.GL", lab.names=c("R1.13C.GL", "R2.13C.GL", "R3.13C.GL"), method="td.pos.resid", unlab.SD.percentile=0.5, lab.SD.percentile=0.5, min.num.nongrowers=10))
      #Write the find.labeled.correction output tables to text files:
      write.table(R1.13C.GL.td.pos.resid.list$putative.nongrower.metrics, "qSIP_output/RM_DIM_R1_13C_GL_td_pos_resid_putative_nongrower_metrics.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R1.13C.GL.td.pos.resid.list$taxa.in, "qSIP_output/RM_DIM_R1_13C_GL_td_pos_resid_taxa_in.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R1.13C.GL.td.pos.resid.list$taxa.out, "qSIP_output/RM_DIM_R1_13C_GL_td_pos_resid_taxa_out.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R1.13C.GL.td.pos.resid.list$data.low.SD, "qSIP_output/RM_DIM_R1_13C_GL_td_pos_resid_data_low_SD.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      #Find the best iteration:
      R1.13C.GL.td.pos.resid.best.iteration.list <- select.best.iteration(putative.nongrower.metrics=R1.13C.GL.td.pos.resid.list$putative.nongrower.metrics, quantile.threshold=0.10)
      find.labeled.correction.plot(find.labeled.correction.list=R1.13C.GL.td.pos.resid.list, filename="RM_DIM_GL_R1_13C_GL_td_pos_resid_plots", path="qSIP_output/Figures/", highlight=R1.13C.GL.td.pos.resid.best.iteration.list$best.iteration, highlight.col="magenta")
    system.time(R2.13C.GL.td.pos.resid.list <- find.labeled.correction(LIST=unlab.WAD.corr.GL.list, lab.replicate="R2.13C.GL", lab.names=c("R1.13C.GL", "R2.13C.GL", "R3.13C.GL"), method="td.pos.resid", unlab.SD.percentile=0.5, lab.SD.percentile=0.5, min.num.nongrowers=10))
      #Write the find.labeled.correction output tables to text files:
      write.table(R2.13C.GL.td.pos.resid.list$putative.nongrower.metrics, "qSIP_output/RM_DIM_R2_13C_GL_td_pos_resid_putative_nongrower_metrics.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R2.13C.GL.td.pos.resid.list$taxa.in, "qSIP_output/RM_DIM_R2_13C_GL_td_pos_resid_taxa_in.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R2.13C.GL.td.pos.resid.list$taxa.out, "qSIP_output/RM_DIM_R2_13C_GL_td_pos_resid_taxa_out.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R2.13C.GL.td.pos.resid.list$data.low.SD, "qSIP_output/RM_DIM_R2_13C_GL_td_pos_resid_data_low_SD.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      #Find the best iteration:
      R2.13C.GL.td.pos.resid.best.iteration.list <- select.best.iteration(putative.nongrower.metrics=R2.13C.GL.td.pos.resid.list$putative.nongrower.metrics, quantile.threshold=0.10)
      find.labeled.correction.plot(find.labeled.correction.list=R2.13C.GL.td.pos.resid.list, filename="RM_DIM_GL_R2_13C_GL_td_pos_resid_plots", path="qSIP_output/Figures/", highlight=R2.13C.GL.td.pos.resid.best.iteration.list$best.iteration, highlight.col="magenta")
    system.time(R3.13C.GL.td.pos.resid.list <- find.labeled.correction(LIST=unlab.WAD.corr.GL.list, lab.replicate="R3.13C.GL", lab.names=c("R1.13C.GL", "R2.13C.GL", "R3.13C.GL"), method="td.pos.resid", unlab.SD.percentile=0.5, lab.SD.percentile=0.5, min.num.nongrowers=10))
      #Write the find.labeled.correction output tables to text files:
      write.table(R3.13C.GL.td.pos.resid.list$putative.nongrower.metrics, "qSIP_output/RM_DIM_R3_13C_GL_td_pos_resid_putative_nongrower_metrics.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R3.13C.GL.td.pos.resid.list$taxa.in, "qSIP_output/RM_DIM_R3_13C_GL_td_pos_resid_taxa_in.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R3.13C.GL.td.pos.resid.list$taxa.out, "qSIP_output/RM_DIM_R3_13C_GL_td_pos_resid_taxa_out.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R3.13C.GL.td.pos.resid.list$data.low.SD, "qSIP_output/RM_DIM_R3_13C_GL_td_pos_resid_data_low_SD.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      #Find the best iteration:
      R3.13C.GL.td.pos.resid.best.iteration.list <- select.best.iteration(putative.nongrower.metrics=R3.13C.GL.td.pos.resid.list$putative.nongrower.metrics, quantile.threshold=0.10)
      find.labeled.correction.plot(find.labeled.correction.list=R3.13C.GL.td.pos.resid.list, filename="RM_DIM_GL_R3_13C_GL_td_pos_resid_plots", path="qSIP_output/Figures/", highlight=R3.13C.GL.td.pos.resid.best.iteration.list$best.iteration, highlight.col="magenta")
  #13C_N.GL:
    system.time(R1.13C_N.GL.td.pos.resid.list <- find.labeled.correction(LIST=unlab.WAD.corr.GL.list, lab.replicate="R1.13C_N.GL", lab.names=c("R1.13C_N.GL", "R2.13C_N.GL", "R3.13C_N.GL"), method="td.pos.resid", unlab.SD.percentile=0.5, lab.SD.percentile=0.5, min.num.nongrowers=10))
      #Write the find.labeled.correction output tables to text files:
      write.table(R1.13C_N.GL.td.pos.resid.list$putative.nongrower.metrics, "qSIP_output/RM_DIM_R1_13C_N_GL_td_pos_resid_putative_nongrower_metrics.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R1.13C_N.GL.td.pos.resid.list$taxa.in, "qSIP_output/RM_DIM_R1_13C_N_GL_td_pos_resid_taxa_in.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R1.13C_N.GL.td.pos.resid.list$taxa.out, "qSIP_output/RM_DIM_R1_13C_N_GL_td_pos_resid_taxa_out.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R1.13C_N.GL.td.pos.resid.list$data.low.SD, "qSIP_output/RM_DIM_R1_13C_N_GL_td_pos_resid_data_low_SD.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      #Find the best iteration:
      R1.13C_N.GL.td.pos.resid.best.iteration.list <- select.best.iteration(putative.nongrower.metrics=R1.13C_N.GL.td.pos.resid.list$putative.nongrower.metrics, quantile.threshold=0.10)
      find.labeled.correction.plot(find.labeled.correction.list=R1.13C_N.GL.td.pos.resid.list, filename="RM_DIM_GL_R1_13C_N_GL_td_pos_resid_plots", path="qSIP_output/Figures/", highlight=R1.13C_N.GL.td.pos.resid.best.iteration.list$best.iteration, highlight.col="magenta")
    system.time(R2.13C_N.GL.td.pos.resid.list <- find.labeled.correction(LIST=unlab.WAD.corr.GL.list, lab.replicate="R2.13C_N.GL", lab.names=c("R1.13C_N.GL", "R2.13C_N.GL", "R3.13C_N.GL"), method="td.pos.resid", unlab.SD.percentile=0.5, lab.SD.percentile=0.5, min.num.nongrowers=10))
      #Write the find.labeled.correction output tables to text files:
      write.table(R2.13C_N.GL.td.pos.resid.list$putative.nongrower.metrics, "qSIP_output/RM_DIM_R2_13C_N_GL_td_pos_resid_putative_nongrower_metrics.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R2.13C_N.GL.td.pos.resid.list$taxa.in, "qSIP_output/RM_DIM_R2_13C_N_GL_td_pos_resid_taxa_in.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R2.13C_N.GL.td.pos.resid.list$taxa.out, "qSIP_output/RM_DIM_R2_13C_N_GL_td_pos_resid_taxa_out.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R2.13C_N.GL.td.pos.resid.list$data.low.SD, "qSIP_output/RM_DIM_R2_13C_N_GL_td_pos_resid_data_low_SD.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      #Find the best iteration:
      R2.13C_N.GL.td.pos.resid.best.iteration.list <- select.best.iteration(putative.nongrower.metrics=R2.13C_N.GL.td.pos.resid.list$putative.nongrower.metrics, quantile.threshold=0.10)
      find.labeled.correction.plot(find.labeled.correction.list=R2.13C_N.GL.td.pos.resid.list, filename="RM_DIM_GL_R2_13C_N_GL_td_pos_resid_plots", path="qSIP_output/Figures/", highlight=R2.13C_N.GL.td.pos.resid.best.iteration.list$best.iteration, highlight.col="magenta")
    system.time(R3.13C_N.GL.td.pos.resid.list <- find.labeled.correction(LIST=unlab.WAD.corr.GL.list, lab.replicate="R3.13C_N.GL", lab.names=c("R1.13C_N.GL", "R2.13C_N.GL", "R3.13C_N.GL"), method="td.pos.resid", unlab.SD.percentile=0.5, lab.SD.percentile=0.5, min.num.nongrowers=10))
      #Write the find.labeled.correction output tables to text files:
      write.table(R3.13C_N.GL.td.pos.resid.list$putative.nongrower.metrics, "qSIP_output/RM_DIM_R3_13C_N_GL_td_pos_resid_putative_nongrower_metrics.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R3.13C_N.GL.td.pos.resid.list$taxa.in, "qSIP_output/RM_DIM_R3_13C_N_GL_td_pos_resid_taxa_in.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R3.13C_N.GL.td.pos.resid.list$taxa.out, "qSIP_output/RM_DIM_R3_13C_N_GL_td_pos_resid_taxa_out.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R3.13C_N.GL.td.pos.resid.list$data.low.SD, "qSIP_output/RM_DIM_R3_13C_N_GL_td_pos_resid_data_low_SD.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      #Find the best iteration:
      R3.13C_N.GL.td.pos.resid.best.iteration.list <- select.best.iteration(putative.nongrower.metrics=R3.13C_N.GL.td.pos.resid.list$putative.nongrower.metrics, quantile.threshold=0.10)
      find.labeled.correction.plot(find.labeled.correction.list=R3.13C_N.GL.td.pos.resid.list, filename="RM_DIM_GL_R3_13C_N_GL_td_pos_resid_plots", path="qSIP_output/Figures/", highlight=R3.13C_N.GL.td.pos.resid.best.iteration.list$best.iteration, highlight.col="magenta")
  #12C_18O.GL:
    system.time(R1.12C_18O.GL.td.pos.resid.list <- find.labeled.correction(LIST=unlab.WAD.corr.GL.list, lab.replicate="R1.12C_18O.GL", lab.names=c("R1.12C_18O.GL", "R2.12C_18O.GL", "R3.12C_18O.GL"), method="td.pos.resid", unlab.SD.percentile=0.5, lab.SD.percentile=0.5, min.num.nongrowers=10))
      #Write the find.labeled.correction output tables to text files:
      write.table(R1.12C_18O.GL.td.pos.resid.list$putative.nongrower.metrics, "qSIP_output/RM_DIM_R1_12C_18O_GL_td_pos_resid_putative_nongrower_metrics.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R1.12C_18O.GL.td.pos.resid.list$taxa.in, "qSIP_output/RM_DIM_R1_12C_18O_GL_td_pos_resid_taxa_in.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R1.12C_18O.GL.td.pos.resid.list$taxa.out, "qSIP_output/RM_DIM_R1_12C_18O_GL_td_pos_resid_taxa_out.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R1.12C_18O.GL.td.pos.resid.list$data.low.SD, "qSIP_output/RM_DIM_R1_12C_18O_GL_td_pos_resid_data_low_SD.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      #Find the best iteration:
      R1.12C_18O.GL.td.pos.resid.best.iteration.list <- select.best.iteration(putative.nongrower.metrics=R1.12C_18O.GL.td.pos.resid.list$putative.nongrower.metrics, quantile.threshold=0.10)
      find.labeled.correction.plot(find.labeled.correction.list=R1.12C_18O.GL.td.pos.resid.list, filename="RM_DIM_GL_R1_12C_18O_GL_td_pos_resid_plots", path="qSIP_output/Figures/", highlight=R1.12C_18O.GL.td.pos.resid.best.iteration.list$best.iteration, highlight.col="magenta")
    system.time(R2.12C_18O.GL.td.pos.resid.list <- find.labeled.correction(LIST=unlab.WAD.corr.GL.list, lab.replicate="R2.12C_18O.GL", lab.names=c("R1.12C_18O.GL", "R2.12C_18O.GL", "R3.12C_18O.GL"), method="td.pos.resid", unlab.SD.percentile=0.5, lab.SD.percentile=0.5, min.num.nongrowers=10))
      #Write the find.labeled.correction output tables to text files:
      write.table(R2.12C_18O.GL.td.pos.resid.list$putative.nongrower.metrics, "qSIP_output/RM_DIM_R2_12C_18O_GL_td_pos_resid_putative_nongrower_metrics.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R2.12C_18O.GL.td.pos.resid.list$taxa.in, "qSIP_output/RM_DIM_R2_12C_18O_GL_td_pos_resid_taxa_in.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R2.12C_18O.GL.td.pos.resid.list$taxa.out, "qSIP_output/RM_DIM_R2_12C_18O_GL_td_pos_resid_taxa_out.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R2.12C_18O.GL.td.pos.resid.list$data.low.SD, "qSIP_output/RM_DIM_R2_12C_18O_GL_td_pos_resid_data_low_SD.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      #Find the best iteration:
      R2.12C_18O.GL.td.pos.resid.best.iteration.list <- select.best.iteration(putative.nongrower.metrics=R2.12C_18O.GL.td.pos.resid.list$putative.nongrower.metrics, quantile.threshold=0.10)
      find.labeled.correction.plot(find.labeled.correction.list=R2.12C_18O.GL.td.pos.resid.list, filename="RM_DIM_GL_R2_12C_18O_GL_td_pos_resid_plots", path="qSIP_output/Figures/", highlight=R2.12C_18O.GL.td.pos.resid.best.iteration.list$best.iteration, highlight.col="magenta")
    system.time(R3.12C_18O.GL.td.pos.resid.list <- find.labeled.correction(LIST=unlab.WAD.corr.GL.list, lab.replicate="R3.12C_18O.GL", lab.names=c("R1.12C_18O.GL", "R2.12C_18O.GL", "R3.12C_18O.GL"), method="td.pos.resid", unlab.SD.percentile=0.5, lab.SD.percentile=0.5, min.num.nongrowers=10))
      #Write the find.labeled.correction output tables to text files:
      write.table(R3.12C_18O.GL.td.pos.resid.list$putative.nongrower.metrics, "qSIP_output/RM_DIM_R3_12C_18O_GL_td_pos_resid_putative_nongrower_metrics.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R3.12C_18O.GL.td.pos.resid.list$taxa.in, "qSIP_output/RM_DIM_R3_12C_18O_GL_td_pos_resid_taxa_in.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R3.12C_18O.GL.td.pos.resid.list$taxa.out, "qSIP_output/RM_DIM_R3_12C_18O_GL_td_pos_resid_taxa_out.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R3.12C_18O.GL.td.pos.resid.list$data.low.SD, "qSIP_output/RM_DIM_R3_12C_18O_GL_td_pos_resid_data_low_SD.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      #Find the best iteration:
      R3.12C_18O.GL.td.pos.resid.best.iteration.list <- select.best.iteration(putative.nongrower.metrics=R3.12C_18O.GL.td.pos.resid.list$putative.nongrower.metrics, quantile.threshold=0.10)
      find.labeled.correction.plot(find.labeled.correction.list=R3.12C_18O.GL.td.pos.resid.list, filename="RM_DIM_GL_R3_12C_18O_GL_td_pos_resid_plots", path="qSIP_output/Figures/", highlight=R3.12C_18O.GL.td.pos.resid.best.iteration.list$best.iteration, highlight.col="magenta")
  #12C_18O_N.GL:
    system.time(R1.12C_18O_N.GL.td.pos.resid.list <- find.labeled.correction(LIST=unlab.WAD.corr.GL.list, lab.replicate="R1.12C_18O_N.GL", lab.names=c("R1.12C_18O_N.GL", "R2.12C_18O_N.GL", "R3.12C_18O_N.GL"), method="td.pos.resid", unlab.SD.percentile=0.5, lab.SD.percentile=0.5, min.num.nongrowers=10))
      #Write the find.labeled.correction output tables to text files:
      write.table(R1.12C_18O_N.GL.td.pos.resid.list$putative.nongrower.metrics, "qSIP_output/RM_DIM_R1_12C_18O_N_GL_td_pos_resid_putative_nongrower_metrics.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R1.12C_18O_N.GL.td.pos.resid.list$taxa.in, "qSIP_output/RM_DIM_R1_12C_18O_N_GL_td_pos_resid_taxa_in.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R1.12C_18O_N.GL.td.pos.resid.list$taxa.out, "qSIP_output/RM_DIM_R1_12C_18O_N_GL_td_pos_resid_taxa_out.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R1.12C_18O_N.GL.td.pos.resid.list$data.low.SD, "qSIP_output/RM_DIM_R1_12C_18O_N_GL_td_pos_resid_data_low_SD.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      #Find the best iteration:
      R1.12C_18O_N.GL.td.pos.resid.best.iteration.list <- select.best.iteration(putative.nongrower.metrics=R1.12C_18O_N.GL.td.pos.resid.list$putative.nongrower.metrics, quantile.threshold=0.10)
      find.labeled.correction.plot(find.labeled.correction.list=R1.12C_18O_N.GL.td.pos.resid.list, filename="RM_DIM_GL_R1_12C_18O_N_GL_td_pos_resid_plots", path="qSIP_output/Figures/", highlight=R1.12C_18O_N.GL.td.pos.resid.best.iteration.list$best.iteration, highlight.col="magenta")
    system.time(R2.12C_18O_N.GL.td.pos.resid.list <- find.labeled.correction(LIST=unlab.WAD.corr.GL.list, lab.replicate="R2.12C_18O_N.GL", lab.names=c("R1.12C_18O_N.GL", "R2.12C_18O_N.GL", "R3.12C_18O_N.GL"), method="td.pos.resid", unlab.SD.percentile=0.5, lab.SD.percentile=0.5, min.num.nongrowers=10))
      #Write the find.labeled.correction output tables to text files:
      write.table(R2.12C_18O_N.GL.td.pos.resid.list$putative.nongrower.metrics, "qSIP_output/RM_DIM_R2_12C_18O_N_GL_td_pos_resid_putative_nongrower_metrics.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R2.12C_18O_N.GL.td.pos.resid.list$taxa.in, "qSIP_output/RM_DIM_R2_12C_18O_N_GL_td_pos_resid_taxa_in.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R2.12C_18O_N.GL.td.pos.resid.list$taxa.out, "qSIP_output/RM_DIM_R2_12C_18O_N_GL_td_pos_resid_taxa_out.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R2.12C_18O_N.GL.td.pos.resid.list$data.low.SD, "qSIP_output/RM_DIM_R2_12C_18O_N_GL_td_pos_resid_data_low_SD.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      #Find the best iteration:
      R2.12C_18O_N.GL.td.pos.resid.best.iteration.list <- select.best.iteration(putative.nongrower.metrics=R2.12C_18O_N.GL.td.pos.resid.list$putative.nongrower.metrics, quantile.threshold=0.10)
      find.labeled.correction.plot(find.labeled.correction.list=R2.12C_18O_N.GL.td.pos.resid.list, filename="RM_DIM_GL_R2_12C_18O_N_GL_td_pos_resid_plots", path="qSIP_output/Figures/", highlight=R2.12C_18O_N.GL.td.pos.resid.best.iteration.list$best.iteration, highlight.col="magenta")
    system.time(R3.12C_18O_N.GL.td.pos.resid.list <- find.labeled.correction(LIST=unlab.WAD.corr.GL.list, lab.replicate="R3.12C_18O_N.GL", lab.names=c("R1.12C_18O_N.GL", "R2.12C_18O_N.GL", "R3.12C_18O_N.GL"), method="td.pos.resid", unlab.SD.percentile=0.5, lab.SD.percentile=0.5, min.num.nongrowers=10))
      #Write the find.labeled.correction output tables to text files:
      write.table(R3.12C_18O_N.GL.td.pos.resid.list$putative.nongrower.metrics, "qSIP_output/RM_DIM_R3_12C_18O_N_GL_td_pos_resid_putative_nongrower_metrics.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R3.12C_18O_N.GL.td.pos.resid.list$taxa.in, "qSIP_output/RM_DIM_R3_12C_18O_N_GL_td_pos_resid_taxa_in.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R3.12C_18O_N.GL.td.pos.resid.list$taxa.out, "qSIP_output/RM_DIM_R3_12C_18O_N_GL_td_pos_resid_taxa_out.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      write.table(R3.12C_18O_N.GL.td.pos.resid.list$data.low.SD, "qSIP_output/RM_DIM_R3_12C_18O_N_GL_td_pos_resid_data_low_SD.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
      #Find the best iteration:
      R3.12C_18O_N.GL.td.pos.resid.best.iteration.list <- select.best.iteration(putative.nongrower.metrics=R3.12C_18O_N.GL.td.pos.resid.list$putative.nongrower.metrics, quantile.threshold=0.10)
      find.labeled.correction.plot(find.labeled.correction.list=R3.12C_18O_N.GL.td.pos.resid.list, filename="RM_DIM_GL_R3_12C_18O_N_GL_td_pos_resid_plots", path="qSIP_output/Figures/", highlight=R3.12C_18O_N.GL.td.pos.resid.best.iteration.list$best.iteration, highlight.col="magenta")


#Correct the summarized and raw data for tube-level shifts in WAD identified above for ALL unlabeled and labeled replicates in the Grassland ecosystem:
  #Apply the tube-level shift to the already-summarized taxon-level labeled replicate WADs to 'correct' them and create a new summary list (analagous to 'unlab.WAD.corr.GL.list') that includes those corrected values:
    WAD.corr.GL.list <- add.lab.WAD.corr.summary(summary.list=unlab.WAD.corr.GL.list, reps.by.trt=WAD.by.taxon.GL$reps.by.trt, lab.reps.to.add=c("R1.18O.GL", "R2.18O.GL", "R3.18O.GL", "R1.13C.GL", "R2.13C.GL", "R3.13C.GL", "R1.13C_N.GL", "R2.13C_N.GL", "R3.13C_N.GL", "R1.12C_18O.GL", "R2.12C_18O.GL", "R3.12C_18O.GL", "R1.12C_18O_N.GL", "R2.12C_18O_N.GL", "R3.12C_18O_N.GL"), lab.shifts=c(R1.18O.GL.td.pos.resid.best.iteration.list$norm.correction, R2.18O.GL.td.pos.resid.best.iteration.list$norm.correction, R3.18O.GL.td.pos.resid.best.iteration.list$norm.correction, R1.13C.GL.td.pos.resid.best.iteration.list$norm.correction, R2.13C.GL.td.pos.resid.best.iteration.list$norm.correction, R3.13C.GL.td.pos.resid.best.iteration.list$norm.correction, R1.13C_N.GL.td.pos.resid.best.iteration.list$norm.correction, R2.13C_N.GL.td.pos.resid.best.iteration.list$norm.correction, R3.13C_N.GL.td.pos.resid.best.iteration.list$norm.correction, R1.12C_18O.GL.td.pos.resid.best.iteration.list$norm.correction, R2.12C_18O.GL.td.pos.resid.best.iteration.list$norm.correction, R3.12C_18O.GL.td.pos.resid.best.iteration.list$norm.correction, R1.12C_18O_N.GL.td.pos.resid.best.iteration.list$norm.correction, R2.12C_18O_N.GL.td.pos.resid.best.iteration.list$norm.correction, R3.12C_18O_N.GL.td.pos.resid.best.iteration.list$norm.correction), CI=0.90)
    #Look at the results:
      names(WAD.corr.GL.list)                           #names of the two data frames in the output list
      WAD.corr.GL.list$WAD.norm.fit.parms               #looking at the first object -- a data frame
      WAD.corr.GL.list[[1]]                             #another way to look at the first object
      WAD.corr.GL.list$corr.names                       #looking at the second object in the list -- a vector of the names of the corrected unlabeled replicates
      WAD.corr.GL.list[[2]]                             #another way to look at the second object
      head(WAD.corr.GL.list$WAD.table.corr)             #looking at the head of the third object -- a data frame
      head(WAD.corr.GL.list[[3]])                       #another way to look at the head of the third object

  #Apply the tube-level shifts to 'correct' the unlabeled WADs in 'data' (corrects all unlabeled replicates at once):
    data.corr <- apply.unlabeled.correction(raw.data=data, correction.table=unlab.WAD.corr.GL.list$WAD.norm.fit.parms, reps.by.trt=WAD.by.taxon.GL$reps.by.trt, vars=c("Density.g.ml", "Sample", "iso.treat.eco"))
  #Apply the tube-level shift to 'correct' the labeled WADs in 'data.corr' (corrects one labeled replicate at a time):
    data.corr <- apply.labeled.correction(raw.data=data, raw.data.corr=data.corr, lab.replicate="R1.18O.GL", correction.value=R1.18O.GL.td.pos.resid.best.iteration.list$norm.correction, reps.by.trt=WAD.by.taxon.GL$reps.by.trt, vars=c("Density.g.ml", "Sample", "iso.treat.eco"))
    data.corr <- apply.labeled.correction(raw.data=data, raw.data.corr=data.corr, lab.replicate="R2.18O.GL", correction.value=R2.18O.GL.td.pos.resid.best.iteration.list$norm.correction, reps.by.trt=WAD.by.taxon.GL$reps.by.trt, vars=c("Density.g.ml", "Sample", "iso.treat.eco"))
    data.corr <- apply.labeled.correction(raw.data=data, raw.data.corr=data.corr, lab.replicate="R3.18O.GL", correction.value=R3.18O.GL.td.pos.resid.best.iteration.list$norm.correction, reps.by.trt=WAD.by.taxon.GL$reps.by.trt, vars=c("Density.g.ml", "Sample", "iso.treat.eco"))
    data.corr <- apply.labeled.correction(raw.data=data, raw.data.corr=data.corr, lab.replicate="R1.13C.GL", correction.value=R1.13C.GL.td.pos.resid.best.iteration.list$norm.correction, reps.by.trt=WAD.by.taxon.GL$reps.by.trt, vars=c("Density.g.ml", "Sample", "iso.treat.eco"))
    data.corr <- apply.labeled.correction(raw.data=data, raw.data.corr=data.corr, lab.replicate="R2.13C.GL", correction.value=R2.13C.GL.td.pos.resid.best.iteration.list$norm.correction, reps.by.trt=WAD.by.taxon.GL$reps.by.trt, vars=c("Density.g.ml", "Sample", "iso.treat.eco"))
    data.corr <- apply.labeled.correction(raw.data=data, raw.data.corr=data.corr, lab.replicate="R3.13C.GL", correction.value=R3.13C.GL.td.pos.resid.best.iteration.list$norm.correction, reps.by.trt=WAD.by.taxon.GL$reps.by.trt, vars=c("Density.g.ml", "Sample", "iso.treat.eco"))
    data.corr <- apply.labeled.correction(raw.data=data, raw.data.corr=data.corr, lab.replicate="R1.13C_N.GL", correction.value=R1.13C_N.GL.td.pos.resid.best.iteration.list$norm.correction, reps.by.trt=WAD.by.taxon.GL$reps.by.trt, vars=c("Density.g.ml", "Sample", "iso.treat.eco"))
    data.corr <- apply.labeled.correction(raw.data=data, raw.data.corr=data.corr, lab.replicate="R2.13C_N.GL", correction.value=R2.13C_N.GL.td.pos.resid.best.iteration.list$norm.correction, reps.by.trt=WAD.by.taxon.GL$reps.by.trt, vars=c("Density.g.ml", "Sample", "iso.treat.eco"))
    data.corr <- apply.labeled.correction(raw.data=data, raw.data.corr=data.corr, lab.replicate="R3.13C_N.GL", correction.value=R3.13C_N.GL.td.pos.resid.best.iteration.list$norm.correction, reps.by.trt=WAD.by.taxon.GL$reps.by.trt, vars=c("Density.g.ml", "Sample", "iso.treat.eco"))
    data.corr <- apply.labeled.correction(raw.data=data, raw.data.corr=data.corr, lab.replicate="R1.12C_18O.GL", correction.value=R1.12C_18O.GL.td.pos.resid.best.iteration.list$norm.correction, reps.by.trt=WAD.by.taxon.GL$reps.by.trt, vars=c("Density.g.ml", "Sample", "iso.treat.eco"))
    data.corr <- apply.labeled.correction(raw.data=data, raw.data.corr=data.corr, lab.replicate="R2.12C_18O.GL", correction.value=R2.12C_18O.GL.td.pos.resid.best.iteration.list$norm.correction, reps.by.trt=WAD.by.taxon.GL$reps.by.trt, vars=c("Density.g.ml", "Sample", "iso.treat.eco"))
    data.corr <- apply.labeled.correction(raw.data=data, raw.data.corr=data.corr, lab.replicate="R3.12C_18O.GL", correction.value=R3.12C_18O.GL.td.pos.resid.best.iteration.list$norm.correction, reps.by.trt=WAD.by.taxon.GL$reps.by.trt, vars=c("Density.g.ml", "Sample", "iso.treat.eco"))
    data.corr <- apply.labeled.correction(raw.data=data, raw.data.corr=data.corr, lab.replicate="R1.12C_18O_N.GL", correction.value=R1.12C_18O_N.GL.td.pos.resid.best.iteration.list$norm.correction, reps.by.trt=WAD.by.taxon.GL$reps.by.trt, vars=c("Density.g.ml", "Sample", "iso.treat.eco"))
    data.corr <- apply.labeled.correction(raw.data=data, raw.data.corr=data.corr, lab.replicate="R2.12C_18O_N.GL", correction.value=R2.12C_18O_N.GL.td.pos.resid.best.iteration.list$norm.correction, reps.by.trt=WAD.by.taxon.GL$reps.by.trt, vars=c("Density.g.ml", "Sample", "iso.treat.eco"))
    data.corr <- apply.labeled.correction(raw.data=data, raw.data.corr=data.corr, lab.replicate="R3.12C_18O_N.GL", correction.value=R3.12C_18O_N.GL.td.pos.resid.best.iteration.list$norm.correction, reps.by.trt=WAD.by.taxon.GL$reps.by.trt, vars=c("Density.g.ml", "Sample", "iso.treat.eco"))

  #Re-calculate number of copies per uL, based on relative abundance and total number of copies per uL:
    ncopies.corr <- data.corr$neat.avg.16S.copies*data.corr[,7:(ncol(data.corr)-1)]
    dim(ncopies.corr)
    ncopies.corr <- cbind(data.corr[,1:6], ncopies.corr)  # add first 6 columns of data.corr to ncopies.corr
    dim(ncopies.corr)
    head(ncopies.corr)

  #Melt data.corr into long format by tube, sample, tmt, rep, fraction, DNA conc, and density;
  #Do this for copies.ul and for relative abundance, which is just our data.corr file. Merge these to into 1 masterfile: data.corr.melted
    ncopies.corr.melted <- melt(ncopies.corr, id=c("Sample", "fraction", "iso.treat.eco", "neat.avg.16S.copies", "Density.g.ml", "DNA.ng.ul"), measure.vars=as.character(1:length(taxa.id$taxon)), variable.name="taxon", value.name="copies.ul")
    rel.abundance.corr.melted <- melt(data.corr, id=c("Sample", "fraction", "iso.treat.eco", "neat.avg.16S.copies", "Density.g.ml", "DNA.ng.ul"), measure.vars=as.character(1:length(taxa.id$taxon)), variable.name="taxon", value.name="rel.abundance")
    data.corr.melted <- merge(ncopies.corr.melted, rel.abundance.corr.melted)
    head(data.corr.melted)

  #Merge taxa data and reorder data frame by taxon and SampleID AND FRACTION
    data.corr.melted <- merge(data.corr.melted, taxa.id)
    data.corr.melted <- data.corr.melted[order(data.corr.melted$taxon, data.corr.melted$Sample, data.corr.melted$fraction),]
    row.names(data.corr.melted) <- 1:dim(data.corr.melted)[1]   #rename observations to be sequential
    head(data.corr.melted)


#Calculate corrected WADs for each tube and taxon (Grassland only):
  data.corr.melted.GL <- data.corr.melted[data.corr.melted$iso.treat.eco %in% c("12C_18O_N.GL", "12C_18O.GL", "12C_N.GL", "12C.GL", "13C_N.GL", "13C.GL", "16O.GL", "18O.GL"),]
  data.corr.melted.GL$taxon <- factor(data.corr.melted.GL$taxon)
  data.corr.melted.GL$Sample <- factor(data.corr.melted.GL$Sample)
  data.corr.melted.GL$iso.treat.eco <- factor(data.corr.melted.GL$iso.treat.eco)
  row.names(data.corr.melted.GL) <- 1:dim(data.corr.melted.GL)[1]
  WAD.by.taxon.GL.corr <- WAD.by.taxon.func(X=data.corr.melted.GL, vars=c("taxon", "Density.g.ml", "copies.ul", "Sample", "iso.treat.eco"))
  #Look at the results:
    names(WAD.by.taxon.GL.corr)          #names of the two data frames in the output list
    head(WAD.by.taxon.GL.corr$obs.wads)  #looking at the head of the first data frame
    head(WAD.by.taxon.GL.corr[[1]])      #another way to look at the head of the first data frame
    WAD.by.taxon.GL.corr$reps.by.trt     #looking at the head of the second data frame
    WAD.by.taxon.GL.corr[[2]]            #another way to look at the head of the second data frame

  #Quick and dirty look at the level of replication among treatments and the occurrences of data for taxa in tubes 
    table(data.corr.melted.GL$iso.treat.eco, data.corr.melted.GL$Sample)
    table(data.corr.melted.GL$taxon, data.corr.melted.GL$Sample)


#Plot standard error of WADs for each taxon and treatment:
  graphics.off()
  pdf(file="qSIP_output/Figures/RM_DIM_GL_SE_WAD_by_taxon_plots.pdf", width=6, height=11)
  WAD.SE.by.taxon.GL.corr <- SE.WAD.by.taxon.plot(LIST=WAD.by.taxon.GL.corr, percentile=0.95)
  dev.off()
  #look at the head of the output data frame:
    head(WAD.SE.by.taxon.GL.corr)


#Spot-checking that SE's match up with observations:
  WAD.by.taxon.GL.corr$reps.by.trt
  head(WAD.by.taxon.GL.corr$obs.wads)
  head(WAD.SE.by.taxon.GL.corr)
  sum(as.numeric(apply(WAD.by.taxon.GL.corr$obs.wads[,16:18], 1, function(x) sd(x)/sqrt(sum(!is.na(x))))) != WAD.SE.by.taxon.GL.corr[,2], na.rm=TRUE)
  sum(as.numeric(apply(WAD.by.taxon.GL.corr$obs.wads[,c(12,14,15)], 1, function(x) sd(x)/sqrt(sum(!is.na(x))))) != WAD.SE.by.taxon.GL.corr[,3], na.rm=TRUE)
  sum(as.numeric(apply(WAD.by.taxon.GL.corr$obs.wads[,c(2,13,19)], 1, function(x) sd(x)/sqrt(sum(!is.na(x))))) != WAD.SE.by.taxon.GL.corr[,8], na.rm=TRUE)


#Save workspace:
  save(list=ls(all.names=TRUE), file="qSIP_workspaces/RM_DIM_01/.RData", envir=.GlobalEnv)



