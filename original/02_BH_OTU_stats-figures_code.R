# This code performs various analyses and produces various figures from the growth, flux, and excess atom fraction results


#Reload the saved workspace (only necessary if starting fresh from this script and not from 01_BH_OTU_calcs_code...):
#(Adjust path to most recently saved workspace as neccesary)
  setwd("/Users/bk/Research/Projects/SIP_Modeling/qSIP")
  load("qSIP_workspaces/2015_07_30_01_BH_OTU_calcs_code_piecemeal/.RData")
  
#Until 01_BH_OTU_calcs_code... is re-run, need to update functions here to reflect change in 'sample':
  source(paste(func.path, "/sample.vec.R", sep=""))   #sample.vec
  source(paste(func.path, "/WAD.func.R", sep=""))   #WAD.func
  source(paste(func.path, "/fit.norm.func.R", sep=""))   #fit.norm.func
  source(paste(func.path, "/boot.WAD.func.R", sep=""))    #boot.WAD.func
  source(paste(func.path, "/diff.wad.calc.R", sep=""))    #diff.wad.calc
  source(paste(func.path, "/boot.diff.wad.R", sep=""))   #boot.diff.wad
  source(paste(func.path, "/MW.calc.R", sep=""))   #MW.calc
  source(paste(func.path, "/MW.calc.Schildkraut.R", sep=""))   #MW.calc.Schildkraut
  source(paste(func.path, "/comparison.message.R", sep=""))    #comparison.message
  source(paste(func.path, "/ape.calc.R", sep=""))    #ape.calc
  source(paste(func.path, "/boot.diff.ape.R", sep=""))    #boot.diff.ape
  source(paste(func.path, "/r.calc.R", sep=""))   #r.calc
  source(paste(func.path, "/boot.diff.r.R", sep=""))   #boot.diff.r
  source(paste(func.path, "/boot.TUBE.func.R", sep=""))    #boot.TUBE.func
  source(paste(func.path, "/f.calc.R", sep=""))   #f.calc
  source(paste(func.path, "/boot.diff.f.R", sep=""))   #boot.diff.f
  source(paste(func.path, "/all.taxa.calcs.R", sep=""))   #all.taxa.calcs


#_0a_Evaluate how the number of replicates affects the estimate of WAD for the unlabeled treatment:
  #{
  #Subset the data to only those OTUs that occured in all 12 unlabeled (16-O, 12-C) tubes:
    #First, subset the data (that was previously filtered to contain only OTUS that occurred in 34 of 36 tubes) to grab only those taxon-fractions associated with one of the 16O,12C tubes:
      data.melted.unlab <- data.melted[grep(pattern="\\d\\_\\w\\_(?!13C).*16O", data.melted$trt.code, perl=TRUE),]
        #Re-convert the factor columns to factor:
          data.melted.unlab$taxon <- factor(data.melted.unlab$taxon)
          data.melted.unlab$SampleID <- factor(data.melted.unlab$SampleID)
          data.melted.unlab$trt.code <- factor(data.melted.unlab$trt.code)
          data.melted.unlab$tube <- factor(data.melted.unlab$tube)
          data.melted.unlab$kingdom <- factor(data.melted.unlab$kingdom)
          data.melted.unlab$phylum <- factor(data.melted.unlab$phylum)
          data.melted.unlab$class <- factor(data.melted.unlab$class)
          data.melted.unlab$order <- factor(data.melted.unlab$order)
          data.melted.unlab$family <- factor(data.melted.unlab$family)
          data.melted.unlab$genus <- factor(data.melted.unlab$genus)
          data.melted.unlab$species <- factor(data.melted.unlab$species)
      dim(data.melted)
      dim(data.melted.unlab)
      levels(data.melted.unlab$trt.code)

    #Subset the unlabeled data created above to grab only those taxon-fractions with copies present:
      data.melted.unlab.occurrences <- data.melted.unlab[data.melted.unlab$copies > 0,]
        #Re-convert the factor columns to factor:
          data.melted.unlab.occurrences$taxon <- factor(data.melted.unlab.occurrences$taxon)
          data.melted.unlab.occurrences$SampleID <- factor(data.melted.unlab.occurrences$SampleID)
          data.melted.unlab.occurrences$trt.code <- factor(data.melted.unlab.occurrences$trt.code)
          data.melted.unlab.occurrences$tube <- factor(data.melted.unlab.occurrences$tube)
          data.melted.unlab.occurrences$kingdom <- factor(data.melted.unlab.occurrences$kingdom)
          data.melted.unlab.occurrences$phylum <- factor(data.melted.unlab.occurrences$phylum)
          data.melted.unlab.occurrences$class <- factor(data.melted.unlab.occurrences$class)
          data.melted.unlab.occurrences$order <- factor(data.melted.unlab.occurrences$order)
          data.melted.unlab.occurrences$family <- factor(data.melted.unlab.occurrences$family)
          data.melted.unlab.occurrences$genus <- factor(data.melted.unlab.occurrences$genus)
          data.melted.unlab.occurrences$species <- factor(data.melted.unlab.occurrences$species)
      dim(data.melted)
      dim(data.melted.unlab)
      dim(data.melted.unlab.occurrences)

    #Calculate the number of unique unlabeled tubes with copies present for each taxon (OTU):
      unlab.tubes.per.OTU <- tapply(data.melted.unlab.occurrences$tube, data.melted.unlab.occurrences$taxon, length.unique)
      unlab.tubes.per.OTU <- sort(unlab.tubes.per.OTU)

    #Keep only those OTUs that occur in all 12 of the 12 unlabeled tubes of the experiment:
      dev.off()
      dev.new(width=3.5, height=7)
      par(mfrow=c(2,1))
      hist(unlab.tubes.per.OTU, main="")
      abline(v=11.5, col="red")
      plot(x=1:length(unlab.tubes.per.OTU), y=as.numeric(unlab.tubes.per.OTU), xlab="Number of OTUs", ylab="Unlabeled tubes per OTU", type="l")
      abline(h=11.5, col="red")
      par(mfrow=c(1,1))
      # unlab.tubes.per.OTU[as.numeric(unlab.tubes.per.OTU) >= 12]
      length(unlab.tubes.per.OTU[as.numeric(unlab.tubes.per.OTU) >= 12])

      dev.off()

    #Subset the data.melted.unlab dataframe so that it only contains taxon-fractions for OTUs occurring in all 12 of the 12 unlabeled tubes:
      data.melted.unlab <- data.melted.unlab[data.melted.unlab$taxon %in% names(unlab.tubes.per.OTU[as.numeric(unlab.tubes.per.OTU) >= 12]),]
      row.names(data.melted.unlab) <- 1:dim(data.melted.unlab)[1]
        #Re-convert the factor columns to factor:
          data.melted.unlab$taxon <- factor(data.melted.unlab$taxon)
          data.melted.unlab$SampleID <- factor(data.melted.unlab$SampleID)
          data.melted.unlab$trt.code <- factor(data.melted.unlab$trt.code)
          data.melted.unlab$tube <- factor(data.melted.unlab$tube)
          data.melted.unlab$kingdom <- factor(data.melted.unlab$kingdom)
          data.melted.unlab$phylum <- factor(data.melted.unlab$phylum)
          data.melted.unlab$class <- factor(data.melted.unlab$class)
          data.melted.unlab$order <- factor(data.melted.unlab$order)
          data.melted.unlab$family <- factor(data.melted.unlab$family)
          data.melted.unlab$genus <- factor(data.melted.unlab$genus)
          data.melted.unlab$species <- factor(data.melted.unlab$species)
      dim(data.melted.all)
      dim(data.melted)
      dim(data.melted.unlab)


  #Calculate WAD for each OTU by bootstrapping each level of possible replication for all treatment-week combinations:
    #Set the number of unlabeled OTUs for which to do the bootstrapping (i.e., do the first X OTUs of the total number of OTUs):
      length(levels(data.melted.unlab$taxon))
      num.unlab.OTUs <- 2
    #Set the number of bootstrapping iterations:
      draws <- 1000
    #Create an empty data frame for all results:
      unlab.wad.boots <- data.frame(taxon=character(), trt.combo=character(), id=character(), num.reps=numeric(), stringsAsFactors=FALSE)
    #Create a template data frame for each OTU's results:
      levels(data.melted.unlab$trt.code)
      trt.combo.names <- c(levels(data.melted.unlab$trt.code), paste(levels(data.melted.unlab$trt.code)[1:2], collapse=" ; "), paste(levels(data.melted.unlab$trt.code)[3:4], collapse=" ; "), paste(levels(data.melted.unlab$trt.code)[c(1,3)], collapse=" ; "), paste(levels(data.melted.unlab$trt.code)[c(2,4)], collapse=" ; "), paste(levels(data.melted.unlab$trt.code), collapse=" ; "))
      trt.combo.ids <- c("wk1noC", "wk1plusC", "wk6noC", "wk6plusC", "wk1", "wk6", "noC", "plusC", "all")
      trt.combo.max.reps <- c(rep(3, 4), rep(6, 4), 12)
      trt.combo.key <- data.frame(names=trt.combo.names, id=trt.combo.ids, max.reps=trt.combo.max.reps, stringsAsFactors=FALSE)
      trt.combo.key
      unlab.wad.boots.OTU.template <- data.frame(taxon=character(sum(trt.combo.max.reps)), trt.combo=character(sum(trt.combo.max.reps)), id=character(sum(trt.combo.max.reps)), num.reps=numeric(sum(trt.combo.max.reps)), stringsAsFactors=FALSE)
      counter <- 1
      for (i in 1:dim(trt.combo.key)[1]){
        unlab.wad.boots.OTU.template$trt.combo[counter:(counter+trt.combo.key$max.reps[i]-1)] <- rep(trt.combo.key$names[i], trt.combo.key$max.reps[i])
        unlab.wad.boots.OTU.template$id[counter:(counter+trt.combo.key$max.reps[i]-1)] <- rep(trt.combo.key$id[i], trt.combo.key$max.reps[i])
        unlab.wad.boots.OTU.template$num.reps[counter:(counter+trt.combo.key$max.reps[i]-1)] <- 1:trt.combo.key$max.reps[i]
        if (i != dim(trt.combo.key)[1]){
          counter <- min(which(unlab.wad.boots.OTU.template$num.reps == 0))
        }
        else  rm(counter)
      }
    #Fill in the data frame for all results with the (selected) OTUs:
      for (i in 1:num.unlab.OTUs){
        unlab.wad.boots.OTU.to.add <- unlab.wad.boots.OTU.template
        unlab.wad.boots.OTU.to.add$taxon[1:dim(unlab.wad.boots.OTU.template)[1]] <- levels(data.melted.unlab$taxon)[i]
        unlab.wad.boots <- rbind(unlab.wad.boots, unlab.wad.boots.OTU.to.add)
      }
      unlab.wad.boots$taxon <- factor(unlab.wad.boots$taxon)
      unlab.wad.boots$trt.combo <- factor(unlab.wad.boots$trt.combo)
      unlab.wad.boots$id <- factor(unlab.wad.boots$id)
      summary(unlab.wad.boots)
      dim(unlab.wad.boots)
    #Add empty columns to stor results of each bootstrapping iteration:
      empty.boots.only <- data.frame(matrix(NA, dim(unlab.wad.boots)[1], draws))
      names(empty.boots.only) <- paste("wad", 1:draws, sep="")
      unlab.wad.boots <- data.frame(unlab.wad.boots, empty.boots.only)
    #Do the bootstrapping WAD calculations for the different treatment combos and levels of replication for the OTUs:
      set.seed(100)
      for (i in 1:dim(unlab.wad.boots)[1]){
        #Define the treatment codes to use in subsetting the data:
          current.trt.codes <- strsplit(x=as.character(unlab.wad.boots$trt.combo[i]), split="\\s;\\s", perl=TRUE)[[1]]
        #Define the subset of data to pass to boot.WAD.func according to those treatment codes (and the number of reps):
          current.data <- data.melted.unlab[data.melted.unlab$taxon==as.character(unlab.wad.boots$taxon[i]) & data.melted.unlab$trt.code %in% current.trt.codes,]
          current.reps.all <- levels(factor(as.character(current.data$tube)))
          # current.reps <- current.reps.all[1:unlab.wad.boots$num.reps[i]]
          current.data <- current.data[current.data$tube %in% current.reps.all,]
        #Calculate the WAD for the current dataset with bootstrapping (using the flexible version of boot.WAD.func to enable variable values for 'size'):
          current.boot.out <- boot.WAD.func(X=current.data, vars=c("density.g.ml", "copies", "tube"), CI=0.90, draws=1000, size=unlab.wad.boots$num.reps[i])
        #Write the bootstrapped WADs to the correct row of the output dataframe:
          unlab.wad.boots[i,5:dim(unlab.wad.boots)[2]] <- current.boot.out$boot.wads
      }


  #Create figures of WAD vs number of reps (for one OTU to start with):
    dev.off()
    # dev.new(width=7, height=7)
    pdf(file="qSIP_output/Figures/BH_OTU_Taxon623634_unlabeled_wad_vs_reps.pdf", width=7, height=7)
    par(mai=c(1.02, 1.02, 0.82, 0.22))
    y.min <- min(unlab.wad.boots[unlab.wad.boots$taxon == 623634, 5:dim(unlab.wad.boots)[2]])
    y.max <- max(unlab.wad.boots[unlab.wad.boots$taxon == 623634, 5:dim(unlab.wad.boots)[2]])
    plot(x=c(0, max(unlab.wad.boots$num.reps)), y=c(min(unlab.wad.boots[,5:dim(unlab.wad.boots)[2]]), max(unlab.wad.boots[,5:dim(unlab.wad.boots)[2]])), bty="l", type="n", ylab=expression(paste("Weighted average density (g mL"^-1, ")", sep="")), xlab="Number of replicates", ylim=c(y.min, y.max))
      arrows(x0=unlab.wad.boots$num.reps[unlab.wad.boots$id == "all" & unlab.wad.boots$taxon == 623634], y0=as.numeric(apply(unlab.wad.boots[unlab.wad.boots$id == "all"  & unlab.wad.boots$taxon == 623634, 5:dim(unlab.wad.boots)[2]], 1, quantile, probs=0.05)), x1=unlab.wad.boots$num.reps[unlab.wad.boots$id == "all" & unlab.wad.boots$taxon == 623634], y1=as.numeric(apply(unlab.wad.boots[unlab.wad.boots$id == "all"  & unlab.wad.boots$taxon == 623634, 5:dim(unlab.wad.boots)[2]], 1, quantile, probs=0.95)), length=0.00, angle=90, code=3, col="red", lwd=2)
      points(x=unlab.wad.boots$num.reps[unlab.wad.boots$id == "all" & unlab.wad.boots$taxon == 623634], y=as.numeric(apply(unlab.wad.boots[unlab.wad.boots$id == "all"  & unlab.wad.boots$taxon == 623634, 5:dim(unlab.wad.boots)[2]], 1, median)), pch=21, col="red", bg="red", cex=1.1)
      mtext(text="OTU 623634", side=3, line=0, cex=1)

    dev.off()


    dev.off()
    # dev.new(width=7, height=7)
    pdf(file="qSIP_output/Figures/BH_OTU_Taxon623634_unlabeled_wad_vs_reps_with_tmts.pdf", width=7, height=7)
    par(mai=c(1.02, 1.02, 0.82, 0.22))
    y.min <- min(unlab.wad.boots[unlab.wad.boots$taxon == 623634, 5:dim(unlab.wad.boots)[2]])
    y.max <- max(unlab.wad.boots[unlab.wad.boots$taxon == 623634, 5:dim(unlab.wad.boots)[2]])
    plot(x=c(0, max(unlab.wad.boots$num.reps)), y=c(min(unlab.wad.boots[,5:dim(unlab.wad.boots)[2]]), max(unlab.wad.boots[,5:dim(unlab.wad.boots)[2]])), bty="l", type="n", ylab=expression(paste("Weighted average density (g mL"^-1, ")", sep="")), xlab="Number of replicates", ylim=c(y.min, y.max))
      T623634.unlab.wad.all.reps <- as.numeric(apply(unlab.wad.boots[unlab.wad.boots$id == "all"  & unlab.wad.boots$taxon == 623634 & unlab.wad.boots$num.reps == 1, 5:dim(unlab.wad.boots)[2]], 1, unique))
      points(x=jitter(rep(1, length(T623634.unlab.wad.all.reps)), factor=4), y=T623634.unlab.wad.all.reps, pch=21, col="black", bg="white", cex=1)
      T623634.unlab.all.tmts.jitter <- jitter(unlab.wad.boots$num.reps[c(3,6,9,12)], factor=3)
      arrows(x0=T623634.unlab.all.tmts.jitter, y0=as.numeric(apply(unlab.wad.boots[c(3,6,9,12), 5:dim(unlab.wad.boots)[2]], 1, quantile, probs=0.05)), x1=T623634.unlab.all.tmts.jitter, y1=as.numeric(apply(unlab.wad.boots[c(3,6,9,12), 5:dim(unlab.wad.boots)[2]], 1, quantile, probs=0.95)), length=0.03, angle=90, code=3, col=c("darkgreen", "blue", "purple", "darkorange4"))
      points(x=T623634.unlab.all.tmts.jitter, y=as.numeric(apply(unlab.wad.boots[c(3,6,9,12), 5:dim(unlab.wad.boots)[2]], 1, median)), pch=21, col="black", bg=c("darkgreen", "blue", "purple", "darkorange4"), cex=1)
      T623634.unlab.paired.tmts.jitter <- jitter(unlab.wad.boots$num.reps[c(18,24,30,36)], factor=2)
      arrows(x0=T623634.unlab.paired.tmts.jitter, y0=as.numeric(apply(unlab.wad.boots[c(18,24,30,36), 5:dim(unlab.wad.boots)[2]], 1, quantile, probs=0.05)), x1=T623634.unlab.paired.tmts.jitter, y1=as.numeric(apply(unlab.wad.boots[c(18,24,30,36), 5:dim(unlab.wad.boots)[2]], 1, quantile, probs=0.95)), length=0.03, angle=90, code=3, col=c("darkolivegreen4", "cornflowerblue", "magenta2", "coral1"))
      points(x=T623634.unlab.paired.tmts.jitter, y=as.numeric(apply(unlab.wad.boots[c(18,24,30,36), 5:dim(unlab.wad.boots)[2]], 1, median)), pch=21, col="black", bg=c("darkolivegreen4", "cornflowerblue", "magenta2", "coral1"), cex=1)
      arrows(x0=unlab.wad.boots$num.reps[unlab.wad.boots$id == "all" & unlab.wad.boots$taxon == 623634], y0=as.numeric(apply(unlab.wad.boots[unlab.wad.boots$id == "all"  & unlab.wad.boots$taxon == 623634, 5:dim(unlab.wad.boots)[2]], 1, quantile, probs=0.05)), x1=unlab.wad.boots$num.reps[unlab.wad.boots$id == "all" & unlab.wad.boots$taxon == 623634], y1=as.numeric(apply(unlab.wad.boots[unlab.wad.boots$id == "all"  & unlab.wad.boots$taxon == 623634, 5:dim(unlab.wad.boots)[2]], 1, quantile, probs=0.95)), length=0.00, angle=90, code=3, col="red", lwd=2)
      points(x=unlab.wad.boots$num.reps[unlab.wad.boots$id == "all" & unlab.wad.boots$taxon == 623634], y=as.numeric(apply(unlab.wad.boots[unlab.wad.boots$id == "all"  & unlab.wad.boots$taxon == 623634, 5:dim(unlab.wad.boots)[2]], 1, median)), pch=21, col="red", bg="red", cex=1.1)
      mtext(text="OTU 623634", side=3, line=0, cex=1)
      legend(x=12, y=y.max, xjust=1, yjust=1, bty="n", legend=c("raw reps", "wk1noC", "wk1plusC", "wk6noC", "wk6plusC", "wk1", "wk6", "noC", "plusC", "combined reps"), pch=21, col=c(rep("black", 9), "red"), pt.bg=c("white", "darkgreen", "blue", "purple", "darkorange4", "darkolivegreen4", "cornflowerblue", "magenta2", "coral1", "red"), cex=0.5, pt.cex=c(rep(1, 9), 1.1))

    dev.off()
  #}
