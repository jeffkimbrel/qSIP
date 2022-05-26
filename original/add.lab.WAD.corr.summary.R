#Function to apply the previously identified tube-level shifts to the already-summarized taxon-level labeled replicate WADs to 'correct' them and create a new summary list (analagous to the list that is output from 'find.unlabeled.correction') that includes those corrected values:
  
  add.lab.WAD.corr.summary <- function(summary.list, reps.by.trt, lab.reps.to.add, lab.shifts, CI=0.90){
    WAD.corr.list <- summary.list
    #Add to the third object in the summary list:
      lab.corr <- data.frame(matrix(NA, nrow=dim(summary.list$WAD.table.corr)[1], ncol=length(lab.reps.to.add)))
      for (i in 1:length(lab.reps.to.add)){
        names(lab.corr)[i] <- paste(lab.reps.to.add[i], "corr", sep=".")
        lab.corr[,i] <- summary.list$WAD.table.corr[, names(summary.list$WAD.table.corr) == lab.reps.to.add[i]] - lab.shifts[i]
      }
      WAD.corr.list$WAD.table.corr <- data.frame(summary.list$WAD.table.corr, lab.corr)
    #Add to the second object in the summary list:
      WAD.corr.list$corr.names <- c(summary.list$corr.names, names(lab.corr))
    #Add to the first object in the summary list:
      reps.to.add.df.1 <- data.frame(matrix(NA, nrow=length(lab.reps.to.add), ncol=dim(summary.list$WAD.norm.fit.parms)[2]))
      names(reps.to.add.df.1) <- names(summary.list$WAD.norm.fit.parms)
      trt.names <- gsub(pattern="^(R\\d+)\\.(.+)$", replacement="\\2", x=lab.reps.to.add, perl=TRUE)   #pull out just the treatment names
      rep.names <- gsub(pattern="^(R\\d+)\\.(.+)$", replacement="\\1", x=lab.reps.to.add, perl=TRUE)   #pull out just the replicate names
      new.order <- order(summary.list$WAD.table.corr$mean.unlab.WAD.order)
      for (i in 1:length(lab.reps.to.add)){
        WADs.norm.fit.temp <- fit.norm.func(x=eval(parse(text=paste("summary.list$WAD.table.corr$", rep.names[i], ".", trt.names[i], sep="")))[new.order], z=summary.list$WAD.table.corr$mean.unlab.WAD.order[new.order]/max(summary.list$WAD.table.corr$mean.unlab.WAD.order), CI=CI)
        reps.to.add.df.1$trt[i] <- trt.names[i]
        reps.to.add.df.1$tube[i] <- as.character(reps.by.trt[as.character(reps.by.trt[,1]) == trt.names[i], names(reps.by.trt) == rep.names[i]])
        reps.to.add.df.1$rep[i] <- rep.names[i]
        reps.to.add.df.1$tcode[i] <- lab.reps.to.add[i]
        reps.to.add.df.1$mean[i] <- WADs.norm.fit.temp$mean
        reps.to.add.df.1$stdev[i] <- WADs.norm.fit.temp$stdev
        reps.to.add.df.1$NLL[i] <- WADs.norm.fit.temp$NLL
        reps.to.add.df.1$shift.mean[i] <- lab.shifts[i]
      }
      reps.to.add.df.1$corr.mean <- reps.to.add.df.1$mean - reps.to.add.df.1$shift.mean
      WAD.corr.list$WAD.norm.fit.parms <- rbind(WAD.corr.list$WAD.norm.fit.parms, reps.to.add.df.1)
    WAD.corr.list
  }



