#Function to create the 'WAD.norm.fit.parms' and 'WAD.table.corr' data frames along with the vector 'corr.names':

#Calculates tube-level WAD correction values and applies them to each replicate of the specified unlabeled treatments using WADS of the taxa common to ALL replicates of ALL treatments (unlabeled and labeled) specificed in the function call
#And summarizes the unlabeled correction (though it does not actually apply it to the raw data)
#input: a list of the same form as the output from 'WAD.by.taxon.func'; there are two data frames: (1) observed WADs by taxa (rows) for each relicate (column) and (2) a table giving the IDs of the replicates associated with each treatment
#treatments must be subsets of those included in LIST
#requires the function 'fit.norm.func'
#MAKE THE DEFAULT TO INCLUDE ALL LABELED AND UNLABELED TREATMENTS (FOR DETERIMINGIN TAXA TO INCLUDE IN CALCS)?
#DEAL WITH THE GRAPHS WITHIN THIS FUNCTION



find.unlabeled.correction <- function(LIST, unlab.tmts, lab.tmts, CI=0.90){
  treatments <- c(unlab.tmts, lab.tmts)
  obs.wads.by.taxon <- LIST[[1]]
  reps.by.trt <- LIST[[2]]
  
  #Throw an error if the treatments aren't included in the data:
  if (sum(!(treatments %in% reps.by.trt[,1])) > 0){
    stop("at least one of the specified treatments does not occur in the data")
  }
  
  #Subset the data to include only those treatments specified:
  reps.by.trt <- reps.by.trt[reps.by.trt[,1] %in% treatments,]
    #Convert the factor columns to factor:
    reps.by.trt <- as.data.frame(lapply(reps.by.trt, function(x) if(is.factor(x)) factor(x) else x))
  curr.data <- cbind(taxon=obs.wads.by.taxon[,1], obs.wads.by.taxon[,names(obs.wads.by.taxon) %in% as.character(unlist(reps.by.trt[,2:dim(reps.by.trt)[2]]))])
  
  #Remove any taxon that does not occur in all replicates of all specified treatments (unlabeled & labeled combined):
  curr.data <- curr.data[apply(is.na(curr.data[2:dim(curr.data)[2]]), 1, sum) == 0, ]
  curr.data$taxon <- factor(curr.data$taxon)
  
  #Verify that there are no longer any missing taxa from tubes in the full dataset (commented-out):
  # apply(is.na(curr.data[2:dim(curr.data)[2]]), 1, sum)   #The number of replicates in curr.data for which a given taxon does not occur
  # sort(apply(is.na(curr.data[2:dim(curr.data)[2]]), 1, sum))   #The number of replicates in curr.data for which a given taxon does not occur
  # curr.data$taxon[apply(is.na(curr.data[2:dim(curr.data)[2]]), 1, sum) == dim(curr.data)[2]-1]   #Those taxa that do not occur in at least one of the replicates of curr.data after filtering missing taxa

  #Find the rank-order of taxa in each tube (ranked by WAD) & add this info to the dataframe of WADs by tube (actually creates a new dataframe with WADs & ranks):
  curr.num.reps <- dim(reps.by.trt)[1] * (dim(reps.by.trt)[2]-1)
  rep.orders <- data.frame(matrix(NA, nrow=dim(curr.data)[1], ncol=curr.num.reps))
  counter <- 1
  for (i in 1:length(treatments)){
    for (j in 1:(dim(reps.by.trt)[2]-1)){
      if (!is.na(as.character(reps.by.trt[reps.by.trt[,1] == treatments[i], 1+j]))){
        d <- curr.data[, names(curr.data) == as.character(reps.by.trt[reps.by.trt[,1] == treatments[i], 1+j])]
        o <- order(d)
        r <- d[o]
        s <- order(r)
        rep.orders[,counter] <- s[order(o)]
      }
      names(rep.orders)[counter] <- paste(names(reps.by.trt)[1+j], as.character(treatments[i]), "order", sep=".")
      counter <- counter + 1
    }
  }
  curr.data.orders <- data.frame(rep.orders, curr.data)
  names(curr.data.orders) <- c(names(rep.orders), names(curr.data))
  for (i in 2:dim(curr.data)[2]){
    rep.name <- as.character(names(reps.by.trt[2:dim(reps.by.trt)[2]])[apply(reps.by.trt[,2:dim(reps.by.trt)[2]] == as.character(names(curr.data)[i]), 2, sum, na.rm=TRUE) == 1])
    tmt.name <- as.character(reps.by.trt[apply(reps.by.trt[,2:dim(reps.by.trt)[2]] == as.character(names(curr.data)[i]), 1, sum, na.rm=TRUE) == 1, 1])
    names(curr.data.orders)[curr.num.reps+i] <- paste(rep.name, tmt.name, sep=".")
  }

  #Select the (unlabeled) treatments to include in a global average:
  tmts.to.average <- unlab.tmts
  #Get the rank-order of taxa according to mean WAD across unlabeled replicates & add this to the dataframe:
  d <- curr.data[, as.character(names(curr.data)) %in% as.character(unlist(reps.by.trt[reps.by.trt[,1] %in% tmts.to.average, 2:dim(reps.by.trt)[2]]))]
  m <- apply(d, 1, mean, na.rm=TRUE)
  o <- order(m)
  r <- m[o]
  s <- order(r)
  mean.order <- s[order(o)]
  assign(x=paste("mean", "unlab", "WAD", "order", sep="."), value=mean.order)
  old.names <- names(curr.data.orders)
  curr.data.orders <- data.frame(mean.unlab.WAD.order, curr.data.orders)
  #Check that illegal names aren't used for treatment IDs:
  new.names <- names(curr.data.orders)[2:dim(curr.data.orders)[2]]
  if (sum(old.names != new.names) > 0){
    stop("treatment names must be syntactically valid variable names in R (see ?make.names)")    
  }

  #Create a dataframe containing only the WADs for all replicates of the unlabeled treatments:
  max.reps.trt <- dim(reps.by.trt)[2] - 1
  rep.nums <- paste("R", 1:max.reps.trt, sep="")
  tmts.to.average.col.names <- as.vector(outer(rep.nums, tmts.to.average, paste, sep="."))
  unlab.WADs <- curr.data.orders[, names(curr.data.orders) %in% tmts.to.average.col.names]

  #Verify that there are no longer any missing taxa from tubes in only the subset of unlabeled treatments (commented-out):
  # sort(apply(unlab.WADs, 1, function(x) sum(is.na(x))))

  #Fit a normal cdf to all unlabeled WADs arranged according to mean WAD across unlabeled reps:
  #There are two alternative ways of doing this; the second one yields crazy-high variance (so it is commented-out):
    # (1) Normal cdf fit to the mean of all unlabeled data (takes several seconds to compute):
    new.order <- order(curr.data.orders$mean.unlab.WAD.order)
    unlab.WADs.norm.fit <- fit.norm.func(x=apply(unlab.WADs, 1, mean, na.rm=TRUE)[new.order], z=curr.data.orders$mean.unlab.WAD.order[new.order]/max(curr.data.orders$mean.unlab.WAD.order), CI=CI)
    # (2) Normal cdf fit to all unlabeled data points (takes several minutes to compute):
    # unlab.WADs.norm.fit.all <- fit.norm.func(x=as.numeric(unlist(unlab.WADs)), z=rep(curr.data.orders$mean.unlab.WAD.order[new.order], dim(unlab.WADs)[2])/max(curr.data.orders$mean.unlab.WAD.order), CI=CI)

  #Fit a normal cdf to the data (arranged according to mean WAD across reps) within each unlabeled tube; also include the global unlabeled mean from above (takes several minutes to compute):
  WAD.norm.fit.parms <- data.frame(matrix(NA, nrow=dim(unlab.WADs)[2]+1, ncol=7))
  names(WAD.norm.fit.parms) <- c("trt", "tube", "rep", "tcode", "mean", "stdev", "NLL")
    #Global mean of all unlabeled tubes (normal cdf fit to the taxon-means of all unlabeled data):
    WAD.norm.fit.parms$trt[1] <- "unlabeled"
    WAD.norm.fit.parms$tube[1] <- "global"
    WAD.norm.fit.parms$rep[1] <- "G"
    WAD.norm.fit.parms$tcode[1] <- "G.unlab"
    WAD.norm.fit.parms$mean[1] <- unlab.WADs.norm.fit$mean
    WAD.norm.fit.parms$stdev[1] <- unlab.WADs.norm.fit$stdev
    WAD.norm.fit.parms$NLL[1] <- unlab.WADs.norm.fit$NLL
    #Calculate parameters for each tube:
    trt.names <- gsub(pattern="^(R\\d+)\\.(.+)$", replacement="\\2", x=names(unlab.WADs), perl=TRUE)   #pull out just the treatment names
    rep.names <- gsub(pattern="^(R\\d+)\\.(.+)$", replacement="\\1", x=names(unlab.WADs), perl=TRUE)   #pull out just the replicate names
    for (i in 2:(dim(unlab.WADs)[2]+1)){
      trt.name <- trt.names[i-1]
      rep.name <- rep.names[i-1]
      WADs.norm.fit.temp <- fit.norm.func(x=eval(parse(text=paste("curr.data.orders$", rep.name, ".", trt.name, sep="")))[new.order], z=curr.data.orders$mean.unlab.WAD.order[new.order]/max(curr.data.orders$mean.unlab.WAD.order), CI=CI)
      WAD.norm.fit.parms$trt[i] <- trt.name
      WAD.norm.fit.parms$tube[i] <- as.character(reps.by.trt[as.character(reps.by.trt[,1]) == trt.name, names(reps.by.trt) == rep.name])
      WAD.norm.fit.parms$rep[i] <- rep.name
      WAD.norm.fit.parms$tcode[i] <- paste(rep.name, trt.name, sep=".")
      WAD.norm.fit.parms$mean[i] <- WADs.norm.fit.temp$mean
      WAD.norm.fit.parms$stdev[i] <- WADs.norm.fit.temp$stdev
      WAD.norm.fit.parms$NLL[i] <- WADs.norm.fit.temp$NLL
      rm(rep.name, trt.name, WADs.norm.fit.temp)
    }
    #Re-order the rows in 'WAD.norm.fit.parms':
    WAD.norm.fit.parms <- WAD.norm.fit.parms[order(WAD.norm.fit.parms$trt, decreasing=TRUE), ]

  #For each unlabeled replicate, calculate the shift from the global unlabeled mean and apply that shift to get the "corrected mean":
  WAD.norm.fit.parms$shift.mean <- WAD.norm.fit.parms$mean - WAD.norm.fit.parms$mean[WAD.norm.fit.parms$tcode == "G.unlab"]
  WAD.norm.fit.parms$corr.mean <- WAD.norm.fit.parms$mean - WAD.norm.fit.parms$shift.mean

  #Correct the summarized data (WADs by taxon) using the shifts identified here:
  #(Note that this DOES NOT correct the raw data)
  WAD.norm.fit.parms.temp <- WAD.norm.fit.parms[WAD.norm.fit.parms$shift.mean != 0 & !is.na(WAD.norm.fit.parms$shift.mean),]
  row.names(WAD.norm.fit.parms.temp) <- 1:dim(WAD.norm.fit.parms.temp)[1]
  unlab.corr <- data.frame(matrix(NA, nrow=dim(curr.data.orders)[1], ncol=dim(WAD.norm.fit.parms.temp)[1]))
  for (i in 1:dim(WAD.norm.fit.parms.temp)[1]){
    names(unlab.corr)[i] <- paste(WAD.norm.fit.parms.temp$tcode[i], "corr", sep=".")
    unlab.corr[,i] <- curr.data.orders[, names(curr.data.orders) == WAD.norm.fit.parms.temp$tcode[i]] - WAD.norm.fit.parms.temp$shift.mean[i]
  }
  curr.data.corr <- data.frame(curr.data.orders, unlab.corr)

  list(WAD.norm.fit.parms=WAD.norm.fit.parms, corr.names=names(unlab.corr), WAD.table.corr=curr.data.corr)
}



