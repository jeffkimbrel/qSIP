### Given a data frame, calculate weighted average density (WAD) and bootstrapped CI based using a subset of all reps that occur in that data frame
### Note that this is a more flexible version of the original version of 'boot.WAD.func', where the number of reps to resample (i.e., the 'size') can be chosen by the user (this is useful in simulating different levels of replication)
#
#     output = boot.WAD.func(X, vars=c("density.g.ml", "copies", "tube"), CI=0.90, draws=1000, size=NULL)
#
#     X: data frame that includes x values (e.g., density of DNA) and y values (e.g., number of copies) and the replicate IDs (e.g., tube number)
#     vars: vector of variable names for the x, y, and replicate ID columns in the data frame; default is c("density.g.ml", "copies" ,"tube")
#     CI: confidence interval (0-1); default is 0.90
#     draws: number of bootstrap iterations; default is 1000
#     size: the number of reps to draw at each bootstrap iteration; default is the same number of total reps that occur in the data (i.e., the correct way to do a true bootstrap)
#     -------------------------------------------------------
#     output: 
#     list of:  boot.wads: bootstrapped weighted-average densities (vector of length=draws)
#               obs.wads: data frame of observed weighted average density for each replicate
#               obs.wad.mean: mean of observed weighted average density for each replicate (single value)
#               boot.wads.mean: mean of bootstrapped weighted average densities (single value)
#               boot.wads.median: median of bootstrapped weighted average densities (single value)
#               boot.wads.CI: confidence interval of bootstrapped weighted average densities (named vector of length=2)
#               message: warning message listing replicates that had no occurrences (i.e., all y values for a rep are 0); value for message is "none" when all replicates are present (single character value)
#
#     notes:    requires the function 'WAD.func'
# 
#     Written by Ben Koch & Natasja van Gestel


  boot.WAD.func <- function(X, vars=c("density.g.ml", "copies", "tube"), CI=0.90, draws=1000, size=NULL){

    # Create a dataframe of only x, y, and rep: 
     test.data <- data.frame(x=X[,vars[1]], y=X[,vars[2]], rep=factor(X[,vars[3]]))

    # Calculate observed weighted average density (WAD) for each rep:
     obs.wads <- data.frame(matrix(nrow=length(levels(test.data$rep)), ncol=2))
     names(obs.wads) <- c("wad", "rep")
     for (r in 1:length(levels(test.data$rep))){
       obs.wads$rep[r] <- levels(test.data$rep)[r]
       obs.wads$wad[r] <- WAD.func(y=test.data$y[test.data$rep == levels(test.data$rep)[r]], x=test.data$x[test.data$rep == levels(test.data$rep)[r]])
        }
      obs.wads$rep <- factor(obs.wads$rep)

    # Bootstrapping: Calculate a bootstrap vector of mean WADs across reps:
     if (is.null(size)){
       size <- dim(obs.wads)[1]
     }
     boot.wads <- numeric()
     for (i in 1:draws){
       boot.wads[i] <- mean(sample.vec(obs.wads$wad, size, replace=TRUE), na.rm=T)
     }
     boot.wads.CI <- quantile(boot.wads, probs=c((1-CI)/2, 1-((1-CI)/2)), na.rm=T)
     reps.NAs <- obs.wads$rep[is.na(obs.wads$wad)]
     if (length(reps.NAs) == 0){
       message <- "none"
     }
     else  message <- paste("Warning: no occurrences in rep ", paste(reps.NAs, collapse=" & "), sep="")

  return(list(boot.wads=boot.wads, obs.wads=obs.wads, obs.wad.mean=mean(obs.wads$wad, na.rm=T), boot.wads.mean=mean(boot.wads, na.rm=T), boot.wads.median=median(boot.wads, na.rm=T), boot.wads.CI=boot.wads.CI, message=message))

  }
