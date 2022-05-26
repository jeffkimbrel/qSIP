### Given a SIP data frame and extraction data, calculate weighted average density (WAD) and bootstrapped CI based on all reps along with the corresponding values for soil mass and total 16S copies from the tubes chosen for each bootstrap sample
#
#     output = boot.TUBE.func(X, M.soil, vars=c("density.g.ml", "copies", "tube", "g.soil"), v.frac, CI=0.90, draws=1000)
#
#     X: data frame that includes x values (e.g., density of DNA) and y values (e.g., number of copies) and the replicate IDs (e.g., tube number)
#     M.soil: data frame with mass of soil for each replicate; must contain two columns, one for replicate ID and one for mass of soil
#     vars: vector of variable names for the x, y, replicate ID, and soil mass columns in the data frame; default is c("density.g.ml", "copies" ,"tube", "g.soil")
#     v.frac: volume of a fraction within a tube (uL); default is 50
#     CI: confidence interval (0-1); default is 0.90
#     draws: number of bootstrap iterations; default is 1000
#     -------------------------------------------------------
#     output: 
#     list of:  boot.tot.copies: matrix (draws x number of replicates) of randomly sampled total copies per tube that comprise each bootstrap sample
#               boot.g.soil: matrix (draws x number of replicates) of randomly sampled grams of soil per tube that comprise each bootstrap sample
#               boot.wads: bootstrapped weighted-average densities (vector of length=draws)
#               obs.wads: data frame of observed weighted average density for each replicate
#               obs.wad.mean: mean of observed weighted average density for each replicate (single value)
#               boot.wads.mean: mean of bootstrapped weighted average densities (single value)
#               boot.wads.median: median of bootstrapped weighted average densities (single value)
#               boot.wads.CI: confidence interval of bootstrapped weighted average densities (named vector of length=2)
#               message: warning message listing replicates that had no occurrences (i.e., all y values for a rep are 0); value for message is "none" when all replicates are present (single character value)
#
#     notes:    requires the function 'WAD.func'
#               names of the replicate ID columns in the SIP data frame (X) must also match the name of the replicate ID column in the M.soil data frame
#               values for replicate ID in the M.soil data frame must be unique among treatments and must match exactly with those in the SIP data frame
#
#     Written by Ben Koch & Natasja van Gestel


boot.TUBE.func <- function(X, M.soil, vars=c("density.g.ml", "copies", "tube", "g.soil"), v.frac, CI=0.90, draws=1000){

  # Create a dataframe of only x, y, and rep: 
  test.data <- data.frame(x=X[,vars[1]], y=X[,vars[2]], rep=factor(X[,vars[3]]))

  # Calculate observed weighted average density (WAD), total 16S copies, and mass of soil for each rep:
  obs.wads <- data.frame(matrix(nrow=length(levels(test.data$rep)), ncol=4))
  names(obs.wads) <- c("wad", "tot.copies", "g.soil", "rep")
  for (r in 1:length(levels(test.data$rep))){
   obs.wads$rep[r] <- levels(test.data$rep)[r]
   obs.wads$wad[r] <- WAD.func(y=test.data$y[test.data$rep == levels(test.data$rep)[r]], x=test.data$x[test.data$rep == levels(test.data$rep)[r]])
   obs.wads$tot.copies[r] <- sum(test.data$y[test.data$rep == levels(test.data$rep)[r]])*v.frac
   obs.wads$g.soil[r] <- M.soil[M.soil[,vars[3]] == levels(test.data$rep)[r], vars[4]]
  }
  obs.wads$rep <- factor(obs.wads$rep)

  # Bootstrapping: Calculate a bootstrap vector of mean WADs across reps along with total 16S copies & mass of soil:
  boot.wads <- numeric()
  boot.tot.copies <- boot.g.soil <- matrix(nrow=draws, ncol=dim(obs.wads)[1])
  for (i in 1:draws){
    indices <- sample.vec(1:dim(obs.wads)[1], dim(obs.wads)[1], replace=TRUE)
    boot.wads[i] <- mean(obs.wads$wad[indices], na.rm=T)
    boot.tot.copies[i,] <- obs.wads$tot.copies[indices]
    boot.g.soil[i,] <- obs.wads$g.soil[indices]
  }
  boot.wads.CI <- quantile(boot.wads, probs=c((1-CI)/2, 1-((1-CI)/2)), na.rm=T)
  
  reps.NAs <- obs.wads$rep[is.na(obs.wads$wad)]
  if (length(reps.NAs) == 0){
   message <- "none"
  }
  else  message <- paste("Warning: no occurrences in rep ", paste(reps.NAs, collapse=" & "), sep="")

  return(list(boot.tot.copies=boot.tot.copies, 
              boot.g.soil=boot.g.soil, 
              boot.wads=boot.wads, 
              obs.wads=obs.wads, 
              obs.wad.mean=mean(obs.wads$wad, na.rm=T), 
              boot.wads.mean=mean(boot.wads, na.rm=T), 
              boot.wads.median=median(boot.wads, na.rm=T), 
              boot.wads.CI=boot.wads.CI, 
              message=message))
}
