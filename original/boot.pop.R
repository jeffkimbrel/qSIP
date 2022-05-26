### Given a data frame, calculate number of copies (per uL) and bootstrapped CI based on all reps
#
#     output = boot.pop(X, vars=c("copies", "tube"), CI=0.90, draws=1000)
#
#     Written by Ben Koch


boot.pop <- function(X, vars=c("copies", "tube"), CI=0.90, draws=1000){

  # Create dataframe of only abundance and rep: 
  test.data <- data.frame(N=X[,vars[1]], rep=factor(X[,vars[2]]))

  # Pull out observed abundance for each rep:
  obs.N <- data.frame(matrix(nrow=length(levels(test.data$rep)), ncol=2))
  names(obs.N) <- c("copies", "rep")
  for (r in 1:length(levels(test.data$rep))){
    obs.N$rep[r] <- levels(test.data$rep)[r]
    obs.N$copies[r] <- test.data$N[test.data$rep == levels(test.data$rep)[r]]
  }
  obs.N$rep <- factor(obs.N$rep)

  # Bootstrapping: Calculate a bootstrap vector of mean abundance across reps:
  boot.N <- numeric()
  for (i in 1:draws){
    boot.N[i] <- mean(sample.vec(obs.N$copies, dim(obs.N)[1], replace=TRUE), na.rm=T)
  }
  boot.N.CI <- quantile(boot.N, probs=c((1-CI)/2, 1-((1-CI)/2)), na.rm=T)

  reps.NAs <- obs.N$rep[is.na(obs.N$copies) | obs.N$copies == 0]
  if (length(reps.NAs) == 0){
    message <- "none"
  }
  else  message <- paste("Warning: no occurrences in rep ", paste(reps.NAs, collapse=" & "), sep="")

  return(list(boot.N=boot.N, obs.N=obs.N, obs.N.mean=mean(obs.N$copies, na.rm=T), boot.N.mean=mean(boot.N, na.rm=T), boot.N.median=median(boot.N, na.rm=T), boot.N.CI=boot.N.CI, message=message))
}
