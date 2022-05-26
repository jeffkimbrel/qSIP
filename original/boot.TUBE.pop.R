
boot.TUBE.pop <- function(X, M.soil, vars=c("copies", "tube", "g.soil"), vol=50, CI=0.90, draws=1000){

  # Create dataframe of only abundance and rep: 
  test.data <- data.frame(N=X[,vars[1]], rep=factor(X[,vars[2]]))

  # Calculate observed abundance, total 16S copies, and mass of soil for each rep:
  obs.N <- data.frame(matrix(nrow=length(levels(factor(test.data$rep))), ncol=4))
  names(obs.N) <- c("copies", "tot.copies", "g.soil", "rep")
  for (r in 1:length(levels(test.data$rep))){
    obs.N$rep[r] <- levels(test.data$rep)[r]
    obs.N$copies[r] <- test.data$N[test.data$rep == levels(test.data$rep)[r]]
    obs.N$tot.copies[r] <- test.data$N[test.data$rep == levels(test.data$rep)[r]]*vol
    obs.N$g.soil[r] <- M.soil[M.soil[,vars[2]] == levels(test.data$rep)[r], vars[3]]
  }
  obs.N$rep <- factor(obs.N$rep)

  # Bootstrapping: Calculate a bootstrap vector of mean WADs across reps along with total 16S copies & mass of soil:
  boot.N <- numeric()
  boot.tot.copies <- boot.g.soil <- matrix(nrow=draws, ncol=dim(obs.N)[1])
  for (i in 1:draws){
    indices <- sample.vec(1:dim(obs.N)[1], dim(obs.N)[1], replace=TRUE)
    boot.N[i] <- mean(obs.N$copies[indices], na.rm=T)
    boot.tot.copies[i,] <- obs.N$tot.copies[indices]
    boot.g.soil[i,] <- obs.N$g.soil[indices]
  }
  boot.N.CI <- quantile(boot.N, probs=c((1-CI)/2, 1-((1-CI)/2)), na.rm=T)

  reps.NAs <- obs.N$rep[is.na(obs.N$copies) | obs.N$copies == 0]
  if (length(reps.NAs) == 0){
    message <- "none"
  }
  else  message <- paste("Warning: no occurrences in rep ", paste(reps.NAs, collapse=" & "), sep="")

  return(list(boot.tot.copies=boot.tot.copies, 
              boot.g.soil=boot.g.soil, 
              boot.N=boot.N, 
              obs.N=obs.N, 
              obs.N.mean=mean(obs.N$copies, na.rm=T), 
              boot.N.mean=mean(boot.N, na.rm=T), 
              boot.N.median=median(boot.N, na.rm=T), 
              boot.N.CI=boot.N.CI, 
              message=message))
}
