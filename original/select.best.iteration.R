#Function to quantitatively assess the fit of each iterated subset of taxa and identify candidate sets of non-growers along with one final set of non-growers:

#quantile.threshold: the threshold percentile for identifying a set of potential non-growers according to a subset of selected metrics
select.best.iteration <- function(putative.nongrower.metrics, quantile.threshold=0.10){
  #Find the rows with the lowest absolute value of difference in slope (corrected-free - fixed):
    abs.diff.slope.threshold <- quantile(x=abs(putative.nongrower.metrics$corr.free.slope - putative.nongrower.metrics$fixed.slope), probs=quantile.threshold)
    row.set1 <- which(abs(putative.nongrower.metrics$corr.free.slope - putative.nongrower.metrics$fixed.slope) <= abs.diff.slope.threshold)
  #Find the rows with the lowest absolute value of difference in normal-fitted SD (labeled - unlabeled):
    abs.diff.norm.SD.threshold <- quantile(x=abs(putative.nongrower.metrics$corr.free.norm.sd.Y-putative.nongrower.metrics$norm.sd.X), probs=quantile.threshold)
    row.set2 <- which(abs(putative.nongrower.metrics$corr.free.norm.sd.Y-putative.nongrower.metrics$norm.sd.X) <= abs.diff.norm.SD.threshold)
  #Find the rows with the lowest absolute value of difference in empirical SD (labeled - unlabeled):
    abs.diff.empir.SD.threshold <- quantile(x=abs(putative.nongrower.metrics$corr.free.sd.Y-putative.nongrower.metrics$sd.X), probs=quantile.threshold)
    row.set3 <- which(abs(putative.nongrower.metrics$corr.free.sd.Y-putative.nongrower.metrics$sd.X) <= abs.diff.empir.SD.threshold)
  #Find the rows with the lowest absolute value of difference in empirical mean (labeled - unlabeled):
    abs.diff.empir.mean.threshold <- quantile(x=abs(putative.nongrower.metrics$corr.free.mean.Y-putative.nongrower.metrics$mean.X), probs=quantile.threshold)
    row.set4 <- which(abs(putative.nongrower.metrics$corr.free.mean.Y-putative.nongrower.metrics$mean.X) <= abs.diff.empir.mean.threshold)
  #Create a data frame of the best candidate sets of non-growers (i.e., rows or iterations) of 'putative.nongrower.metrics':
    cand.iterations <- as.data.frame(table(c(row.set1, row.set2, row.set3, row.set4)))
    names(cand.iterations) <- c("iteration", "freq")
    cand.iterations$iteration <- as.numeric(as.character(cand.iterations$iteration))
    cand.iterations$abs.diff.slope <- abs(putative.nongrower.metrics$corr.free.slope[cand.iterations$iteration] - putative.nongrower.metrics$fixed.slope[cand.iterations$iteration])
    cand.iterations$abs.diff.norm.SD <- abs(putative.nongrower.metrics$corr.free.norm.sd.Y[cand.iterations$iteration]-putative.nongrower.metrics$norm.sd.X[cand.iterations$iteration])
    cand.iterations$abs.diff.empir.SD <- abs(putative.nongrower.metrics$corr.free.sd.Y[cand.iterations$iteration]-putative.nongrower.metrics$sd.X[cand.iterations$iteration])
    cand.iterations$abs.diff.empir.mean <- abs(putative.nongrower.metrics$corr.free.mean.Y[cand.iterations$iteration]-putative.nongrower.metrics$mean.X[cand.iterations$iteration])
    cand.iterations$norm.correction <- putative.nongrower.metrics$fixed.norm.mean.Y[cand.iterations$iteration] - putative.nongrower.metrics$norm.mean.X[cand.iterations$iteration]
  #Iteration(s) occuring within the specified quantile most frequently across all metrics:
    freq.iterations <- cand.iterations$iteration[cand.iterations$freq == max(cand.iterations$freq)]
  #Iterations(s) with the minimum absolute difference in slope (from the free-regression of corrected labeled WADs vs unlabeled WADs):
    min.abs.diff.slope.iteration <- cand.iterations$iteration[cand.iterations$abs.diff.slope == min(cand.iterations$abs.diff.slope)]
  #Iterations(s) with the minimum absolute difference in normal-fitted SD (between the normal-fitted SD of corrected labeled WADs vs normal-fitted SD of unlabeled WADs):
    min.abs.diff.norm.SD.iteration <- cand.iterations$iteration[cand.iterations$abs.diff.norm.SD == min(cand.iterations$abs.diff.norm.SD)]
  #Iterations(s) with the minimum absolute difference in empirical SD (between the empirical SD of corrected labeled WADs vs empirical SD of unlabeled WADs):
    min.abs.diff.empir.SD.iteration <- cand.iterations$iteration[cand.iterations$abs.diff.empir.SD == min(cand.iterations$abs.diff.empir.SD)]
  #Iterations(s) with the minimum absolute difference in empirical mean (between the empirical mean of corrected labeled WADs vs empirical mean of unlabeled WADs):
    min.abs.diff.empir.mean.iteration <- cand.iterations$iteration[cand.iterations$abs.diff.empir.mean == min(cand.iterations$abs.diff.empir.mean)]
  #The single best iteration; chosen by selecting the the most frequently occuring iteration within the specified quantile most frequently across all metrics,
    #and then if there are more than one of those, breaking ties using 'abs.diff.norm.SD' and then 'abs.diff.slope' and then 'abs.diff.empir.mean' and then 'abs.diff.empir.SD':
    best.iteration <- freq.iterations
    if(length(best.iteration) > 1){
      best.iteration <- best.iteration[which(cand.iterations$abs.diff.norm.SD[cand.iterations$iteration %in% best.iteration] == min(cand.iterations$abs.diff.norm.SD[cand.iterations$iteration %in% best.iteration]))]
    }
    if(length(best.iteration) > 1){
      best.iteration <- best.iteration[which(cand.iterations$abs.diff.slope[cand.iterations$iteration %in% best.iteration] == min(cand.iterations$abs.diff.slope[cand.iterations$iteration %in% best.iteration]))]
    }
    if(length(best.iteration) > 1){
      best.iteration <- best.iteration[which(cand.iterations$abs.diff.empir.mean[cand.iterations$iteration %in% best.iteration] == min(cand.iterations$abs.diff.empir.mean[cand.iterations$iteration %in% best.iteration]))]
    }
    if(length(best.iteration) > 1){
      best.iteration <- best.iteration[which(cand.iterations$abs.diff.empir.SD[cand.iterations$iteration %in% best.iteration] == min(cand.iterations$abs.diff.empir.SD[cand.iterations$iteration %in% best.iteration]))]
    }
  #Shift in density associated with the best iteration:
    norm.correction <- cand.iterations$norm.correction[cand.iterations$iteration %in% best.iteration]
  #Return a list of the best iteration along with the other candidate iterations:
    list(best.iteration=best.iteration, norm.correction=norm.correction, freq.iterations=freq.iterations, min.abs.diff.slope.iteration=min.abs.diff.slope.iteration, min.abs.diff.norm.SD.iteration=min.abs.diff.norm.SD.iteration, min.abs.diff.empir.SD.iteration=min.abs.diff.empir.SD.iteration, min.abs.diff.empir.mean.iteration=min.abs.diff.empir.mean.iteration, candidate.iterations=cand.iterations)
}



