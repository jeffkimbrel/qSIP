#Apply the shift to 'correct' the unlabeled data in data (corrects all unlabeled replicates at once):

apply.unlabeled.correction <- function(raw.data, correction.table, reps.by.trt, vars=c("density.g.ml", "tube", "trt.code")){
  raw.data.corr <- raw.data
  unlab.WAD.norm.fit.parms <- correction.table
  for (i in 1:dim(unlab.WAD.norm.fit.parms)[1]){
    indices1 <- as.character(raw.data.corr[,vars[3]]) == unlab.WAD.norm.fit.parms$trt[i] & as.character(raw.data.corr[,vars[2]]) == as.character(reps.by.trt[reps.by.trt[,vars[3]] == unlab.WAD.norm.fit.parms$trt[i], names(reps.by.trt) == unlab.WAD.norm.fit.parms$rep[i]])
    indices2 <- as.character(raw.data[,vars[3]]) == unlab.WAD.norm.fit.parms$trt[i] & as.character(raw.data[,vars[2]]) == as.character(reps.by.trt[reps.by.trt[,vars[3]] == unlab.WAD.norm.fit.parms$trt[i], names(reps.by.trt) == unlab.WAD.norm.fit.parms$rep[i]])
    raw.data.corr[indices1, vars[1]] <- raw.data[indices2, vars[1]] - unlab.WAD.norm.fit.parms$shift.mean[i]
  }
  raw.data.corr
}



