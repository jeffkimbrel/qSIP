#Apply the shift to 'correct' the labeled data in data.corr (corrects one labeled replicate at a time):


apply.labeled.correction <- function(raw.data, raw.data.corr, lab.replicate, correction.value, reps.by.trt, vars=c("density.g.ml", "tube", "trt.code")){
  trt.name <- gsub(pattern="^(R\\d+)\\.(.+)$", replacement="\\2", x=lab.replicate, perl=TRUE)   #pull out just the treatment name
  rep.name <- gsub(pattern="^(R\\d+)\\.(.+)$", replacement="\\1", x=lab.replicate, perl=TRUE)   #pull out just the replicate name
  indices1 <- as.character(raw.data.corr[,vars[3]]) == trt.name & as.character(raw.data.corr[,vars[2]]) == as.character(reps.by.trt[reps.by.trt[,vars[3]] == trt.name, names(reps.by.trt) == rep.name])
  indices2 <- as.character(raw.data[,vars[3]]) == trt.name & as.character(raw.data[,vars[2]]) == as.character(reps.by.trt[reps.by.trt[,vars[3]] == trt.name, names(reps.by.trt) == rep.name])
  raw.data.corr[indices1, vars[1]] <- raw.data[indices2, vars[1]] - correction.value
  raw.data.corr
}



