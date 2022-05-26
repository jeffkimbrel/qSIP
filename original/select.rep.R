#Function to isolate the density and number of 16S copies associated  with a given replicate of a given treatment:

select.rep <- function(DATA, focal.tmt, replicate.index, vars=c("density.g.ml", "copies", "tube", "trt.code")){
  focal.reps <- id.reps(DATA=DATA, focal.tmt=focal.tmt, vars=c(vars[3], vars[4]))
  indices <- DATA[,vars[3]] %in% focal.reps[replicate.index]  
  X <- DATA[indices, vars[1]]
  Y <- DATA[indices, vars[2]]
  list(focal.tmt=focal.tmt, focal.rep=focal.reps[replicate.index], X=X, Y=Y)
}



