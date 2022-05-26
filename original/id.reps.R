#Function to identify the names of all replicates associated with a given treatment:

id.reps <- function(DATA, focal.tmt, vars=c("tube", "trt.code")){
  focal.reps <- as.character(unique(DATA[DATA[,vars[2]] == focal.tmt, vars[1]]))
  focal.reps  
}



