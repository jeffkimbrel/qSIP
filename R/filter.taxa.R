#' Explore Filter Taxa
#'
#' Define a function to use to filter the data so that only taxa with an appropriate level of occurrence and replication among the tubes and treatments of the experiment are retained for further analysis:
#'
#' @param DATA DATA
#' @param trt.code.1 trt.code.1
#' @param trt.code.2 trt.code.2
#' @param trt.refs trt.refs
#' @param vars vars
#' @param min.reps min.reps
#'
#' @export

explore.filter.taxa <- function(DATA, trt.code.1=NULL, trt.code.2=NULL, trt.refs=NULL, vars=c("taxon", "copies", "tube", "trt.code"), min.reps){

  #First, determine the specified treatments according to those specified from the treatment comparisons data.frame:
  pattern <- "([[:space:]]*(,|;)[[:space:]]*)|([[:space:]]+)"
  m1 <- gregexpr(pattern, trt.code.1, perl=TRUE)
  m2 <- gregexpr(pattern, trt.code.2, perl=TRUE)
  m3 <- gregexpr(pattern, trt.refs, perl=TRUE)
  if(length(m1) == 0){
    trts1 <- NULL
  }  else  trts1 <- sort(unlist(regmatches(trt.code.1, m1, invert=TRUE)[[1]]))
  if(length(m2) == 0){
    trts2 <- NULL
  }  else  trts2 <- sort(unlist(regmatches(trt.code.2, m2, invert=TRUE)[[1]]))
  if(length(m3) == 0){
    trts3 <- NULL
  }  else  trts3 <- sort(unlist(regmatches(trt.refs, m3, invert=TRUE)[[1]]))
  trts.to.filter <- unique(c(trts1, trts2, trts3))

  #Subset the data into only those taxon-reps with copies present in the specified treatments:
  DATA.occurrences <- DATA[!is.na(DATA[,vars[2]]) & DATA[,vars[2]] > 0 & DATA[,vars[4]] %in% trts.to.filter,]
  #Renumber row.names:
  row.names(DATA.occurrences) <- 1:dim(DATA.occurrences)[1]
  #Convert the factor columns to factor:
  DATA.occurrences <- as.data.frame(lapply(DATA.occurrences, function(x) if(is.factor(x)) factor(x) else x))

  #Calculate the number of unique tubes in the specified treatments with copies present for each taxon:
  tubes.per.taxon <- tapply(DATA.occurrences[,vars[3]], DATA.occurrences[,vars[1]], function(x) length(unique(x)))
  tubes.per.taxon <- sort(tubes.per.taxon)

  #Visualize keeping only those taxa that occur in at least 'min.reps' tubes across all the specified treatments:
  dev.new(width=3.5, height=7)
  par(mfrow=c(2,1))
  hist(tubes.per.taxon, xlim=c(0,max(as.numeric(tubes.per.taxon))), xlab="# of tubes per taxon", ylab="# of taxa", main="")
  abline(v=min.reps-0.1, col="red")
  plot(x=1:length(tubes.per.taxon), y=as.numeric(tubes.per.taxon), ylim=c(0,max(as.numeric(tubes.per.taxon))), xlab="# of taxa", ylab="# of tubes per taxon", type="l")
  abline(h=min.reps-0.1, col="red")
  par(mfrow=c(1,1))

  #Taxa filtered:
  tot.starting.taxa <- length(levels(factor(DATA[,vars[1]])))
  print("Taxa occurring in (all) the specified treatment(s): ", sep="")
  print(tubes.per.taxon[as.numeric(tubes.per.taxon)])
  print(paste("Taxa occurring in ≥", min.reps, " total replicates of the specified treatment(s): ", sep=""))
  print(tubes.per.taxon[as.numeric(tubes.per.taxon) >= min.reps])
  print(paste("Total number of taxa occuring in the specified data: ", tot.starting.taxa, sep=""))
  print(paste("Number of taxa occuring in (all) the specified treatment(s): ", length(tubes.per.taxon), sep=""))
  print(paste("Number of taxa that occurred in ≥ ", min.reps, " total replicates of the specified treatment(s): ", length(tubes.per.taxon[as.numeric(tubes.per.taxon) >= min.reps]), sep=""))

  #Now, explore what subsetting the DATA dataframe would mean so that it only contains taxon-tubes for taxa occurring in at least 'min.reps' tubes across the specified treatment(s):
  print(paste("Dimensions of the specified data frame (before filtering): ", paste(dim(DATA), collapse="   "), sep=""))
  print(paste("Dimensions of the specified data frame including only those taxa that do not occur in the specified treatment(s): ", paste(dim(DATA[DATA[,vars[1]] %in% as.character(rep(1:tot.starting.taxa))[!(as.character(rep(1:tot.starting.taxa)) %in% names(tubes.per.taxon))],]), collapse="   "), sep=""))
  print(paste("Dimensions of the specified data frame including only those taxa that do not occur in ≥", min.reps, " total replicates of the specified treatment(s): ", paste(dim(DATA[DATA[,vars[1]] %in% names(tubes.per.taxon[as.numeric(tubes.per.taxon) < min.reps]),]), collapse="   "), sep=""))
  DATA <- DATA[DATA[,vars[1]] %in% names(tubes.per.taxon[as.numeric(tubes.per.taxon) >= min.reps]),]
  #Renumber row.names:
  row.names(DATA) <- 1:dim(DATA)[1]
  #Convert the factor columns to factor:
  DATA <- as.data.frame(lapply(DATA, function(x) if(is.factor(x)) factor(x) else x))
  print(paste("Dimensions of the specified data frame (after filtering): ", paste(dim(DATA), collapse="   "), sep=""))
}
