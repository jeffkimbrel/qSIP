#Create a function to extract the sequence of taxa added or removed at each iteration:
  get.seq.taxa.nums <- function(DATA){
    seq.nums <- gsub(pattern="taxon\\.(.*)", replacement="\\1", x=names(DATA), perl=TRUE)
    seq.nums[seq.nums == "none"] <- "0"
    as.numeric(seq.nums)
  }



