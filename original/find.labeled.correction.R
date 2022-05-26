#Function to produce a list of three data frames useful for finding the optimal correction to apply to labeled replicates:
  #(1) putative.nongrower.metrics
  #(2) taxa.in
  #(3) taxa.out
  
  #NOTE: 'lab.names' indicates the names of all replicates in the same labeled treatment as the target replicate specified in 'lab.replicate'
  
find.labeled.correction <- function(LIST, lab.replicate, lab.names, method=c("td.pos.resid", "td.abs.resid", "bu.abs.resid"), unlab.SD.percentile=0.5, lab.SD.percentile=0.5, min.num.nongrowers=10){
  curr.data.corr <- LIST$WAD.table.corr
  unlab.corr.names <- LIST$corr.names
  #Filter the taxa to be included in the algorithm (method) for identifying the set of nongrowers for determining the correction factor for the labeled replicate:
    #Calculate the mean and SD of WAD for each taxon across all unlabeled tubes:
    curr.data.corr$unlab.mean <- apply(curr.data.corr[, names(curr.data.corr) %in% unlab.corr.names], 1, mean)
    curr.data.corr$unlab.SD <- apply(curr.data.corr[, names(curr.data.corr) %in% unlab.corr.names], 1, sd)
    #Calculate the mean and SD of WAD for each taxon across all labeled tubes:
    curr.data.corr$lab.mean <- apply(curr.data.corr[, names(curr.data.corr) %in% lab.names], 1, mean)
    curr.data.corr$lab.SD <- apply(curr.data.corr[, names(curr.data.corr) %in% lab.names], 1, sd)
    #Apply the specified unlabeled percentile to set the threshold SD:
    unlab.SD.threshold <- quantile(curr.data.corr$unlab.SD, unlab.SD.percentile)
    #Exclude all taxa that have SD for WADs across all unlabeled tubes that exceeds the threshold SD (i.e., keep only taxa with relatively low variability in WAD):
    curr.data.lowSD <- curr.data.corr[curr.data.corr$unlab.SD <= unlab.SD.threshold,]
    #Apply the specified labeled percentile to set the threshold SD:
    lab.SD.threshold <- quantile(curr.data.corr$lab.SD, lab.SD.percentile)
    #Exclude all taxa that have SD for WADs across all labeled tubes that exceeds the threshold SD (i.e., keep only taxa with relatively low variability in WAD):
    curr.data.lowSD <- curr.data.lowSD[curr.data.lowSD$lab.SD <= lab.SD.threshold,]
  #Use the specified method to iterate through the remaining taxa to create a data frame of metrics associated with different subsets of putative nongrowing taxa in the labeled treatment for use in assessing a 'match' with unlabeled taxa.
  #And simultaneously create two more data frames: one that lists the taxa included and one the taxa excluded for each iteration
    if(method == "td.pos.resid"){
      LIST <- td.pos.resid(DATA=curr.data.lowSD, lab.replicate=lab.replicate, min.num.nongrowers=min.num.nongrowers)
    } else if(method == "td.abs.resid"){
      LIST <- td.abs.resid(DATA=curr.data.lowSD, lab.replicate=lab.replicate, min.num.nongrowers=min.num.nongrowers)
    } else if(method == "bu.abs.resid"){
      LIST <- bu.abs.resid(DATA=curr.data.lowSD, lab.replicate=lab.replicate, min.num.nongrowers=min.num.nongrowers)
    } else {
      stop("specified method does not exist")
    }
  list(putative.nongrower.metrics=LIST$putative.nongrower.metrics, taxa.in=LIST$taxa.in, taxa.out=LIST$taxa.out, data.low.SD=LIST$data.low.SD, lab.replicate=lab.replicate, unlab.corr.names=unlab.corr.names, lab.names=lab.names, method=method, unlab.SD.percentile=unlab.SD.percentile, lab.SD.percentile=lab.SD.percentile, min.num.nongrowers=min.num.nongrowers)
}



