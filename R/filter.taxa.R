#' Filter Taxa
#'
#' Define a function to use to filter the data so that only taxa with an appropriate level of occurrence and replication among the tubes and treatments of the experiment are retained for further analysis:
#'
#' @param DATA DATA
#' @param trt.code.1 trt.code.1
#' @param trt.code.2 trt.code.2
#' @param trt.refs trt.refs
#' @param min.reps min.reps
#'
#' @export

filter.taxa = function(DATA, trt.code.1=NULL, trt.code.2=NULL, trt.refs=NULL, min.reps, ...) {
  require(tidyverse)
  require(patchwork)

  #First, determine the specified treatments according to those specified from the treatment comparisons data.frame:
  trts1 = unlist(str_split(trt.code.1, ";"))
  trts2 = unlist(str_split(trt.code.2, ";"))
  trts3 = unlist(str_split(trt.refs, ";"))
  trts.to.filter <- str_trim(unique(c(trts1, trts2, trts3)))

  #Subset the data into only those taxon-reps with copies present in the specified treatments:
  DATA.occurrences = DATA %>%
    filter(unique.tmt %in% trts.to.filter) %>%
    filter(!is.na(t.DNA.ng.fraction)) %>%
    filter(t.DNA.ng.fraction > 0)

  #Calculate the number of unique tubes in the specified treatments with copies present for each taxon:
  tubes.per.taxon = DATA.occurrences %>%
    group_by(taxon) %>%
    summarize(tubes.per.taxon = n_distinct(unique.tube)) %>%
    arrange(tubes.per.taxon) %>%
    deframe()

  #Visualize keeping only those taxa that occur in at least 'min.reps' tubes across all the specified treatments:
  p1 = tubes.per.taxon %>%
    enframe(name = "taxon", value = "tubes.per.taxon") %>%
    ggplot(aes(x = tubes.per.taxon)) +
    theme_bw() +
    geom_histogram(binwidth = 1, color = "black", aes(fill = ifelse(tubes.per.taxon >= min.reps, "PASS", "FAIL"))) +
    scale_fill_manual(values = c("PASS" = "#037bcf", "FAIL" = "gray80")) +
    theme(legend.position = "bottom") +
    expand_limits(x = 0) +
    labs(x = "# of tubes per taxon", y = "# of taxa", fill = "Meets min.rep requirements")

  p2 = tubes.per.taxon %>%
    enframe(name = "taxon", value = "tubes.per.taxon") %>%
    ggplot(aes(x = reorder(taxon, tubes.per.taxon), y = tubes.per.taxon)) +
    theme_bw() +
    geom_point(pch = 21, aes(fill = ifelse(tubes.per.taxon >= min.reps, "PASS", "FAIL"))) +
    scale_fill_manual(values = c("PASS" = "#037bcf", "FAIL" = "gray80")) +
    theme(legend.position = "none",
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    expand_limits(y = 0) +
    labs(y = "# of tubes per taxon", x = "# of taxa", fill = "Meets min.rep requirements")

  p = p1 / p2

  #Taxa filtered:
  tot.starting.taxa = DATA %>%
    select(taxon) %>%
    unique() %>%
    pull() %>%
    length()

  print("Taxa occurring in (all) the specified treatment(s): ", sep="")
  tubes.per.taxon = tubes.per.taxon
  print(tubes.per.taxon)

  print(paste("Taxa occurring in ≥", min.reps, " total replicates of the specified treatment(s): ", sep=""))
  print(tubes.per.taxon[as.numeric(tubes.per.taxon) >= min.reps])
  print(paste("Total number of taxa occuring in the specified data: ", tot.starting.taxa, sep=""))
  print(paste("Number of taxa occuring in (all) the specified treatment(s): ", length(tubes.per.taxon), sep=""))
  print(paste("Number of taxa that occurred in ≥ ", min.reps, " total replicates of the specified treatment(s): ", length(tubes.per.taxon[as.numeric(tubes.per.taxon) >= min.reps]), sep=""))

  #Now, explore what subsetting the DATA dataframe would mean so that it only contains taxon-tubes for taxa occurring in at least 'min.reps' tubes across the specified treatment(s):
  print(paste("Dimensions of the specified data frame (before filtering): ", paste(dim(DATA), collapse="   "), sep=""))
  print(paste("Dimensions of the specified data frame including only those taxa that do not occur in the specified treatment(s): ", paste(dim(DATA[DATA[,"taxon"] %in% as.character(rep(1:tot.starting.taxa))[!(as.character(rep(1:tot.starting.taxa)) %in% names(tubes.per.taxon))],]), collapse="   "), sep=""))
  print(paste("Dimensions of the specified data frame including only those taxa that do not occur in ≥", min.reps, " total replicates of the specified treatment(s): ", paste(dim(DATA[DATA[,"taxon"] %in% names(tubes.per.taxon[as.numeric(tubes.per.taxon) < min.reps]),]), collapse="   "), sep=""))

  DATA = DATA %>%
    filter(taxon %in% names(tubes.per.taxon[as.numeric(tubes.per.taxon) >= min.reps])) %>%
    mutate(across(where(is.character),as_factor))

  print(paste("Dimensions of the specified data frame (after filtering): ", paste(dim(DATA), collapse="   "), sep=""))

  return(DATA)
}

#' Explore Filter Taxa
#'
#' Define a function to use to filter the data so that only taxa with an appropriate level of occurrence and replication among the tubes and treatments of the experiment are retained for further analysis:
#'
#' @param DATA DATA
#' @param trt.code.1 trt.code.1
#' @param trt.code.2 trt.code.2
#' @param trt.refs trt.refs
#' @param min.reps min.reps
#'
#' @export

explore.filter.taxa = function(DATA, trt.code.1=NULL, trt.code.2=NULL, trt.refs=NULL, min.reps, ...) {
  .Deprecated("filter.taxa()")
}
