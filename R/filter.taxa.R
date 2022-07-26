#' Filter Taxa
#'
#' A function to filter the data so that only taxa with an appropriate level of
#' occurrence and replication among the tubes and treatments of the experiment
#' are retained for further analysis.
#'
#' @param df Dataframe
#' @param treatment1 trt.code.1
#' @param treatment2 trt.code.2
#' @param treatment_refs trt.refs
#' @param min_reps min.reps
#' @param taxon_column Column header with taxon identifier (formerly vars[1])
#' @param filter_column Column header with abundance or copy numbers (formerly vars[2])
#' @param tube_column Column header with tube or sample identifier (formerly vars[3])
#' @param treatment_column Column header with treatment identifier (formerly vars[4])
#'
#' @return A list with filtered data and extras
#' \itemize{
#'   \item df: The filtered dataframe
#'   \item plot: A summary ggplot object
#'   \item treatment1: The argument passed to treatment1
#'   \item treatment2: The argument passed to treatment2
#'   \item treatment_refs: The argument passed to treatment_refs
#'   \item min_reps: The argument passed to min_reps
#'   \item taxon_column: The argument passed to taxon_column
#'   \item treatment_column: The argument passed to treatment_column
#'   \item filter_column: The argument passed to filter_column
#'   \item tube_column: The argument passed to tube_column
#' }
#'
#' @export

filter_taxa = function(df, treatment1=NULL, treatment2=NULL, treatment_refs=NULL, min_reps = 3, taxon_column = "taxon", treatment_column = "unique.tmt", filter_column = "t.copies.ul", tube_column = "unique.tube") {
  require(tidyverse)

  #First, determine the specified treatments according to those specified from the treatment comparisons data.frame:
  trts1 = parse_treatments(treatment1)
  trts2 = parse_treatments(treatment2)
  trts3 = parse_treatments(treatment_refs)
  treatments_to_filter <- c(trts1, trts2, trts3)

  #Subset the data into only those taxon-reps with copies present in the specified treatments:
  df_subset = df %>%
    filter(!!as.name(treatment_column) %in% treatments_to_filter) %>%
    filter(!is.na(!!as.name(filter_column))) %>%
    filter(!!as.name(filter_column) > 0)

  #Calculate the number of unique tubes in the specified treatments with copies present for each taxon:
  tubes_per_taxon = df_subset %>%
    group_by(!!as.name(taxon_column)) %>%
    summarize(tubes_per_taxon = n_distinct(!!as.name(tube_column))) %>%
    arrange(tubes_per_taxon) %>%
    deframe()

  #Visualize keeping only those taxa that occur in at least 'min_reps' tubes across all the specified treatments:
  p = tubes_per_taxon %>%
    enframe(name = "taxon", value = "tubes_per_taxon") %>%
    ggplot(aes(x = tubes_per_taxon)) +
    theme_bw() +
    geom_histogram(binwidth = 1, color = "black", aes(fill = ifelse(tubes_per_taxon >= min_reps, "PASS", "FAIL"))) +
    scale_fill_manual(values = c("PASS" = "#037bcf", "FAIL" = "gray80")) +
    theme(legend.position = "bottom") +
    expand_limits(x = 0) +
    labs(title = paste("Filtered Treatments:", treatments_to_filter, ", Min Reps:", min_reps ), x = "# of tubes per taxon", y = "# of taxa", fill = "Meets min.rep requirements")

  #Taxa filtered:
  tot_starting_taxa = df %>%
    select(taxon) %>%
    unique() %>%
    pull() %>%
    length()

  taxa_passing_filter = tubes_per_taxon %>%
    enframe(name = "taxon", value = "tubes_per_taxon") %>%
    filter(tubes_per_taxon >= min_reps) %>%
    pull(taxon)

  message(paste("Total number of taxa occuring in the specified data: ", tot_starting_taxa, sep=""))
  message(paste("Number of taxa occuring in (all) the specified treatment(s): ", length(tubes_per_taxon), sep=""))
  message(paste("Number of taxa that occurred in ≥ ", min_reps, " total replicates of the specified treatment(s): ", length(taxa_passing_filter), sep=""))

  df_filtered = df %>%
    filter(taxon %in% taxa_passing_filter) %>%
    mutate(across(where(is.character),as_factor))

  message(paste("Dimensions of the specified data frame (before filtering): ", paste(dim(df), collapse="   "), sep=""))
  message(paste("Dimensions of the specified data frame (after filtering): ", paste(dim(df_filtered), collapse="   "), sep=""))

  l = list("df" = df_filtered,
           "plot" = p,
           "taxa" = taxa_passing_filter,
           "treatment1" = treatment1,
           "treatment2" = treatment2,
           "treatment_refs" = treatment_refs,
           "min_reps" = min_reps,
           "taxon_column" = taxon_column,
           "treatment_column" = treatment_column,
           "filter_column" = filter_column,
           "tube_column" = tube_column)

  return(l)
}

#' Explore Filter Taxa (Deprecated)
#'
#' @keywords internal
#' @export

explore.filter.taxa = function(...) {
  .Defunct("filter_taxa()")
}

#' Filter Taxa (Deprecated)
#'
#' @keywords internal
#' @export

filter.taxa = function(...) {
  .Defunct("filter_taxa()")
}
