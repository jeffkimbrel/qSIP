#' Calculates light isotope tube shift
#'
#' @param data A dataframe containing taxon wad values per tube
#' @param light_isotopes List of light isotopes
#' @export

tube_shift_light = function(data, light_isotopes = c("12C")) {
  require(tidyverse)
  data = data %>%
    filter(Isotope %in% light_isotopes & wads != "NA" & wads != "NaN") %>%
    group_by(taxon) %>%
    mutate(taxon_median_light = median(wads),
           diff_from_median_light = taxon_median_light - wads,
           taxon_uncorrected_stress = diff_from_median_light^2) %>% # mean or median?
    ungroup() %>%
    group_by(SampleID) %>%
    mutate(tube_adjustment = median(diff_from_median_light),
           wads_tube_corrected = wads + tube_adjustment,
           taxon_corrected_stress = (taxon_median_light - wads_tube_corrected)^2,
           tube_rank_wads = percent_rank(wads),
           tube_rank_diff = percent_rank(-diff_from_median_light),
           inactive = TRUE) %>%
    ungroup()

  return(data)
}
