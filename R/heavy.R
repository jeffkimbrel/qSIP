#' Calculates heavy isotope tube shift
#' Adds light and heavy tube shifts together
#'
#' @param data A dataframe containing taxon wad values per tube
#' @param light A dataframe containing light isotope adjusted taxon wad values averaged across all tubes
#' @param heavy_isotopes List of heavy isotopes
#' @param low Low
#' @param high High
#'
#' @export

tube_shift_heavy = function(data, light, heavy_isotopes = c("13C"), low = 0, high = 0.1) {
  require(tidyverse)

  taxon_median_light = light %>%
    select(taxon, taxon_median_light) %>%
    unique()

  rowcount = nrow(data)

  data = data %>%
    filter(Isotope %in% heavy_isotopes & wads != "NA") %>%
    left_join(taxon_median_light, by = "taxon") %>%
    filter(taxon_median_light != "NA") %>%
    group_by(SampleID) %>%
    mutate(diff_from_median_light = taxon_median_light - wads, # if wad density is below standard, diff will be positive number
           taxon_uncorrected_stress = diff_from_median_light^2,
           tube_rank_wads = percent_rank(wads),
           tube_rank_diff = percent_rank(-diff_from_median_light),
           inactive = case_when(tube_rank_wads < low ~ FALSE,
                                tube_rank_wads > high ~ FALSE,
                                tube_rank_wads >= low & tube_rank_wads <= high ~ TRUE)) %>%
    # inactive = case_when(tube_rank_diff < low ~ FALSE,
    #               tube_rank_diff > high ~ FALSE,
    #               tube_rank_diff >= low & tube_rank_diff <= high ~ TRUE)) %>%
    # inactive = ifelse(tube_rank_wads <= fraction, TRUE, FALSE)) %>%
    ungroup()

  inactive = data %>%
    filter(inactive == TRUE) %>%
    group_by(SampleID) %>%
    summarize(tube_adjustment = median(diff_from_median_light)) %>%
  # summarize(tube_adjustment = median(taxon_median_light - wads))
  ungroup()

  data = data %>%
    left_join(inactive, by = "SampleID") %>%
    mutate(wads_tube_corrected = wads + tube_adjustment,
           taxon_corrected_stress = (taxon_median_light - wads_tube_corrected)^2) %>% # mean or median?
    ungroup()

  data = bind_rows(light, data)
  message(rowcount - nrow(data), " ASVs were removed from original dataframe")
  return(data)
}
