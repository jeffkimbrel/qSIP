#' Calculates and plots stress
#' Window and step size can be specified by user
#'
#' @param data A dataframe containing taxon wad values per tube
#' @param light A dataframe containing light isotope adjusted taxon wad values averaged across all tubes
#' @param window Size of the sliding window
#' @param step Step size for the sliding window
#' @export

stress_test = function(data, light, window = 0.1, step = 0.1){

  f = seq(0, 1 - window, by = step)
  df = data.frame()

  for (i in f) {
    j = i + window
    S = tube_shift_heavy(data = data, light = light, high = j, low = i) %>%
      filter(inactive == TRUE) %>%
      summarize(S = mean(taxon_corrected_stress)) %>%
      pull()
    df = rbind(df, data.frame("low" = i, "high" = j, "stress" = S))

  }
  df = arrange(df, stress)
  message("lowest stress window ", df[1,1], "-", df[1,2])
  return(df)
}

#' Stress test plot
#' @param stress_test_results Dataframe output of stress_test function
#' @export

stress_test_plot = function(stress_test_results){
    stress_test_results %>%
    ggplot(aes(x = low, y = stress)) +
      geom_point()
}
