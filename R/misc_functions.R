#' Parse Treatments
#'
#' A function to parse the treatment string in a qSIP dataset to a vector of
#' treatments with whitespace removed.
#'
#' @param treatment The treatment string
#' @param sep The delimiter
#'
#' @export

parse_treatments = function(treatment, sep = ";") {

  if (is.null(treatment)) {
    return(NULL)
  } else {
    treatments = str_trim(unlist(str_split(treatment, sep)))

    if (is.na(treatments[1])) {
      treatments = NULL
    }
    return(treatments)
  }
}

#' Calculate Tiny Abundance
#'
#' A function to return the adjustments to zero-value data in the specified
#' column.
#'
#' @param df The starting qSIP dataframe
#' @param column The column header to adjust
#' @param multiplier The adjustment multiplier
#'
#' @export

calc_tiny_abundance = function(df, column = "t.copies.ul", multiplier = 0.5) {

  if (multiplier <= 0) {
    stop("multiplier must be positive")
  }

  tryCatch(df %>% filter(!!as.name(column) > 0) %>% pull(column) %>% min() * multiplier,
           error = function(e)
             message(paste(e, "Column either doesn't exist or is not numeric")))

}

#' Create Sample Vector
#'
#' Create a version of 'sample' that avoids the problem when length(x)==1
#'
#' Written by Ben Koch & Natasja van Gestel
#'
#' @param x vector of any length
#'
#' @return vector of length 'size' with elements drawn from x
#'
#' @note forces 'sample' to always treat x as a vector
#' @note (see examples on the ?sample page and also here: http://stackoverflow.com/questions/13990125/sampling-in-r-from-vector-of-varying-length)
#'
#' @export


sample_vec <- function(x, ...){
  x[sample(length(x), ...)]
}

#' Create Sample Vector (Deprecated)
#'
#' @keywords internal
#' @export

sample.vec = function(...) {
  .Defunct("sample_vec()")
}
