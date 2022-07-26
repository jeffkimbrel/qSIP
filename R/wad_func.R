#' Calculate Weighted Average Density (WAD)
#'
#' Given vectors of x values and y values, calculate the weighted-average of
#' the x-values (e.g., the weighted average density (WAD))
#'
#' Written by Ben Koch & Natasja van Gestel
#'
#' @param y vector of y-values (e.g., number of 16S copies)
#' @param x vector of x-values (e.g., density of DNA)
#'
#' @return weighted-average of the x-values (single value)
#'
#' @export

WAD_func <- function(y, x){
  sum(x*(y/sum(y)))
}

#' Calculate Weighted Average Density (WAD) (Deprecated)
#'
#' @keywords internal
#' @export

WAD.func = function(...) {
  .Defunct("WAD_func()")
}
