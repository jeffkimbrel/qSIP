#' Calculate Weighted Average Density (WAD)
#'
#' Given vectors of x values and y values, calculate the weighted-average of
#' the x-values (e.g., the weighted average density (WAD))
#'
#'
#' @param y vector of y-values (e.g., number of 16S copies)
#' @param x vector of x-values (e.g., density of DNA)
#'
#' @return weighted-average of the x-values (single value)
#'
#' @export
#' @examples
#' WAD_func(c(1, 10, 100), c(1, 2, 3))

WAD_func <- function(y, x) {
  if (length(x) != length(y)) {
    stop("x and y are different lengths")
  }

  sum(x * (y / sum(y)))
}

#' Calculate Weighted Average Density (WAD) (Deprecated)
#'
#' @keywords internal
#' @export

WAD.func = function(...) {
  .Defunct("WAD_func()")
}
