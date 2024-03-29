% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wad_func.R
\name{WAD_func}
\alias{WAD_func}
\title{Calculate Weighted Average Density (WAD)}
\usage{
WAD_func(y, x)
}
\arguments{
\item{y}{vector of y-values (e.g., number of 16S copies)}

\item{x}{vector of x-values (e.g., density of DNA)}
}
\value{
weighted-average of the x-values (single value)
}
\description{
Given vectors of x values and y values, calculate the weighted-average of
the x-values (e.g., the weighted average density (WAD))
}
\examples{
WAD_func(c(1, 10, 100), c(1, 2, 3))
}
