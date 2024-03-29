% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stress_test.R
\name{stress_test}
\alias{stress_test}
\title{Calculates and plots stress
Window and step size can be specified by user}
\usage{
stress_test(data, light, window = 0.1, step = 0.1)
}
\arguments{
\item{data}{A dataframe containing taxon wad values per tube}

\item{light}{A dataframe containing light isotope adjusted taxon wad values averaged across all tubes}

\item{window}{Size of the sliding window}

\item{step}{Step size for the sliding window}
}
\description{
Calculates and plots stress
Window and step size can be specified by user
}
