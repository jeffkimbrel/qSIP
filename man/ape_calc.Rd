% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ape_calc.R
\name{ape_calc}
\alias{ape_calc}
\title{Calculate APE}
\usage{
ape_calc(
  X.light,
  X.heavy,
  boot.out.light,
  boot.out.heavy,
  MW.out,
  iso.compare,
  var = "trt.code",
  CI = 0.9
)
}
\arguments{
\item{X.light}{data frame with data for treatment 1, the unlabeled, or 'light', treatment}

\item{X.heavy}{data frame with data for treatment 2, the labeled, or 'heavy', treatment}

\item{boot.out.light}{list containing output from boot.WAD.func or boot.TUBE.func for the 'light' treatment (corresponding to X.light)}

\item{boot.out.heavy}{list containing output from boot.WAD.func or boot.TUBE.func for the 'heavy' treatment (corresponding to X.heavy)}

\item{MW.out}{data frame containing output from MW.calc for the unlabeled treatment(s) for a taxon}

\item{iso.compare}{code specifying which isotopes are being compared; "13C" specifies the comparison between 13C (heavy) vs 12C (light) treatments; "18O" specifies the comparison between 18O (heavy) vs 16O (light) treatments; "15N" specifies the comparison between 15N (heavy) vs 14N (light) treatments; any other value returns an error}

\item{var}{variable name for the treatment ID columns in the data frames}

\item{CI}{confidence interval (0-1)}
}
\value{
A list with
\itemize{
  \item boot.apes: vector of bootstrapped estimates of excess atom fraction (unitless proportion) based on the two treatments
  \item obs.ape: observed excess atom fraction based on the two treatments
  \item boot.apes.mean: mean of the bootstrapped excess atom fraction values
  \item boot.apes.median: median of the bootstrapped excess atom fraction values
  \item boot.apes.CI: upper and lower confidence intervals (defined by argument 'CI') of the bootstrapped excess atom fraction values
  \item message: warning message listing replicates in either of the 2 treatments (does not include the reference group) that had no occurrences (i.e., all y values for a rep are 0); value for message is "none" when all replicates are present (single character value)
}
}
\description{
Given two treatments, returns an estimate (with uncertainty) of excess atom
fraction (or atom percent excess) of the heavy isotope.
}
\details{
Assumptions:
* assumes a fixed GC content, molecular weight, carbon, and nitrogen content for each taxon based on the mean observed WAD
* oxygen content is fixed for DNA regardless of GC content
* the weighted average density of unlabeled DNA is a good proxy for GC content of the taxon, according to the well-known relationship between GC content and the density of DNA (Schildkraut et al. 1962)
* (NOTE: This is not caused by the effect of GC content on molecular weight, which is negligible. It is caused by differences in the relationship between base composition and binding with water that occurs in CsCl.)
* any change in weighted average density of the 16S DNA molecule is caused by incorporation of stable isotopes into the DNA
* the relative change in weighted average density is equal to the relative change in molecular weight
* (Incorporation of stable isotopes will not alter the binding properties between DNA and water in the CsCl gradient. Those properties are really important for affecting density, but they're a function of GC content, and should not change just because of heavy isotope incorporation.)

Written by Ben Koch & Natasja van Gestel
}
\note{
requires the functions 'WAD.func', 'boot.WAD.func', 'MW.calc', & 'comparison.message'

number of reps need not be equal between the two groups
}
