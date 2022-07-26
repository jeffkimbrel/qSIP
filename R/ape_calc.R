#' Calculate APE
#'
#' Given two treatments, returns an estimate (with uncertainty) of excess atom
#' fraction (or atom percent excess) of the heavy isotope.
#'
#' @details Assumptions:
#' * assumes a fixed GC content, molecular weight, carbon, and nitrogen content for each taxon based on the mean observed WAD
#' * oxygen content is fixed for DNA regardless of GC content
#' * the weighted average density of unlabeled DNA is a good proxy for GC content of the taxon, according to the well-known relationship between GC content and the density of DNA (Schildkraut et al. 1962)
#' * (NOTE: This is not caused by the effect of GC content on molecular weight, which is negligible. It is caused by differences in the relationship between base composition and binding with water that occurs in CsCl.)
#' * any change in weighted average density of the 16S DNA molecule is caused by incorporation of stable isotopes into the DNA
#' * the relative change in weighted average density is equal to the relative change in molecular weight
#' * (Incorporation of stable isotopes will not alter the binding properties between DNA and water in the CsCl gradient. Those properties are really important for affecting density, but they're a function of GC content, and should not change just because of heavy isotope incorporation.)
#'
#' Written by Ben Koch & Natasja van Gestel
#'
#' @param X.light data frame with data for treatment 1, the unlabeled, or 'light', treatment
#' @param X.heavy data frame with data for treatment 2, the labeled, or 'heavy', treatment
#' @param boot.out.light list containing output from boot.WAD.func or boot.TUBE.func for the 'light' treatment (corresponding to X.light)
#' @param boot.out.heavy list containing output from boot.WAD.func or boot.TUBE.func for the 'heavy' treatment (corresponding to X.heavy)
#' @param MW.out data frame containing output from MW.calc for the unlabeled treatment(s) for a taxon
#' @param iso.compare code specifying which isotopes are being compared; "13C" specifies the comparison between 13C (heavy) vs 12C (light) treatments; "18O" specifies the comparison between 18O (heavy) vs 16O (light) treatments; "15N" specifies the comparison between 15N (heavy) vs 14N (light) treatments; any other value returns an error
#' @param var variable name for the treatment ID columns in the data frames
#' @param CI confidence interval (0-1)
#'
#' @return A list with
#' \itemize{
#'   \item boot.apes: vector of bootstrapped estimates of excess atom fraction (unitless proportion) based on the two treatments
#'   \item obs.ape: observed excess atom fraction based on the two treatments
#'   \item boot.apes.mean: mean of the bootstrapped excess atom fraction values
#'   \item boot.apes.median: median of the bootstrapped excess atom fraction values
#'   \item boot.apes.CI: upper and lower confidence intervals (defined by argument 'CI') of the bootstrapped excess atom fraction values
#'   \item message: warning message listing replicates in either of the 2 treatments (does not include the reference group) that had no occurrences (i.e., all y values for a rep are 0); value for message is "none" when all replicates are present (single character value)
#' }
#'
#' @note requires the functions 'WAD.func', 'boot.WAD.func', 'MW.calc', & 'comparison.message'
#' @note number of reps need not be equal between the two groups
#'
#' @export

ape_calc <- function(X.light, X.heavy, boot.out.light, boot.out.heavy, MW.out, iso.compare, var="trt.code", CI=0.90){
  obs.labeled.MW <- MW.out$MW * (boot.out.heavy$obs.wad.mean/boot.out.light$obs.wad.mean)
  boot.labeled.MW <- MW.out$MW * (boot.out.heavy$boot.wads/boot.out.light$boot.wads)

  if (is.na(iso.compare)){                                            #if the iso.compare code is 'NA', print an error message
    stop(paste("Error: unable to calculate ape for the specified isotope comparison '", iso.compare, "'", sep=""))
  } else if (iso.compare == "13C"){                                     #if the treatments being compared are 13C (heavy) vs 12C (light) do the 13-carbon calculation...
    ape.of.fully.labeled.13C <- 1-(11237.2/(1000000+11237.2))         #atom percent excess of a substance with 100% 13C atoms (relative to VPDB; see IAEA 1995, Werner & Brand 2001)
    heavy.max.MW <- MW.out$MW + (MW.out$Catoms * (1.008665 * (1000000/(1000000+11237.2))))    #assumes unlabeled DNA already contains a minute amount of 13C (at the natural abundance level of VPDB)
    obs.ape <- (obs.labeled.MW-MW.out$MW)/(heavy.max.MW-MW.out$MW) * ape.of.fully.labeled.13C
    boot.apes <- (boot.labeled.MW-MW.out$MW)/(heavy.max.MW-MW.out$MW) * ape.of.fully.labeled.13C
  } else if (iso.compare == "18O"){                                     #if the treatments being compared are 18O (heavy) vs 16O (light) do the 18-oxygen calculation...
    ape.of.fully.labeled.18O <- 1-(2005.20/(1000000+379.9+2005.20))   #atom percent excess of a substance with 100% 18O atoms (relative to VSMOW; see IAEA 1995, Werner & Brand 2001)
    heavy.max.MW <- MW.out$MW + (MW.out$Oatoms * ((1.008665 * 2 * (1000000/(1000000+379.9+2005.20))) + (1.008665 * 1 * (379.9/(1000000+379.9+2005.20)))))    #assumes unlabeled DNA already contains minute amounts of 18O and 17O (at the natural abundance levels of those isotopes in VSMOW)
    obs.ape <- (obs.labeled.MW-MW.out$MW)/(heavy.max.MW-MW.out$MW) * ape.of.fully.labeled.18O
    boot.apes <- (boot.labeled.MW-MW.out$MW)/(heavy.max.MW-MW.out$MW) * ape.of.fully.labeled.18O
  } else if (iso.compare == "15N"){                                                                 #if the treatments being compared are 15N (heavy) vs 14N (light) do the 15-nitrogen calculation...
    ape.of.fully.labeled.15N <- 1-((1000000/272)/(1000000+(1000000/272)))                         #atom percent excess of a substance with 100% 15N atoms (relative to AIR-N2; see IAEA 1995, de Laeter et al. 2003; accepted 14N/15N ratio = 272)
    heavy.max.MW <- MW.out$MW + (MW.out$Natoms * (1.008665 * (1000000/(1000000+(1000000/272)))))  #assumes unlabeled DNA already contains minute amounts of 15N (at the natural abundance level of AIR-N2)
    obs.ape <- (obs.labeled.MW-MW.out$MW)/(heavy.max.MW-MW.out$MW) * ape.of.fully.labeled.15N
    boot.apes <- (boot.labeled.MW-MW.out$MW)/(heavy.max.MW-MW.out$MW) * ape.of.fully.labeled.15N
  } else {                                                              #if the iso.compare code is anything else, print an error message
    stop(paste("Error: unable to calculate ape for the specified isotope comparison '", iso.compare, "'", sep=""))
  }

  boot.apes.CI <- quantile(boot.apes, probs=c((1-CI)/2, 1-((1-CI)/2)), na.rm=T)

  message <- comparison_message(X.light=X.light, X.heavy=X.heavy, boot.out.light=boot.out.light, boot.out.heavy=boot.out.heavy, var=var)

  # Collect all output into a list:
  return(list(boot.apes=boot.apes,
              obs.ape=obs.ape,
              boot.apes.mean=mean(boot.apes, na.rm=T),
              boot.apes.median=median(boot.apes, na.rm=T),
              boot.apes.CI=boot.apes.CI,
              message=message))
}



#' Calculate APE (Deprecated)
#'
#' @keywords internal
#' @export

ape.calc = function(...) {
  .Defunct("ape_calc()")
}

