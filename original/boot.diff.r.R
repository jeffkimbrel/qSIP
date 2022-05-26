### Given 2 data frames (of 2 treatments) and the incubation period, return an estimate (with uncertainty) of the intrinsic rate of increase (r, units for exponential model: per day; units for linear model: copies per day per g soil) based on the degree of labeling by the heavy isotope
#
#     output = boot.diff.r(X.light, X.heavy, X.reference, M.soil=NULL, iso.compare, days, vars=c("density.g.ml", "copies", "tube", "trt.code", "g.soil"), growth.model="exponential", prop.O.from.water=0.50, v.frac=NULL, CI=0.90, draws=1000)
#
#     X.light: data frame with data for treatment 1, the unlabeled, or 'light', treatment
#     X.heavy: data frame with data for treatment 2, the labeled, or 'heavy', treatment
#     X.reference: data frame with data for the taxon from all treatments without heavy isotopes
#     M.soil: data frame with mass of soil for each replicate; must contain two columns, one for replicate ID and one for mass of soil; default is NULL; M.soil is only required for the linear growth model
#     iso.compare: code specifying which isotopes are being compared; "13C" specifies the comparison between 13C (heavy) vs 12C (light) treatments; "18O" specifies the comparison between 18O (heavy) vs 16O (light) treatments; "15N" specifies the comparison between 15N (heavy) vs 14N (light) treatments; any other value returns an error
#     days: duration of the incubation period in days
#     vars: vector of variable names for the x, y, replicate ID, treatment ID, and soil mass columns in the data frames; default is c("density.g.ml", "copies" ,"tube", "trt.code", "g.soil"); the soil mass column ('g.soil') is only used for the linear growth model
#     vars: vector of variable names for the x, y, replicate ID, and treatment ID columns in the data frames; default is c("density.g.ml", "copies" ,"tube", "trt.code")
#     growth.model: model to use for calculating growth ("linear" or "exponential"); default is "exponential"
#     prop.O.from.water: proportion of oxygen atoms in DNA that come from environmental water; default is 0.50
#     v.frac: volume of a fraction within a tube (uL); default is NULL; only required for linear growth model
#     CI: confidence interval (0-1); default is 0.90
#     draws: number of bootstrap iterations; default is 1000
#     -------------------------------------------------------
#     output: 
#     list of:        boot.r: vector of bootstrapped estimates of intrinsic rate of increase (per day OR copies per day per g soil) based on the two treatments
#                     obs.r: observed intrinsic rate of increase based on the two treatments
#                     boot.r.mean: mean of the bootstrapped intrinsic rate of increase values
#                     boot.r.median: median of the bootstrapped intrinsic rate of increase values
#                     boot.r.CI: upper and lower confidence intervals (defined by argument 'CI') of the bootstrapped intrinsic rate of increase values
#                     message: warning message listing replicates in either of the 2 treatments (does not include the reference group) that had no occurrences (i.e., all y values for a rep are 0); value for message is "none" when all replicates are present (single character value)
#
#     notes:          requires the functions 'WAD.func', 'boot.WAD.func', 'boot.TUBE.func', 'MW.calc', 'r.calc', & 'comparison.message'
#                     names of the x, y, replicate ID, and treatment ID columns in the three data frames (two treatments plus reference) must match exactly
#                     names of the replicate ID columns in the three data frames must also match the name of the replicate ID column in the M.soil data frame
#                     values for replicate ID in the three data frames must be unique among treatments and must match exactly with those in the M.soil data frame
#                     number of reps need not be equal between the two groups
#     assumptions:    assumes either a linear growth model ( Nt=N0+(rt) ) or an exponential growth model ( Nt=N0*e^(rt) )
#                     note that the units of r differ between the linear model (copies * day-1 * g soil-1) and the exponential model (1/day)
#                     assumes a fixed GC content, molecular weight, carbon, and nitrogen content for each taxon based on the mean observed WAD
#                     oxygen content is fixed for DNA regardless of GC content
#                     the weighted average density of unlabeled DNA is a good proxy for GC content of the taxon, according to the well-known relationship between GC content and the density of DNA (Schildkraut et al. 1962)
#                       (NOTE: This is not caused by the effect of GC content on molecular weight, which is negligible. It is caused by differences in the relationship between base composition and binding with water that occurs in CsCl.)
#                     any change in weighted average density of the 16S DNA molecule is caused by incorporation of stable isotopes into the DNA
#                     the relative change in weighted average density is equal to the relative change in molecular weight
#                       (Incorporation of stable isotopes will not alter the binding properties between DNA and water in the CsCl gradient. Those properties are really important for affecting density, but they're a function of GC content, and should not change just because of heavy isotope incorporation.)
#
#     Written by Ben Koch & Natasja van Gestel


boot.diff.r <- function(X.light, X.heavy, X.reference, M.soil=NULL, iso.compare, days, vars=c("density.g.ml", "copies", "tube", "trt.code", "g.soil"), growth.model="exponential", prop.O.from.water=0.50, v.frac=NULL, CI=0.90, draws=1000){
  if (growth.model == "exponential"){
    group1 <- boot.WAD.func(X=X.light, vars=vars[1:3], CI=CI, draws=draws)
    group2 <- boot.WAD.func(X=X.heavy, vars=vars[1:3], CI=CI, draws=draws)
  }
  
  else if (growth.model == "linear" & (is.null(M.soil) | is.null(v.frac))){
    stop("Error: unable to calculate r for a linear growth model without data on total number of copies and volume of fractions\n  Specify a data frame for the 'M.soil' argument and a value for the 'v.frac' argument")
  }

  else if (growth.model == "linear" & !is.null(M.soil) & !is.null(v.frac)){
    group1 <- boot.TUBE.func(X=X.light, M.soil=M.soil, vars=vars[c(1:3,5)], v.frac=v.frac, CI=CI, draws=draws)
    group2 <- boot.TUBE.func(X=X.heavy, M.soil=M.soil, vars=vars[c(1:3,5)], v.frac=v.frac, CI=CI, draws=draws)
  }

  MW.output <- MW.calc(X.reference, vars=c(vars[1], vars[2], vars[3]))

  r.out <- r.calc(X.light=X.light, X.heavy=X.heavy, boot.out.light=group1, boot.out.heavy=group2, MW.out=MW.output, iso.compare=iso.compare, days=days, var=vars[4], growth.model=growth.model, prop.O.from.water=prop.O.from.water, CI=CI)
  r.out
}
