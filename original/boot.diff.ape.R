### Given 2 data frames (of 2 treatments), return an estimate (with uncertainty) of excess atom fraction (or atom percent excess) of the heavy isotope
#
#     output = boot.diff.ape(X.light, X.heavy, X.reference, iso.compare, vars=c("density.g.ml", "copies", "tube", "trt.code"), CI=0.90, draws=1000)
#
#     X.light: data frame with data for treatment 1, the unlabeled, or 'light', treatment
#     X.heavy: data frame with data for treatment 2, the labeled, or 'heavy', treatment
#     X.reference: data frame with data for the taxon from all treatments without heavy isotopes
#     iso.compare: code specifying which isotopes are being compared; "13C" specifies the comparison between 13C (heavy) vs 12C (light) treatments; "18O" specifies the comparison between 18O (heavy) vs 16O (light) treatments; "15N" specifies the comparison between 15N (heavy) vs 14N (light) treatments; any other value returns an error
#     vars: vector of variable names for the x, y, replicate ID, and treatment ID columns in the data frames; default is c("density.g.ml", "copies" ,"tube", "trt.code")
#     CI: confidence interval (0-1); default is 0.90
#     draws: number of bootstrap iterations; default is 1000
#     -------------------------------------------------------
#     output: 
#     list of:        boot.apes: vector of bootstrapped estimates of excess atom fraction (unitless proportion) based on the two treatments
#                     obs.ape: observed excess atom fraction based on the two treatments
#                     boot.apes.mean: mean of the bootstrapped excess atom fraction values
#                     boot.apes.median: median of the bootstrapped excess atom fraction values
#                     boot.apes.CI: upper and lower confidence intervals (defined by argument 'CI') of the bootstrapped excess atom fraction values
#                     message: warning message listing replicates in either of the 2 treatments (does not include the reference group) that had no occurrences (i.e., all y values for a rep are 0); value for message is "none" when all replicates are present (single character value)
#
#     notes:          requires the functions 'WAD.func', 'boot.WAD.func', 'MW.calc', 'ape.calc', & 'comparison.message'
#                     names of the x, y, replicate ID, and treatment ID columns in the three data frames (two treatments plus reference) must match exactly
#                     number of reps need not be equal between the two groups
#     assumptions:    assumes a fixed GC content, molecular weight, carbon, and nitrogen content for each taxon based on the mean observed WAD
#                     oxygen content is fixed for DNA regardless of GC content
#                     the weighted average density of unlabeled DNA is a good proxy for GC content of the taxon, according to the well-known relationship between GC content and the density of DNA (Schildkraut et al. 1962)
#                       (NOTE: This is not caused by the effect of GC content on molecular weight, which is negligible. It is caused by differences in the relationship between base composition and binding with water that occurs in CsCl.)
#                     any change in weighted average density of the 16S DNA molecule is caused by incorporation of stable isotopes into the DNA
#                     the relative change in weighted average density is equal to the relative change in molecular weight
#                       (Incorporation of stable isotopes will not alter the binding properties between DNA and water in the CsCl gradient. Those properties are really important for affecting density, but they're a function of GC content, and should not change just because of heavy isotope incorporation.)
#
#     Written by Ben Koch & Natasja van Gestel


boot.diff.ape <- function(X.light, X.heavy, X.reference, iso.compare, vars=c("density.g.ml", "copies", "tube", "trt.code"), CI=0.90, draws=1000){
  group1 <- boot.WAD.func(X=X.light, vars=vars[1:3], CI=CI, draws=draws)
  group2 <- boot.WAD.func(X=X.heavy, vars=vars[1:3], CI=CI, draws=draws)
  MW.output <- MW.calc(X.reference, vars=c(vars[1], vars[2], vars[3]))
  
  ape.out <- ape.calc(X.light=X.light, X.heavy=X.heavy, boot.out.light=group1, boot.out.heavy=group2, MW.out=MW.output, iso.compare=iso.compare, var=vars[4], CI=CI)
  ape.out
}
