### Returns a warning message about which replicate(s) are missing data when comparing two treatments
#
#     output = comparison.message(X.light, X.heavy, boot.out.light, boot.out.heavy, var="trt.code")
#
#     X.light: data frame with data for treatment 1, the unlabeled, or 'light', treatment
#     X.heavy: data frame with data for treatment 2, the labeled, or 'heavy', treatment
#     boot.out.light: list containing output from boot.WAD.func or boot.TUBE.func for the 'light' treatment (corresponding to X.light)
#     boot.out.heavy: list containing output from boot.WAD.func or boot.TUBE.func for the 'heavy' treatment (corresponding to X.heavy)
#     var: variable name for the treatment ID columns in the data frames; default is "trt.code"
#     -------------------------------------------------------
#     output: 
#               message: warning message listing replicates in either of the 2 treatments (does not include the reference group) that had no occurrences (i.e., all y values for a rep are 0); value for message is "none" when all replicates are present (single character value)
#
#     notes:    number of reps need not be equal between the two groups
#
#     Written by Ben Koch & Natasja van Gestel


comparison.message <- function(X.light, X.heavy, boot.out.light, boot.out.heavy, var="trt.code"){
  g1.reps.NAs <- boot.out.light$obs.wads$rep[is.na(boot.out.light$obs.wads$wad)]
  g2.reps.NAs <- boot.out.heavy$obs.wads$rep[is.na(boot.out.heavy$obs.wads$wad)]
  
  if (length(c(g1.reps.NAs, g2.reps.NAs)) == 0){
    message <- "none"
  }
  else if (length(g1.reps.NAs) != 0 & length(g2.reps.NAs) == 0){
    message <- paste("Warning: no occurrences in rep ", paste(g1.reps.NAs, collapse=" & "), " (", paste(unique(X.light[,var]), collapse="/"), ")", sep="")
  }
  else if (length(g1.reps.NAs) == 0 & length(g2.reps.NAs) != 0){
    message <- paste("Warning: no occurrences in rep ", paste(g2.reps.NAs, collapse=" & "), " (", paste(unique(X.heavy[,var]), collapse="/"), ")", sep="")
  }
  else  message <- paste("Warning: no occurrences in rep ", paste(g1.reps.NAs, collapse=" & "), " (", paste(unique(X.light[,var]), collapse="/"), ") and rep ", paste(g2.reps.NAs, collapse=" & "), " (", paste(unique(X.heavy[,var]), collapse="/"), ")", sep="")
  message
}
