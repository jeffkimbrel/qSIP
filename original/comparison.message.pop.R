
comparison.message.pop <- function(T0, Tt, boot.out.T0, boot.out.Tt, var="trt.code"){
  g1.reps.NAs <- boot.out.T0$obs.N$rep[is.na(boot.out.T0$obs.N$copies) | boot.out.T0$obs.N$copies == 0]
  g2.reps.NAs <- boot.out.Tt$obs.N$rep[is.na(boot.out.Tt$obs.N$copies) | boot.out.Tt$obs.N$copies == 0]

  if (length(c(g1.reps.NAs, g2.reps.NAs)) == 0){
    message <- "none"
  }
  else if (length(g1.reps.NAs) != 0 & length(g2.reps.NAs) == 0){
    message <- paste("Warning: no occurrences in rep ", paste(g1.reps.NAs, collapse=" & "), " (", paste(unique(T0[,var]), collapse="/"), ")", sep="")
  }
  else if (length(g1.reps.NAs) == 0 & length(g2.reps.NAs) != 0){
    message <- paste("Warning: no occurrences in rep ", paste(g2.reps.NAs, collapse=" & "), " (", paste(unique(Tt[,var]), collapse="/"), ")", sep="")
  }
  else  message <- paste("Warning: no occurrences in rep ", paste(g1.reps.NAs, collapse=" & "), " (", paste(unique(T0[,var]), collapse="/"), ") and rep ", paste(g2.reps.NAs, collapse=" & "), " (", paste(unique(Tt[,var]), collapse="/"), ")", sep="")
  message
}
