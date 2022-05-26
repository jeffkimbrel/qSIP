### Given two treatments, returns statistics to assess whether the means of weighted average density (WAD) differ significantly between treatments
#
#     output = diff.wad.calc(X1, X2, boot.out1, boot.out2, var="trt.code", CI=0.90, draws=1000, tailed.test=2)
#
#     X1: data frame with data for treatment 1
#     X2: data frame with data for treatment 2
#     boot.out1: list containing output from boot.WAD.func or boot.TUBE.func for treatment 1 (corresponding to X1; the 'light' treatment if a 1-tailed permutation test is specified)
#     boot.out2: list containing output from boot.WAD.func or boot.TUBE.func for treatment 2 (corresponding to X2; the 'heavy' treatment if a 1-tailed permutation test is specified)
#     var: variable name for the treatment ID columns in the data frames; default is "trt.code"
#     CI: confidence interval (0-1); default is 0.90
#     draws: number of iterations for the permutation test; default is 1000
#     tailed.test: specifies whether a 2-tailed or 1-tailed permutation test should be performed in calculating the p-value; default is 2; note: if 1-tailed is specified, then the function assumes group 2 (X2) is the 'heavy' treatment and group 1 (X1) is the light treatment
#     -------------------------------------------------------
#     output: 
#     list of:  boot.wads1: vector of bootstrapped WADs for treatment group 1 (X1)
#               boot.wads2: vector of bootstrapped WADs for treatment group 2 (X2)
#               boot.diffs: vector of bootstrapped difference in WAD between group 2 and group 1 (calculated as group2 - group1)
#               obs.wads.by.group: data frame of mean (across reps) of the observed WAD for each group
#               boot.wads.by.group: data frame of mean, median, CI-lower, and CI-upper of the bootstrapped WADs for each group
#               obs.diff: observed difference in WADs between treatments (i.e. subtracting mean WAD -- across reps -- for group 1 from that for group 2)
#               boot.diffs.mean: mean of the bootstrapped differences in WAD between group 2 and group 1
#               boot.diffs.median: median of the bootstrapped differences in WAD between group 2 and group 1
#               boot.diffs.CI: upper and lower confidence intervals (defined by argument 'CI') of the bootstrapped differences in WAD between group 2 and group 1
#               p.value: P-value (obtained by permutation test); the probability of observing a difference between treatments that is as large or larger than the observed difference, if the null hypothesis (that the treatments have equal WADs) is true
#               message: warning message listing replicates in either group that had no occurrences (i.e., all y values for a rep are 0); value for message is "none" when all replicates are present (single character value)
#
#     notes:    requires the functions 'WAD.func', 'boot.WAD.func', & 'comparison.message'
#               names of the x, y, replicate ID, and treatment ID columns in the two data frames must match exactly
#               number of reps need not be equal between the two groups
#               the test statistic for the P-value calculation is the difference in means of the 2 groups (calculated as group2 - group1)
#               the P-value estimates the probability that the observed difference (group2-group1) is as extreme or more extreme than would be expected if group identities were randomized
#
#     Written by Ben Koch & Natasja van Gestel


diff.wad.calc <- function(X1, X2, boot.out1, boot.out2, var="trt.code", CI=0.90, draws=1000, tailed.test=2){

  # Calculate quantities for output:
   obs.diff <- boot.out2$obs.wad.mean - boot.out1$obs.wad.mean
   boot.diffs <- boot.out2$boot.wads - boot.out1$boot.wads
   boot.diffs.CI <- quantile(boot.diffs, probs=c((1-CI)/2, 1-((1-CI)/2)), na.rm=T)

  # Create obs.wads.by.group dataframe for output:
   group1.names <- sort(levels(factor(X1[,var])))
   group2.names <- sort(levels(factor(X2[,var])))
   if (length(group1.names) > 1 & length(group2.names) <= 1) warning("treatment 1 in 'diff.wad.calc' contains multiple levels")
   if (length(group1.names) <= 1 & length(group2.names) > 1) warning("treatment 2 in 'diff.wad.calc' contains multiple levels")
   if (length(group1.names) > 1 & length(group2.names) > 1) warning("treatments 1 & 2 in 'diff.wad.calc' both contain multiple levels")
   group.names <- c(paste(sort(levels(factor(X1[,var]))), collapse="_"), paste(sort(levels(factor(X2[,var]))), collapse="_"))
   obs.wads.by.group <- data.frame(matrix(nrow=2, ncol=2))
   names(obs.wads.by.group) <- c("obs.wad.mean", "group")
   obs.wads.by.group$group <- group.names
   obs.wads.by.group$group <- factor(obs.wads.by.group$group)
   obs.wads.by.group$obs.wad.mean[1] <- boot.out1$obs.wad.mean
   obs.wads.by.group$obs.wad.mean[2] <- boot.out2$obs.wad.mean

  # Conduct the permutation test and calculate the p-value:
   obs.wads.all <- data.frame(wad=c(boot.out1$obs.wads$wad, boot.out2$obs.wads$wad), group=c(rep("group1", length(boot.out1$obs.wads$wad)), rep("group2", length(boot.out2$obs.wads$wad))))
   N <- length(obs.wads.all[,1])
   perm.diff <- numeric(draws)  
   for (p in 1:draws){
     perm <- sample.vec(1:N)
     perm.diff[p] <- diff(by(obs.wads.all[perm,1], obs.wads.all[,2], mean, na.rm=T))
   }
   if (tailed.test == 2){
     p.value <- (sum(abs(perm.diff)>=abs(obs.diff), na.rm=T))/length(perm.diff[!is.na(perm.diff)])      #proportion of abs(Permuted diff) values that exceed (or equal) abs(obs diff) (equivalent to a p-value); note that this is a two-tailed test
   }
   else if (tailed.test == 1){
     p.value <- sum(perm.diff>=obs.diff, na.rm=T)/length(perm.diff[!is.na(perm.diff)])				          #proportion of permuted diff values that exceed (or equal) obs diff (equivalent to a p-value); note that this is a one-tailed test and assumes that the WAD for group 2 is expected to be greater than that for group 1
   }
   
   message <- comparison.message(X.light=X1, X.heavy=X2, boot.out.light=boot.out1, boot.out.heavy=boot.out2, var=var)

  # Create boot.wads.by.group dataframe for output:
   boot.wads.by.group <- data.frame(matrix(nrow=2, ncol=5))
   names(boot.wads.by.group) <- c("boot.wads.mean", "boot.wads.median", "boot.wads.CI.L", "boot.wads.CI.U", "group")
   boot.wads.by.group$group <- group.names
   boot.wads.by.group$group <- factor(boot.wads.by.group$group)
   boot.wads.by.group$boot.wads.mean[1] <- boot.out1$boot.wads.mean
   boot.wads.by.group$boot.wads.mean[2] <- boot.out2$boot.wads.mean
   boot.wads.by.group$boot.wads.median[1] <- boot.out1$boot.wads.median
   boot.wads.by.group$boot.wads.median[2] <- boot.out2$boot.wads.median
   boot.wads.by.group[1,c(3:4)] <- boot.out1$boot.wads.CI
   boot.wads.by.group[2,c(3:4)] <- boot.out2$boot.wads.CI

  # Collect all output into a list:
   return(list(boot.wads1=boot.out1$boot.wads, 
               boot.wads2=boot.out2$boot.wads, 
               boot.diffs=boot.diffs, 
               obs.wads.by.group=obs.wads.by.group, 
               boot.wads.by.group=boot.wads.by.group, 
               obs.diff=obs.diff, 
               boot.diffs.mean=mean(boot.diffs, na.rm=T), 
               boot.diffs.median=median(boot.diffs, na.rm=T), 
               boot.diffs.CI=boot.diffs.CI, 
               p.value=p.value, 
               message=message))
}
