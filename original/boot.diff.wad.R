### Given 2 data frames (of 2 treatments), return statistics to assess whether the means of weighted average density (WAD) differ significantly between treatments
#
#     output = boot.diff.wad(X1, X2, vars=c("density.g.ml", "copies", "tube", "trt.code"), CI=0.90, draws=1000, tailed.test=2)
#
#     X1: data frame with data for treatment 1
#     X2: data frame with data for treatment 2
#     vars: vector of variable names for the x, y, replicate ID, and treatment ID columns in the data frames; default is c("density.g.ml", "copies" ,"tube", "trt.code")
#     CI: confidence interval (0-1); default is 0.90
#     draws: number of bootstrap iterations and number of iterations for the permutation test; default is 1000
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
#     notes:    requires the functions 'WAD.func', 'boot.WAD.func', 'diff.wad.calc', & 'comparison.message'
#               names of the x, y, replicate ID, and treatment ID columns in the two data frames must match exactly
#               number of reps need not be equal between the two groups
#               the test statistic for the P-value calculation is the difference in means of the 2 groups (calculated as group2 - group1)
#               the P-value estimates the probability that the observed difference (group2-group1) is as extreme or more extreme than would be expected if group identities were randomized
#
#     Written by Ben Koch & Natasja van Gestel


boot.diff.wad <- function(X1, X2, vars=c("density.g.ml", "copies", "tube", "trt.code"), CI=0.90, draws=1000, tailed.test=2){
   group1 <- boot.WAD.func(X=X1, vars=vars[1:3], CI=CI, draws=draws)
   group2 <- boot.WAD.func(X=X2, vars=vars[1:3], CI=CI, draws=draws)

   diff.wad.out <- diff.wad.calc(X1=X1, X2=X2, boot.out1=group1, boot.out2=group2, var=vars[4], CI=CI, draws=draws, tailed.test=tailed.test)
   diff.wad.out
 }
