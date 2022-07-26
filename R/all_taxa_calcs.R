#' All Taxa Calcs
#'
#' Given a data frame (of all treatments and all taxa), along with a data frame
#' specifying treatment comparisons to make.
#'
#' For all taxa, calculate growth rates, carbon fluxes, abundances, excess atom
#' fraction, and return statistics to assess whether the means of weighted
#' average density (WAD) differ significantly between treatments
#'
#' For the written files, values are the quantities specified in the name of
#' the data frame, columns are individual bootstrap samples (i.e., mean of
#' resampled data for a single iteration) and rows correspond to different taxa
#' and treatment comparisons
#'
#' Written by Ben Koch & Natasja van Gestel
#'
#' @param X.all data frame with data for all taxa and treatments
#' @param comparisons data frame specifying which treatments to compare; must contain these columns: "comparisonID", "iso.compare", "trt.code.1", "trt.code.2", "trt.refs", "days"
#' @param M.soil data frame with mass of soil for each replicate; must contain two columns, one for replicate ID and one for mass of soil; default is NULL meaning that is M.soil is not provided, only WAD and ape results are returned
#' @param taxon_column column with taxon IDs (formerly vars[1])
#' @param density_column column with density data (formerly #vars[2])
#' @param copies_ul_column column with abundance or copy number data (formerly vars[3])
#' @param tube_column column with tube or sample IDs (formerly vars[4])
#' @param treatment_column column with treatment IDs (formerly vars[5])
#' @param soil_g_column column with soil weight (formerly vars[6])
#' @param growth.model model to use for calculating growth ("linear" or "exponential")
#' @param prop.O.from.water proportion of oxygen atoms in DNA that come from environmental water
#' @param v.frac volume of a fraction within a tube (uL)
#' @param copies.cell assumed number of 16S copies per cell
#' @param pgC.cell assumed mass of carbon per cell (units: picograms)
#' @param CI confidence interval (0-1)
#' @param draws number of bootstrap iterations
#' @param tailed.test specifies whether a 2-tailed or 1-tailed permutation test should be performed in calculating the p-value; default is 2; note: if 1-tailed is specified, then the function assumes group 2 (trt.code.2) is the 'heavy' treatment and group 1 (trt.code.1) is the 'light' treatment
#'
#' @return A list with
#' \itemize{
#'   \item taxonID: identifier for the taxon (factor)
#'   \item comparisonID: identifier for the treatment comparison
#'   \item trt.code.1: ID code for treatment 1 (factor)
#'   \item trt.code.2:  ID code for treatment 1 (factor)
#'   \item wad1.obs.mean: mean of observed weighted average density values for each tube (rep) in treatment 1
#'   \item wad2.obs.mean: mean of observed weighted average density values for each tube (rep) in treatment 2
#'   \item wad1.boot.mean: mean of bootstrapped weighted average densities treatment 1
#'   \item wad1.boot.median: median of bootstrapped weighted average densities treatment 1
#'   \item wad1.boot.CI.L: lower end of confidence interval of bootstrapped weighted average densities for treatment 1
#'   \item wad1.boot.CI.U: upper end of confidence interval of bootstrapped weighted average densities for treatment 1
#'   \item wad2.boot.mean: mean of bootstrapped weighted average densities treatment 2
#'   \item wad2.boot.median: median of bootstrapped weighted average densities treatment 2
#'   \item wad2.boot.CI.L: lower end of confidence interval of bootstrapped weighted average densities for treatment 2
#'   \item wad2.boot.CI.U: upper end of confidence interval of bootstrapped weighted average densities for treatment 2
#'   \item wad.diff.obs: observed difference in WADs between treatments (i.e. subtracting mean WAD -- across reps -- for group 1 from that for group 2)
#'   \item wad.diff.boot.mean: mean of the bootstrapped differences in WAD between group 2 and group 1
#'   \item wad.diff.boot.median: median of the bootstrapped differences in WAD between group 2 and group 1
#'   \item wad.diff.boot.CI.L: lower end of confidence interval of the bootstrapped differences in WAD between group 2 and group 1
#'   \item wad.diff.boot.CI.U: upper end of confidence interval of the bootstrapped differences in WAD between group 2 and group 1
#'   \item wad.diff.p.value: P-value (obtained by permutation test); the probability of observing a difference between treatments that is as large or larger than the observed difference, if the null hypothesis (that the treatments have equal WADs) is true
#'   \item ape.obs: observed excess atom fraction based on the two treatments in the comparison
#'   \item ape.boot.mean: mean of the bootstrapped excess atom fraction values
#'   \item ape.boot.median: median of the bootstrapped excess atom fraction values
#'   \item ape.boot.CI.L: lower confidence intervals of the bootstrapped excess atom fraction values
#'   \item ape.boot.CI.U: upper confidence intervals of the bootstrapped excess atom fraction values
#'   \item r.obs: observed intrinsic rate of increase based on the two treatments in the comparison
#'   \item r.boot.mean: mean of the bootstrapped intrinsic rate of increase values
#'   \item r.boot.median: median of the bootstrapped intrinsic rate of increase values
#'   \item r.boot.CI.L: lower confidence intervals of the bootstrapped intrinsic rate of increase values
#'   \item r.boot.CI.U: upper confidence intervals of the bootstrapped intrinsic rate of increase values
#'   \item f.obs: observed carbon flux based on the two treatments in the comparison
#'   \item f.boot.mean: mean of the bootstrapped carbon flux values
#'   \item f.boot.median: median of the bootstrapped carbon flux values
#'   \item f.boot.CI.L: lower confidence intervals of the bootstrapped carbon flux values
#'   \item f.boot.CI.U: upper confidence intervals of the bootstrapped carbon flux values
#'   \item N1.obs.mean: mean of observed abundances for each tube (rep) in treatment 1
#'   \item N2.obs.mean: mean of observed abundances for each tube (rep) in treatment 2
#'   \item N1.boot.mean: mean of bootstrapped abundances treatment 1
#'   \item N1.boot.median: median of bootstrapped abundances treatment 1
#'   \item N1.boot.CI.L: lower end of confidence interval of bootstrapped abundances for treatment 1
#'   \item N1.boot.CI.U: upper end of confidence interval of bootstrapped abundances for treatment 1
#'   \item N2.boot.mean: mean of bootstrapped abundances treatment 2
#'   \item N2.boot.median: median of bootstrapped abundances treatment 2
#'   \item N2.boot.CI.L: lower end of confidence interval of bootstrapped abundances for treatment 2
#'   \item N2.boot.CI.U: upper end of confidence interval of bootstrapped abundances for treatment 2
#'   \item N.obs.mean: mean of observed abundances for each tube (rep) in treatments 1 and 2
#'   \item N.boot.mean: mean of bootstrapped abundances treatments 1 and 2
#'   \item N.boot.median: median of bootstrapped abundances treatments 1 and 2
#'   \item N.boot.CI.L: lower end of confidence interval of bootstrapped abundances for treatments 1 and 2
#'   \item N.boot.CI.U: upper end of confidence interval of bootstrapped abundances for treatments 1 and 2
#'   \item message: warning message listing replicates in either of the 2 treatments (does not include the reference group) that had no occurrences (i.e., all y values for a rep are 0); value for message is "none" when all replicates are present (single character value)
#' }
#'
#' @return The following files will be written
#' \itemize{
#'   \item bootstrapped_wad1.txt
#'   \item bootstrapped_wad2.txt
#'   \item bootstrapped_wad_diff.txt
#'   \item bootstrapped_ape.txt
#'   \item bootstrapped_r.txt
#'   \item bootstrapped_C_fluxes.txt
#'   \item bootstrapped_N1.txt
#'   \item bootstrapped_N2.txt
#'   \item bootstrapped_N.txt
#' }
#'
#' @note requires the functions 'WAD.func', 'boot.TUBE.func', 'diff.wad.calc', 'MW.calc', 'ape.calc', 'r.calc', 'f.calc', & 'comparison.message'
#' @note not implemented yet: include an argument to specify which quantities the user wants calculated e.g., ('wad.diff', 'ape', 'f', 'r', 'all')
#' @note (not all inputs would be needed for certain quantities)
#' @note If M.soil is not provided, then only WAD and ape results are returned; if M.soil is provided then all results (WAD, ape, r, f, N) are returned
#' @note treatment code names cannot contain spaces
#' @note the treatments in the comparisons data frame (columns trt.code.1, trt.code.2, trt.refs) can be separated by semicolons, commas, or spaces (or a combination of semicolons and spaces or commas and spaces)
#' @note if the comparisons data frame is imported as a csv file, using commas to separate reference treatments may be problematic unless quoting is used (recommend using tab-delimited text files for import)
#' @note this function (along with ape, r, f) cannot handle comparing doubly-labeled treatments (yet)
#' @note 'copies' is really 'copies.ul'
#' @note units of abundance (N) are in copies/unit medium, where the unit medium is specified by vars[6] (i.e., mass of soil: 'g.soil' by default)
#'
#' @export
#'


all_taxa_calcs <- function(X.all,
                           comparisons,
                           M.soil=Sdat,
                           taxon_column = "taxon",            #vars[1]
                           density_column = "Density",        #vars[2]
                           copies_ul_column = "t.copies.ul",  #vars[3]
                           tube_column = "unique.tube",       #vars[4]
                           treatment_column = "unique.tmt",   #vars[5]
                           soil_g_column = "g.dry.soil.tube", #vars[6]
                           #vars=c("taxon", "Density", "t.copies.ul", "unique.tube", "unique.tmt", "g.dry.soil.tube"),
                           growth.model="linear",
                           prop.O.from.water=0.60,
                           v.frac=30,
                           copies.cell=6,
                           pgC.cell=0.1,
                           CI=0.90,
                           draws=1000,
                           tailed.test=1) {

  require(progress)

  #If M.soil is unspecified, then create a sham M.soil data.frame for the sake of doing the calculations, but r, f, and N results will not be included in output:
  sham <- FALSE
  if (is.null(M.soil)){
    sham <- TRUE
    M.soil = tibble("A" = levels(X.all$df$unique.tube), "B" = 1) %>%
      select(!!as.name(tube_column) := A, !!as.name(soil_g_column) := B)
  }

  # (L108) ------

  #Create an empty data frame for the output:
  info <- data.frame(matrix(nrow=0, ncol=51))

  #Create empty data frames for the bootstrapped abundance (N1, N2, N), C flux, and bootstrapped growth rate outputs, and also the bootstrapped ape, wad.diff, wad2, and wad1 values:
  N.boots <- N2.boots <- N1.boots <- C.flux.boots <- r.boots <- ape.boots <- wad.diff.boots <- wad2.boots <- wad1.boots <- data.frame(matrix(nrow=0, ncol=draws+4))

  #Establish the number of comparisons and the number of taxa:
  N.comparisons <- length(levels(factor(comparisons$comparisonID)))
  N.taxa = length(X.all$taxa)

  # (L118) ------
  #Create a key relating treatment codes in the raw data to the "effective" treatments specified by (potentially) multiple trt codes in the comparisons data frame:

  trt.names = tibble("trt" = c(1, 2),
                     "eff.trt.name" = c(paste(parse_treatments(comparisons$trt.code.1), collapse = '_'),
                                        paste(parse_treatments(comparisons$trt.code.2), collapse = '_')),
                     "trt.name" = c(comparisons$trt.code.1,
                                    comparisons$trt.code.2))

  # (L136) ------

  #Establish the number of unique treatments and 'effective' treatments specified in the 'comparisons' data.frame:
  # LOOKS LIKE THIS WILL ALWAYS BE TWO?
  all.eff.treatments <- sort(unique(as.character(trt.names$eff.trt.name)))
  N.eff.treatments <- length(all.eff.treatments)

  #Create an empty dataframe for building a key that relates treatment codes to the names of treatments used in the output lists:
  TRTID <- data.frame(trt.mat.name=character(N.eff.treatments), trt.code=character(N.eff.treatments), stringsAsFactors=FALSE)

  #Generate bootstrap samples for all taxon-treatment combinations:
  for (p in 1:N.eff.treatments){
    #Calculate the number of reps over which to resample - for each 'effective' treatment:
    curr.trts.incl.NAs <- as.character(unique(trt.names[trt.names$eff.trt.name == all.eff.treatments[p], 3:dim(trt.names)[2]]))
    test.data <- X.all$df %>%
      filter(!!as.name(treatment_column) %in% parse_treatments(curr.trts.incl.NAs))
    N.reps <- test.data %>% select(!!as.name(tube_column)) %>% unique() %>% pull() %>% length()

    #Calculate bootstrap vectors of the resampled indices by which to calculate mean WADs, mean of total 16S copies, and mass of soil across reps:
    trt.indices <- matrix(nrow=draws, ncol=N.reps)
    for (h in 1:draws){
      trt.indices[h,] <- sample_vec(1:N.reps, N.reps, replace=TRUE)
    }

    #Name the matrix for the current (effective) treatment:
    assign(paste("TRT_", as.character(all.eff.treatments[p]), sep=""), trt.indices)
    #Write the current (effective) treatment code and the name of the (effective) treatment used in the output matrix to the 'key' dataframe:
    TRTID$trt.mat.name[p] <- paste("TRT_", as.character(all.eff.treatments[p]), sep="")
    TRTID$trt.code[p] <- as.character(all.eff.treatments[p])
  }

  #Convert columns in the TRTID 'key' dataframe to factor:
  TRTID$trt.mat.name <- factor(TRTID$trt.mat.name)
  TRTID$trt.code <- factor(TRTID$trt.code)

  # (L169) ------

  for (c in 1:N.comparisons){  #for each comparison...
    #Get the code specifying which isotopes are being compared in the current comparison:
    iso.compare <- as.character(comparisons$iso.compare[c])

    if (is.na(iso.compare)) {
    } else {
      #Get a vector of codes for treatment 1 that are to be used as unlabeled 'light' treatments (without heavy isotopes) in the current comparison:
      trt.1 = parse_treatments(comparisons$trt.code.1[c])

      #Get a vector of codes for treatment 2 that are to be used as labeled 'heavy' treatments in the current comparison:
      trt.2 = parse_treatments(comparisons$trt.code.2[c])


      #Get a vector of codes for the treatments that are to be used as 'reference' treatments (without heavy isotopes) in the current comparison:
      trt.refs <- parse_treatments(comparisons$trt.refs[c])

      #Get the duration of the incubation period for the current comparison:
      days <- comparisons$days[c]

      #Get the names of the matrices of resampled indices for the 'light' and 'heavy' treatments:
      mat.name.light <- as.character(TRTID$trt.mat.name[1])
      mat.name.heavy <- as.character(TRTID$trt.mat.name[2])

      message(paste("Running calculations for", length(X.all$taxa), "taxa..."))
      pb <- progress_bar$new(total = length(X.all$taxa))
      for (taxa in X.all$taxa){  #for each taxon...
        pb$tick()

        #Subset the data by taxon and comparison and calculate the necessary values from the 'reference' treatments:
        #X.light.orig <- X.all.orig[X.all.orig[,vars[1]]==levels(X.all.orig[,vars[1]])[i] & X.all.orig[,vars[5]] %in% trt.1,]
        X.light = X.all$df %>%
          filter(!!as.name(taxon_column) == taxa) %>%
          filter(!!as.name(treatment_column) %in% trt.1)
        #dim(X.light.orig) == dim(X.light)

        #X.heavy.orig <- X.all.orig[X.all.orig[,vars[1]]==levels(X.all.orig[,vars[1]])[i] & X.all.orig[,vars[5]] %in% trt.2,]
        X.heavy = X.all$df %>%
          filter(!!as.name(taxon_column) == taxa) %>%
          filter(!!as.name(treatment_column) %in% trt.2)
        #dim(X.heavy.orig) == dim(X.heavy)

        #X.reference.orig <- X.all.orig[X.all.orig[,vars[1]]==levels(X.all.orig[,vars[1]])[i] & X.all.orig[,vars[5]] %in% trt.refs,]
        X.reference = X.all$df %>%
          filter(!!as.name(taxon_column) == taxa) %>%
          filter(!!as.name(treatment_column) %in% trt.refs)
        #dim(X.reference.orig) == dim(X.reference)

        MW.out = MW_calc(df = X.reference, density_column = "Density", copies_ul_column = "t.copies.ul", tube_column = "unique.tube")

        #204

        #First, for the 'light' treatment:
        #Calculate observed weighted average density (WAD), total 16S copies, and mass of soil for each rep:
        obs.wads = X.light %>%
          group_by(!!as.name(tube_column)) %>%
          summarize(wad = WAD_func(y=!!as.name(copies_ul_column), x = !!as.name(density_column)),
                    tot.copies = sum(!!as.name(copies_ul_column))*v.frac) %>%
          left_join(M.soil, by = tube_column) %>%
          select(wad, tot.copies, g.soil = !!as.name(soil_g_column), rep = !!as.name(tube_column))

        #216

        #Calculate a bootstrap vector of mean WADs across reps along with total 16S copies & mass of soil:
        resampled.tot.copies <- matrix(obs.wads$tot.copies[eval(parse(text=mat.name.light))], nrow=draws, ncol=dim(obs.wads)[1])
        resampled.g.soil <- matrix(obs.wads$g.soil[eval(parse(text=mat.name.light))], nrow=draws, ncol=dim(obs.wads)[1])
        resampled.wads <- matrix(obs.wads$wad[eval(parse(text=mat.name.light))], nrow=draws, ncol=dim(obs.wads)[1])
        boot.wads <- apply(resampled.wads, 1, mean, na.rm=T)

        boot.out.light <- list(
          boot.tot.copies = resampled.tot.copies,
          boot.g.soil = resampled.g.soil,
          boot.wads = boot.wads,
          obs.wads = obs.wads,
          obs.wad.mean = mean(obs.wads$wad, na.rm=T),
          boot.wads.mean=mean(boot.wads, na.rm=T),
          boot.wads.median=median(boot.wads, na.rm=T),
          boot.wads.CI=quantile(boot.wads, probs=c((1-CI)/2, 1-((1-CI)/2)), na.rm=T)
        )

        #233

        #Next, for the 'heavy' treatment:
        #Calculate observed weighted average density (WAD), total 16S copies, and mass of soil for each rep:
        obs.wads = X.heavy %>%
          group_by(!!as.name(tube_column)) %>%
          summarize(wad = WAD_func(y=!!as.name(copies_ul_column), x = !!as.name(density_column)),
                    tot.copies = sum(!!as.name(copies_ul_column))*v.frac) %>%
          left_join(M.soil, by = tube_column) %>%
          select(wad, tot.copies, g.soil = !!as.name(soil_g_column), rep = !!as.name(tube_column))

        #Calculate a bootstrap vector of mean WADs across reps along with total 16S copies & mass of soil:
        resampled.tot.copies <- matrix(obs.wads$tot.copies[eval(parse(text=mat.name.heavy))], nrow=draws, ncol=dim(obs.wads)[1])
        resampled.g.soil <- matrix(obs.wads$g.soil[eval(parse(text=mat.name.heavy))], nrow=draws, ncol=dim(obs.wads)[1])
        resampled.wads <- matrix(obs.wads$wad[eval(parse(text=mat.name.heavy))], nrow=draws, ncol=dim(obs.wads)[1])
        boot.wads <- apply(resampled.wads, 1, mean, na.rm=T)

        boot.out.heavy <- list(
          boot.tot.copies = resampled.tot.copies,
          boot.g.soil = resampled.g.soil,
          boot.wads = boot.wads,
          obs.wads = obs.wads,
          obs.wad.mean = mean(obs.wads$wad, na.rm=T),
          boot.wads.mean=mean(boot.wads, na.rm=T),
          boot.wads.median=median(boot.wads, na.rm=T),
          boot.wads.CI=quantile(boot.wads, probs=c((1-CI)/2, 1-((1-CI)/2)), na.rm=T)
        )

        #262

        #Calculate the appropriate output values (suppress warnings in diff.wad.calc; these warn about multiple treatment levels in 'X.light' and 'X.heavy' which may be commonly (and correctly) specified in the 'comparisons' dataframe):
        diff.wad.out <- suppressWarnings(diff_wad_calc(X1=X.light, X2=X.heavy, boot.out1=boot.out.light, boot.out2=boot.out.heavy, var=treatment_column, CI=CI, draws=draws, tailed.test=tailed.test))
        ape.out <- ape_calc(X.light=X.light, X.heavy=X.heavy, boot.out.light=boot.out.light, boot.out.heavy=boot.out.heavy, MW.out=MW.out, iso.compare=iso.compare, var=treatment_column, CI=CI)
        r.out <- r_calc(X.light=X.light, X.heavy=X.heavy, boot.out.light=boot.out.light, boot.out.heavy=boot.out.heavy, MW.out=MW.out, iso.compare=iso.compare, days=days, var=treatment_column, growth.model=growth.model, prop.O.from.water=prop.O.from.water, CI=CI)
        f.out <- f_calc(X.light=X.light, X.heavy=X.heavy, boot.out.light=boot.out.light, boot.out.heavy=boot.out.heavy, MW.out=MW.out, iso.compare=iso.compare, days=days, var=treatment_column, growth.model=growth.model, prop.O.from.water=prop.O.from.water, copies.cell=copies.cell, pgC.cell=pgC.cell, CI=CI)


        #Write the vector of bootstrapped abundances (N1, N2, N), C fluxes, and bootstrapped growth rates (and bootstrapped ape, wad.diff, wad2, & wad1) to output dataframes to enable downstream calculations (such as assemblage-level C flux):
        N.to.add <- N2.to.add <- N1.to.add <- f.to.add <- r.to.add <- ape.to.add <- wad.diff.to.add <- wad2.to.add <- wad1.to.add <- data.frame(matrix(nrow=1, ncol=draws+4))
        N.to.add[1] <- N2.to.add[1] <- N1.to.add[1] <- f.to.add[1] <- r.to.add[1] <- ape.to.add[1] <- wad.diff.to.add[1] <- wad2.to.add[1] <- wad1.to.add[1] <- taxa
        N.to.add[2] <- N2.to.add[2] <- N1.to.add[2] <- f.to.add[2] <- r.to.add[2] <- ape.to.add[2] <- wad.diff.to.add[2] <- wad2.to.add[2] <- wad1.to.add[2] <- as.character(comparisons$comparisonID[c])
        N.to.add[3] <- N2.to.add[3] <- N1.to.add[3] <- f.to.add[3] <- r.to.add[3] <- ape.to.add[3] <- wad.diff.to.add[3] <- wad2.to.add[3] <- wad1.to.add[3] <- as.character(diff.wad.out$obs.wads.by.group$group[1])
        N.to.add[4] <- N2.to.add[4] <- N1.to.add[4] <- f.to.add[4] <- r.to.add[4] <- ape.to.add[4] <- wad.diff.to.add[4] <- wad2.to.add[4] <- wad1.to.add[4] <- as.character(diff.wad.out$obs.wads.by.group$group[2])

        wad1.to.add[5:dim(wad1.to.add)[2]] <- diff.wad.out$boot.wads1
        wad2.to.add[5:dim(wad2.to.add)[2]] <- diff.wad.out$boot.wads2
        wad.diff.to.add[5:dim(wad.diff.to.add)[2]] <- diff.wad.out$boot.diffs
        ape.to.add[5:dim(ape.to.add)[2]] <- ape.out$boot.apes
        r.to.add[5:dim(r.to.add)[2]] <- r.out$boot.r
        f.to.add[5:dim(f.to.add)[2]] <- f.out$boot.f
        N1.to.add[5:dim(N1.to.add)[2]] <- apply(boot.out.light$boot.tot.copies/boot.out.light$boot.g.soil, 1, mean, na.rm=T)
        N2.to.add[5:dim(N2.to.add)[2]] <- apply(boot.out.heavy$boot.tot.copies/boot.out.heavy$boot.g.soil, 1, mean, na.rm=T)
        N.to.add[5:dim(N.to.add)[2]] <- apply(cbind(boot.out.light$boot.tot.copies/boot.out.light$boot.g.soil, boot.out.heavy$boot.tot.copies/boot.out.heavy$boot.g.soil), 1, mean, na.rm=T)
        names(wad1.to.add) <- c("taxonID", "comparisonID", "trt.code.1", "trt.code.2", paste("wad1.", 1:draws, sep=""))
        names(wad2.to.add) <- c("taxonID", "comparisonID", "trt.code.1", "trt.code.2", paste("wad2.", 1:draws, sep=""))
        names(wad.diff.to.add) <- c("taxonID", "comparisonID", "trt.code.1", "trt.code.2", paste("wad.diff", 1:draws, sep=""))
        names(ape.to.add) <- c("taxonID", "comparisonID", "trt.code.1", "trt.code.2", paste("ape", 1:draws, sep=""))
        names(r.to.add) <- c("taxonID", "comparisonID", "trt.code.1", "trt.code.2", paste("r", 1:draws, sep=""))
        names(f.to.add) <- c("taxonID", "comparisonID", "trt.code.1", "trt.code.2", paste("f", 1:draws, sep=""))
        names(N1.to.add) <- c("taxonID", "comparisonID", "trt.code.1", "trt.code.2", paste("N1.", 1:draws, sep=""))
        names(N2.to.add) <- c("taxonID", "comparisonID", "trt.code.1", "trt.code.2", paste("N2.", 1:draws, sep=""))
        names(N.to.add) <- c("taxonID", "comparisonID", "trt.code.1", "trt.code.2", paste("N.", 1:draws, sep=""))

        wad1.boots <- rbind(wad1.boots, wad1.to.add)
        wad2.boots <- rbind(wad2.boots, wad2.to.add)
        wad.diff.boots <- rbind(wad.diff.boots, wad.diff.to.add)
        ape.boots <- rbind(ape.boots, ape.to.add)
        r.boots <- rbind(r.boots, r.to.add)
        C.flux.boots <- rbind(C.flux.boots, f.to.add)
        N1.boots <- rbind(N1.boots, N1.to.add)
        N2.boots <- rbind(N2.boots, N2.to.add)
        N.boots <- rbind(N.boots, N.to.add)

        #Write the appropriate output values for the current taxon-comparison to a data frame:
        info.to.add <- data.frame(
          taxonID=taxa,
          comparisonID=as.character(comparisons$comparisonID[c]),
          trt.code.1=as.character(diff.wad.out$obs.wads.by.group$group[1]),
          trt.code.2=as.character(diff.wad.out$obs.wads.by.group$group[2]),

          wad1.obs.mean=diff.wad.out$obs.wads.by.group$obs.wad.mean[1],
          wad2.obs.mean=diff.wad.out$obs.wads.by.group$obs.wad.mean[2],
          wad1.boot.mean=diff.wad.out$boot.wads.by.group$boot.wads.mean[1],
          wad1.boot.median=diff.wad.out$boot.wads.by.group$boot.wads.median[1],
          wad1.boot.CI.L=diff.wad.out$boot.wads.by.group$boot.wads.CI.L[1],
          wad1.boot.CI.U=diff.wad.out$boot.wads.by.group$boot.wads.CI.U[1],
          wad2.boot.mean=diff.wad.out$boot.wads.by.group$boot.wads.mean[2],
          wad2.boot.median=diff.wad.out$boot.wads.by.group$boot.wads.median[2],
          wad2.boot.CI.L=diff.wad.out$boot.wads.by.group$boot.wads.CI.L[2],
          wad2.boot.CI.U=diff.wad.out$boot.wads.by.group$boot.wads.CI.U[2],
          wad.diff.obs=diff.wad.out$obs.diff,
          wad.diff.boot.mean=diff.wad.out$boot.diffs.mean,
          wad.diff.boot.median=diff.wad.out$boot.diffs.median,
          wad.diff.boot.CI.L=as.numeric(diff.wad.out$boot.diffs.CI[1]),
          wad.diff.boot.CI.U=as.numeric(diff.wad.out$boot.diffs.CI[2]),
          wad.diff.p.value=diff.wad.out$p.value,

          ape.obs=ape.out$obs.ape,
          ape.boot.mean=ape.out$boot.apes.mean,
          ape.boot.median=ape.out$boot.apes.median,
          ape.boot.CI.L=as.numeric(ape.out$boot.apes.CI[1]),
          ape.boot.CI.U=as.numeric(ape.out$boot.apes.CI[2]),

          r.obs=r.out$obs.r,
          r.boot.mean=r.out$boot.r.mean,
          r.boot.median=r.out$boot.r.median,
          r.boot.CI.L=as.numeric(r.out$boot.r.CI[1]),
          r.boot.CI.U=as.numeric(r.out$boot.r.CI[2]),

          f.obs=f.out$obs.f,
          f.boot.mean=f.out$boot.f.mean,
          f.boot.median=f.out$boot.f.median,
          f.boot.CI.L=as.numeric(f.out$boot.f.CI[1]),
          f.boot.CI.U=as.numeric(f.out$boot.f.CI[2]),

          N1.obs.mean=mean(boot.out.light$obs.wads$tot.copies/boot.out.light$obs.wads$g.soil, na.rm=T),
          N2.obs.mean=mean(boot.out.heavy$obs.wads$tot.copies/boot.out.heavy$obs.wads$g.soil, na.rm=T),
          N1.boot.mean=mean(apply(boot.out.light$boot.tot.copies/boot.out.light$boot.g.soil, 1, mean, na.rm=T), na.rm=T),
          N1.boot.median=median(apply(boot.out.light$boot.tot.copies/boot.out.light$boot.g.soil, 1, mean, na.rm=T), na.rm=T),
          N1.boot.CI.L=quantile(apply(boot.out.light$boot.tot.copies/boot.out.light$boot.g.soil, 1, mean, na.rm=T), probs=(1-CI)/2, na.rm=T),
          N1.boot.CI.U=quantile(apply(boot.out.light$boot.tot.copies/boot.out.light$boot.g.soil, 1, mean, na.rm=T), probs=1-((1-CI)/2), na.rm=T),
          N2.boot.mean=mean(apply(boot.out.heavy$boot.tot.copies/boot.out.heavy$boot.g.soil, 1, mean, na.rm=T), na.rm=T),
          N2.boot.median=median(apply(boot.out.heavy$boot.tot.copies/boot.out.heavy$boot.g.soil, 1, mean, na.rm=T), na.rm=T),
          N2.boot.CI.L=quantile(apply(boot.out.heavy$boot.tot.copies/boot.out.heavy$boot.g.soil, 1, mean, na.rm=T), probs=(1-CI)/2, na.rm=T),
          N2.boot.CI.U=quantile(apply(boot.out.heavy$boot.tot.copies/boot.out.heavy$boot.g.soil, 1, mean, na.rm=T), probs=1-((1-CI)/2), na.rm=T),
          N.obs.mean=mean(c(boot.out.light$obs.wads$tot.copies/boot.out.light$obs.wads$g.soil, boot.out.heavy$obs.wads$tot.copies/boot.out.heavy$obs.wads$g.soil), na.rm=T),
          N.boot.mean=mean(apply(cbind(boot.out.light$boot.tot.copies/boot.out.light$boot.g.soil, boot.out.heavy$boot.tot.copies/boot.out.heavy$boot.g.soil), 1, mean, na.rm=T), na.rm=T),
          N.boot.median=median(apply(cbind(boot.out.light$boot.tot.copies/boot.out.light$boot.g.soil, boot.out.heavy$boot.tot.copies/boot.out.heavy$boot.g.soil), 1, mean, na.rm=T), na.rm=T),
          N.boot.CI.L=quantile(apply(cbind(boot.out.light$boot.tot.copies/boot.out.light$boot.g.soil, boot.out.heavy$boot.tot.copies/boot.out.heavy$boot.g.soil), 1, mean, na.rm=T), probs=(1-CI)/2, na.rm=T),
          N.boot.CI.U=quantile(apply(cbind(boot.out.light$boot.tot.copies/boot.out.light$boot.g.soil, boot.out.heavy$boot.tot.copies/boot.out.heavy$boot.g.soil), 1, mean, na.rm=T), probs=1-((1-CI)/2), na.rm=T),

          message=diff.wad.out$message
        )
        info <- rbind(info, info.to.add)
      }
    }
  }
  #Convert appropriate variables in 'N.boots', 'N2.boots', 'N1.boots', 'C.flux.boots', and 'r.boots' (and 'ape.boots', 'wad.diff.boots', 'wad2.boots', 'wad1.boots') to factors:
  wad1.boots$taxonID <- factor(wad1.boots$taxonID)
  wad1.boots$comparisonID <- factor(wad1.boots$comparisonID)
  wad1.boots$trt.code.1 <- factor(wad1.boots$trt.code.1)
  wad1.boots$trt.code.2 <- factor(wad1.boots$trt.code.2)
  wad2.boots$taxonID <- factor(wad2.boots$taxonID)
  wad2.boots$comparisonID <- factor(wad2.boots$comparisonID)
  wad2.boots$trt.code.1 <- factor(wad2.boots$trt.code.1)
  wad2.boots$trt.code.2 <- factor(wad2.boots$trt.code.2)
  wad.diff.boots$taxonID <- factor(wad.diff.boots$taxonID)
  wad.diff.boots$comparisonID <- factor(wad.diff.boots$comparisonID)
  wad.diff.boots$trt.code.1 <- factor(wad.diff.boots$trt.code.1)
  wad.diff.boots$trt.code.2 <- factor(wad.diff.boots$trt.code.2)
  ape.boots$taxonID <- factor(ape.boots$taxonID)
  ape.boots$comparisonID <- factor(ape.boots$comparisonID)
  ape.boots$trt.code.1 <- factor(ape.boots$trt.code.1)
  ape.boots$trt.code.2 <- factor(ape.boots$trt.code.2)
  r.boots$taxonID <- factor(r.boots$taxonID)
  r.boots$comparisonID <- factor(r.boots$comparisonID)
  r.boots$trt.code.1 <- factor(r.boots$trt.code.1)
  r.boots$trt.code.2 <- factor(r.boots$trt.code.2)
  C.flux.boots$taxonID <- factor(C.flux.boots$taxonID)
  C.flux.boots$comparisonID <- factor(C.flux.boots$comparisonID)
  C.flux.boots$trt.code.1 <- factor(C.flux.boots$trt.code.1)
  C.flux.boots$trt.code.2 <- factor(C.flux.boots$trt.code.2)
  N1.boots$taxonID <- factor(N1.boots$taxonID)
  N1.boots$comparisonID <- factor(N1.boots$comparisonID)
  N1.boots$trt.code.1 <- factor(N1.boots$trt.code.1)
  N1.boots$trt.code.2 <- factor(N1.boots$trt.code.2)
  N2.boots$taxonID <- factor(N2.boots$taxonID)
  N2.boots$comparisonID <- factor(N2.boots$comparisonID)
  N2.boots$trt.code.1 <- factor(N2.boots$trt.code.1)
  N2.boots$trt.code.2 <- factor(N2.boots$trt.code.2)
  N.boots$taxonID <- factor(N.boots$taxonID)
  N.boots$comparisonID <- factor(N.boots$comparisonID)
  N.boots$trt.code.1 <- factor(N.boots$trt.code.1)
  N.boots$trt.code.2 <- factor(N.boots$trt.code.2)

  #Export the abundance bootstrapped estimates, C flux bootstrapped estimates, and growth rate bootstrapped estimates (and ape, wad.diff, wad2, & wad1 bootstrapped estimates) to a text file:
  dir.create(path=paste(getwd(), "/qSIP_output", sep=""), showWarnings=FALSE)
  write.table(wad1.boots, "qSIP_output/bootstrapped_wad1.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
  write.table(wad2.boots, "qSIP_output/bootstrapped_wad2.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
  write.table(wad.diff.boots, "qSIP_output/bootstrapped_wad_diff.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
  write.table(ape.boots, "qSIP_output/bootstrapped_ape.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
  if (!sham){   #If M.soil was not a sham, then provide this output:
    write.table(r.boots, "qSIP_output/bootstrapped_r.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
    write.table(C.flux.boots, "qSIP_output/bootstrapped_C_fluxes.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
    write.table(N1.boots, "qSIP_output/bootstrapped_N1.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
    write.table(N2.boots, "qSIP_output/bootstrapped_N2.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
    write.table(N.boots, "qSIP_output/bootstrapped_N.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
  }

  #Convert appropriate variables in 'info' to factors:
  info$taxonID <- factor(info$taxonID)
  info$comparisonID <- factor(info$comparisonID)
  info$trt.code.1 <- factor(info$trt.code.1)
  info$trt.code.2 <- factor(info$trt.code.2)

  info %>%
    rownames_to_column("PERCENT") %>% pull(PERCENT) %>% table()

  if (!sham) {   #If M.soil was not a sham, then provide the complete output:
    return(info)
  } else {        #If M.soil was a sham, then do not provide the r, f, and N output:
    cols.to.keep <- names(info)[!(names(info) %in% c("r.obs", "r.boot.mean", "r.boot.median", "r.boot.CI.L", "r.boot.CI.U", "f.obs", "f.boot.mean", "f.boot.median", "f.boot.CI.L", "f.boot.CI.U", "N1.obs.mean", "N2.obs.mean", "N1.boot.mean", "N1.boot.median", "N1.boot.CI.L", "N1.boot.CI.U", "N2.boot.mean", "N2.boot.median", "N2.boot.CI.L", "N2.boot.CI.U", "N.obs.mean", "N.boot.mean", "N.boot.median", "N.boot.CI.L", "N.boot.CI.U"))]
    info <- info[,cols.to.keep]
    return(info)
  }

}




#' All Taxa Calcs (Deprecated)
#'
#' @keywords internal
#' @export

all.taxa.calcs = function(...) {
  .Defunct("all_taxa_calcs()")
}

