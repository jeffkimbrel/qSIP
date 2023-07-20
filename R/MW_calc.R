#' Calculate Molecular Weight
#'
#' Given a data frame for unlabeled treatment(s) for a specific taxon (i.e., no
#' added heavy isotopes), calculate the weighted average density (WAD), GC
#' content, molecular weight, and carbon content
#'
#' assumes a fixed GC content, molecular weight, carbon, and nitrogen content
#' for each taxon based on the mean observed WAD
#' oxygen content is fixed for DNA regardless of GC content
#' the weighted average density of unlabeled DNA is a good proxy for GC content of the taxon, according to the well-known relationship between GC content and the density of DNA (uses equation from McHugh & Morrissey unpublished data)
#' (NOTE: This is not caused by the effect of GC content on molecular weight, which is negligible. It is caused by differences in the relationship between base composition and binding with water that occurs in CsCl.)
#'
#' Written by Ben Koch & Natasja van Gestel
#'
#' @param df data frame that includes x values (e.g., density of DNA) and y values (e.g., number of copies) and the replicate IDs (e.g., tube number)
#' @param density_column Column header with fraction density (formerly vars[1])
#' @param copies_ul_column Column header with abundance or copy numbers (formerly vars[2])
#' @param tube_column Column header with tube or sample identifier (formerly vars[3])
#'
#' @return A data.frame with the following columns
#' \itemize{
#'   \item WAD: mean of observed weighted average density for each replicate
#'   \item GC: guanine + cytosine content of the taxon's DNA
#'   \item MW: molecular weight (g/mol) of the taxon's DNA
#'   \item Catoms: mean number of carbon atoms per nucleotide for the taxon
#'   \item Oatoms: mean number of carbon atoms per nucleotide for the taxon
#'   \item Natoms: mean number of nitrogen atoms per nucleotide for the taxon
#'   \item message: warning message listing replicates that had no occurrences (i.e., all y values for a rep are 0); value for message is "none" when all replicates are present
#' }
#'
#' @note requires the function 'WAD_func'
#'
#' @export

MW_calc <- function(df, density_column = "Density", copies_ul_column = "t.copies.ul", tube_column = "unique.tube"){

  # Create a dataframe of only x, y, and rep:
  test.data = df %>%
    select(x = !!as.name(density_column),
           y = !!as.name(copies_ul_column),
           rep = !!as.name(tube_column)) %>%
    mutate(rep = as.factor(rep))

  test.data$rep = droplevels(test.data$rep)

  # Calculate observed weighted average density (WAD) for each rep:
  obs.wads <- data.frame(matrix(nrow=length(levels(test.data$rep)), ncol=2))
  names(obs.wads) <- c("wad", "rep")

  for (r in 1:length(levels(test.data$rep))){
    obs.wads$rep[r] <- levels(test.data$rep)[r]
    obs.wads$wad[r] <- WAD_func(y=test.data$y[test.data$rep == levels(test.data$rep)[r]],
                                x=test.data$x[test.data$rep == levels(test.data$rep)[r]])
  }

  obs.wads$rep <- factor(obs.wads$rep)

  reps.NAs <- obs.wads$rep[is.na(obs.wads$wad)]
  if (length(reps.NAs) == 0){
    message <- "none"
  } else  {
    message <- paste("Warning: no occurrences in rep ", paste(reps.NAs, collapse=" & "), sep="")
  }

  # Calculate the mean of observed weighted average density across reps:
  density <- mean(obs.wads$wad, na.rm=T)

  # Calculate GC content, molecular weight, and C content:
  # GC <- (1/0.098)*(density-1.66)  #From Schildkraut et al. 1962
  GC <- (1/0.0835059954345993)*(density-1.64605745338531)  #From McHugh & Morrissey unpublished data
  MW <- (GC*0.496) + 307.691      #Molecular weight (g/mol) varies with GC content
  Catoms <- (-0.5*GC) + 10        #Carbon content varies with GC content
  Oatoms <- 6                     #Oxygen content is constant for DNA regardless of GC content
  Natoms <- (0.5*GC) + 3.5        #Nitrogen content varies with GC content
  data.frame(WAD=density,
             GC=GC,
             MW=MW,
             Catoms=Catoms,
             Oatoms=Oatoms,
             Natoms=Natoms,
             message=message)
}

#' Calculate Molecular Weight (Deprecated)
#'
#' @keywords internal
#' @export

MW.calc = function(...) {
  .Defunct("MW_calc()")
}
