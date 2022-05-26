### Given a data frame for unlabeled treatment(s) for a specific taxon (i.e., no added heavy isotopes), calculate the weighted average density (WAD), GC content, molecular weight, and carbon content

#     output = MW.calc(X, vars=c("density.g.ml", "copies", "tube"))
#
#     X: data frame that includes x values (e.g., density of DNA) and y values (e.g., number of copies) and the replicate IDs (e.g., tube number)
#     vars: vector of variable names for the x, y, and replicate ID columns in the data frame; default is c("density.g.ml", "copies" ,"tube")
#     -------------------------------------------------------
#     output: 
#     data frame of:  WAD: mean of observed weighted average density for each replicate
#                     GC: guanine + cytosine content of the taxon's DNA
#                     MW: molecular weight (g/mol) of the taxon's DNA
#                     Catoms: mean number of carbon atoms per nucleotide for the taxon
#                     Oatoms: mean number of carbon atoms per nucleotide for the taxon
#                     Natoms: mean number of nitrogen atoms per nucleotide for the taxon
#                     message: warning message listing replicates that had no occurrences (i.e., all y values for a rep are 0); value for message is "none" when all replicates are present
#
#     notes:          requires the function 'WAD.func'
#     assumptions:    assumes a fixed GC content, molecular weight, carbon, and nitrogen content for each taxon based on the mean observed WAD
#                     oxygen content is fixed for DNA regardless of GC content
#                     the weighted average density of unlabeled DNA is a good proxy for GC content of the taxon, according to the well-known relationship between GC content and the density of DNA (uses equation from McHugh & Morrissey unpublished data)
#                       (NOTE: This is not caused by the effect of GC content on molecular weight, which is negligible. It is caused by differences in the relationship between base composition and binding with water that occurs in CsCl.)
# 
#     Written by Ben Koch & Natasja van Gestel


MW.calc <- function(X, vars=c("density.g.ml", "copies", "tube")){

  # Create a dataframe of only x, y, and rep: 
   test.data <- data.frame(x=X[,vars[1]], y=X[,vars[2]], rep=factor(X[,vars[3]]))

  # Calculate observed weighted average density (WAD) for each rep:
   obs.wads <- data.frame(matrix(nrow=length(levels(test.data$rep)), ncol=2))
   names(obs.wads) <- c("wad", "rep")
   for (r in 1:length(levels(test.data$rep))){
     obs.wads$rep[r] <- levels(test.data$rep)[r]
     obs.wads$wad[r] <- WAD.func(y=test.data$y[test.data$rep == levels(test.data$rep)[r]], x=test.data$x[test.data$rep == levels(test.data$rep)[r]])
      }
   obs.wads$rep <- factor(obs.wads$rep)
    
   reps.NAs <- obs.wads$rep[is.na(obs.wads$wad)]
   if (length(reps.NAs) == 0){
     message <- "none"
   }
   else  message <- paste("Warning: no occurrences in rep ", paste(reps.NAs, collapse=" & "), sep="")
  
  # Calculate the mean of observed weighted average density across reps:
   density <- mean(obs.wads$wad, na.rm=T)
  
  # Calculate GC content, molecular weight, and C content:
   # GC <- (1/0.098)*(density-1.66)  #From Schildkraut et al. 1962
   GC <- (1/0.0835059954345993)*(density-1.64605745338531)  #From McHugh & Morrissey unpublished data
   MW <- (GC*0.496) + 307.691      #Molecular weight (g/mol) varies with GC content
   Catoms <- (-0.5*GC) + 10        #Carbon content varies with GC content
   Oatoms <- 6                     #Oxygen content is constant for DNA regardless of GC content
   Natoms <- (0.5*GC) + 3.5        #Nitrogen content varies with GC content
   data.frame(WAD=density, GC=GC, MW=MW, Catoms=Catoms, Oatoms=Oatoms, Natoms=Natoms, message=message)
}
