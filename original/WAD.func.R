### Given vectors of x values and y values, calculate the weighted-average of the x-values (e.g., the weighted average density (WAD))
#
#     output = WAD.func(y, x)
#
#     y: vector of y-values (e.g., number of 16S copies)
#     x: vector of x-values (e.g., density of DNA)
#     -------------------------------------------------------
#     output: weighted-average of the x-values (single value)
#     Written by Ben Koch & Natasja van Gestel


  WAD.func <- function(y, x){
    WAD <- sum(x*(y/sum(y)))
    WAD
  }
