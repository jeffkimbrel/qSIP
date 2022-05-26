### Given a data frame, calculate weighted average density (WAD) of each tube (replicate) for each taxon and also specify the tubes that compose each treatment
#
#     output = WAD.by.taxon.func(X, vars=c("taxon", "density.g.ml", "copies", "tube", "trt.code"))
#
#     X: data frame that includes taxa IDs, x values (e.g., density of DNA), y values (e.g., number of copies), and the replicate IDs (e.g., tube number)
#     vars: vector of variable names for the taxon, x, y, replicate ID, and treatment ID columns in the data frame; default is c("taxon", "density.g.ml", "copies", "tube", "trt.code")
#     -------------------------------------------------------
#     output: 
#     list of:  obs.wads: a data.frame of WAD values organized by taxon (rows) and by replicate ID (i.e., tube; the columns)
#               reps.by.trt: a data.frame listing the replicate IDs associated with each treatment
#
#     notes:    requires the function 'WAD.func'
# 
#     Written by Ben Koch


  WAD.by.taxon.func <- function(X, vars=c("taxon", "density.g.ml", "copies", "tube", "trt.code")){

    # Create a dataframe of only taxon, x, y, rep, and treatment: 
     test.data <- data.frame(taxon=factor(as.character(X[,vars[1]])), x=X[,vars[2]], y=X[,vars[3]], rep=factor(as.character(X[,vars[4]])), trt=factor(as.character(X[,vars[5]])))

    # Establish the number of taxa and the number of tubes (reps):
     num.taxa <- length(levels(test.data$taxon))
     num.reps <- length(levels(test.data$rep))

    # Calculate observed weighted average density (WAD) for each rep:
     obs.wads.by.taxon <- data.frame(matrix(NA, nrow=num.taxa, ncol=num.reps+1))
     names(obs.wads.by.taxon) <-c(vars[1], levels(test.data$rep))
     for (i in 1:num.taxa){
       obs.wads.by.taxon$taxon[i] <- levels(test.data$taxon)[i]
       for (r in 1:num.reps){
         obs.wads.by.taxon[i,(r+1)] <- WAD.func(y=test.data$y[test.data$taxon == levels(test.data$taxon)[i] & test.data$rep == levels(test.data$rep)[r]], x=test.data$x[test.data$taxon == levels(test.data$taxon)[i] & test.data$rep == levels(test.data$rep)[r]])
       }
     } 
     obs.wads.by.taxon$taxon <- factor(obs.wads.by.taxon$taxon)
     obs.wads.by.taxon <- obs.wads.by.taxon[order(as.numeric(as.character(obs.wads.by.taxon$taxon))),]
     row.names(obs.wads.by.taxon) <- 1:dim(obs.wads.by.taxon)[1]

    # Create a data.frame of replicate IDs for each treatment:
     num.trts <- length(levels(test.data$trt))
     reps.by.trt.list <- tapply(test.data$rep, test.data$trt, function(x) unique(as.character(x)))
     max.reps.trt <- max(sapply(reps.by.trt.list, length), na.rm=TRUE)
     reps.by.trt <- data.frame(matrix(NA, nrow=num.trts, ncol=max.reps.trt+1))
     names(reps.by.trt) <- c(vars[5], paste("R", 1:max.reps.trt, sep=""))
     reps.by.trt[,1] <- names(reps.by.trt.list)
     for (m in 1:num.trts){
       reps.by.trt[m,2:(2+length(reps.by.trt.list[[m]])-1)] <- reps.by.trt.list[[m]]
     }
     reps.by.trt <- data.frame(lapply(reps.by.trt, function(x) factor(as.character(x))), check.names=FALSE)  #convert columns to factors

  return(list(obs.wads=obs.wads.by.taxon, reps.by.trt=reps.by.trt))

  }
