bu.abs.resid <- function(DATA, lab.replicate=lab.replicate, min.num.nongrowers=min.num.nongrowers){
  curr.data.lowSD <- DATA
  lab.replicate <- lab.replicate
  min.num.nongrowers <- min.num.nongrowers
  #Create three empty data.frames for results:
  putative.nongrower.metrics <- data.frame(matrix(NA, 1 + dim(curr.data.lowSD)[1] - min.num.nongrowers, 29))
  names(putative.nongrower.metrics) <- c("next.taxon", "N", "mean.X", "sd.X", "norm.mean.X", "norm.sd.X", "norm.NLL.X", "fixed.intercept", "fixed.slope", "fixed.correl", "fixed.r2", "fixed.r2.adj", "fixed.resid.intercept", "fixed.resid.slope", "fixed.mean.Y", "fixed.sd.Y", "fixed.norm.mean.Y", "fixed.norm.sd.Y", "fixed.norm.NLL.Y", "corr.free.intercept", "corr.free.slope", "corr.free.correl", "corr.free.r2", "corr.free.r2.adj", "corr.free.mean.Y", "corr.free.sd.Y", "corr.free.norm.mean.Y", "corr.free.norm.sd.Y", "corr.free.norm.NLL.Y")
  taxa.in <- taxa.out <- data.frame(matrix(NA, dim(curr.data.lowSD)[1], 1 + dim(curr.data.lowSD)[1] - min.num.nongrowers))
  #Isolate the data and create the appropriate initial regression model:
  X.all <- curr.data.lowSD$unlab.mean
  Y.all <- curr.data.lowSD[, names(curr.data.lowSD) == lab.replicate]
  TAX.all <- curr.data.lowSD$taxon
  delta.WAD <- Y.all - X.all
  X <- X.all[delta.WAD %in% sort(delta.WAD)[1:min.num.nongrowers]]
  Y <- Y.all[delta.WAD %in% sort(delta.WAD)[1:min.num.nongrowers]]
  TAX <- curr.data.lowSD$taxon[delta.WAD %in% sort(delta.WAD)[1:min.num.nongrowers]]
  mod1 <- lm(I(Y-1*X)~1)         #fixed slope of 1; note that the function "I()" is used to inhibit the interpretation of operators such as "+", "-", "*" and "^" as formula operators, so they are instead used as arithmetical operators.
  #Iterate through the taxa:
  for (j in 0:(dim(curr.data.lowSD)[1] - min.num.nongrowers)){
    if (j == 0){
      nxt.index.to.add <- 0
      indices.to.keep <- which(Y.all %in% Y)
      putative.nongrower.metrics$next.taxon[j+1] <- names(taxa.out)[j+1] <- names(taxa.in)[j+1] <- "none"
      X <- X.all[indices.to.keep]
      Y <- Y.all[indices.to.keep]
      TAX <- TAX.all[indices.to.keep]
      taxa.in[1:length(TAX),j+1] <- as.character(TAX)
      taxa.out[1:(dim(taxa.out)[1]-min.num.nongrowers-j),j+1] <- as.character(TAX.all[which(!(Y.all %in% Y))])
      taxa.in[,j+1] <- factor(as.character(taxa.in[,j+1]))
      taxa.out[,j+1] <- factor(as.character(taxa.out[,j+1]))
    } else if (j != 0){
      Y.hats.all <- X.all + as.numeric(coef(mod1))
      Y.resids.all <- Y.all - Y.hats.all
      nxt.index.to.add <- which(abs(Y.resids.all) == min(abs(Y.resids.all[!(Y.all %in% Y)])))   #get the index of the taxon having an absolute value residual that matches the smallest absolute value residual of the taxa not already in the included taxa
      indices.to.keep <- sort(c(which(Y.all %in% Y), nxt.index.to.add))
      putative.nongrower.metrics$next.taxon[j+1] <- as.character(TAX.all[nxt.index.to.add])
      names(taxa.out)[j+1] <- names(taxa.in)[j+1] <- paste("taxon.", as.character(TAX.all[nxt.index.to.add]), sep="")
      X <- X.all[indices.to.keep]
      Y <- Y.all[indices.to.keep]
      TAX <- TAX.all[indices.to.keep]
      taxa.in[1:length(TAX),j+1] <- as.character(TAX)
      if (length(TAX.all[which(!(Y.all %in% Y))]) != 0){   #if there are still taxa left out, write them; if there are no taxa still left out, don't write anything (leave the column full of NAs)
        taxa.out[1:(dim(taxa.out)[1]-min.num.nongrowers-j),j+1] <- as.character(TAX.all[which(!(Y.all %in% Y))])
      }
      taxa.in[,j+1] <- factor(as.character(taxa.in[,j+1]))
      taxa.out[,j+1] <- factor(as.character(taxa.out[,j+1]))
    }
    #Fill the data frame of metrics:
    norm.fit.X <- fit.norm.func(x=sort(X), z=(1:length(X))/length(X), CI=0.90)                        #fit normal cdf
    norm.fit.Y <- fit.norm.func(x=sort(Y), z=(1:length(Y))/length(Y), CI=0.90)                        #fit normal cdf
    mod1 <- lm(I(Y-1*X)~1)                                                                            #fixed slope of 1; note that the function "I()" is used to inhibit the interpretation of operators such as "+", "-", "*" and "^" as formula operators, so they are instead used as arithmetical operators.
    Y.corr <- Y - (norm.fit.Y$mean - norm.fit.X$mean)                                                 #correct the labeled WADs using the difference in normal-fitted means of putative nongrowing taxa in labeled tubes and unlabeled tubes
    # Y.corr <- Y - as.numeric(coef(mod1))                                                              #correct the labeled WADs using the intercept term from the regression with fixed slope of 1
    norm.fit.Y.corr <- fit.norm.func(x=sort(Y.corr), z=(1:length(Y.corr))/length(Y.corr), CI=0.90)    #fit normal cdf
    mod1.resid <- lm(as.numeric(resid(mod1))~X)                                                       #regression of residuals on x-values
    mod2 <- lm(Y.corr~X)                                                                              #free intercept and slope of 'corrected' data
    putative.nongrower.metrics$N[j+1] <- length(Y)
    putative.nongrower.metrics$mean.X[j+1] <- mean(X)                                                 #empirical mean
    putative.nongrower.metrics$sd.X[j+1] <- sd(X)                                                     #empirical sd
    putative.nongrower.metrics$norm.mean.X[j+1] <- norm.fit.X$mean                                    #fitted normal distribution mean
    putative.nongrower.metrics$norm.sd.X[j+1] <- norm.fit.X$stdev                                     #fitted normal distribution sd
    putative.nongrower.metrics$norm.NLL.X[j+1] <- norm.fit.X$NLL                                      #NLL of fitted normal distribution
    putative.nongrower.metrics$fixed.intercept[j+1] <- as.numeric(coef(mod1))
    putative.nongrower.metrics$fixed.slope[j+1] <- 1
    putative.nongrower.metrics$fixed.correl[j+1] <- cor(X, Y)
    Y.hats <- X+as.numeric(coef(mod1))
    fixed.r2 <- 1-(sum((Y-Y.hats)^2)/sum((Y-mean(Y))^2))
    putative.nongrower.metrics$fixed.r2[j+1] <- fixed.r2
    putative.nongrower.metrics$fixed.r2.adj[j+1] <- 1-(((1-fixed.r2)*(length(Y)-1))/(length(Y)-length(coef(mod1))-1))
    putative.nongrower.metrics$fixed.resid.intercept[j+1] <- as.numeric(coef(mod1.resid)[1])          #intercept of residuals vs x-values of the regression with a fixed slope of 1
    putative.nongrower.metrics$fixed.resid.slope[j+1] <-  as.numeric(coef(mod1.resid)[2])             #slope of residuals vs x-values of the regression with a fixed slope of 1
    putative.nongrower.metrics$fixed.mean.Y[j+1] <- mean(Y)                                           #empirical mean
    putative.nongrower.metrics$fixed.sd.Y[j+1] <- sd(Y)                                               #empirical sd
    putative.nongrower.metrics$fixed.norm.mean.Y[j+1] <- norm.fit.Y$mean                              #fitted normal distribution mean
    putative.nongrower.metrics$fixed.norm.sd.Y[j+1] <- norm.fit.Y$stdev                               #fitted normal distribution sd
    putative.nongrower.metrics$fixed.norm.NLL.Y[j+1] <- norm.fit.Y$NLL                                #NLL of fitted normal distribution
    putative.nongrower.metrics$corr.free.intercept[j+1] <- as.numeric(coef(mod2))[1]
    putative.nongrower.metrics$corr.free.slope[j+1] <- as.numeric(coef(mod2))[2]
    putative.nongrower.metrics$corr.free.correl[j+1] <- cor(X, Y.corr)
    putative.nongrower.metrics$corr.free.r2[j+1] <- summary(mod2)$r.squared
    putative.nongrower.metrics$corr.free.r2.adj[j+1] <- summary(mod2)$adj.r.squared
    putative.nongrower.metrics$corr.free.mean.Y[j+1] <- mean(Y.corr)                                  #empirical mean
    putative.nongrower.metrics$corr.free.sd.Y[j+1] <- sd(Y.corr)                                      #empirical sd
    putative.nongrower.metrics$corr.free.norm.mean.Y[j+1] <- norm.fit.Y.corr$mean                     #fitted normal distribution mean
    putative.nongrower.metrics$corr.free.norm.sd.Y[j+1] <- norm.fit.Y.corr$stdev                      #fitted normal distribution sd
    putative.nongrower.metrics$corr.free.norm.NLL.Y[j+1] <- norm.fit.Y.corr$NLL                       #NLL of fitted normal distribution
  }
  list(putative.nongrower.metrics=putative.nongrower.metrics, taxa.in=taxa.in, taxa.out=taxa.out, data.low.SD=curr.data.lowSD)
}



