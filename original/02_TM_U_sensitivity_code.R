# This code analyzes the sensitivity of the U parameter, constrains U, and compares growth models for the TM qSIP data


graphics.off()	#close all graphics windows


#Set working directory:
  #ALREADY DONE FOR THIS WORKSPACE IN PREVIOUSLY RUN CODE; only reset it here if loading previously saved workspace (see below)


#Reload the saved workspace resulting from the previous script:
  setwd("/Users/bk/Research/Projects/SIP_Modeling/qSIP")
  load("qSIP_workspaces/TM_01/.RData")


#Load libraries & scripts:
  source("qSIP_repo/U.sens.func.R")         #U.sens.func
  source("qSIP_repo/U.sens.plot.func.R")    #U.sens.plot.func


#Import colors for phyla for use below:
  cols.for.phyla <- read.table("qSIP_data/TM_PhylaColors.txt", header=T, sep="\t", comment.char="")
  phyla.cols <- cols.for.phyla$col


#Quick look at the sensitivity of the results to the parameter U (taxon 190):__________________________________________________________________________________
  #{
    #EXPONENTIAL GROWTH MODEL:
    #Calculate sensitivity to U (observed) EXPONENTIAL GROWTH MODEL:
      T190T1 <- data.melted[data.melted$taxon==190 & data.melted$tmt=="16O",]
      T190T2 <- data.melted[data.melted$taxon==190 & data.melted$tmt=="18O",]
      T190ref <- data.melted[data.melted$taxon==190 & (data.melted$tmt=="16O"),]
      T190T0.tube <- ncopies.tube.melted[ncopies.tube.melted$taxon==190 & ncopies.tube.melted$tmt=="Time0",]
      T190T1.tube <- ncopies.tube.melted[ncopies.tube.melted$taxon==190 & ncopies.tube.melted$tmt=="16O",]
      T190T2.tube <- ncopies.tube.melted[ncopies.tube.melted$taxon==190 & ncopies.tube.melted$tmt=="18O",]
      T190T0v12.r.out <- boot.r.pop(T0=T190T0.tube, Tt=rbind(T190T1.tube, T190T2.tube), M.soil=Sdat, days=10, vars=c("copies", "tube", "tmt", "g.soil"), growth.model="exponential", vol=c(100, 50), CI=0.90, draws=1000)
      T190T0.tube.copies.out <- boot.TUBE.pop(X=T190T0.tube, M.soil=Sdat, vars=c("copies", "tube", "g.soil"), vol=100, CI=0.90, draws=1000)
      T190T12.tube.copies.out <- boot.TUBE.pop(X=rbind(T190T1.tube, T190T2.tube), M.soil=Sdat, vars=c("copies", "tube", "g.soil"), vol=50, CI=0.90, draws=1000)
      U.r.sens.T190.obs <- data.frame(U=seq(0.01, 1, 0.001))
      U.r.sens.T190.obs$tot.copies.g.soil.t0 <- mean(T190T0.tube.copies.out$obs.N$tot.copies / T190T0.tube.copies.out$obs.N$g.soil)
      U.r.sens.T190.obs$tot.copies.g.soil.t10 <- mean(T190T12.tube.copies.out$obs.N$tot.copies / T190T12.tube.copies.out$obs.N$g.soil)
      U.r.sens.T190.obs$light.copies.g.soil.t10 <- rep(NA, length(U.r.sens.T190.obs$U))
      U.r.sens.T190.obs$r.net <- T190T0v12.r.out$obs.r
      for (i in 1:length(U.r.sens.T190.obs$U)){
        T190T1v2.r.out <- boot.diff.r(X.light=T190T1, X.heavy=T190T2, X.reference=T190ref, M.soil=Sdat, iso.compare="18O", days=10, vars=c("density.g.ml", "copies", "tube", "tmt", "g.soil"), growth.model="exponential", prop.O.from.water=U.r.sens.T190.obs$U[i], v.frac=50, CI=0.90, draws=1000)
        U.r.sens.T190.obs$light.copies.g.soil.t10[i] <- U.r.sens.T190.obs$tot.copies.g.soil.t10[i] * (1 / exp(T190T1v2.r.out$obs.r * 10))    #transformation if exponential growth model was used
        U.r.sens.T190.obs$r.gross[i] <- T190T1v2.r.out$obs.r
        U.r.sens.T190.obs$d[i] <- T190T0v12.r.out$obs.r - T190T1v2.r.out$obs.r
      }

    #Graph sensitivity of abundances to U (observed) EXPONENTIAL GROWTH MODEL:
      dev.off()    
      # dev.new(width=7.5, height=7.5)
      pdf(file="qSIP_output/Figures/TM_abundance_sensitivity_to_U_taxon190_obs_exponential.pdf", width=7.5, height=7.5)
      par(mai=c(0.44,0.63,0.2,0.05), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
      y.values <- c(U.r.sens.T190.obs$tot.copies.g.soil.t0, U.r.sens.T190.obs$tot.copies.g.soil.t10, U.r.sens.T190.obs$light.copies.g.soil.t10)
      y.values <- y.values[y.values != -Inf & y.values != Inf & !is.na(y.values)]
      y.min <- min(y.values)
      y.max <- max(y.values)*1.1

      plot(y=U.r.sens.T190.obs$tot.copies.g.soil.t10, x=U.r.sens.T190.obs$U, type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(0, 1), ylim=c(y.min, y.max), main="Taxon 190")
        points(y=U.r.sens.T190.obs$tot.copies.g.soil.t10, x=U.r.sens.T190.obs$U, pch=21, cex=0.6, col=as.character(phyla.cols[12]), bg=as.character(phyla.cols[12]))
        points(y=U.r.sens.T190.obs$tot.copies.g.soil.t0, x=U.r.sens.T190.obs$U, pch=21, cex=0.6, col="green", bg="green")
        points(y=U.r.sens.T190.obs$light.copies.g.soil.t10, x=U.r.sens.T190.obs$U, pch=21, cex=0.6, col="gray30", bg="gray30")
        abline(h=0, col="black")
        par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
        axis(side=1, tck=-0.008, las=1, cex.axis=0.6)
        mtext("U", side=1, line=1.35, cex=0.75)
        par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
        axis(side=2, tck=-0.008, las=1, cex.axis=0.6)
        mtext(expression(paste("Abundance (16S copies g soil"^-1, ")", sep="")), side=2, line=2.4, cex=0.75)
        legend(x=1, y=y.max, legend=c("total, time 10", "total (unlabeled), time 0", "unlabeled, time 10"), pch=21, col=c(as.character(phyla.cols[12]), "green", "gray30"), pt.bg=c(as.character(phyla.cols[12]), "green", "gray30"), bty="o", cex=0.6, xjust=1, yjust=1)

      dev.off()
  
    #Graph sensitivity of growth rates to U (observed) EXPONENTIAL GROWTH MODEL:
      dev.off()    
      # dev.new(width=7.5, height=7.5)
      pdf(file="qSIP_output/Figures/TM_growth_sensitivity_to_U_taxon190_obs_exponential.pdf", width=7.5, height=7.5)
      par(mai=c(0.44,0.63,0.2,0.05), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
      y.values <- c(U.r.sens.T190.obs$r.net, U.r.sens.T190.obs$r.gross, U.r.sens.T190.obs$d)
      y.values <- y.values[y.values != -Inf & y.values != Inf & !is.na(y.values)]
      # y.min <- min(y.values)
      # y.max <- max(y.values)
      y.min <- -1
      y.max <- 1

      plot(y=U.r.sens.T190.obs$r.net, x=U.r.sens.T190.obs$U, type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(0, 1), ylim=c(y.min, y.max), main="Taxon 190")
        points(y=U.r.sens.T190.obs$r.net, x=U.r.sens.T190.obs$U, pch=21, cex=0.6, col=as.character(phyla.cols[12]), bg=as.character(phyla.cols[12]))
        points(y=U.r.sens.T190.obs$r.gross, x=U.r.sens.T190.obs$U, pch=21, cex=0.6, col=as.character(phyla.cols[12]), bg="white")
        points(y=U.r.sens.T190.obs$d, x=U.r.sens.T190.obs$U, pch=21, cex=0.6, col="gray30", bg=as.character(phyla.cols[12]))
        abline(h=0, col="black")
        par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
        axis(side=1, tck=-0.008, las=1, cex.axis=0.6)
        mtext("U", side=1, line=1.35, cex=0.75)
        par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
        axis(side=2, tck=-0.008, las=1, cex.axis=0.6)
        mtext(expression(paste("day"^-1, sep="")), side=2, line=2.4, cex=0.75)                         #units if exponential growth model was used
        legend(x=1, y=y.max, legend=c("net", "gross", "death"), pch=21, col=c(as.character(phyla.cols[12]), as.character(phyla.cols[12]), "gray30"), pt.bg=c(as.character(phyla.cols[12]), "white", as.character(phyla.cols[12])), bty="o", cex=0.6, xjust=1, yjust=1)

      dev.off()
  
    #Calculate sensitivity to U (bootstrapping) EXPONENTIAL GROWTH MODEL:
      T190T1 <- data.melted[data.melted$taxon==190 & data.melted$tmt=="16O",]
      T190T2 <- data.melted[data.melted$taxon==190 & data.melted$tmt=="18O",]
      T190ref <- data.melted[data.melted$taxon==190 & (data.melted$tmt=="16O"),]
      T190T0.tube <- ncopies.tube.melted[ncopies.tube.melted$taxon==190 & ncopies.tube.melted$tmt=="Time0",]
      T190T1.tube <- ncopies.tube.melted[ncopies.tube.melted$taxon==190 & ncopies.tube.melted$tmt=="16O",]
      T190T2.tube <- ncopies.tube.melted[ncopies.tube.melted$taxon==190 & ncopies.tube.melted$tmt=="18O",]
      T190T0v12.r.out <- boot.r.pop(T0=T190T0.tube, Tt=rbind(T190T1.tube, T190T2.tube), M.soil=Sdat, days=10, vars=c("copies", "tube", "tmt", "g.soil"), growth.model="exponential", vol=c(100, 50), CI=0.90, draws=1000)
      U.r.sens.T190.boot <- data.frame(U=seq(0.01, 1, 0.001))
      U.r.sens.T190.boot$r.net <- T190T0v12.r.out$boot.r.median
      for (i in 1:length(U.r.sens.T190.boot$U)){
        T190T1v2.r.out <- boot.diff.r(X.light=T190T1, X.heavy=T190T2, X.reference=T190ref, M.soil=Sdat, iso.compare="18O", days=10, vars=c("density.g.ml", "copies", "tube", "tmt", "g.soil"), growth.model="exponential", prop.O.from.water=U.r.sens.T190.boot$U[i], v.frac=50, CI=0.90, draws=1000)
        U.r.sens.T190.boot$r.gross[i] <- T190T1v2.r.out$boot.r.median
        U.r.sens.T190.boot$d[i] <- median(T190T0v12.r.out$boot.r - T190T1v2.r.out$boot.r, na.rm=TRUE)
        U.r.sens.T190.boot$d.CI.L[i] <- quantile(T190T0v12.r.out$boot.r - T190T1v2.r.out$boot.r, probs=0.05, na.rm=TRUE)
        U.r.sens.T190.boot$d.CI.U[i] <- quantile(T190T0v12.r.out$boot.r - T190T1v2.r.out$boot.r, probs=0.95, na.rm=TRUE)
      }

    #Graph sensitivity of growth rates to U (bootstrapping) EXPONENTIAL GROWTH MODEL:
      dev.off()    
      # dev.new(width=7.5, height=7.5)
      pdf(file="qSIP_output/Figures/TM_growth_sensitivity_to_U_taxon190_boots_exponential.pdf", width=7.5, height=7.5)
      par(mai=c(0.44,0.63,0.2,0.05), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
      y.values <- c(U.r.sens.T190.boot$r.net, U.r.sens.T190.boot$r.gross, U.r.sens.T190.boot$d, U.r.sens.T190.boot$d.CI.L, U.r.sens.T190.boot$d.CI.U)
      y.values <- y.values[y.values != -Inf & y.values != Inf & !is.na(y.values)]
      # y.min <- min(y.values)    #if linear growth model was used
      # y.max <- max(y.values)    #if linear growth model was used
      y.min <- -1    #if exponential growth model was used
      y.max <- 1     #if exponential growth model was used

      plot(y=U.r.sens.T190.boot$r.net, x=U.r.sens.T190.boot$U, type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(0, 1), ylim=c(y.min, y.max), main="Taxon 190")
        arrows(x0=U.r.sens.T190.boot$U, y0=U.r.sens.T190.boot$d.CI.L, x1=U.r.sens.T190.boot$U, y1=U.r.sens.T190.boot$d.CI.U, length=0, angle=90, code=3, col="gray30", lwd=0.5)
        points(y=U.r.sens.T190.boot$d, x=U.r.sens.T190.boot$U, pch=21, cex=0.6, col="gray30", bg=as.character(phyla.cols[12]))
        points(y=U.r.sens.T190.boot$r.net, x=U.r.sens.T190.boot$U, pch=21, cex=0.6, col=as.character(phyla.cols[12]), bg=as.character(phyla.cols[12]))
        points(y=U.r.sens.T190.boot$r.gross, x=U.r.sens.T190.boot$U, pch=21, cex=0.6, col=as.character(phyla.cols[12]), bg="white")
        abline(h=0, col="black")
        par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
        axis(side=1, tck=-0.008, las=1, cex.axis=0.6)
        mtext("U", side=1, line=1.35, cex=0.75)
        par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
        axis(side=2, tck=-0.008, las=1, cex.axis=0.6)
        mtext(expression(paste("day"^-1, sep="")), side=2, line=2.4, cex=0.75)                         #if exponential growth model was used
        # mtext(expression(paste("copies g soil"^-1, " day"^-1, sep="")), side=2, line=2.4, cex=0.75)    #if linear growth model was used
        legend(x=1, y=y.max, legend=c("net", "gross", "death"), pch=21, col=c(as.character(phyla.cols[12]), as.character(phyla.cols[12]), "gray30"), pt.bg=c(as.character(phyla.cols[12]), "white", as.character(phyla.cols[12])), bty="o", cex=0.6, xjust=1, yjust=1)

      dev.off()


    #LINEAR GROWTH MODEL:
    #Calculate sensitivity to U (observed) LINEAR GROWTH MODEL:
      T190T1 <- data.melted[data.melted$taxon==190 & data.melted$tmt=="16O",]
      T190T2 <- data.melted[data.melted$taxon==190 & data.melted$tmt=="18O",]
      T190ref <- data.melted[data.melted$taxon==190 & (data.melted$tmt=="16O"),]
      T190T0.tube <- ncopies.tube.melted[ncopies.tube.melted$taxon==190 & ncopies.tube.melted$tmt=="Time0",]
      T190T1.tube <- ncopies.tube.melted[ncopies.tube.melted$taxon==190 & ncopies.tube.melted$tmt=="16O",]
      T190T2.tube <- ncopies.tube.melted[ncopies.tube.melted$taxon==190 & ncopies.tube.melted$tmt=="18O",]
      T190T0v12.r.out <- boot.r.pop(T0=T190T0.tube, Tt=rbind(T190T1.tube, T190T2.tube), M.soil=Sdat, days=10, vars=c("copies", "tube", "tmt", "g.soil"), growth.model="linear", vol=c(100, 50), CI=0.90, draws=1000)
      T190T0.tube.copies.out <- boot.TUBE.pop(X=T190T0.tube, M.soil=Sdat, vars=c("copies", "tube", "g.soil"), vol=100, CI=0.90, draws=1000)
      T190T12.tube.copies.out <- boot.TUBE.pop(X=rbind(T190T1.tube, T190T2.tube), M.soil=Sdat, vars=c("copies", "tube", "g.soil"), vol=50, CI=0.90, draws=1000)
      U.r.sens.T190.obs <- data.frame(U=seq(0.01, 1, 0.001))
      U.r.sens.T190.obs$tot.copies.g.soil.t0 <- mean(T190T0.tube.copies.out$obs.N$tot.copies / T190T0.tube.copies.out$obs.N$g.soil)
      U.r.sens.T190.obs$tot.copies.g.soil.t10 <- mean(T190T12.tube.copies.out$obs.N$tot.copies / T190T12.tube.copies.out$obs.N$g.soil)
      U.r.sens.T190.obs$light.copies.g.soil.t10 <- rep(NA, length(U.r.sens.T190.obs$U))
      U.r.sens.T190.obs$r.net <- T190T0v12.r.out$obs.r
      for (i in 1:length(U.r.sens.T190.obs$U)){
        T190T1v2.r.out <- boot.diff.r(X.light=T190T1, X.heavy=T190T2, X.reference=T190ref, M.soil=Sdat, iso.compare="18O", days=10, vars=c("density.g.ml", "copies", "tube", "tmt", "g.soil"), growth.model="linear", prop.O.from.water=U.r.sens.T190.obs$U[i], v.frac=50, CI=0.90, draws=1000)
        U.r.sens.T190.obs$light.copies.g.soil.t10[i] <- U.r.sens.T190.obs$tot.copies.g.soil.t10[i] - (T190T1v2.r.out$obs.r * 10)             #transformation if linear growth model was used
        U.r.sens.T190.obs$r.gross[i] <- T190T1v2.r.out$obs.r
        U.r.sens.T190.obs$d[i] <- T190T0v12.r.out$obs.r - T190T1v2.r.out$obs.r
      }

    #Graph sensitivity of abundances to U (observed) LINEAR GROWTH MODEL:
      dev.off()    
      # dev.new(width=7.5, height=7.5)
      pdf(file="qSIP_output/Figures/TM_abundance_sensitivity_to_U_taxon190_obs_linear.pdf", width=7.5, height=7.5)
      par(mai=c(0.44,0.63,0.2,0.05), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
      y.values <- c(U.r.sens.T190.obs$tot.copies.g.soil.t0, U.r.sens.T190.obs$tot.copies.g.soil.t10, U.r.sens.T190.obs$light.copies.g.soil.t10)
      y.values <- y.values[y.values != -Inf & y.values != Inf & !is.na(y.values)]
      y.min <- min(y.values)
      y.max <- max(y.values)*1.1

      plot(y=U.r.sens.T190.obs$tot.copies.g.soil.t10, x=U.r.sens.T190.obs$U, type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(0, 1), ylim=c(y.min, y.max), main="Taxon 190")
        points(y=U.r.sens.T190.obs$tot.copies.g.soil.t10, x=U.r.sens.T190.obs$U, pch=21, cex=0.6, col=as.character(phyla.cols[12]), bg=as.character(phyla.cols[12]))
        points(y=U.r.sens.T190.obs$tot.copies.g.soil.t0, x=U.r.sens.T190.obs$U, pch=21, cex=0.6, col="green", bg="green")
        points(y=U.r.sens.T190.obs$light.copies.g.soil.t10, x=U.r.sens.T190.obs$U, pch=21, cex=0.6, col="gray30", bg="gray30")
        abline(h=0, col="black")
        par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
        axis(side=1, tck=-0.008, las=1, cex.axis=0.6)
        mtext("U", side=1, line=1.35, cex=0.75)
        par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
        axis(side=2, tck=-0.008, las=1, cex.axis=0.6)
        mtext(expression(paste("Abundance (16S copies g soil"^-1, ")", sep="")), side=2, line=2.4, cex=0.75)
        legend(x=1, y=y.max, legend=c("total, time 10", "total (unlabeled), time 0", "unlabeled, time 10"), pch=21, col=c(as.character(phyla.cols[12]), "green", "gray30"), pt.bg=c(as.character(phyla.cols[12]), "green", "gray30"), bty="o", cex=0.6, xjust=1, yjust=1)

      dev.off()
  
    #Graph sensitivity of growth rates to U (observed) LINEAR GROWTH MODEL:
      dev.off()    
      # dev.new(width=7.5, height=7.5)
      pdf(file="qSIP_output/Figures/TM_growth_sensitivity_to_U_taxon190_obs_linear.pdf", width=7.5, height=7.5)
      par(mai=c(0.44,0.63,0.2,0.05), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
      y.values <- c(U.r.sens.T190.obs$r.net, U.r.sens.T190.obs$r.gross, U.r.sens.T190.obs$d)
      y.values <- y.values[y.values != -Inf & y.values != Inf & !is.na(y.values)]
      y.min <- min(y.values)
      y.max <- max(y.values)

      plot(y=U.r.sens.T190.obs$r.net, x=U.r.sens.T190.obs$U, type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(0, 1), ylim=c(y.min, y.max), main="Taxon 190")
        points(y=U.r.sens.T190.obs$r.net, x=U.r.sens.T190.obs$U, pch=21, cex=0.6, col=as.character(phyla.cols[12]), bg=as.character(phyla.cols[12]))
        points(y=U.r.sens.T190.obs$r.gross, x=U.r.sens.T190.obs$U, pch=21, cex=0.6, col=as.character(phyla.cols[12]), bg="white")
        points(y=U.r.sens.T190.obs$d, x=U.r.sens.T190.obs$U, pch=21, cex=0.6, col="gray30", bg=as.character(phyla.cols[12]))
        abline(h=0, col="black")
        par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
        axis(side=1, tck=-0.008, las=1, cex.axis=0.6)
        mtext("U", side=1, line=1.35, cex=0.75)
        par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
        axis(side=2, tck=-0.008, las=1, cex.axis=0.6)
        mtext(expression(paste("copies g soil"^-1, " day"^-1, sep="")), side=2, line=2.4, cex=0.75)    #units if linear growth model was used
        legend(x=1, y=y.max, legend=c("net", "gross", "death"), pch=21, col=c(as.character(phyla.cols[12]), as.character(phyla.cols[12]), "gray30"), pt.bg=c(as.character(phyla.cols[12]), "white", as.character(phyla.cols[12])), bty="o", cex=0.6, xjust=1, yjust=1)

      dev.off()
  
    #Calculate sensitivity to U (bootstrapping) LINEAR GROWTH MODEL:
      T190T1 <- data.melted[data.melted$taxon==190 & data.melted$tmt=="16O",]
      T190T2 <- data.melted[data.melted$taxon==190 & data.melted$tmt=="18O",]
      T190ref <- data.melted[data.melted$taxon==190 & (data.melted$tmt=="16O"),]
      T190T0.tube <- ncopies.tube.melted[ncopies.tube.melted$taxon==190 & ncopies.tube.melted$tmt=="Time0",]
      T190T1.tube <- ncopies.tube.melted[ncopies.tube.melted$taxon==190 & ncopies.tube.melted$tmt=="16O",]
      T190T2.tube <- ncopies.tube.melted[ncopies.tube.melted$taxon==190 & ncopies.tube.melted$tmt=="18O",]
      T190T0v12.r.out <- boot.r.pop(T0=T190T0.tube, Tt=rbind(T190T1.tube, T190T2.tube), M.soil=Sdat, days=10, vars=c("copies", "tube", "tmt", "g.soil"), growth.model="linear", vol=c(100, 50), CI=0.90, draws=1000)
      U.r.sens.T190.boot <- data.frame(U=seq(0.01, 1, 0.001))
      U.r.sens.T190.boot$r.net <- T190T0v12.r.out$boot.r.median
      for (i in 1:length(U.r.sens.T190.boot$U)){
        T190T1v2.r.out <- boot.diff.r(X.light=T190T1, X.heavy=T190T2, X.reference=T190ref, M.soil=Sdat, iso.compare="18O", days=10, vars=c("density.g.ml", "copies", "tube", "tmt", "g.soil"), growth.model="linear", prop.O.from.water=U.r.sens.T190.boot$U[i], v.frac=50, CI=0.90, draws=1000)
        U.r.sens.T190.boot$r.gross[i] <- T190T1v2.r.out$boot.r.median
        U.r.sens.T190.boot$d[i] <- median(T190T0v12.r.out$boot.r - T190T1v2.r.out$boot.r, na.rm=TRUE)
        U.r.sens.T190.boot$d.CI.L[i] <- quantile(T190T0v12.r.out$boot.r - T190T1v2.r.out$boot.r, probs=0.05, na.rm=TRUE)
        U.r.sens.T190.boot$d.CI.U[i] <- quantile(T190T0v12.r.out$boot.r - T190T1v2.r.out$boot.r, probs=0.95, na.rm=TRUE)
      }

    #Graph sensitivity of growth rates to U (bootstrapping) LINEAR GROWTH MODEL:
      dev.off()    
      # dev.new(width=7.5, height=7.5)
      pdf(file="qSIP_output/Figures/TM_growth_sensitivity_to_U_taxon190_boots_linear.pdf", width=7.5, height=7.5)
      par(mai=c(0.44,0.63,0.2,0.05), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
      y.values <- c(U.r.sens.T190.boot$r.net, U.r.sens.T190.boot$r.gross, U.r.sens.T190.boot$d, U.r.sens.T190.boot$d.CI.L, U.r.sens.T190.boot$d.CI.U)
      y.values <- y.values[y.values != -Inf & y.values != Inf & !is.na(y.values)]
      y.min <- min(y.values)    #if linear growth model was used
      y.max <- max(y.values)    #if linear growth model was used
      # y.min <- -1    #if exponential growth model was used
      # y.max <- 1     #if exponential growth model was used

      plot(y=U.r.sens.T190.boot$r.net, x=U.r.sens.T190.boot$U, type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(0, 1), ylim=c(y.min, y.max), main="Taxon 190")
        arrows(x0=U.r.sens.T190.boot$U, y0=U.r.sens.T190.boot$d.CI.L, x1=U.r.sens.T190.boot$U, y1=U.r.sens.T190.boot$d.CI.U, length=0, angle=90, code=3, col="gray30", lwd=0.5)
        points(y=U.r.sens.T190.boot$d, x=U.r.sens.T190.boot$U, pch=21, cex=0.6, col="gray30", bg=as.character(phyla.cols[12]))
        points(y=U.r.sens.T190.boot$r.net, x=U.r.sens.T190.boot$U, pch=21, cex=0.6, col=as.character(phyla.cols[12]), bg=as.character(phyla.cols[12]))
        points(y=U.r.sens.T190.boot$r.gross, x=U.r.sens.T190.boot$U, pch=21, cex=0.6, col=as.character(phyla.cols[12]), bg="white")
        abline(h=0, col="black")
        par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
        axis(side=1, tck=-0.008, las=1, cex.axis=0.6)
        mtext("U", side=1, line=1.35, cex=0.75)
        par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
        axis(side=2, tck=-0.008, las=1, cex.axis=0.6)
        # mtext(expression(paste("day"^-1, sep="")), side=2, line=2.4, cex=0.75)                         #if exponential growth model was used
        mtext(expression(paste("copies g soil"^-1, " day"^-1, sep="")), side=2, line=2.4, cex=0.75)    #if linear growth model was used
        legend(x=1, y=y.max, legend=c("net", "gross", "death"), pch=21, col=c(as.character(phyla.cols[12]), as.character(phyla.cols[12]), "gray30"), pt.bg=c(as.character(phyla.cols[12]), "white", as.character(phyla.cols[12])), bty="o", cex=0.6, xjust=1, yjust=1)

      dev.off()
  #}
  

#More in-depth analysis of the sensitivity to U (all taxa)...USING AN EXPONENTIAL GROWTH MODEL:________________________________________________________________
  #{
    #Create dataframes (two per taxon) of relevant variables across a range of values for U:
      #Example for taxon 190:
      set.seed(100)
      system.time(U.sens.T190.exp <- U.sens.func(DATA=data.melted, DATA.POP=ncopies.tube.melted, taxonID=190, time0="Time0", tmt1="16O", tmt2="18O", M.soil=Sdat, iso.compare="18O", days=10, vars=c("taxon", "density.g.ml", "copies", "tube", "tmt", "g.soil"), growth.model="exponential", vol=c(100, 50), CI=0.90, draws=1000, sens.seq=seq(0.01, 1, 0.01)))

    #Create graphs of relevant variables across a range of values for U (for a specific taxon):
      #Example for taxon 190:
      #note: set a priori value of U=0.33 according to McHugh et al. unpublished data from E. coli pure cultures incubated with 18O-H2O
      graphics.off()    
      # dev.new(width=8.5, height=11)
      pdf(file="qSIP_output/Figures/TM_U_sensitivity_taxon190_exponential.pdf", width=8.5, height=11)
      U.sens.plot.func(LIST=U.sens.T190.exp, taxonID=190, U=0.33, growth.model="exponential", selected.Us=seq(0.1,1,length.out=10))

      dev.off()

    #Conduct sensitivity analysis for all taxa; write output dataframes to files (two for each taxon); create a summary table of sensitivity results for all taxa:
      #Create a directory (if it doesn't exist already) for writing the taxon-specific data frames with sensitivity results:
        dir.create(path=paste(getwd(), "/qSIP_output/TM_Sensitivity_Exponential", sep=""), showWarnings=FALSE)
      #Set desired target U for analysis:
        #note: this is based on U=0.33 according to McHugh et al. unpublished data from E. coli pure cultures incubated with 18O-H2O
        U.target <- 0.33
      #Create a dataframe containing summary info from sensitivity analysis for all taxa:
        U.sens.all.taxa <- data.frame(taxon=levels(data.melted$taxon))
        U.sens.all.taxa$tot.copies.g.soil.t0.boot.CI.U <- U.sens.all.taxa$tot.copies.g.soil.t0.boot.CI.L <- U.sens.all.taxa$tot.copies.g.soil.t0.boot.median <- U.sens.all.taxa$tot.copies.g.soil.t0.obs <- rep(NA, length(U.sens.all.taxa$taxon))
        U.sens.all.taxa$light.copies.g.soil.tT.boot.CI.U.at.U.target <- U.sens.all.taxa$light.copies.g.soil.tT.boot.CI.L.at.U.target <- U.sens.all.taxa$light.copies.g.soil.tT.boot.median.at.U.target <- U.sens.all.taxa$light.copies.g.soil.tT.obs.at.U.target <- rep(NA, length(U.sens.all.taxa$taxon))
        U.sens.all.taxa$Tabun.less.equal.0abun.at.U.target.boot.median <- U.sens.all.taxa$Tabun.less.equal.0abun.at.U.target.obs <- rep(FALSE, length(U.sens.all.taxa$taxon))
        U.sens.all.taxa$U.max.Tabun.less.equal.0abun.at.U.target.obs <- U.sens.all.taxa$U.min.Tabun.less.equal.0abun.at.U.target.obs <- rep(NA, length(U.sens.all.taxa$taxon))
        U.sens.all.taxa$U.max.Tabun.less.equal.0abun.at.U.target.boot.median <- U.sens.all.taxa$U.min.Tabun.less.equal.0abun.at.U.target.boot.median <- rep(NA, length(U.sens.all.taxa$taxon))
      #Cycle through all taxa to create summary info from sensitivity analysis:
        set.seed(100)
        for(i in 1:length(U.sens.all.taxa$taxon)){
          curr.taxon <- as.character(U.sens.all.taxa$taxon[i])
          curr.U.sens.list <- U.sens.func(DATA=data.melted, DATA.POP=ncopies.tube.melted, taxonID=as.numeric(curr.taxon), time0="Time0", tmt1="16O", tmt2="18O", M.soil=Sdat, iso.compare="18O", days=10, vars=c("taxon", "density.g.ml", "copies", "tube", "tmt", "g.soil"), growth.model="exponential", vol=c(100, 50), CI=0.90, draws=1000, sens.seq=seq(0.01, 1, 0.01))
          U.target.index <- which(abs(curr.U.sens.list$U.sens$U-U.target) == min(abs(curr.U.sens.list$U.sens$U-U.target)))
          U.sens.all.taxa$tot.copies.g.soil.t0.obs[i] <- unique(curr.U.sens.list$U.sens$tot.copies.g.soil.t0.obs)
          U.sens.all.taxa$tot.copies.g.soil.t0.boot.median[i] <- unique(curr.U.sens.list$U.sens$tot.copies.g.soil.t0.boot.median)
          U.sens.all.taxa$tot.copies.g.soil.t0.boot.CI.L[i] <- unique(curr.U.sens.list$U.sens$tot.copies.g.soil.t0.boot.CI.L)
          U.sens.all.taxa$tot.copies.g.soil.t0.boot.CI.U[i] <- unique(curr.U.sens.list$U.sens$tot.copies.g.soil.t0.boot.CI.U)
          U.sens.all.taxa$light.copies.g.soil.tT.obs.at.U.target[i] <- curr.U.sens.list$U.sens$light.copies.g.soil.tT.obs[U.target.index]
          U.sens.all.taxa$light.copies.g.soil.tT.boot.median.at.U.target[i] <- curr.U.sens.list$U.sens$light.copies.g.soil.tT.boot.median[U.target.index]
          U.sens.all.taxa$light.copies.g.soil.tT.boot.CI.L.at.U.target[i] <- curr.U.sens.list$U.sens$light.copies.g.soil.tT.boot.CI.L[U.target.index]
          U.sens.all.taxa$light.copies.g.soil.tT.boot.CI.U.at.U.target[i] <- curr.U.sens.list$U.sens$light.copies.g.soil.tT.boot.CI.U[U.target.index]
          U.sens.all.taxa$Tabun.less.equal.0abun.at.U.target.obs[i] <- U.sens.all.taxa$light.copies.g.soil.tT.obs.at.U.target[i] <= U.sens.all.taxa$tot.copies.g.soil.t0.obs[i]
          U.sens.all.taxa$Tabun.less.equal.0abun.at.U.target.boot.median[i] <- U.sens.all.taxa$light.copies.g.soil.tT.boot.median.at.U.target[i] <= U.sens.all.taxa$tot.copies.g.soil.t0.boot.median[i]
          U.range.obs.indices <- which(curr.U.sens.list$U.sens$light.copies.g.soil.tT.obs <= curr.U.sens.list$U.sens$tot.copies.g.soil.t0.obs)
          if(length(U.range.obs.indices) > 0){
            U.sens.all.taxa$U.min.Tabun.less.equal.0abun.at.U.target.obs[i] <- min(curr.U.sens.list$U.sens$U[U.range.obs.indices])
            U.sens.all.taxa$U.max.Tabun.less.equal.0abun.at.U.target.obs[i] <- max(curr.U.sens.list$U.sens$U[U.range.obs.indices])
          }
          U.range.boot.indices <- which(curr.U.sens.list$U.sens$light.copies.g.soil.tT.boot.median <= curr.U.sens.list$U.sens$tot.copies.g.soil.t0.boot.median)
          if(length(U.range.boot.indices) > 0){
            U.sens.all.taxa$U.min.Tabun.less.equal.0abun.at.U.target.boot.median[i] <- min(curr.U.sens.list$U.sens$U[U.range.boot.indices])
            U.sens.all.taxa$U.max.Tabun.less.equal.0abun.at.U.target.boot.median[i] <- max(curr.U.sens.list$U.sens$U[U.range.boot.indices])
          }
          #Write the results (curr.U.sens.list$U.sens & curr.U.sens.list$U.sens.MW) to a text file:
          write.table(curr.U.sens.list$U.sens, paste("qSIP_output/TM_Sensitivity_Exponential/TM_U_sens_taxon_", sprintf("%03d", as.numeric(curr.taxon)), ".txt", sep=""), append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
          write.table(curr.U.sens.list$U.sens.MW, paste("qSIP_output/TM_Sensitivity_Exponential/TM_U_sens_MW_taxon_", sprintf("%03d", as.numeric(curr.taxon)), ".txt", sep=""), append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
        }
      #Description of the variables in the summary table (U.sens.all.taxa):
        # taxon                                                 ...  taxon
        # tot.copies.g.soil.t0.obs                              ...  Q.light.0.obs
        # tot.copies.g.soil.t0.boot.median                      ...  Q.light.0.boot.median
        # tot.copies.g.soil.t0.boot.CI.L                        ...  Q.light.0.boot.CI.L
        # tot.copies.g.soil.t0.boot.CI.U                        ...  Q.light.0.boot.CI.U
        # light.copies.g.soil.tT.obs.at.U.target                ...  Q.light.10.obs (calculated at the target value of U)
        # light.copies.g.soil.tT.boot.median.at.U.target        ...  Q.light.10.boot.median (calculated at the target value of U)
        # light.copies.g.soil.tT.boot.CI.L.at.U.target          ...  Q.light.10.boot.CI.L (calculated at the target value of U)
        # light.copies.g.soil.tT.boot.CI.U.at.U.target          ...  Q.light.10.boot.CI.U (calculated at the target value of U)
        # Tabun.less.equal.0abun.at.U.target.obs                ...  Q.light.10.obs <= Q.light.0.obs (at the target value of U)?
        # Tabun.less.equal.0abun.at.U.target.boot.median        ...  Q.light.10.boot.median <= Q.light.0.boot.median (at the target value of U)?
        # U.min.Tabun.less.equal.0abun.at.U.target.obs          ...  min of range in U for which Q.light.10.obs <= Q.light.0.obs
        # U.max.Tabun.less.equal.0abun.at.U.target.obs          ...  max of range in U for which Q.light.10.obs <= Q.light.0.obs
        # U.min.Tabun.less.equal.0abun.at.U.target.boot.median  ...  min of range in U for which Q.light.10.boot.median <= Q.light.0.boot.median
        # U.max.Tabun.less.equal.0abun.at.U.target.boot.median  ...  max of range in U for which Q.light.10.boot.median <= Q.light.0.boot.median

    #Create graphs of relevant variables across a range of values for U (for all taxa):
      #note: set a priori value of U=0.33 according to McHugh et al. unpublished data from E. coli pure cultures incubated with 18O-H2O
      graphics.off()
      #Create a vector of the taxon-specific sensitivity file names to use in plotting:
      all.sens.files <-list.files("qSIP_output/TM_Sensitivity_Exponential")
      MW.sens.files <- sort(all.sens.files[1:length(levels(data.melted$taxon))])
      sens.files <- sort(all.sens.files[(length(levels(data.melted$taxon))+1):length(all.sens.files)])
      # dev.new(width=8.5, height=11)
      pdf(file="qSIP_output/Figures/TM_U_sensitivity_all_taxa_exponential.pdf", width=8.5, height=11)
      for (i in 1:length(levels(data.melted$taxon))){
        curr.U.sens <- read.table(paste("qSIP_output/TM_Sensitivity_Exponential/", sens.files[i], sep=""), header=TRUE, sep="\t", colClasses="numeric")
        curr.U.sens.MW <- read.table(paste("qSIP_output/TM_Sensitivity_Exponential/", MW.sens.files[i], sep=""), header=TRUE, sep="\t", colClasses="numeric")
        curr.taxon <- gsub(pattern=".+\\_(\\d+)\\.txt", replacement="\\1", x=sens.files[i], perl=TRUE)   #get the ID (name) of the current taxon
        curr.U.sens.list <- list(U.sens=curr.U.sens, U.sens.MW=curr.U.sens.MW)
        U.sens.plot.func(LIST=curr.U.sens.list, taxonID=as.numeric(curr.taxon), U=0.33, growth.model="exponential", selected.Us=seq(0.1,1,length.out=10))
      }
      dev.off()


    #Create a summary graph of the abundance of unlabeled 'light' copies (t=0, t=10) at U=0.33 (all taxa):
      #note: this is based on U=0.33 according to McHugh et al. unpublished data from E. coli pure cultures incubated with 18O-H2O
      graphics.off()
      # dev.new(width=9.7, height=7.5)
      pdf(file="qSIP_output/Figures/TM_Taxa_same_order_PhylumGroups_unlabeled_abundance_exponential.pdf", width=9.7, height=7.5)
      par(mfcol=c(1,2))
      par(mai=c(0.44,0.42,0.2,0.05), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
      phyla <- levels(data.melted$phylum)
      phyla.cols <- cols.for.phyla$col
      ranks <- order(U.sens.all.taxa$tot.copies.g.soil.t0.boot.median)
      tax.order <- data.frame(axis.loc=seq(1, dim(U.sens.all.taxa)[1]), ranks)
      x.mins.N <- c(U.sens.all.taxa$tot.copies.g.soil.t0.obs, U.sens.all.taxa$tot.copies.g.soil.t0.boot.median, U.sens.all.taxa$tot.copies.g.soil.t0.boot.CI.L, U.sens.all.taxa$light.copies.g.soil.tT.obs.at.U.target, U.sens.all.taxa$light.copies.g.soil.tT.boot.median.at.U.target, U.sens.all.taxa$light.copies.g.soil.tT.boot.CI.L.at.U.target)
      x.maxs.N <- c(U.sens.all.taxa$tot.copies.g.soil.t0.obs, U.sens.all.taxa$tot.copies.g.soil.t0.boot.median, U.sens.all.taxa$tot.copies.g.soil.t0.boot.CI.U, U.sens.all.taxa$light.copies.g.soil.tT.obs.at.U.target, U.sens.all.taxa$light.copies.g.soil.tT.boot.median.at.U.target, U.sens.all.taxa$light.copies.g.soil.tT.boot.CI.U.at.U.target)
      x.min.N <- min(x.mins.N[x.mins.N != -Inf], na.rm=TRUE)
      x.max.N <- max(x.maxs.N[x.maxs.N != Inf], na.rm=TRUE)

      N.add.panel <- function(DATA, first){
        plot(y=1:dim(DATA)[1], x=DATA$tot.copies.g.soil.t0.boot.median[tax.order$ranks], type="n", bty="n", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min.N, x.max.N), main="")
        counter <- 1
        for (p in 1:length(phyla)){
          curr.comp.ranked <- DATA[tax.order$ranks,]
          mids.0 <- curr.comp.ranked$tot.copies.g.soil.t0.boot.median[as.character(curr.comp.ranked$taxon) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
          lowers.0 <- curr.comp.ranked$tot.copies.g.soil.t0.boot.CI.L[as.character(curr.comp.ranked$taxon) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
          uppers.0 <- curr.comp.ranked$tot.copies.g.soil.t0.boot.CI.U[as.character(curr.comp.ranked$taxon) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
          tax.nums <- counter:(counter+length(mids.0)-1)
          counter <- counter+length(mids.0)
          points(x=mids.0, y=tax.nums, pch=21, cex=0.6, col=as.character(phyla.cols[p]), bg=as.character(phyla.cols[p]))
          arrows(x0=lowers.0, y0=tax.nums, x1=uppers.0, y1=tax.nums, length=0, angle=90, code=3, col=as.character(phyla.cols[p]))
          if (first){
            par(xpd=NA)
            if(!is.element(p, c(5,6,9,10,16,17,20,21))){
              text(x=x.max.N*0.8, y=mean(tax.nums), labels=phyla[p], adj=c(0,0.5), pos=4, offset=-1.5, cex=0.6, col=as.character(phyla.cols)[p])
            }
            else if (is.element(p, c(5,10))){
              text(x=x.max.N*0.8, y=mean(tax.nums), labels=phyla[p], adj=c(0,0.5), pos=4, offset=0, cex=0.6, col=as.character(phyla.cols)[p])
            }
            else if (is.element(p, c(16,20))){
              text(x=x.max.N*0.8, y=mean(tax.nums), labels=phyla[p], adj=c(0,0.5), pos=4, offset=1.5, cex=0.6, col=as.character(phyla.cols)[p])
            }
            else if (is.element(p, c(6,9,17,21))){
              text(x=x.max.N*0.8, y=mean(tax.nums), labels=phyla[p], adj=c(0,0.5), pos=4, offset=3, cex=0.6, col=as.character(phyla.cols)[p])
            }
            par(xpd=FALSE)
          }
          mids.T <- curr.comp.ranked$light.copies.g.soil.tT.boot.median.at.U.target[as.character(curr.comp.ranked$taxon) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
          lowers.T <- curr.comp.ranked$light.copies.g.soil.tT.boot.CI.L.at.U.target[as.character(curr.comp.ranked$taxon) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
          uppers.T <- curr.comp.ranked$light.copies.g.soil.tT.boot.CI.U.at.U.target[as.character(curr.comp.ranked$taxon) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
          arrows(x0=lowers.T, y0=tax.nums, x1=uppers.T, y1=tax.nums, length=0, angle=90, code=3, col=as.character(phyla.cols[p]), lwd=0.5)
          points(x=mids.T, y=tax.nums, pch=21, cex=0.6, col=as.character(phyla.cols[p]), bg="white")
          mids.0.obs <- curr.comp.ranked$tot.copies.g.soil.t0.obs[as.character(curr.comp.ranked$taxon) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
          points(x=mids.0.obs, y=tax.nums, pch=24, cex=0.3, col="gray30", bg=as.character(phyla.cols[p]))
          mids.T.obs <- curr.comp.ranked$light.copies.g.soil.tT.obs.at.U.target[as.character(curr.comp.ranked$taxon) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
          points(x=mids.T.obs, y=tax.nums, pch=24, cex=0.3, col=as.character(phyla.cols[p]), bg="white")
        }
        abline(v=0, col="black")
        par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
        axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
        mtext(expression(paste("Abundance (16S copies g soil"^-1, ")", sep="")), side=1, line=1.4, cex=0.75)
        if (first){
          axis(side=1, at=c(x.min.N-((x.max.N-x.min.N)*0.04), x.max.N+((x.max.N-x.min.N)*0.04)), labels=FALSE, tck=0, las=1, cex.axis=0.6)
          par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
          axis(side=2, at=c(1-((dim(DATA)[1]-1)*0.04), dim(DATA)[1]+((dim(DATA)[1]-1)*0.04)), labels=FALSE, tck=0, las=1, cex.axis=0.6)
          axis(side=2, at=1:dim(DATA)[1], labels=rep(NA, 103), tck=-0.015, las=1, cex.axis=0.6)
          mtext("Taxon", side=2, line=1.35, cex=0.75)
        }
        else {
          axis(side=1, at=c(x.min.N-((x.max.N-x.min.N)*0.04), x.max.N+((x.max.N-x.min.N)*0.04)), labels=FALSE, tck=0, las=1, cex.axis=0.6)
        }
      }

      N.add.panel(DATA=U.sens.all.taxa, first=TRUE)
      mtext("(full scale)", side=3, line=0, cex=0.85)
      par(xpd=NA)
      legend(x=x.max.N, y=dim(U.sens.all.taxa)[1]-20, legend=c("time 0, boot", "time 10, boot", "time 0, obs", "time 10, obs"), pch=c(21,21,24,24), pt.cex=c(0.6, 0.6, 0.3, 0.3), col=c(as.character(phyla.cols[12]), as.character(phyla.cols[12]), "gray30", as.character(phyla.cols[12])), pt.bg=c(as.character(phyla.cols[12]), "white", as.character(phyla.cols[12]), "white"), bty="o", cex=0.6, xjust=0.5, yjust=1)
      par(xpd=FALSE)
      x.max.N <- 1000000000
      N.add.panel(DATA=U.sens.all.taxa, first=FALSE)
      mtext("(zoomed in)", side=3, line=0, cex=0.85)
      # legend(x=0.5*x.max.N, y=0.6*dim(U.sens.all.taxa)[1], legend=phyla, bty="n", text.col=as.character(phyla.cols), cex=0.6)
      rm(N.add.panel)
      par(mfcol=c(1,1))

      dev.off()  


    #Plot of ranked taxa vs range in U (obs, boot) for which Q.light.10 <= Q.light.0
      graphics.off()
      # dev.new(width=3.75+1.5, height=7.5)
      pdf(file="qSIP_output/Figures/TM_Taxa_same_order_PhylumGroups_U_ranges_exponential.pdf", width=3.75+1.5, height=7.5)
      par(mfcol=c(1,1))
      par(mai=c(0.44,0.42,0.2,1.55), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
      phyla <- levels(data.melted$phylum)
      phyla.cols <- cols.for.phyla$col
      ranks <- order(U.sens.all.taxa$U.max.Tabun.less.equal.0abun.at.U.target.obs)
      tax.order <- data.frame(axis.loc=seq(1, dim(U.sens.all.taxa)[1]), ranks)
      x.mins.U <- c(U.sens.all.taxa$U.min.Tabun.less.equal.0abun.at.U.target.obs, U.sens.all.taxa$U.min.Tabun.less.equal.0abun.at.U.target.boot.median)
      x.maxs.U <- c(U.sens.all.taxa$U.max.Tabun.less.equal.0abun.at.U.target.obs, U.sens.all.taxa$U.max.Tabun.less.equal.0abun.at.U.target.boot.median)
      x.min.U <- min(x.mins.U[x.mins.U != -Inf], na.rm=TRUE)
      x.max.U <- max(x.maxs.U[x.maxs.U != Inf], na.rm=TRUE)

      U.add.panel <- function(DATA){
        plot(y=1:dim(DATA)[1], x=DATA$U.max.Tabun.less.equal.0abun.at.U.target.obs[tax.order$ranks], type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min.U, x.max.U), main="")
        abline(v=0.33, col="black", lwd=0.5)     #note: this is based on U=0.33 according to McHugh et al. unpublished data from E. coli pure cultures incubated with 18O-H2O
        counter <- 1
        for (p in 1:length(phyla)){
          curr.comp.ranked <- DATA[tax.order$ranks,]
          lowers.o <- curr.comp.ranked$U.min.Tabun.less.equal.0abun.at.U.target.obs[as.character(curr.comp.ranked$taxon) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
          uppers.o <- curr.comp.ranked$U.max.Tabun.less.equal.0abun.at.U.target.obs[as.character(curr.comp.ranked$taxon) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
          tax.nums <- counter:(counter+length(lowers.o)-1)
          counter <- counter+length(lowers.o)
          arrows(x0=lowers.o, y0=tax.nums, x1=uppers.o, y1=tax.nums, length=0, angle=90, code=3, col=as.character(phyla.cols[p]), lwd=2.5)
          par(xpd=NA)
          if(!is.element(p, c(5,6,9,10,16,17,20))){
            text(x=x.max.U, y=mean(tax.nums), labels=phyla[p], adj=c(0,0.5), pos=4, offset=0.2, cex=0.6, col=as.character(phyla.cols)[p])
          }
          else if (is.element(p, c(5,10,20))){
            text(x=x.max.U, y=mean(tax.nums), labels=phyla[p], adj=c(0,0.5), pos=4, offset=1.7, cex=0.6, col=as.character(phyla.cols)[p])
          }
          else if (is.element(p, c(16))){
            text(x=x.max.U, y=mean(tax.nums), labels=phyla[p], adj=c(0,0.5), pos=4, offset=3.2, cex=0.6, col=as.character(phyla.cols)[p])
          }
          else if (is.element(p, c(6,9,17))){
            text(x=x.max.U, y=mean(tax.nums), labels=phyla[p], adj=c(0,0.5), pos=4, offset=4.7, cex=0.6, col=as.character(phyla.cols)[p])
          }
          par(xpd=FALSE)
          lowers.b <- curr.comp.ranked$U.min.Tabun.less.equal.0abun.at.U.target.boot.median[as.character(curr.comp.ranked$taxon) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
          uppers.b <- curr.comp.ranked$U.max.Tabun.less.equal.0abun.at.U.target.boot.median[as.character(curr.comp.ranked$taxon) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
          arrows(x0=lowers.b, y0=tax.nums, x1=uppers.b, y1=tax.nums, length=0, angle=90, code=3, col="gray30", lwd=1.0)
          arrows(x0=lowers.b, y0=tax.nums, x1=uppers.b, y1=tax.nums, length=0, angle=90, code=3, col=as.character(phyla.cols[p]), lwd=0.5)
        }
        par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
        axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
        mtext("U", side=1, line=1.4, cex=0.75)
        par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
        axis(side=2, at=1:dim(DATA)[1], labels=rep(NA, 103), tck=-0.015, las=1, cex.axis=0.6)
        mtext("Taxon", side=2, line=1.35, cex=0.75)
      }

      U.add.panel(DATA=U.sens.all.taxa)
      mtext("Range in U for which Q.light.10 <= Q.light.0", side=3, line=0, cex=0.85)
      par(xpd=NA)
      legend(x=1.537, y=dim(U.sens.all.taxa)[1]+((dim(U.sens.all.taxa)[1]-1)*0.07), legend=c("for obs Q.light", "for boot Q.light"), lwd=c(2.5, 1.0), col=c(as.character(phyla.cols[12]), "gray30"), bty="o", cex=0.6, xjust=1, yjust=1)
      legend(x=1.537, y=dim(U.sens.all.taxa)[1]+((dim(U.sens.all.taxa)[1]-1)*0.07), legend=c("for obs Q.light", "for boot Q.light"), lwd=c(2.5, 0.5), col=c(as.character(phyla.cols[12]), as.character(phyla.cols[12])), bty="o", cex=0.6, xjust=1, yjust=1)
      par(xpd=FALSE)
      # legend(x=0.5*x.max.U, y=0.6*dim(U.sens.all.taxa)[1], legend=phyla, bty="n", text.col=as.character(phyla.cols), cex=0.6)
      rm(U.add.panel)
      par(mfcol=c(1,1))

      dev.off()  


    #Rename U.sens.all.taxa to denote that exponential growth model was used:
      U.sens.all.taxa.exp <- U.sens.all.taxa
      rm(U.sens.all.taxa)


    #Write the sensitivity summary results for all taxa (U.sens.all.taxa.exp) to a text file:
      write.table(U.sens.all.taxa.exp, "qSIP_output/TM_U_sens_all_taxa_exponential.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)


    #Constrain U by finding the maximum range of plausbile U values common to all taxa and taking the mean of that range:
      #Range in plausible U values based on observed data & corresponding consensus value of U:
        U.min.obs <- max(U.sens.all.taxa.exp$U.min.Tabun.less.equal.0abun.at.U.target.obs)   #across all taxa, this is the highest value of the lower bound for the range in plausible U values
        U.max.obs <- min(U.sens.all.taxa.exp$U.max.Tabun.less.equal.0abun.at.U.target.obs)   #across all taxa, this is the lowest value of the upper bound for the range in plausible U values
        U.consensus.obs <- mean(c(U.max.obs, U.min.obs))
        U.min.obs
        U.max.obs
        U.consensus.obs
      #Range in plausible U values based on bootstrapped data & corresponding consensus value of U:
        U.min.boot <- max(U.sens.all.taxa.exp$U.min.Tabun.less.equal.0abun.at.U.target.boot.median)   #across all taxa, this is the highest value of the lower bound for the range in plausible U values
        U.max.boot <- min(U.sens.all.taxa.exp$U.max.Tabun.less.equal.0abun.at.U.target.boot.median)   #across all taxa, this is the lowest value of the upper bound for the range in plausible U values
        U.consensus.boot <- mean(c(U.max.boot, U.min.boot))
        U.min.boot
        U.max.boot
        U.consensus.boot
      #Consensus value of U:
        U.consensus.obs.exp <- U.consensus.obs
        rm(U.consensus.obs, U.min.obs, U.max.obs, U.consensus.boot, U.min.boot, U.max.boot)
        U.consensus.obs.exp     #USE THIS ONE (OBSERVED VALUES HAVE LESS POTENTIAL FOR BIAS)


    #Get all growth and abundance estimates for all taxa at the consensus value of U:
      #Create a vector of the taxon-specific sensitivity file names to use in plotting:
        all.sens.files.exp <-list.files("qSIP_output/TM_Sensitivity_Exponential")
        sens.files.exp <- sort(all.sens.files.exp[(length(levels(data.melted$taxon))+1):length(all.sens.files.exp)])
      #Create a data frame containing all growth and abundance estimates for all taxa at the consensus value of U:
        example.sens.frame <- curr.U.sens.list$U.sens     #example sensitivity data.frame (from above) for setting proper names and dimensions
        all.growth.abund.exp <- data.frame(matrix(NA, nrow=length(levels(data.melted$taxon)), ncol=1+dim(example.sens.frame)[2]))
        names(all.growth.abund.exp) <- c("taxon", names(example.sens.frame))
        for (i in 1:length(levels(data.melted$taxon))){
          curr.U.sens <- read.table(paste("qSIP_output/TM_Sensitivity_Exponential/", sens.files.exp[i], sep=""), header=TRUE, sep="\t", colClasses="numeric")
          curr.taxon <- gsub(pattern=".+\\_(\\d+)\\.txt", replacement="\\1", x=sens.files.exp[i], perl=TRUE)   #get the ID (name) of the current taxon
          all.growth.abund.exp$taxon[i] <- as.numeric(curr.taxon)
          all.growth.abund.exp[i,2:dim(all.growth.abund.exp)[2]] <- curr.U.sens[which(abs(curr.U.sens$U-round(U.consensus.obs.exp,2)) == min(abs(curr.U.sens$U-round(U.consensus.obs.exp,2)))),]
        }
        summary(all.growth.abund.exp)
  #}


#More in-depth analysis of the sensitivity to U (all taxa)...USING A LINEAR GROWTH MODEL:______________________________________________________________________
  #{
    #Create dataframes (two per taxon) of relevant variables across a range of values for U:
      #Example for taxon 190:
      set.seed(100)
      system.time(U.sens.T190.lin <- U.sens.func(DATA=data.melted, DATA.POP=ncopies.tube.melted, taxonID=190, time0="Time0", tmt1="16O", tmt2="18O", M.soil=Sdat, iso.compare="18O", days=10, vars=c("taxon", "density.g.ml", "copies", "tube", "tmt", "g.soil"), growth.model="linear", vol=c(100, 50), CI=0.90, draws=1000, sens.seq=seq(0.01, 1, 0.01)))

    #Create graphs of relevant variables across a range of values for U (for a specific taxon):
      #Example for taxon 190:
      #note: set a priori value of U=0.33 according to McHugh et al. unpublished data from E. coli pure cultures incubated with 18O-H2O
      graphics.off()    
      # dev.new(width=8.5, height=11)
      pdf(file="qSIP_output/Figures/TM_U_sensitivity_taxon190_linear.pdf", width=8.5, height=11)
      U.sens.plot.func(LIST=U.sens.T190.lin, taxonID=190, U=0.33, growth.model="linear", selected.Us=seq(0.1,1,length.out=10))

      dev.off()

    #Conduct sensitivity analysis for all taxa; write output dataframes to files (two for each taxon); create a summary table of sensitivity results for all taxa:
      #Create a directory (if it doesn't exist already) for writing the taxon-specific data frames with sensitivity results:
        dir.create(path=paste(getwd(), "/qSIP_output/TM_Sensitivity_Linear", sep=""), showWarnings=FALSE)
      #Set desired target U for analysis:
        #note: this is based on U=0.33 according to McHugh et al. unpublished data from E. coli pure cultures incubated with 18O-H2O
        U.target <- 0.33
      #Create a dataframe containing summary info from sensitivity analysis for all taxa:
        U.sens.all.taxa <- data.frame(taxon=levels(data.melted$taxon))
        U.sens.all.taxa$tot.copies.g.soil.t0.boot.CI.U <- U.sens.all.taxa$tot.copies.g.soil.t0.boot.CI.L <- U.sens.all.taxa$tot.copies.g.soil.t0.boot.median <- U.sens.all.taxa$tot.copies.g.soil.t0.obs <- rep(NA, length(U.sens.all.taxa$taxon))
        U.sens.all.taxa$light.copies.g.soil.tT.boot.CI.U.at.U.target <- U.sens.all.taxa$light.copies.g.soil.tT.boot.CI.L.at.U.target <- U.sens.all.taxa$light.copies.g.soil.tT.boot.median.at.U.target <- U.sens.all.taxa$light.copies.g.soil.tT.obs.at.U.target <- rep(NA, length(U.sens.all.taxa$taxon))
        U.sens.all.taxa$Tabun.less.equal.0abun.at.U.target.boot.median <- U.sens.all.taxa$Tabun.less.equal.0abun.at.U.target.obs <- rep(FALSE, length(U.sens.all.taxa$taxon))
        U.sens.all.taxa$U.max.Tabun.less.equal.0abun.at.U.target.obs <- U.sens.all.taxa$U.min.Tabun.less.equal.0abun.at.U.target.obs <- rep(NA, length(U.sens.all.taxa$taxon))
        U.sens.all.taxa$U.max.Tabun.less.equal.0abun.at.U.target.boot.median <- U.sens.all.taxa$U.min.Tabun.less.equal.0abun.at.U.target.boot.median <- rep(NA, length(U.sens.all.taxa$taxon))
      #Cycle through all taxa to create summary info from sensitivity analysis:
        set.seed(100)
        for(i in 1:length(U.sens.all.taxa$taxon)){
          curr.taxon <- as.character(U.sens.all.taxa$taxon[i])
          curr.U.sens.list <- U.sens.func(DATA=data.melted, DATA.POP=ncopies.tube.melted, taxonID=as.numeric(curr.taxon), time0="Time0", tmt1="16O", tmt2="18O", M.soil=Sdat, iso.compare="18O", days=10, vars=c("taxon", "density.g.ml", "copies", "tube", "tmt", "g.soil"), growth.model="linear", vol=c(100, 50), CI=0.90, draws=1000, sens.seq=seq(0.01, 1, 0.01))
          U.target.index <- which(abs(curr.U.sens.list$U.sens$U-U.target) == min(abs(curr.U.sens.list$U.sens$U-U.target)))
          U.sens.all.taxa$tot.copies.g.soil.t0.obs[i] <- unique(curr.U.sens.list$U.sens$tot.copies.g.soil.t0.obs)
          U.sens.all.taxa$tot.copies.g.soil.t0.boot.median[i] <- unique(curr.U.sens.list$U.sens$tot.copies.g.soil.t0.boot.median)
          U.sens.all.taxa$tot.copies.g.soil.t0.boot.CI.L[i] <- unique(curr.U.sens.list$U.sens$tot.copies.g.soil.t0.boot.CI.L)
          U.sens.all.taxa$tot.copies.g.soil.t0.boot.CI.U[i] <- unique(curr.U.sens.list$U.sens$tot.copies.g.soil.t0.boot.CI.U)
          U.sens.all.taxa$light.copies.g.soil.tT.obs.at.U.target[i] <- curr.U.sens.list$U.sens$light.copies.g.soil.tT.obs[U.target.index]
          U.sens.all.taxa$light.copies.g.soil.tT.boot.median.at.U.target[i] <- curr.U.sens.list$U.sens$light.copies.g.soil.tT.boot.median[U.target.index]
          U.sens.all.taxa$light.copies.g.soil.tT.boot.CI.L.at.U.target[i] <- curr.U.sens.list$U.sens$light.copies.g.soil.tT.boot.CI.L[U.target.index]
          U.sens.all.taxa$light.copies.g.soil.tT.boot.CI.U.at.U.target[i] <- curr.U.sens.list$U.sens$light.copies.g.soil.tT.boot.CI.U[U.target.index]
          U.sens.all.taxa$Tabun.less.equal.0abun.at.U.target.obs[i] <- U.sens.all.taxa$light.copies.g.soil.tT.obs.at.U.target[i] <= U.sens.all.taxa$tot.copies.g.soil.t0.obs[i]
          U.sens.all.taxa$Tabun.less.equal.0abun.at.U.target.boot.median[i] <- U.sens.all.taxa$light.copies.g.soil.tT.boot.median.at.U.target[i] <= U.sens.all.taxa$tot.copies.g.soil.t0.boot.median[i]
          U.range.obs.indices <- which(curr.U.sens.list$U.sens$light.copies.g.soil.tT.obs <= curr.U.sens.list$U.sens$tot.copies.g.soil.t0.obs)
          if(length(U.range.obs.indices) > 0){
            U.sens.all.taxa$U.min.Tabun.less.equal.0abun.at.U.target.obs[i] <- min(curr.U.sens.list$U.sens$U[U.range.obs.indices])
            U.sens.all.taxa$U.max.Tabun.less.equal.0abun.at.U.target.obs[i] <- max(curr.U.sens.list$U.sens$U[U.range.obs.indices])
          }
          U.range.boot.indices <- which(curr.U.sens.list$U.sens$light.copies.g.soil.tT.boot.median <= curr.U.sens.list$U.sens$tot.copies.g.soil.t0.boot.median)
          if(length(U.range.boot.indices) > 0){
            U.sens.all.taxa$U.min.Tabun.less.equal.0abun.at.U.target.boot.median[i] <- min(curr.U.sens.list$U.sens$U[U.range.boot.indices])
            U.sens.all.taxa$U.max.Tabun.less.equal.0abun.at.U.target.boot.median[i] <- max(curr.U.sens.list$U.sens$U[U.range.boot.indices])
          }
          #Write the results (curr.U.sens.list$U.sens & curr.U.sens.list$U.sens.MW) to a text file:
          write.table(curr.U.sens.list$U.sens, paste("qSIP_output/TM_Sensitivity_Linear/TM_U_sens_taxon_", sprintf("%03d", as.numeric(curr.taxon)), ".txt", sep=""), append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
          write.table(curr.U.sens.list$U.sens.MW, paste("qSIP_output/TM_Sensitivity_Linear/TM_U_sens_MW_taxon_", sprintf("%03d", as.numeric(curr.taxon)), ".txt", sep=""), append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
        }
      #Description of the variables in the summary table (U.sens.all.taxa):
        # taxon                                                 ...  taxon
        # tot.copies.g.soil.t0.obs                              ...  Q.light.0.obs
        # tot.copies.g.soil.t0.boot.median                      ...  Q.light.0.boot.median
        # tot.copies.g.soil.t0.boot.CI.L                        ...  Q.light.0.boot.CI.L
        # tot.copies.g.soil.t0.boot.CI.U                        ...  Q.light.0.boot.CI.U
        # light.copies.g.soil.tT.obs.at.U.target                ...  Q.light.10.obs (calculated at the target value of U)
        # light.copies.g.soil.tT.boot.median.at.U.target        ...  Q.light.10.boot.median (calculated at the target value of U)
        # light.copies.g.soil.tT.boot.CI.L.at.U.target          ...  Q.light.10.boot.CI.L (calculated at the target value of U)
        # light.copies.g.soil.tT.boot.CI.U.at.U.target          ...  Q.light.10.boot.CI.U (calculated at the target value of U)
        # Tabun.less.equal.0abun.at.U.target.obs                ...  Q.light.10.obs <= Q.light.0.obs (at the target value of U)?
        # Tabun.less.equal.0abun.at.U.target.boot.median        ...  Q.light.10.boot.median <= Q.light.0.boot.median (at the target value of U)?
        # U.min.Tabun.less.equal.0abun.at.U.target.obs          ...  min of range in U for which Q.light.10.obs <= Q.light.0.obs
        # U.max.Tabun.less.equal.0abun.at.U.target.obs          ...  max of range in U for which Q.light.10.obs <= Q.light.0.obs
        # U.min.Tabun.less.equal.0abun.at.U.target.boot.median  ...  min of range in U for which Q.light.10.boot.median <= Q.light.0.boot.median
        # U.max.Tabun.less.equal.0abun.at.U.target.boot.median  ...  max of range in U for which Q.light.10.boot.median <= Q.light.0.boot.median

    #Create graphs of relevant variables across a range of values for U (for all taxa):
      #note: set a priori value of U=0.33 according to McHugh et al. unpublished data from E. coli pure cultures incubated with 18O-H2O
      graphics.off()
      #Create a vector of the taxon-specific sensitivity file names to use in plotting:
      all.sens.files <-list.files("qSIP_output/TM_Sensitivity_Linear")
      MW.sens.files <- sort(all.sens.files[1:length(levels(data.melted$taxon))])
      sens.files <- sort(all.sens.files[(length(levels(data.melted$taxon))+1):length(all.sens.files)])
      # dev.new(width=8.5, height=11)
      pdf(file="qSIP_output/Figures/TM_U_sensitivity_all_taxa_linear.pdf", width=8.5, height=11)
      for (i in 1:length(levels(data.melted$taxon))){
        curr.U.sens <- read.table(paste("qSIP_output/TM_Sensitivity_Linear/", sens.files[i], sep=""), header=TRUE, sep="\t", colClasses="numeric")
        curr.U.sens.MW <- read.table(paste("qSIP_output/TM_Sensitivity_Linear/", MW.sens.files[i], sep=""), header=TRUE, sep="\t", colClasses="numeric")
        curr.taxon <- gsub(pattern=".+\\_(\\d+)\\.txt", replacement="\\1", x=sens.files[i], perl=TRUE)   #get the ID (name) of the current taxon
        curr.U.sens.list <- list(U.sens=curr.U.sens, U.sens.MW=curr.U.sens.MW)
        U.sens.plot.func(LIST=curr.U.sens.list, taxonID=as.numeric(curr.taxon), U=0.33, growth.model="linear", selected.Us=seq(0.1,1,length.out=10))
      }
      dev.off()


    #Create a summary graph of the abundance of unlabeled 'light' copies (t=0, t=10) at U=0.33 (all taxa):
      #note: this is based on U=0.33 according to McHugh et al. unpublished data from E. coli pure cultures incubated with 18O-H2O
      graphics.off()
      # dev.new(width=9.7, height=7.5)
      pdf(file="qSIP_output/Figures/TM_Taxa_same_order_PhylumGroups_unlabeled_abundance_linear.pdf", width=9.7, height=7.5)
      par(mfcol=c(1,2))
      par(mai=c(0.44,0.42,0.2,0.05), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
      phyla <- levels(data.melted$phylum)
      phyla.cols <- cols.for.phyla$col
      ranks <- order(U.sens.all.taxa$tot.copies.g.soil.t0.boot.median)
      tax.order <- data.frame(axis.loc=seq(1, dim(U.sens.all.taxa)[1]), ranks)
      x.mins.N <- c(U.sens.all.taxa$tot.copies.g.soil.t0.obs, U.sens.all.taxa$tot.copies.g.soil.t0.boot.median, U.sens.all.taxa$tot.copies.g.soil.t0.boot.CI.L, U.sens.all.taxa$light.copies.g.soil.tT.obs.at.U.target, U.sens.all.taxa$light.copies.g.soil.tT.boot.median.at.U.target, U.sens.all.taxa$light.copies.g.soil.tT.boot.CI.L.at.U.target)
      x.maxs.N <- c(U.sens.all.taxa$tot.copies.g.soil.t0.obs, U.sens.all.taxa$tot.copies.g.soil.t0.boot.median, U.sens.all.taxa$tot.copies.g.soil.t0.boot.CI.U, U.sens.all.taxa$light.copies.g.soil.tT.obs.at.U.target, U.sens.all.taxa$light.copies.g.soil.tT.boot.median.at.U.target, U.sens.all.taxa$light.copies.g.soil.tT.boot.CI.U.at.U.target)
      x.min.N <- min(x.mins.N[x.mins.N != -Inf], na.rm=TRUE)
      x.max.N <- max(x.maxs.N[x.maxs.N != Inf], na.rm=TRUE)

      N.add.panel <- function(DATA, first){
        plot(y=1:dim(DATA)[1], x=DATA$tot.copies.g.soil.t0.boot.median[tax.order$ranks], type="n", bty="n", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min.N, x.max.N), main="")
        counter <- 1
        for (p in 1:length(phyla)){
          curr.comp.ranked <- DATA[tax.order$ranks,]
          mids.0 <- curr.comp.ranked$tot.copies.g.soil.t0.boot.median[as.character(curr.comp.ranked$taxon) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
          lowers.0 <- curr.comp.ranked$tot.copies.g.soil.t0.boot.CI.L[as.character(curr.comp.ranked$taxon) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
          uppers.0 <- curr.comp.ranked$tot.copies.g.soil.t0.boot.CI.U[as.character(curr.comp.ranked$taxon) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
          tax.nums <- counter:(counter+length(mids.0)-1)
          counter <- counter+length(mids.0)
          points(x=mids.0, y=tax.nums, pch=21, cex=0.6, col=as.character(phyla.cols[p]), bg=as.character(phyla.cols[p]))
          arrows(x0=lowers.0, y0=tax.nums, x1=uppers.0, y1=tax.nums, length=0, angle=90, code=3, col=as.character(phyla.cols[p]))
          if (first){
            par(xpd=NA)
            if(!is.element(p, c(5,6,9,10,16,17,20,21))){
              text(x=x.max.N*0.8, y=mean(tax.nums), labels=phyla[p], adj=c(0,0.5), pos=4, offset=-1.5, cex=0.6, col=as.character(phyla.cols)[p])
            }
            else if (is.element(p, c(5,10))){
              text(x=x.max.N*0.8, y=mean(tax.nums), labels=phyla[p], adj=c(0,0.5), pos=4, offset=0, cex=0.6, col=as.character(phyla.cols)[p])
            }
            else if (is.element(p, c(16,20))){
              text(x=x.max.N*0.8, y=mean(tax.nums), labels=phyla[p], adj=c(0,0.5), pos=4, offset=1.5, cex=0.6, col=as.character(phyla.cols)[p])
            }
            else if (is.element(p, c(6,9,17,21))){
              text(x=x.max.N*0.8, y=mean(tax.nums), labels=phyla[p], adj=c(0,0.5), pos=4, offset=3, cex=0.6, col=as.character(phyla.cols)[p])
            }
            par(xpd=FALSE)
          }
          mids.T <- curr.comp.ranked$light.copies.g.soil.tT.boot.median.at.U.target[as.character(curr.comp.ranked$taxon) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
          lowers.T <- curr.comp.ranked$light.copies.g.soil.tT.boot.CI.L.at.U.target[as.character(curr.comp.ranked$taxon) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
          uppers.T <- curr.comp.ranked$light.copies.g.soil.tT.boot.CI.U.at.U.target[as.character(curr.comp.ranked$taxon) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
          arrows(x0=lowers.T, y0=tax.nums, x1=uppers.T, y1=tax.nums, length=0, angle=90, code=3, col=as.character(phyla.cols[p]), lwd=0.5)
          points(x=mids.T, y=tax.nums, pch=21, cex=0.6, col=as.character(phyla.cols[p]), bg="white")
          mids.0.obs <- curr.comp.ranked$tot.copies.g.soil.t0.obs[as.character(curr.comp.ranked$taxon) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
          points(x=mids.0.obs, y=tax.nums, pch=24, cex=0.3, col="gray30", bg=as.character(phyla.cols[p]))
          mids.T.obs <- curr.comp.ranked$light.copies.g.soil.tT.obs.at.U.target[as.character(curr.comp.ranked$taxon) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
          points(x=mids.T.obs, y=tax.nums, pch=24, cex=0.3, col=as.character(phyla.cols[p]), bg="white")
        }
        abline(v=0, col="black")
        par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
        axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
        mtext(expression(paste("Abundance (16S copies g soil"^-1, ")", sep="")), side=1, line=1.4, cex=0.75)
        if (first){
          axis(side=1, at=c(x.min.N-((x.max.N-x.min.N)*0.04), x.max.N+((x.max.N-x.min.N)*0.04)), labels=FALSE, tck=0, las=1, cex.axis=0.6)
          par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
          axis(side=2, at=c(1-((dim(DATA)[1]-1)*0.04), dim(DATA)[1]+((dim(DATA)[1]-1)*0.04)), labels=FALSE, tck=0, las=1, cex.axis=0.6)
          axis(side=2, at=1:dim(DATA)[1], labels=rep(NA, 103), tck=-0.015, las=1, cex.axis=0.6)
          mtext("Taxon", side=2, line=1.35, cex=0.75)
        }
        else {
          axis(side=1, at=c(x.min.N-((x.max.N-x.min.N)*0.04), x.max.N+((x.max.N-x.min.N)*0.04)), labels=FALSE, tck=0, las=1, cex.axis=0.6)
        }
      }

      N.add.panel(DATA=U.sens.all.taxa, first=TRUE)
      mtext("(full scale)", side=3, line=0, cex=0.85)
      par(xpd=NA)
      legend(x=x.max.N, y=dim(U.sens.all.taxa)[1]-20, legend=c("time 0, boot", "time 10, boot", "time 0, obs", "time 10, obs"), pch=c(21,21,24,24), pt.cex=c(0.6, 0.6, 0.3, 0.3), col=c(as.character(phyla.cols[12]), as.character(phyla.cols[12]), "gray30", as.character(phyla.cols[12])), pt.bg=c(as.character(phyla.cols[12]), "white", as.character(phyla.cols[12]), "white"), bty="o", cex=0.6, xjust=0.5, yjust=1)
      par(xpd=FALSE)
      x.max.N <- 1000000000
      N.add.panel(DATA=U.sens.all.taxa, first=FALSE)
      mtext("(zoomed in)", side=3, line=0, cex=0.85)
      # legend(x=0.5*x.max.N, y=0.6*dim(U.sens.all.taxa)[1], legend=phyla, bty="n", text.col=as.character(phyla.cols), cex=0.6)
      rm(N.add.panel)
      par(mfcol=c(1,1))

      dev.off()  


    #Plot of ranked taxa vs range in U (obs, boot) for which Q.light.10 <= Q.light.0
      graphics.off()
      # dev.new(width=3.75+1.5, height=7.5)
      pdf(file="qSIP_output/Figures/TM_Taxa_same_order_PhylumGroups_U_ranges_linear.pdf", width=3.75+1.5, height=7.5)
      par(mfcol=c(1,1))
      par(mai=c(0.44,0.42,0.2,1.55), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
      phyla <- levels(data.melted$phylum)
      phyla.cols <- cols.for.phyla$col
      ranks <- order(U.sens.all.taxa$U.max.Tabun.less.equal.0abun.at.U.target.obs)
      tax.order <- data.frame(axis.loc=seq(1, dim(U.sens.all.taxa)[1]), ranks)
      x.mins.U <- c(U.sens.all.taxa$U.min.Tabun.less.equal.0abun.at.U.target.obs, U.sens.all.taxa$U.min.Tabun.less.equal.0abun.at.U.target.boot.median)
      x.maxs.U <- c(U.sens.all.taxa$U.max.Tabun.less.equal.0abun.at.U.target.obs, U.sens.all.taxa$U.max.Tabun.less.equal.0abun.at.U.target.boot.median)
      x.min.U <- min(x.mins.U[x.mins.U != -Inf], na.rm=TRUE)
      x.max.U <- max(x.maxs.U[x.maxs.U != Inf], na.rm=TRUE)

      U.add.panel <- function(DATA){
        plot(y=1:dim(DATA)[1], x=DATA$U.max.Tabun.less.equal.0abun.at.U.target.obs[tax.order$ranks], type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min.U, x.max.U), main="")
        abline(v=0.33, col="black", lwd=0.5)     #note: this is based on U=0.33 according to McHugh et al. unpublished data from E. coli pure cultures incubated with 18O-H2O
        counter <- 1
        for (p in 1:length(phyla)){
          curr.comp.ranked <- DATA[tax.order$ranks,]
          lowers.o <- curr.comp.ranked$U.min.Tabun.less.equal.0abun.at.U.target.obs[as.character(curr.comp.ranked$taxon) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
          uppers.o <- curr.comp.ranked$U.max.Tabun.less.equal.0abun.at.U.target.obs[as.character(curr.comp.ranked$taxon) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
          tax.nums <- counter:(counter+length(lowers.o)-1)
          counter <- counter+length(lowers.o)
          arrows(x0=lowers.o, y0=tax.nums, x1=uppers.o, y1=tax.nums, length=0, angle=90, code=3, col=as.character(phyla.cols[p]), lwd=2.5)
          par(xpd=NA)
          if(!is.element(p, c(5,6,9,10,16,17,20))){
            text(x=x.max.U, y=mean(tax.nums), labels=phyla[p], adj=c(0,0.5), pos=4, offset=0.2, cex=0.6, col=as.character(phyla.cols)[p])
          }
          else if (is.element(p, c(5,10,20))){
            text(x=x.max.U, y=mean(tax.nums), labels=phyla[p], adj=c(0,0.5), pos=4, offset=1.7, cex=0.6, col=as.character(phyla.cols)[p])
          }
          else if (is.element(p, c(16))){
            text(x=x.max.U, y=mean(tax.nums), labels=phyla[p], adj=c(0,0.5), pos=4, offset=3.2, cex=0.6, col=as.character(phyla.cols)[p])
          }
          else if (is.element(p, c(6,9,17))){
            text(x=x.max.U, y=mean(tax.nums), labels=phyla[p], adj=c(0,0.5), pos=4, offset=4.7, cex=0.6, col=as.character(phyla.cols)[p])
          }
          par(xpd=FALSE)
          lowers.b <- curr.comp.ranked$U.min.Tabun.less.equal.0abun.at.U.target.boot.median[as.character(curr.comp.ranked$taxon) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
          uppers.b <- curr.comp.ranked$U.max.Tabun.less.equal.0abun.at.U.target.boot.median[as.character(curr.comp.ranked$taxon) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
          arrows(x0=lowers.b, y0=tax.nums, x1=uppers.b, y1=tax.nums, length=0, angle=90, code=3, col="gray30", lwd=1.0)
          arrows(x0=lowers.b, y0=tax.nums, x1=uppers.b, y1=tax.nums, length=0, angle=90, code=3, col=as.character(phyla.cols[p]), lwd=0.5)
        }
        par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
        axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
        mtext("U", side=1, line=1.4, cex=0.75)
        par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
        axis(side=2, at=1:dim(DATA)[1], labels=rep(NA, 103), tck=-0.015, las=1, cex.axis=0.6)
        mtext("Taxon", side=2, line=1.35, cex=0.75)
      }

      U.add.panel(DATA=U.sens.all.taxa)
      mtext("Range in U for which Q.light.10 <= Q.light.0", side=3, line=0, cex=0.85)
      par(xpd=NA)
      legend(x=1.537, y=dim(U.sens.all.taxa)[1]+((dim(U.sens.all.taxa)[1]-1)*0.07), legend=c("for obs Q.light", "for boot Q.light"), lwd=c(2.5, 1.0), col=c(as.character(phyla.cols[12]), "gray30"), bty="o", cex=0.6, xjust=1, yjust=1)
      legend(x=1.537, y=dim(U.sens.all.taxa)[1]+((dim(U.sens.all.taxa)[1]-1)*0.07), legend=c("for obs Q.light", "for boot Q.light"), lwd=c(2.5, 0.5), col=c(as.character(phyla.cols[12]), as.character(phyla.cols[12])), bty="o", cex=0.6, xjust=1, yjust=1)
      par(xpd=FALSE)
      # legend(x=0.5*x.max.U, y=0.6*dim(U.sens.all.taxa)[1], legend=phyla, bty="n", text.col=as.character(phyla.cols), cex=0.6)
      rm(U.add.panel)
      par(mfcol=c(1,1))

      dev.off()  


    #Rename U.sens.all.taxa to denote that linear growth model was used:
      U.sens.all.taxa.lin <- U.sens.all.taxa
      rm(U.sens.all.taxa)


    #Write the sensitivity summary results for all taxa (U.sens.all.taxa.lin) to a text file:
      write.table(U.sens.all.taxa.lin, "qSIP_output/TM_U_sens_all_taxa_linear.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)


    #Constrain U by finding the maximum range of plausbile U values common to all taxa and taking the mean of that range:
      #Range in plausible U values based on observed data & corresponding consensus value of U:
        U.min.obs <- max(U.sens.all.taxa.lin$U.min.Tabun.less.equal.0abun.at.U.target.obs)   #across all taxa, this is the highest value of the lower bound for the range in plausible U values
        U.max.obs <- min(U.sens.all.taxa.lin$U.max.Tabun.less.equal.0abun.at.U.target.obs)   #across all taxa, this is the lowest value of the upper bound for the range in plausible U values
        U.consensus.obs <- mean(c(U.max.obs, U.min.obs))
        U.min.obs
        U.max.obs
        U.consensus.obs
      #Range in plausible U values based on bootstrapped data & corresponding consensus value of U:
        U.min.boot <- max(U.sens.all.taxa.lin$U.min.Tabun.less.equal.0abun.at.U.target.boot.median)   #across all taxa, this is the highest value of the lower bound for the range in plausible U values
        U.max.boot <- min(U.sens.all.taxa.lin$U.max.Tabun.less.equal.0abun.at.U.target.boot.median)   #across all taxa, this is the lowest value of the upper bound for the range in plausible U values
        U.consensus.boot <- mean(c(U.max.boot, U.min.boot))
        U.min.boot
        U.max.boot
        U.consensus.boot
      #Consensus value of U:
        U.consensus.obs.lin <- U.consensus.obs
        rm(U.consensus.obs, U.min.obs, U.max.obs, U.consensus.boot, U.min.boot, U.max.boot)
        U.consensus.obs.lin     #USE THIS ONE (OBSERVED VALUES HAVE LESS POTENTIAL FOR BIAS)


    #Get all growth and abundance estimates for all taxa at the consensus value of U:
      #Create a vector of the taxon-specific sensitivity file names to use in plotting:
        all.sens.files.lin <-list.files("qSIP_output/TM_Sensitivity_Linear")
        sens.files.lin <- sort(all.sens.files.lin[(length(levels(data.melted$taxon))+1):length(all.sens.files.lin)])
      #Create a data frame containing all growth and abundance estimates for all taxa at the consensus value of U:
        example.sens.frame <- curr.U.sens.list$U.sens     #example sensitivity data.frame (from above) for setting proper names and dimensions
        all.growth.abund.lin <- data.frame(matrix(NA, nrow=length(levels(data.melted$taxon)), ncol=1+dim(example.sens.frame)[2]))
        names(all.growth.abund.lin) <- c("taxon", names(example.sens.frame))
        for (i in 1:length(levels(data.melted$taxon))){
          curr.U.sens <- read.table(paste("qSIP_output/TM_Sensitivity_Linear/", sens.files.lin[i], sep=""), header=TRUE, sep="\t", colClasses="numeric")
          curr.taxon <- gsub(pattern=".+\\_(\\d+)\\.txt", replacement="\\1", x=sens.files.lin[i], perl=TRUE)   #get the ID (name) of the current taxon
          all.growth.abund.lin$taxon[i] <- as.numeric(curr.taxon)
          all.growth.abund.lin[i,2:dim(all.growth.abund.lin)[2]] <- curr.U.sens[which(abs(curr.U.sens$U-round(U.consensus.obs.lin,2)) == min(abs(curr.U.sens$U-round(U.consensus.obs.lin,2)))),]
        }
        summary(all.growth.abund.lin)
  #}


#Estimate assemblage-level turnover rates and compare to literature values using radiolabeled nucleotides/amino acids:_________________________________________
  #{
    #EXPONENTIAL GROWTH MODEL:
      #Calculate assemblage-level turnover using death rates (convert death rates to positive values) and initial abundances:
      #note: this is an approximate estimate of assemblage-level turnover; strictly speaking, the calculation should be done with complete bootstrap vectors (this is done later in the code)
        assemblage.turnover.d.exp <- sum(-all.growth.abund.exp$d.boot.median * all.growth.abund.exp$tot.copies.g.soil.t0.boot.median) / sum(all.growth.abund.exp$tot.copies.g.soil.t0.boot.median)
        assemblage.turnover.d.exp
      #Calculate assemblage-level turnover using birth rates and initial abundances (may be better to use time t abundances; but death rate is more appropriate anyway):
      #note: this is an approximate estimate of assemblage-level turnover; strictly speaking, the calculation should be done with complete bootstrap vectors (this is done later in the code)
        assemblage.turnover.b.exp <- sum(all.growth.abund.exp$r.gross.boot.median * all.growth.abund.exp$tot.copies.g.soil.t0.boot.median) / sum(all.growth.abund.exp$tot.copies.g.soil.t0.boot.median)
        assemblage.turnover.b.exp

    #LINEAR GROWTH MODEL:
      #Calculate assemblage-level turnover using death rates (convert death rates to positive values) and initial abundances:
      #note: this is an approximate estimate of assemblage-level turnover; strictly speaking, the calculation should be done with complete bootstrap vectors (this is done later in the code)
        assemblage.turnover.d.lin <- sum(-all.growth.abund.lin$d.boot.median) / sum(all.growth.abund.lin$tot.copies.g.soil.t0.boot.median)
        assemblage.turnover.d.lin
      #Calculate assemblage-level turnover using birth rates and initial abundances:
      #note: this is an approximate estimate of assemblage-level turnover; strictly speaking, the calculation should be done with complete bootstrap vectors (this is done later in the code)
        assemblage.turnover.b.lin <- sum(all.growth.abund.lin$r.gross.boot.median) / sum(all.growth.abund.lin$tot.copies.g.soil.t0.boot.median)
        assemblage.turnover.b.lin

    #Import quick & dirty literature values for bacterial assemblage-level turnover from Rousk & Baath 2011 (Table 1):
      Lit.turnover <- read.table("qSIP_data/Lit_turnover_for_R.txt", header=TRUE, sep="\t")
      summary(Lit.turnover)
      
    #Graph quick & dirty literature turnover values to compare with those from these data (rough graph):
      graphics.off()
      x.values.lit <- (1:dim(Lit.turnover)[1])-1
      plot(y=Lit.turnover$mean.turnover.rate.days, x=x.values.lit, bty="l", type="p", pch=21, col="black", bg="black", xlim=c(0,100), ylim=c(0,max(Lit.turnover[,4:6], na.rm=TRUE)), xaxt="n", xlab="", ylab="Assemblage turnover rate (day-1)")
      arrows(x0=x.values.lit, y0=Lit.turnover$min.turnover.rate.days, x1=x.values.lit, y1=Lit.turnover$max.turnover.rate.days, length=0, angle=90, code=3)
      points(x=15, y=assemblage.turnover.d.exp, pch=21, col="black", bg="red")    #death rate is more reflective of true turnover than is birth rate
      points(x=15, y=assemblage.turnover.b.exp, pch=21, col="black", bg="green")
      points(x=20, y=assemblage.turnover.d.lin, pch=22, col="black", bg="red")    #death rate is more reflective of true turnover than is birth rate
      points(x=20, y=assemblage.turnover.b.lin, pch=22, col="black", bg="green")
      legend(x=100, y=1.04*max(Lit.turnover[,4:6], na.rm=TRUE), legend=c("exponential growth model (d)", "linear growth model (d)", "exponential growth model (b)", "linear growth model (b)"), pch=c(21,22,21,22), col="black", pt.bg=c("red", "red", "green", "green"), bty="o", cex=1.0, xjust=1, yjust=1)

      dev.off()
  #}




#Save workspace:
  save(list=ls(all.names=TRUE), file="qSIP_workspaces/TM_01-02/.RData", envir=.GlobalEnv)



