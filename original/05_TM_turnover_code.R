# This code performs the analysis of the literature soil microbial turnover values (including those from the TM qSIP data)


graphics.off()	#close all graphics windows


#Set working directory:
  #ALREADY DONE FOR THIS WORKSPACE IN PREVIOUSLY RUN CODE; only reset it here if loading previously saved workspace (see below) or starting clean (then do not load previously saved workspace and be sure to un-comment the data importing steps below):


#Reload the saved workspace resulting from the previous script:
  setwd("/Users/bk/Research/Projects/SIP_Modeling/qSIP")
  load("qSIP_workspaces/TM_01-02-03-04/.RData")


#Load libraries & scripts:
  library(car)
  library(smatr)

#Import the complete set of literature values for bacterial assemblage-level turnover:
  Lit.turnover3 <- read.table("qSIP_data/Lit_turnover_for_R_complete.txt", header=TRUE, sep="\t")
  summary(Lit.turnover3)


#These data importing and formatting/calculation steps were already done in the TM analysis of the imported workspace._________________________________________
#They do not need to be re-run here unless I am starting this analysis from scratch (i.e., from a brand new workspace).  Thus they are commented out:
  #{
    # #Import quick & dirty literature values for bacterial assemblage-level turnover from Rousk & Baath 2011 (Table 1):
    #   Lit.turnover <- read.table("qSIP_data/Lit_turnover_for_R.txt", header=TRUE, sep="\t")
    #   summary(Lit.turnover)
    #
    # #Import a more comprehensive set of selected literature values for bacterial assemblage-level turnover:
    #   Lit.turnover2 <- read.table("qSIP_data/Lit_turnover_for_R_selected.txt", header=TRUE, sep="\t")
    #   summary(Lit.turnover2)
    #
    # #Read in the necessary bootstrapped growth, C flux, and abundance estimates from the all.taxa.calcs & all.taxa.calcs.pop output of the TM analysis:
    #   all.rates <- read.table("qSIP_output/TM_all_rates.txt", header=TRUE, sep="")
    #   r.gross.boots <- read.table("qSIP_output/TM_bootstrapped_r.txt", header=TRUE, sep="")
    #   r.net.boots <- read.table("qSIP_output/TM_bootstrapped_r_pop.txt", header=TRUE, sep="")
    #   N.Tt.boots <- read.table("qSIP_output/TM_bootstrapped_N_Tt_pop.txt", header=TRUE, sep="")
    #   N.T0.boots <- read.table("qSIP_output/TM_bootstrapped_N_T0_pop.txt", header=TRUE, sep="")
    #
    # #Get indices for the comparisons of interest (there is only one comparison for this dataset)
    #   inds1 <- r.gross.boots$comparisonID == 1            #growth with added water
    #
    # #Calculate death rate as d = r - b (convention is that death rates are negative)...bootstrapped:
    #   r.net.boots.only <- r.net.boots[inds1, 5:dim(r.net.boots)[2]]
    #   b.boots.only <- r.gross.boots[inds1, 5:dim(r.gross.boots)[2]]
    #   d.boots.only <- r.net.boots.only - b.boots.only
    #
    # #Abundance of copies at times 0 & t:
    #   N.Tt.boots.only <- N.Tt.boots[inds1, 5:dim(N.Tt.boots)[2]]
    #   N.T0.boots.only <- N.T0.boots[inds1, 5:dim(N.T0.boots)[2]]
    #
    # #Calculate assemblage-level turnover using death rates and initial abundances (convert death rates to positive values):
    #   assemblage.turnover.d.obs <- sum(-all.rates$d.obs * all.rates$N.tot.T0.obs) / sum(all.rates$N.tot.T0.obs)
    #   assemblage.turnover.d.boots <- apply(-d.boots.only * N.T0.boots.only, 2, sum) / apply(N.T0.boots.only, 2, sum)
    #   assemblage.turnover.d.median <- median(assemblage.turnover.d.boots, na.rm=TRUE)
    #   assemblage.turnover.d.CI.L <- quantile(assemblage.turnover.d.boots, probs=0.05, na.rm=TRUE)
    #   assemblage.turnover.d.CI.U <- quantile(assemblage.turnover.d.boots, probs=0.95, na.rm=TRUE)
    # #Calculate assemblage-level turnover using birth rates and initial abundances (may be better to use time t abundances; but death rate is more appropriate anyway):
    #   assemblage.turnover.b.obs <- sum(all.rates$b.obs * all.rates$N.tot.T0.obs) / sum(all.rates$N.tot.T0.obs)
    #   assemblage.turnover.b.boots <- apply(b.boots.only * N.T0.boots.only, 2, sum) / apply(N.T0.boots.only, 2, sum)
    #   assemblage.turnover.b.median <- median(assemblage.turnover.b.boots, na.rm=TRUE)
    #   assemblage.turnover.b.CI.L <- quantile(assemblage.turnover.b.boots, probs=0.05, na.rm=TRUE)
    #   assemblage.turnover.b.CI.U <- quantile(assemblage.turnover.b.boots, probs=0.95, na.rm=TRUE)
    #
    # #Turnover values to report in the text:
    #   #Theresa's soil rewetting data (qSIP):
    #     assemblage.turnover.d.median
    #     assemblage.turnover.d.CI.L
    #     assemblage.turnover.d.CI.U
  #}


#Estimate assemblage-level turnover rates and compare to the complete collection of literature values (all values) using various methods:______________________
  #{
    #Use the complete set of literature values for bacterial assemblage-level turnover:
      # Lit.turnover3
    #Calculate median and 90%CIs of turnover rate for all literature values:
      Lit.turnover3$median.turnover.rate.d <- NA
      Lit.turnover3$CI.L.turnover.rate.d <- NA
      Lit.turnover3$CI.U.turnover.rate.d <- NA
      set.seed(100)
      for (i in 1:dim(Lit.turnover3)[1]){
        if (!is.na(Lit.turnover3$mean.turnover.time.d[i]) & Lit.turnover3$author[i] != "Baath" & Lit.turnover3$year[i] != 1998){   #Baath 1998 has turnover times and rates already provided
          temp.times.d <- rnorm(1000, mean=Lit.turnover3$mean.turnover.time.d[i], sd=Lit.turnover3$sd.turnover.time.d[i])
          temp.rates.d <- 1/temp.times.d
          Lit.turnover3$mean.turnover.rate.d[i] <- mean(temp.rates.d, na.rm=TRUE)
          Lit.turnover3$sd.turnover.rate.d[i] <- sd(temp.rates.d, na.rm=TRUE)
          Lit.turnover3$median.turnover.rate.d[i] <- median(temp.rates.d, na.rm=TRUE)
          Lit.turnover3$CI.L.turnover.rate.d[i] <- quantile(temp.rates.d, probs=0.05, na.rm=TRUE)
          Lit.turnover3$CI.U.turnover.rate.d[i] <- quantile(temp.rates.d, probs=0.95, na.rm=TRUE)
        } else if (!is.na(Lit.turnover3$mean.turnover.rate.d[i])){
          if (Lit.turnover3$author[i] == "Hunt et al." | (Lit.turnover3$author[i] == "Baath" & Lit.turnover3$year[i] == 1998)){   #Hunt et al. & Baath 1998 only have one estimate -- set this as the median
            temp.rates.d <- rnorm(1000, mean=Lit.turnover3$mean.turnover.rate.d[i], sd=Lit.turnover3$sd.turnover.rate.d[i])
            Lit.turnover3$median.turnover.rate.d[i] <- Lit.turnover3$mean.turnover.rate.d[i]
          } else if (Lit.turnover3$author[i] != "Hunt et al."){
            temp.rates.d <- rnorm(1000, mean=Lit.turnover3$mean.turnover.rate.d[i], sd=Lit.turnover3$sd.turnover.rate.d[i])
            Lit.turnover3$median.turnover.rate.d[i] <- median(temp.rates.d, na.rm=TRUE)
            Lit.turnover3$CI.L.turnover.rate.d[i] <- quantile(temp.rates.d, probs=0.05, na.rm=TRUE)
            Lit.turnover3$CI.U.turnover.rate.d[i] <- quantile(temp.rates.d, probs=0.95, na.rm=TRUE)
          }
        }
      }

    #Graph the complete set of literature bacterial assemblage turnover rates in comparison to those estimated in TM analysis:
        #Turnover rates - in units of day-1:
          graphics.off()
          # dev.new(width=3.0, height=3.0)
          pdf(file="qSIP_output/Figures/TM_assemblage_turnover_comparison_complete.pdf", width=3.0, height=3.0)
          par(mai=c(0.29,0.47,0.02,0.02), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
          y.mins <- c(0, as.numeric(unlist(Lit.turnover3[,26:28])), assemblage.turnover.d.CI.L)
          y.maxs <- c(0, as.numeric(unlist(Lit.turnover3[,26:28])), assemblage.turnover.d.CI.U)
          y.min <- min(y.mins[y.mins != -Inf], na.rm=TRUE)
          y.max <- max(y.maxs[y.maxs != Inf], na.rm=TRUE)
          x.values.lit <- (1:dim(Lit.turnover3)[1])-1
          y.medians.lit <- Lit.turnover3$median.turnover.rate.d[order(Lit.turnover3$median.turnover.rate.d, decreasing=FALSE)]
          y.CI.L.lit <- Lit.turnover3$CI.L.turnover.rate.d[order(Lit.turnover3$median.turnover.rate.d, decreasing=FALSE)]
          y.CI.U.lit <- Lit.turnover3$CI.U.turnover.rate.d[order(Lit.turnover3$median.turnover.rate.d, decreasing=FALSE)]

          plot(y=c(1,1), x=c(0, max(x.values.lit)+3), type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", ylim=c(y.min, y.max), main="")
          arrows(x0=x.values.lit, y0=y.CI.L.lit, x1=x.values.lit, y1=y.CI.U.lit, length=0, angle=90, code=3, col="black", lwd=1.5)
          points(x=x.values.lit, y=y.medians.lit, pch=21, cex=1.0, col="black", bg="black", lwd=1.5)
          arrows(x0=150.48, y0=assemblage.turnover.d.CI.L, x1=150.48, y1=assemblage.turnover.d.CI.U, length=0, angle=90, code=3, col="black", lwd=1.5)
          points(x=150.48, y=assemblage.turnover.d.median, pch=21, cex=1.0, col="black", bg="white", lwd=1.5)    #death rate is more reflective of true turnover than is birth rate

          par(mgp=c(3,0,0))			#set spacing for the x-axis labels (the second value)
          axis(side=1, at=0:15, labels=NA, tck=0, las=1, cex.axis=1.0)
          mtext("Literature values", side=1, line=0.4, at=146/2, cex=0.75)
          mtext("This study", side=1, line=0.4, at=150.48, cex=0.75)
          par(mgp=c(3,0.45,0))				  #set spacing for the y-axis labels (the second value)
          axis(side=2, tck=-0.015, las=1, cex.axis=0.6)
          mtext(expression(paste("Assemblage turnover rate (day"^-1, ")", sep="")), side=2, line=1.7, cex=0.75)

          dev.off()
  #}


#Turnover meta-analysis code:__________________________________________________________________________________________________________________________________
  #{
    #Graph the quick & dirty literature bacterial assemblage turnover rates in comparison to those estimated in TM analysis): SORTED
      # Lit.turnover
        #Turnover rates - in units of day-1:
          graphics.off()
          # dev.new(width=3.0, height=3.0)
          pdf(file="qSIP_output/Figures/TM_assemblage_turnover_comparison_sorted.pdf", width=3.0, height=3.0)
          par(mai=c(0.29,0.47,0.02,0.02), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
          y.mins <- c(0, as.numeric(unlist(Lit.turnover[,4:6])), assemblage.turnover.d.CI.L)
          y.maxs <- c(0, as.numeric(unlist(Lit.turnover[,4:6])), assemblage.turnover.d.CI.U)
          y.min <- min(y.mins[y.mins != -Inf], na.rm=TRUE)
          y.max <- max(y.maxs[y.maxs != Inf], na.rm=TRUE)
          x.values.lit <- (1:dim(Lit.turnover)[1])-1

          plot(y=c(1,1), x=c(0,15), type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", ylim=c(y.min, y.max), main="")
          arrows(x0=x.values.lit, y0=Lit.turnover$min.turnover.rate.days[order(Lit.turnover$mean.turnover.rate.days, decreasing=FALSE)], x1=x.values.lit, y1=Lit.turnover$max.turnover.rate.days[order(Lit.turnover$mean.turnover.rate.days, decreasing=FALSE)], length=0, angle=90, code=3, col="black", lwd=1.5)
          points(x=x.values.lit, y=Lit.turnover$mean.turnover.rate.days[order(Lit.turnover$mean.turnover.rate.days, decreasing=FALSE)], pch=21, cex=1.0, col="black", bg="black", lwd=1.5)
          arrows(x0=13.3, y0=assemblage.turnover.d.CI.L, x1=13.3, y1=assemblage.turnover.d.CI.U, length=0, angle=90, code=3, col="black", lwd=1.5)
          points(x=13.3, y=assemblage.turnover.d.median, pch=21, cex=1.0, col="black", bg="white", lwd=1.5)    #death rate is more reflective of true turnover than is birth rate

          par(mgp=c(3,0,0))			#set spacing for the x-axis labels (the second value)
          axis(side=1, at=0:15, labels=NA, tck=0, las=1, cex.axis=1.0)
          mtext("Literature values", side=1, line=0.4, at=5.5, cex=0.75)
          mtext("This study", side=1, line=0.4, at=13.3, cex=0.75)
          par(mgp=c(3,0.45,0))				  #set spacing for the y-axis labels (the second value)
          axis(side=2, tck=-0.015, las=1, cex.axis=0.6)
          mtext(expression(paste("Assemblage turnover rate (day"^-1, ")", sep="")), side=2, line=1.7, cex=0.75)

          dev.off()

        #Calculate lower and upper range (excluding the 'outlier'):
          Lit.turnover$mean.turnover.rate.days[order(Lit.turnover$mean.turnover.rate.days, decreasing=FALSE)]
          rect.low <- Lit.turnover$mean.turnover.rate.days[order(Lit.turnover$mean.turnover.rate.days, decreasing=FALSE)][1]
          rect.high <- Lit.turnover$mean.turnover.rate.days[order(Lit.turnover$mean.turnover.rate.days, decreasing=FALSE)][11]

    #Graph the quick & dirty literature bacterial assemblage turnover rates in comparison to those estimated in TM analysis): SORTED WITH RECTANGLE
      # Lit.turnover
        #Turnover rates - in units of day-1:
          graphics.off()
          # dev.new(width=3.0, height=3.0)
          pdf(file="qSIP_output/Figures/TM_assemblage_turnover_comparison_sorted_rect.pdf", width=3.0, height=3.0)
          par(mai=c(0.29,0.47,0.02,0.02), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
          y.mins <- c(0, as.numeric(unlist(Lit.turnover[,4:6])), assemblage.turnover.d.CI.L)
          y.maxs <- c(0, as.numeric(unlist(Lit.turnover[,4:6])), assemblage.turnover.d.CI.U)
          y.min <- min(y.mins[y.mins != -Inf], na.rm=TRUE)
          y.max <- max(y.maxs[y.maxs != Inf], na.rm=TRUE)
          x.values.lit <- (1:dim(Lit.turnover)[1])-1

          plot(y=c(1,1), x=c(0,15), type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", ylim=c(y.min, y.max), main="")
          arrows(x0=x.values.lit, y0=Lit.turnover$min.turnover.rate.days[order(Lit.turnover$mean.turnover.rate.days, decreasing=FALSE)], x1=x.values.lit, y1=Lit.turnover$max.turnover.rate.days[order(Lit.turnover$mean.turnover.rate.days, decreasing=FALSE)], length=0, angle=90, code=3, col="black", lwd=1.5)
          points(x=x.values.lit, y=Lit.turnover$mean.turnover.rate.days[order(Lit.turnover$mean.turnover.rate.days, decreasing=FALSE)], pch=21, cex=1.0, col="black", bg="black", lwd=1.5)
          arrows(x0=13.3, y0=assemblage.turnover.d.CI.L, x1=13.3, y1=assemblage.turnover.d.CI.U, length=0, angle=90, code=3, col="black", lwd=1.5)
          points(x=13.3, y=assemblage.turnover.d.median, pch=21, cex=1.0, col="black", bg="white", lwd=1.5)    #death rate is more reflective of true turnover than is birth rate
          rect(xleft=-1000, ybottom=rect.low, xright=1000, ytop=rect.high, density=NA, border=NA, col=rgb(red=0, green=0, blue=0, alpha=0.25))

          par(mgp=c(3,0,0))			#set spacing for the x-axis labels (the second value)
          axis(side=1, at=0:15, labels=NA, tck=0, las=1, cex.axis=1.0)
          mtext("Literature values", side=1, line=0.4, at=5.5, cex=0.75)
          mtext("This study", side=1, line=0.4, at=13.3, cex=0.75)
          par(mgp=c(3,0.45,0))				  #set spacing for the y-axis labels (the second value)
          axis(side=2, tck=-0.015, las=1, cex.axis=0.6)
          mtext(expression(paste("Assemblage turnover rate (day"^-1, ")", sep="")), side=2, line=1.7, cex=0.75)

          dev.off()

    #Graph the complete set of literature bacterial assemblage turnover rates (the qSIP-derived estiate from the TM analysis is not yet included):
        #Turnover rates - in units of day-1:
          graphics.off()
          # dev.new(width=3.0, height=3.0)
          pdf(file="qSIP_output/Figures/TM_assemblage_turnover_complete.pdf", width=3.0, height=3.0)
          par(mai=c(0.29,0.47,0.02,0.02), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
          y.mins <- c(0, as.numeric(unlist(Lit.turnover3[,26:28])))
          y.maxs <- c(0, as.numeric(unlist(Lit.turnover3[,26:28])))
          y.min <- min(y.mins[y.mins != -Inf], na.rm=TRUE)
          y.max <- max(y.maxs[y.maxs != Inf], na.rm=TRUE)
          x.values.lit <- (1:dim(Lit.turnover3)[1])-1
          y.medians.lit <- Lit.turnover3$median.turnover.rate.d[order(Lit.turnover3$median.turnover.rate.d, decreasing=FALSE)]
          y.CI.L.lit <- Lit.turnover3$CI.L.turnover.rate.d[order(Lit.turnover3$median.turnover.rate.d, decreasing=FALSE)]
          y.CI.U.lit <- Lit.turnover3$CI.U.turnover.rate.d[order(Lit.turnover3$median.turnover.rate.d, decreasing=FALSE)]

          plot(y=c(1,1), x=c(0, max(x.values.lit)), type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", ylim=c(y.min, y.max), main="")
          arrows(x0=x.values.lit, y0=y.CI.L.lit, x1=x.values.lit, y1=y.CI.U.lit, length=0, angle=90, code=3, col="black", lwd=1.5)
          points(x=x.values.lit, y=y.medians.lit, pch=21, cex=1.0, col="black", bg="black", lwd=1.5)

          par(mgp=c(3,0,0))			#set spacing for the x-axis labels (the second value)
          axis(side=1, at=0:15, labels=NA, tck=0, las=1, cex.axis=1.0)
          mtext("Literature values", side=1, line=0.4, cex=0.75)
          par(mgp=c(3,0.45,0))				  #set spacing for the y-axis labels (the second value)
          axis(side=2, tck=-0.015, las=1, cex.axis=0.6)
          mtext(expression(paste("Assemblage turnover rate (day"^-1, ")", sep="")), side=2, line=1.7, cex=0.75)

          dev.off()

    #Graph the complete set of literature bacterial assemblage turnover rates WITH RECTANGLE AND LINE:
        #Turnover rates - in units of day-1:
          graphics.off()
          # dev.new(width=3.0, height=3.0)
          pdf(file="qSIP_output/Figures/TM_assemblage_turnover_complete_rect.pdf", width=3.0, height=3.0)
          par(mai=c(0.29,0.47,0.02,0.02), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
          y.mins <- c(0, as.numeric(unlist(Lit.turnover3[,26:28])))
          y.maxs <- c(0, as.numeric(unlist(Lit.turnover3[,26:28])))
          y.min <- min(y.mins[y.mins != -Inf], na.rm=TRUE)
          y.max <- max(y.maxs[y.maxs != Inf], na.rm=TRUE)
          x.values.lit <- (1:dim(Lit.turnover3)[1])-1
          y.medians.lit <- Lit.turnover3$median.turnover.rate.d[order(Lit.turnover3$median.turnover.rate.d, decreasing=FALSE)]
          y.CI.L.lit <- Lit.turnover3$CI.L.turnover.rate.d[order(Lit.turnover3$median.turnover.rate.d, decreasing=FALSE)]
          y.CI.U.lit <- Lit.turnover3$CI.U.turnover.rate.d[order(Lit.turnover3$median.turnover.rate.d, decreasing=FALSE)]

          plot(y=c(1,1), x=c(0, max(x.values.lit)), type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", ylim=c(y.min, y.max), main="")
          arrows(x0=x.values.lit, y0=y.CI.L.lit, x1=x.values.lit, y1=y.CI.U.lit, length=0, angle=90, code=3, col="black", lwd=1.5)
          points(x=x.values.lit, y=y.medians.lit, pch=21, cex=1.0, col="black", bg="black", lwd=1.5)
          rect(xleft=-1000, ybottom=rect.low, xright=1000, ytop=rect.high, density=NA, border=NA, col=rgb(red=0, green=0, blue=0, alpha=0.25))
          abline(h=assemblage.turnover.d.median, col="green", lwd=2)
      
          par(mgp=c(3,0,0))			#set spacing for the x-axis labels (the second value)
          axis(side=1, at=0:15, labels=NA, tck=0, las=1, cex.axis=1.0)
          mtext("Literature values", side=1, line=0.4, cex=0.75)
          par(mgp=c(3,0.45,0))				  #set spacing for the y-axis labels (the second value)
          axis(side=2, tck=-0.015, las=1, cex.axis=0.6)
          mtext(expression(paste("Assemblage turnover rate (day"^-1, ")", sep="")), side=2, line=1.7, cex=0.75)

          dev.off()

    #Graph the complete set of literature bacterial assemblage turnover rates (LOG=SCALE):
        #Turnover rates - in units of day-1:
          graphics.off()
          # dev.new(width=3.0, height=3.0)
          pdf(file="qSIP_output/Figures/TM_assemblage_turnover_complete_log10.pdf", width=3.0, height=3.0)
          par(mai=c(0.29,0.47,0.02,0.02), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
          y.mins <- log10(c(0, as.numeric(unlist(Lit.turnover3[,26:28]))))
          y.maxs <- log10(c(0, as.numeric(unlist(Lit.turnover3[,26:28]))))
          y.min <- min(y.mins[y.mins != -Inf], na.rm=TRUE)
          y.max <- max(y.maxs[y.maxs != Inf], na.rm=TRUE)
          x.values.lit <- (1:dim(Lit.turnover3)[1])-1
          y.medians.lit <- log10(Lit.turnover3$median.turnover.rate.d[order(Lit.turnover3$median.turnover.rate.d, decreasing=FALSE)])
          y.CI.L.lit <- log10(Lit.turnover3$CI.L.turnover.rate.d[order(Lit.turnover3$median.turnover.rate.d, decreasing=FALSE)])        #WARNING: SOME NEGATIVE VALUES EXIST & BECOME UNDEFINED WHEN LOGGED
          y.CI.U.lit <- log10(Lit.turnover3$CI.U.turnover.rate.d[order(Lit.turnover3$median.turnover.rate.d, decreasing=FALSE)])
          AT <- c(0.001, 0.01, 0.1, 1)

          plot(y=c(1,1), x=c(0, max(x.values.lit)), type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", ylim=c(y.min, y.max), main="")
          arrows(x0=x.values.lit, y0=y.CI.L.lit, x1=x.values.lit, y1=y.CI.U.lit, length=0, angle=90, code=3, col="black", lwd=1.5)
          points(x=x.values.lit, y=y.medians.lit, pch=21, cex=1.0, col="black", bg="black", lwd=1.5)

          par(mgp=c(3,0,0))			#set spacing for the x-axis labels (the second value)
          axis(side=1, at=0:15, labels=NA, tck=0, las=1, cex.axis=1.0)
          mtext("Literature values", side=1, line=0.4, cex=0.75)
          par(mgp=c(3,0.45,0))				  #set spacing for the y-axis labels (the second value)
          axis(side=2, at=log10(AT), labels=AT, tck=-0.015, las=1, cex.axis=0.6)
          mtext(expression(paste("Assemblage turnover rate (day"^-1, ")", sep="")), side=2, line=1.7, cex=0.75)

          dev.off()
      
    #Graph the complete set of literature bacterial assemblage turnover rates (LOG=SCALE) with RECTANGLE AND LINE:
        #Turnover rates - in units of day-1:
          graphics.off()
          # dev.new(width=3.0, height=3.0)
          pdf(file="qSIP_output/Figures/TM_assemblage_turnover_complete_log10_rect.pdf", width=3.0, height=3.0)
          par(mai=c(0.29,0.47,0.02,0.02), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
          y.mins <- log10(c(0, as.numeric(unlist(Lit.turnover3[,26:28]))))
          y.maxs <- log10(c(0, as.numeric(unlist(Lit.turnover3[,26:28]))))
          y.min <- min(y.mins[y.mins != -Inf], na.rm=TRUE)
          y.max <- max(y.maxs[y.maxs != Inf], na.rm=TRUE)
          x.values.lit <- (1:dim(Lit.turnover3)[1])-1
          y.medians.lit <- log10(Lit.turnover3$median.turnover.rate.d[order(Lit.turnover3$median.turnover.rate.d, decreasing=FALSE)])
          y.CI.L.lit <- log10(Lit.turnover3$CI.L.turnover.rate.d[order(Lit.turnover3$median.turnover.rate.d, decreasing=FALSE)])        #WARNING: SOME NEGATIVE VALUES EXIST & BECOME UNDEFINED WHEN LOGGED
          y.CI.U.lit <- log10(Lit.turnover3$CI.U.turnover.rate.d[order(Lit.turnover3$median.turnover.rate.d, decreasing=FALSE)])
          AT <- c(0.001, 0.01, 0.1, 1)

          plot(y=c(1,1), x=c(0, max(x.values.lit)), type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", ylim=c(y.min, y.max), main="")
          arrows(x0=x.values.lit, y0=y.CI.L.lit, x1=x.values.lit, y1=y.CI.U.lit, length=0, angle=90, code=3, col="black", lwd=1.5)
          points(x=x.values.lit, y=y.medians.lit, pch=21, cex=1.0, col="black", bg="black", lwd=1.5)
          rect(xleft=-1000, ybottom=log10(rect.low), xright=1000, ytop=log10(rect.high), density=NA, border=NA, col=rgb(red=0, green=0, blue=0, alpha=0.25))
          abline(h=log10(assemblage.turnover.d.median), col="green", lwd=2)

          par(mgp=c(3,0,0))			#set spacing for the x-axis labels (the second value)
          axis(side=1, at=0:15, labels=NA, tck=0, las=1, cex.axis=1.0)
          mtext("Literature values", side=1, line=0.4, cex=0.75)
          par(mgp=c(3,0.45,0))				  #set spacing for the y-axis labels (the second value)
          axis(side=2, at=log10(AT), labels=AT, tck=-0.015, las=1, cex.axis=0.6)
          mtext(expression(paste("Assemblage turnover rate (day"^-1, ")", sep="")), side=2, line=1.7, cex=0.75)

          dev.off()

    #Graph the complete set of literature bacterial assemblage turnover rates (LOG=SCALE) vs temperature (°C):
        #Turnover rates - in units of day-1:
          graphics.off()
          # dev.new(width=3.0, height=3.0)
          pdf(file="qSIP_output/Figures/TM_assemblage_turnover_complete_tempC.pdf", width=3.0, height=3.0)
          par(mai=c(0.34,0.47,0.02,0.02), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
          y.mins <- log10(c(0, as.numeric(unlist(Lit.turnover3[,26:28]))))
          y.maxs <- log10(c(0, as.numeric(unlist(Lit.turnover3[,26:28]))))
          y.min <- min(y.mins[y.mins != -Inf], na.rm=TRUE)
          y.max <- max(y.maxs[y.maxs != Inf], na.rm=TRUE)
          x.values <- Lit.turnover3$incub.temp.C
          y.medians.lit <- log10(Lit.turnover3$median.turnover.rate.d)
          y.CI.L.lit <- log10(Lit.turnover3$CI.L.turnover.rate.d)        #WARNING: SOME NEGATIVE VALUES EXIST & BECOME UNDEFINED WHEN LOGGED
          y.CI.U.lit <- log10(Lit.turnover3$CI.U.turnover.rate.d)
          AT <- c(0.001, 0.01, 0.1, 1)

          lm.temp.C <- lm(y.medians.lit~x.values)                        #WARNING: POINTS ARE NOT YET WEIGHTED BY CI's OR N

          plot(y=c(1,1), x=c(0, max(x.values, na.rm=TRUE)), type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", ylim=c(y.min, y.max), main="")
          arrows(x0=x.values, y0=y.CI.L.lit, x1=x.values, y1=y.CI.U.lit, length=0, angle=90, code=3, col="black", lwd=1.5)
          points(x=x.values, y=y.medians.lit, pch=21, cex=1.0, col="black", bg="black", lwd=1.5)
          regLine(lm.temp.C, col="blue", lwd=2)

          par(mgp=c(3,-0.1,0))			#set spacing for the x-axis labels (the second value)
          axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
          mtext("Temperature (°C)", side=1, line=0.90, cex=0.75)
          par(mgp=c(3,0.45,0))				  #set spacing for the y-axis labels (the second value)
          axis(side=2, at=log10(AT), labels=AT, tck=-0.015, las=1, cex.axis=0.6)
          mtext(expression(paste("Assemblage turnover rate (day"^-1, ")", sep="")), side=2, line=1.7, cex=0.75)

          summary(lm.temp.C)

          dev.off()

    #Graph the complete set of literature bacterial assemblage turnover rates (LOG=SCALE) vs temperature (1/kT):
        #Turnover rates - in units of day-1:
          graphics.off()
          # dev.new(width=3.0, height=3.0)
          pdf(file="qSIP_output/Figures/TM_assemblage_turnover_complete_temp_kT.pdf", width=3.0, height=3.0)
          par(mai=c(0.34,0.47,0.02,0.02), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
          y.mins <- log10(c(0, as.numeric(unlist(Lit.turnover3[,26:28]))))
          y.maxs <- log10(c(0, as.numeric(unlist(Lit.turnover3[,26:28]))))
          y.min <- min(y.mins[y.mins != -Inf], na.rm=TRUE)
          y.max <- max(y.maxs[y.maxs != Inf], na.rm=TRUE)
          x.values <- 1/(((Lit.turnover3$incub.temp.C) + 273.15)*8.617330350e-05)     #where k= Boltzmann's constant (8.6173303(50)×10−5 eV K-1)
          y.medians.lit <- log10(Lit.turnover3$median.turnover.rate.d)
          y.CI.L.lit <- log10(Lit.turnover3$CI.L.turnover.rate.d)         #WARNING: SOME NEGATIVE VALUES EXIST & BECOME UNDEFINED WHEN LOGGED
          y.CI.U.lit <- log10(Lit.turnover3$CI.U.turnover.rate.d)
          AT <- c(0.001, 0.01, 0.1, 1)

          set.seed(100)
          lm.temp.kT <- lm(y.medians.lit~x.values)                        #WARNING: POINTS ARE NOT YET WEIGHTED BY CI's OR N   /   WARNING: THIS IS OLS REGRESSION
          sma.temp.kT <- sma(y.medians.lit~x.values, method="SMA")        #SMA regression

          plot(y=c(1,1), x=c(min(x.values, na.rm=TRUE), max(x.values, na.rm=TRUE)), type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", ylim=c(y.min, y.max), main="")
          arrows(x0=x.values, y0=y.CI.L.lit, x1=x.values, y1=y.CI.U.lit, length=0, angle=90, code=3, col="black", lwd=1.5)
          points(x=x.values, y=y.medians.lit, pch=21, cex=1.0, col="black", bg="black", lwd=1.5)
          regLine(lm.temp.kT, col="blue", lwd=2)                  #regression line is for OLS
          segments(x0=min(x.values, na.rm=TRUE), y0=(min(x.values, na.rm=TRUE)*coef(sma.temp.kT)[2])+coef(sma.temp.kT)[1], x1=max(x.values, na.rm=TRUE), y1=(max(x.values, na.rm=TRUE)*coef(sma.temp.kT)[2])+coef(sma.temp.kT)[1], col="red", lwd=2, lty=2)     #regression line is for SMA

          par(mgp=c(3,-0.1,0))			#set spacing for the x-axis labels (the second value)
          axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
          mtext("Temperature (1/kT)", side=1, line=0.90, cex=0.75)
          par(mgp=c(3,0.45,0))				  #set spacing for the y-axis labels (the second value)
          axis(side=2, at=log10(AT), labels=AT, tck=-0.015, las=1, cex.axis=0.6)
          mtext(expression(paste("Assemblage turnover rate (day"^-1, ")", sep="")), side=2, line=1.7, cex=0.75)

          summary(lm.temp.kT)
          sma.temp.kT
          set.seed(100)
          sma(y.medians.lit~x.values, method="SMA", slope.test=0.63)      #SMA regression, testing if slope = predicted activation energy of 0.63eV (Brown et al. 2004)

          dev.off()

    #Graph the complete set of literature bacterial assemblage turnover rates (LOG=SCALE) vs incubation period:
        #Turnover rates - in units of day-1:
          graphics.off()
          # dev.new(width=3.0, height=3.0)
          pdf(file="qSIP_output/Figures/TM_assemblage_turnover_complete_incubation_period.pdf", width=3.0, height=3.0)
          par(mai=c(0.34,0.47,0.02,0.02), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
          y.mins <- log10(c(0, as.numeric(unlist(Lit.turnover3[,26:28]))))
          y.maxs <- log10(c(0, as.numeric(unlist(Lit.turnover3[,26:28]))))
          y.min <- min(y.mins[y.mins != -Inf], na.rm=TRUE)
          y.max <- max(y.maxs[y.maxs != Inf], na.rm=TRUE)
          x.values <- log10(Lit.turnover3$incub.period.d)
          y.medians.lit <- log10(Lit.turnover3$median.turnover.rate.d)
          y.CI.L.lit <- log10(Lit.turnover3$CI.L.turnover.rate.d)         #WARNING: SOME NEGATIVE VALUES EXIST & BECOME UNDEFINED WHEN LOGGED
          y.CI.U.lit <- log10(Lit.turnover3$CI.U.turnover.rate.d)
          AT <- c(0.001, 0.01, 0.1, 1)
          AT.X <- c(0.01, 0.1, 1, 10, 100, 1000)

          lm.incub.period.d <- lm(y.medians.lit~x.values)                 #WARNING: POINTS ARE NOT YET WEIGHTED BY CI's OR N

          plot(y=c(1,1), x=c(min(x.values, na.rm=TRUE), max(x.values, na.rm=TRUE)), type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", ylim=c(y.min, y.max), main="")
          arrows(x0=x.values, y0=y.CI.L.lit, x1=x.values, y1=y.CI.U.lit, length=0, angle=90, code=3, col="black", lwd=1.5)
          points(x=x.values, y=y.medians.lit, pch=21, cex=1.0, col="black", bg="black", lwd=1.5)
          regLine(lm.incub.period.d, col="blue", lwd=2)

          par(mgp=c(3,-0.1,0))			#set spacing for the x-axis labels (the second value)
          axis(side=1, at=log10(AT.X), labels=AT.X, tck=-0.015, las=1, cex.axis=0.6)
          mtext("Incubation period (days)", side=1, line=0.90, cex=0.75)
          par(mgp=c(3,0.45,0))				  #set spacing for the y-axis labels (the second value)
          axis(side=2, at=log10(AT), labels=AT, tck=-0.015, las=1, cex.axis=0.6)
          mtext(expression(paste("Assemblage turnover rate (day"^-1, ")", sep="")), side=2, line=1.7, cex=0.75)

          summary(lm.incub.period.d)

          dev.off()

# still to look at: water, methods, seiving, prokaryotic vs microbial (+fungi)

  #}




#Save workspace:
  save(list=ls(all.names=TRUE), file="qSIP_workspaces/TM_01-02-03-04-05/.RData", envir=.GlobalEnv)



