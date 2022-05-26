# This code produces figures and performs relevant statistical tests for the TM qSIP data


graphics.off()	#close all graphics windows


#Set working directory:
  #ALREADY DONE FOR THIS WORKSPACE IN PREVIOUSLY RUN CODE; only reset it here if loading previously saved workspace (see below)


#Reload the saved workspace resulting from the previous script:
  setwd("/Users/bk/Research/Projects/SIP_Modeling/qSIP")
  load("qSIP_workspaces/TM_01-02-03/.RData")


#Load libraries & scripts:
  library(car)


#Calculate assemblage-wide 16S abundance (copies g soil^-1) at the beginning and end of the incubation:________________________________________________________
  #{
    #Bootstrapped adundance at Time 0:
      quantile(as.numeric(apply(N.T0.boots.only, 2, sum)), probs=c(0.05, 0.50, 0.95))
    #Bootstrapped abundance at Time 0 (scientific notation):
      quantile(as.numeric(apply(N.T0.boots.only, 2, sum)), probs=c(0.05, 0.50, 0.95))/(10^10)
    #Bootstrapped adundance at Time t:
      quantile(as.numeric(apply(N.Tt.boots.only, 2, sum)), probs=c(0.05, 0.50, 0.95))
    #Bootstrapped abundance at Time t (scientific notation):
      quantile(as.numeric(apply(N.Tt.boots.only, 2, sum)), probs=c(0.05, 0.50, 0.95))/(10^9)
  #}


#Calculate the proportion of genera that did & did not become enriched in 18O at the level of 90%CIs:__________________________________________________________
  #{
    #Taxa that became enriched in 18O (lower CI of ape greater than zero):
      all.rates$taxonID[all.rates$ape.boot.CI.L > 0]
      length(all.rates$taxonID[all.rates$ape.boot.CI.L > 0])
      length(all.rates$taxonID[all.rates$ape.boot.CI.L > 0])/dim(all.rates)[1]    #as a proportion of all genera
    #Taxa that did not become enriched in 18O (CI of ape overlaps zero):
      all.rates$taxonID[all.rates$ape.boot.CI.U >= 0 & all.rates$ape.boot.CI.L <= 0]
      length(all.rates$taxonID[all.rates$ape.boot.CI.U >= 0 & all.rates$ape.boot.CI.L <= 0])
      length(all.rates$taxonID[all.rates$ape.boot.CI.U >= 0 & all.rates$ape.boot.CI.L <= 0])/dim(all.rates)[1]    #as a proportion of all genera
    #Taxa showing 'negative' enrichment in 18O (upper CI of ape less than zero):
      all.rates$taxonID[all.rates$ape.boot.CI.U < 0]
      length(all.rates$taxonID[all.rates$ape.boot.CI.U < 0])
      length(all.rates$taxonID[all.rates$ape.boot.CI.U < 0])/dim(all.rates)[1]    #as a proportion of all genera
  #}


#Calculate the proportion of genera with positive and negative rates of net population growth:_________________________________________________________________
  #{
    #Taxa with negative rates of net population growth (upper CI of r.net less than zero):
      all.rates$taxonID[all.rates$r.net.boot.CI.U < 0]
      length(all.rates$taxonID[all.rates$r.net.boot.CI.U < 0])
      length(all.rates$taxonID[all.rates$r.net.boot.CI.U < 0])/dim(all.rates)[1]    #as a proportion of all genera
    #Taxa with rates of net population growth no different than zero (CI of r.net overlaps zero):
      all.rates$taxonID[all.rates$r.net.boot.CI.U >= 0 & all.rates$r.net.boot.CI.L <= 0]
      length(all.rates$taxonID[all.rates$r.net.boot.CI.U >= 0 & all.rates$r.net.boot.CI.L <= 0])
      length(all.rates$taxonID[all.rates$r.net.boot.CI.U >= 0 & all.rates$r.net.boot.CI.L <= 0])/dim(all.rates)[1]    #as a proportion of all genera
    #Taxa with positive rates of net population growth (lower CI of r.net greater than zero):
      all.rates$taxonID[all.rates$r.net.boot.CI.L > 0]
      length(all.rates$taxonID[all.rates$r.net.boot.CI.L > 0])
      length(all.rates$taxonID[all.rates$r.net.boot.CI.L > 0])/dim(all.rates)[1]    #as a proportion of all genera
      taxa.id[taxa.id$taxon == all.rates$taxonID[all.rates$r.net.boot.CI.L > 0],]   #phylogenetic information for this taxon
  #}


#Calculate the proportion of genera with positive and negative rates of gross population growth (birth):_______________________________________________________
  #{
    #Taxa with negative rates of gross population growth (birth) (upper CI of b less than zero):
      all.rates$taxonID[all.rates$b.boot.CI.U < 0]
      length(all.rates$taxonID[all.rates$b.boot.CI.U < 0])
      length(all.rates$taxonID[all.rates$b.boot.CI.U < 0])/dim(all.rates)[1]    #as a proportion of all genera
    #Taxa with rates of gross population growth (birth) no different than zero (CI of b overlaps zero):
      all.rates$taxonID[all.rates$b.boot.CI.U >= 0 & all.rates$b.boot.CI.L <= 0]
      length(all.rates$taxonID[all.rates$b.boot.CI.U >= 0 & all.rates$b.boot.CI.L <= 0])
      length(all.rates$taxonID[all.rates$b.boot.CI.U >= 0 & all.rates$b.boot.CI.L <= 0])/dim(all.rates)[1]    #as a proportion of all genera
    #Taxa with positive rates of gross population growth (birth) (lower CI of b greater than zero):
      all.rates$taxonID[all.rates$b.boot.CI.L > 0]
      length(all.rates$taxonID[all.rates$b.boot.CI.L > 0])
      length(all.rates$taxonID[all.rates$b.boot.CI.L > 0])/dim(all.rates)[1]    #as a proportion of all genera
  #}


#Calculate the proportion of genera with positive and negative rates of mortality (death):_____________________________________________________________________
  #{
    #Taxa with negative rates of mortality (i.e., high mortality) (upper CI of d less than zero):
      all.rates$taxonID[all.rates$d.boot.CI.U < 0]
      length(all.rates$taxonID[all.rates$d.boot.CI.U < 0])
      length(all.rates$taxonID[all.rates$d.boot.CI.U < 0])/dim(all.rates)[1]    #as a proportion of all genera
    #Taxa with rates of mortality no different than zero (CI of d overlaps zero):
      all.rates$taxonID[all.rates$d.boot.CI.U >= 0 & all.rates$d.boot.CI.L <= 0]
      length(all.rates$taxonID[all.rates$d.boot.CI.U >= 0 & all.rates$d.boot.CI.L <= 0])
      length(all.rates$taxonID[all.rates$d.boot.CI.U >= 0 & all.rates$d.boot.CI.L <= 0])/dim(all.rates)[1]    #as a proportion of all genera
    #Taxa with positive rates of mortality (i.e., mortality rate is inverted to 'produce' cells) (lower CI of d greater than zero):
      all.rates$taxonID[all.rates$d.boot.CI.L > 0]
      length(all.rates$taxonID[all.rates$d.boot.CI.L > 0])
      length(all.rates$taxonID[all.rates$d.boot.CI.L > 0])/dim(all.rates)[1]    #as a proportion of all genera
  #}


#Calculate various carbon flux estimates for the assemblage and particular taxa:_______________________________________________________________________________
  #{
    #Calculate assemblage-level production (C flux) in micrograms C (originally was in picograms C):
      assemblage.f.gross.boots <- apply(f.gross.boots.only, 2, sum)/1000000
      assemblage.f.gross.median <- median(assemblage.f.gross.boots, na.rm=TRUE)
      assemblage.f.gross.CI.L <- quantile(assemblage.f.gross.boots, probs=0.05, na.rm=TRUE)
      assemblage.f.gross.CI.U <- quantile(assemblage.f.gross.boots, probs=0.95, na.rm=TRUE)
      assemblage.f.gross.median
      c(assemblage.f.gross.CI.L, assemblage.f.gross.CI.U)
    #Calculate assemblage-level loss of C to cell death in micrograms C (originally was in picograms C):
      assemblage.f.death.boots <- apply(f.death.boots.only, 2, sum)/1000000
      assemblage.f.death.median <- median(assemblage.f.death.boots, na.rm=TRUE)
      assemblage.f.death.CI.L <- quantile(assemblage.f.death.boots, probs=0.05, na.rm=TRUE)
      assemblage.f.death.CI.U <- quantile(assemblage.f.death.boots, probs=0.95, na.rm=TRUE)
      assemblage.f.death.median
      c(assemblage.f.death.CI.L, assemblage.f.death.CI.U)      
    #Assemblage-level loss of C to cell death was X times higher than asemblage-level prokaryotic production:
      -assemblage.f.death.median/assemblage.f.gross.median
    #Taxon with the highest production:
      #Production of that taxon:
        max(all.rates$f.gross.boot.median, na.rm=TRUE)
      #Identification of that taxon:
        all.rates$taxonID[all.rates$f.gross.boot.median == max(all.rates$f.gross.boot.median, na.rm=TRUE)]
        taxa.id[taxa.id$taxon == all.rates$taxonID[all.rates$f.gross.boot.median == max(all.rates$f.gross.boot.median, na.rm=TRUE)],]
      #Abundance (Tt) of that taxon:
        all.rates$N.tot.Tt.boot.median[all.rates$f.gross.boot.median == max(all.rates$f.gross.boot.median, na.rm=TRUE)]
      #Abundance (T0) of that taxon:
        all.rates$N.tot.T0.boot.median[all.rates$f.gross.boot.median == max(all.rates$f.gross.boot.median, na.rm=TRUE)]
      #Birth rate (b) of that taxon:
        all.rates$b.boot.median[all.rates$f.gross.boot.median == max(all.rates$f.gross.boot.median, na.rm=TRUE)]
      #Visualize how this taxon ranks in terms of production, abundance, and birth rate:
        graphics.off()
        par(mfrow=c(2,2))
        hist(all.rates$f.gross.boot.median, xlab="Production (pg C g soil-1 day-1)", main="")
        abline(v=max(all.rates$f.gross.boot.median, na.rm=TRUE), col="red", lwd=3)
        hist(all.rates$N.tot.Tt.boot.median, xlab="Day 10 abundance (16S copies g soil-1)", main="")
        abline(v=all.rates$N.tot.Tt.boot.median[all.rates$f.gross.boot.median == max(all.rates$f.gross.boot.median, na.rm=TRUE)], col="red", lwd=3)
        hist(all.rates$N.tot.T0.boot.median, xlab="Day 0 abundance (16S copies g soil-1)", main="")
        abline(v=all.rates$N.tot.T0.boot.median[all.rates$f.gross.boot.median == max(all.rates$f.gross.boot.median, na.rm=TRUE)], col="red", lwd=3)
        hist(all.rates$b.boot.median, xlab="b (day-1)", main="")
        abline(v=all.rates$b.boot.median[all.rates$f.gross.boot.median == max(all.rates$f.gross.boot.median, na.rm=TRUE)], col="red", lwd=3)
        par(mfrow=c(1,1))

        graphics.off()

    #A taxon with high production and high growth rate:
      #First identify those taxa with the hishest production:
        top.productive.taxa <- all.rates$taxonID[all.rates$f.gross.boot.median >= quantile(all.rates$f.gross.boot.median, probs=0.95, na.rm=TRUE)]
      #Birth rates (b) of those taxa:
        all.rates$b.boot.median[all.rates$taxonID %in% top.productive.taxa]
      #Abundances (Tt) of those taxa:
        all.rates$N.tot.Tt.boot.median[all.rates$taxonID %in% top.productive.taxa]
      #Table of all this info together:
        all.rates[all.rates$taxonID %in% top.productive.taxa, c(1,17,41)]
      #Visualizing how those taxa rank in terms of abundance, and birth rate:
        graphics.off()
        par(mfrow=c(2,2))
        plot(x=all.rates$N.tot.Tt.boot.median[all.rates$taxonID %in% top.productive.taxa], y=all.rates$b.boot.median[all.rates$taxonID %in% top.productive.taxa], bty="l", type="p", pch=21, col="black", bg="black", xlab="Day 10 abundance (16S copies g soil-1)", ylab="b (day-1)")
        text(x=all.rates$N.tot.Tt.boot.median[all.rates$taxonID %in% top.productive.taxa], y=all.rates$b.boot.median[all.rates$taxonID %in% top.productive.taxa], labels=top.productive.taxa, pos=4, cex=0.4)
        hist(all.rates$N.tot.Tt.boot.median, xlab="Day 10 abundance (16S copies g soil-1)", main="")
        abline(v=all.rates$N.tot.Tt.boot.median[all.rates$taxonID %in% top.productive.taxa], col="blue", lwd=1)
        par(xpd=NA)
        text(x=all.rates$N.tot.Tt.boot.median[all.rates$taxonID %in% top.productive.taxa], y=max(hist(all.rates$N.tot.Tt.boot.median, plot=FALSE)$counts), labels=top.productive.taxa, pos=3, cex=0.4)
        par(xpd=FALSE)
        hist(log10(all.rates$N.tot.Tt.boot.median), xlab="log10 day 10 abundance (16S copies g soil-1)", main="")
        abline(v=log10(all.rates$N.tot.Tt.boot.median[all.rates$taxonID %in% top.productive.taxa]), col="blue", lwd=1)
        par(xpd=NA)
        text(x=log10(all.rates$N.tot.Tt.boot.median[all.rates$taxonID %in% top.productive.taxa]), y=max(hist(log10(all.rates$N.tot.Tt.boot.median), plot=FALSE)$counts), labels=top.productive.taxa, pos=3, cex=0.4)
        par(xpd=FALSE)
        hist(all.rates$b.boot.median, xlab="b (day-1)", main="")
        abline(v=all.rates$b.boot.median[all.rates$taxonID %in% top.productive.taxa], col="blue", lwd=1)
        par(xpd=NA)
        text(x=all.rates$b.boot.median[all.rates$taxonID %in% top.productive.taxa], y=max(hist(all.rates$b.boot.median, plot=FALSE)$counts), labels=top.productive.taxa, pos=3, cex=0.4)
        par(mfrow=c(1,1))
        par(xpd=FALSE)
        par(mfrow=c(1,1))

      #Select taxon 120, according to the above plots:
        #Identification of that taxon:
          taxa.id[taxa.id$taxon == 120,]
        #Abundance (Tt) of that taxon:
          all.rates$N.tot.Tt.boot.median[all.rates$taxonID == 120]
        #Abundance (T0) of that taxon:
          all.rates$N.tot.T0.boot.median[all.rates$taxonID == 120]
        #Birth rate (b) of that taxon:
          all.rates$b.boot.median[all.rates$taxonID == 120]

    graphics.off()
  #}


#Look at birth rates of the phylum Firmicutes:_________________________________________________________________________________________________________________
  #{
    all.rates[all.rates$taxonID %in% taxa.id$taxon[taxa.id$phylum == "Firmicutes"], c(1,18)]
    taxa.id[taxa.id$phylum == "Firmicutes" & !is.na(taxa.id$phylum), c(1,3:8)]
  #}


#Look at the 'precision' of the birth rate estimates:__________________________________________________________________________________________________________
  #{
    apply(b.boots.only, 1, function(x) sd(x, na.rm=TRUE)/mean(x, na.rm=TRUE))
    hist(apply(b.boots.only, 1, function(x) sd(x, na.rm=TRUE)/mean(x, na.rm=TRUE)))
    abline(v=median(apply(b.boots.only, 1, function(x) sd(x, na.rm=TRUE)/mean(x, na.rm=TRUE))), col="blue")
    abline(v=mean(apply(b.boots.only, 1, function(x) sd(x, na.rm=TRUE)/mean(x, na.rm=TRUE))), col="red")
    median(apply(b.boots.only, 1, function(x) sd(x, na.rm=TRUE)/mean(x, na.rm=TRUE)))
    mean(apply(b.boots.only, 1, function(x) sd(x, na.rm=TRUE)/mean(x, na.rm=TRUE)))
  #}


#Graph the relationships between birth rates (new growth), death rates, net growth, and excess atom fracion 18O (rough graphs):________________________________
  #{
    plot(y=all.rates$b.obs, x=all.rates$ape.obs, , type="p", bty="l")
    dev.new()
    plot(y=all.rates$d.obs, x=all.rates$ape.obs, , type="p", bty="l")
    dev.new()
    plot(y=all.rates$r.net.obs, x=all.rates$ape.obs, , type="p", bty="l")
    
    graphics.off()
  #}


#Graph the relationship between birth rates (new growth) and excess atom fracion 18O:__________________________________________________________________________
  #{
    graphics.off()
    # dev.new(width=3.0, height=3.0)
    pdf(file="qSIP_output/Figures/TM_Taxa_birth_vs_eaf18O.pdf", width=3.0, height=3.0)
    #Birth vs. EAF 18O:
      par(mai=c(0.42,0.47,0.02,0.02), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
      x.min <- min(all.rates$ape.obs[all.rates$ape.obs != -Inf], na.rm=TRUE)
      x.max <- max(all.rates$ape.obs[all.rates$ape.obs != Inf], na.rm=TRUE)
      y.min <- min(all.rates$b.obs[all.rates$b.obs != -Inf], na.rm=TRUE)
      y.max <- max(all.rates$b.obs[all.rates$b.obs != Inf], na.rm=TRUE)

      b.add.panel <- function(DATA){
        plot(y=DATA$b.obs, x=DATA$ape.obs, type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min, x.max), ylim=c(y.min, y.max), main="")
        points(x=DATA$ape.obs, y=DATA$b.obs, pch=21, cex=0.6, col="black", bg="black")
        par(mgp=c(3,0,0))			        #set spacing for the x-axis labels (the second value)
        axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
        mtext(expression(paste("Excess atom fraction "^18, "O", sep="")), side=1, line=1.4, cex=0.75)
        par(mgp=c(3,0.45,0))				  #set spacing for the y-axis labels (the second value)
        axis(side=2, tck=-0.015, las=1, cex.axis=0.6)
        mtext(expression(paste(italic("b")~" (day"^-1, ")", sep="")), side=2, line=1.7, at=NA, cex=0.75)
      }

      b.add.panel(DATA=all.rates)

    rm(b.add.panel)

    dev.off()
  #}


#Graph the relationship between birth rates (new growth) and excess atom fracion 18O - formatted for powerpoint:_______________________________________________
  #{
    graphics.off()
    # dev.new(width=5.0, height=5.0)
    pdf(file="qSIP_output/Figures/TM_Taxa_birth_vs_eaf18O_slide.pdf", width=5.0, height=5.0)
    #Birth vs. EAF 18O:
      par(mai=c(0.58,0.80,0.03,0.03), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
      x.min <- min(all.rates$ape.obs[all.rates$ape.obs != -Inf], na.rm=TRUE)
      x.max <- max(all.rates$ape.obs[all.rates$ape.obs != Inf], na.rm=TRUE)
      y.min <- min(all.rates$b.obs[all.rates$b.obs != -Inf], na.rm=TRUE)
      y.max <- max(all.rates$b.obs[all.rates$b.obs != Inf], na.rm=TRUE)

      b.add.panel <- function(DATA){
        plot(y=DATA$b.obs, x=DATA$ape.obs, type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min, x.max), ylim=c(y.min, y.max), main="")
        points(x=DATA$ape.obs, y=DATA$b.obs, pch=21, cex=1.5, col="black", bg="black", lwd=1.5)
        par(mgp=c(3,0.2,0))			        #set spacing for the x-axis labels (the second value)
        axis(side=1, tck=-0.010, las=1, cex.axis=1.0)
        mtext(expression(paste("Excess atom fraction "^18, "O", sep="")), side=1, line=2.6, cex=1.6)
        par(mgp=c(3,0.45,0))				  #set spacing for the y-axis labels (the second value)
        axis(side=2, tck=-0.010, las=1, cex.axis=1.0)
        mtext(expression(paste(italic("b")~" (day"^-1, ")", sep="")), side=2, line=2.8, at=NA, cex=1.6)
      }

      b.add.panel(DATA=all.rates)

    rm(b.add.panel)

    dev.off()
  #}


#Graph the bootstrapped excess atom fraction of 18O (horizontal layout):_______________________________________________________________________________________
  #{
    #Growth rates - in units of r:
      graphics.off()
      # dev.new(width=7.5, height=4)
      pdf(file="qSIP_output/Figures/TM_Taxa_PhylumGroups_eaf18O_horizontal.pdf", width=7.5, height=4)
      par(mai=c(0.29,0.47,0.60,0.03), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
      #Re-code NA's in the taxa.id data.frame as "Other" so that taxa with a phylum of 'NA' is included in the phylum-specific plot:
        taxa.id.temp <- taxa.id
        taxa.id.temp[is.na(taxa.id.temp)] <- "Other"
      phyla <- levels(data.melted$phylum)
      phyla.cols <- cols.for.phyla$col
      ranks <- order(all.rates$ape.boot.median, decreasing=TRUE)
      tax.order <- data.frame(axis.loc=seq(1, dim(all.rates)[1]), ranks)
      y.mins.r <- c(all.rates$ape.boot.median, all.rates$ape.boot.CI.L)
      y.maxs.r <- c(all.rates$ape.boot.median, all.rates$ape.boot.CI.U)
      y.min.r <- min(y.mins.r[y.mins.r != -Inf], na.rm=TRUE)
      y.max.r <- max(y.maxs.r[y.maxs.r != Inf], na.rm=TRUE)

      eaf.add.panel <- function(DATA){
        plot(y=DATA$ape.boot.median[tax.order$ranks], x=1:dim(DATA)[1], type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", ylim=c(y.min.r, y.max.r), main="")
        counter <- 1
        for (p in 1:length(phyla)){
          curr.comp.ranked <- DATA[tax.order$ranks,]
          mids <- curr.comp.ranked$ape.boot.median[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
          lowers <- curr.comp.ranked$ape.boot.CI.L[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
          uppers <- curr.comp.ranked$ape.boot.CI.U[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
          tax.nums <- counter:(counter+length(mids)-1)
          counter <- counter+length(mids)
          arrows(x0=tax.nums, y0=lowers, x1=tax.nums, y1=uppers, length=0, angle=90, code=3, col=as.character(phyla.cols[p]))
          # points(x=tax.nums, y=mids, pch=21, cex=0.6, col=as.character(phyla.cols[p]), bg=as.character(phyla.cols[p]))
          points(x=tax.nums, y=mids, pch=21, cex=0.6, col="gray30", bg=as.character(phyla.cols[p]))
          par(xpd=NA)
          if(!is.element(p, c(3,5,6,8,9,10,11,13,14,16,17,18,21,22))){
            text(x=mean(tax.nums), y=y.max.r, labels=phyla[p], adj=c(0,0.5), pos=3, offset=1.5, cex=0.6, col=as.character(phyla.cols)[p])
          }
          else if (is.element(p, c(3,5,8,13,21))){
            text(x=mean(tax.nums), y=y.max.r, labels=phyla[p], adj=c(0,0.5), pos=3, offset=2.0, cex=0.6, col=as.character(phyla.cols)[p])
          }
          else if (is.element(p, c(6,9,14,22))){
            text(x=mean(tax.nums), y=y.max.r, labels=phyla[p], adj=c(0,0.5), pos=3, offset=1.0, cex=0.6, col=as.character(phyla.cols)[p])
          }
          else if (is.element(p, c(10,16))){
            text(x=mean(tax.nums), y=y.max.r, labels=phyla[p], adj=c(0,0.5), pos=3, offset=2.5, cex=0.6, col=as.character(phyla.cols)[p])
          }
          else if (is.element(p, c(11,17))){
            text(x=mean(tax.nums), y=y.max.r, labels=phyla[p], adj=c(0,0.5), pos=3, offset=0.5, cex=0.6, col=as.character(phyla.cols)[p])
          }
          else if (is.element(p, c(18))){
            text(x=mean(tax.nums), y=y.max.r, labels=phyla[p], adj=c(0,0.5), pos=3, offset=3.0, cex=0.6, col=as.character(phyla.cols)[p])
          }
          par(xpd=FALSE)
        }
        abline(h=0, col="black", lty=2)
        par(mgp=c(3,0,0))			        #set spacing for the x-axis labels (the second value)
        axis(side=1, at=1:dim(DATA)[1], labels=rep(NA, 103), tck=-0.015, las=1, cex.axis=0.6)
        mtext("Genus", side=1, line=0.4, cex=0.75)
        par(mgp=c(3,0.45,0))				  #set spacing for the y-axis labels (the second value)
        axis(side=2, tck=-0.015, las=1, cex.axis=0.6)
        mtext(expression(paste("Excess atom fraction "^18, "O", sep="")), side=2, line=1.7, at=NA, cex=0.75)
      }

      eaf.add.panel(DATA=all.rates)
      rm(eaf.add.panel, taxa.id.temp)
      par(mfrow=c(1,1))

      dev.off()
      
      #Phyla represented in the data ('NA' is included in 'Other' in the graph):
      phyla
      length(phyla)
  #}


#Graph the bootstrapped excess atom fraction of 18O (horizontal layout - formatted for powerpoint):____________________________________________________________
  #{
    #Growth rates - in units of r:
      graphics.off()
      # dev.new(width=10, height=6)
      pdf(file="qSIP_output/Figures/TM_Taxa_PhylumGroups_eaf18O_horizontal_slide.pdf", width=10, height=6)
      par(mai=c(0.49,0.80,1.03,0.03), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
      #Re-code NA's in the taxa.id data.frame as "Other" so that taxa with a phylum of 'NA' is included in the phylum-specific plot:
        taxa.id.temp <- taxa.id
        taxa.id.temp[is.na(taxa.id.temp)] <- "Other"
      phyla <- levels(data.melted$phylum)
      phyla.cols <- cols.for.phyla$col
      ranks <- order(all.rates$ape.boot.median, decreasing=TRUE)
      tax.order <- data.frame(axis.loc=seq(1, dim(all.rates)[1]), ranks)
      y.mins.r <- c(all.rates$ape.boot.median, all.rates$ape.boot.CI.L)
      y.maxs.r <- c(all.rates$ape.boot.median, all.rates$ape.boot.CI.U)
      y.min.r <- min(y.mins.r[y.mins.r != -Inf], na.rm=TRUE)
      y.max.r <- max(y.maxs.r[y.maxs.r != Inf], na.rm=TRUE)

      eaf.add.panel <- function(DATA){
        plot(y=DATA$ape.boot.median[tax.order$ranks], x=1:dim(DATA)[1], type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", ylim=c(y.min.r, y.max.r), main="")
        counter <- 1
        for (p in 1:length(phyla)){
          curr.comp.ranked <- DATA[tax.order$ranks,]
          mids <- curr.comp.ranked$ape.boot.median[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
          lowers <- curr.comp.ranked$ape.boot.CI.L[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
          uppers <- curr.comp.ranked$ape.boot.CI.U[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
          tax.nums <- counter:(counter+length(mids)-1)
          counter <- counter+length(mids)
          arrows(x0=tax.nums, y0=lowers, x1=tax.nums, y1=uppers, length=0, angle=90, code=3, col=as.character(phyla.cols[p]), lwd=1.5)
          # points(x=tax.nums, y=mids, pch=21, cex=1.5, col=as.character(phyla.cols[p]), bg=as.character(phyla.cols[p]))
          points(x=tax.nums, y=mids, pch=21, cex=1.5, col="gray30", bg=as.character(phyla.cols[p]))
          par(xpd=NA)
          if(!is.element(p, c(2,3,5,6,7,8,10,11,12,13,14,16,17,18,21,22))){
            text(x=mean(tax.nums), y=y.max.r, labels=phyla[p], adj=c(0,0.5), pos=3, offset=2.10, cex=1.1, col=as.character(phyla.cols)[p])
          }
          else if (is.element(p, c(3,7,13,22))){
            text(x=mean(tax.nums), y=y.max.r, labels=phyla[p], adj=c(0,0.5), pos=3, offset=2.90, cex=1.1, col=as.character(phyla.cols)[p])
          }
          else if (is.element(p, c(2,5,11,14))){
            text(x=mean(tax.nums), y=y.max.r, labels=phyla[p], adj=c(0,0.5), pos=3, offset=1.30, cex=1.1, col=as.character(phyla.cols)[p])
          }
          else if (is.element(p, c(21))){
            text(x=mean(tax.nums)-5, y=y.max.r, labels=phyla[p], adj=c(0,0.5), pos=3, offset=1.30, cex=1.1, col=as.character(phyla.cols)[p])
          }
          else if (is.element(p, c(8,16))){
            text(x=mean(tax.nums), y=y.max.r, labels=phyla[p], adj=c(0,0.5), pos=3, offset=3.70, cex=1.1, col=as.character(phyla.cols)[p])
          }
          else if (is.element(p, c(6,12,17))){
            text(x=mean(tax.nums), y=y.max.r, labels=phyla[p], adj=c(0,0.5), pos=3, offset=0.50, cex=1.1, col=as.character(phyla.cols)[p])
          }
          else if (is.element(p, c(10,18))){
            text(x=mean(tax.nums), y=y.max.r, labels=phyla[p], adj=c(0,0.5), pos=3, offset=4.50, cex=1.1, col=as.character(phyla.cols)[p])
          }
          par(xpd=FALSE)
        }
        abline(h=0, col="black", lty=2, lwd=1.5)
        par(mgp=c(3,0,0))			        #set spacing for the x-axis labels (the second value)
        axis(side=1, at=1:dim(DATA)[1], labels=rep(NA, 103), tck=-0.010, las=1, cex.axis=1.0)
        mtext("Genus", side=1, line=1.4, cex=1.6)
        par(mgp=c(3,0.45,0))				  #set spacing for the y-axis labels (the second value)
        axis(side=2, tck=-0.010, las=1, cex.axis=1.0)
        mtext(expression(paste("Excess atom fraction "^18, "O", sep="")), side=2, line=2.6, at=NA, cex=1.6)
      }

      eaf.add.panel(DATA=all.rates)
      rm(eaf.add.panel, taxa.id.temp)
      par(mfrow=c(1,1))

      dev.off()
    
      #Phyla represented in the data ('NA' is included in 'Other' in the graph):
      phyla
      length(phyla)
  #}


#Graph the bootstrapped excess atom fraction of 18O (horizontal layout, publication width):____________________________________________________________________
  #{
    #Growth rates - in units of r:
      graphics.off()
      # dev.new(width=6, height=4)
      pdf(file="qSIP_output/Figures/TM_Taxa_PhylumGroups_eaf18O_horizontal_pub_width.pdf", width=6, height=4)
      par(mai=c(0.29,0.47,0.60,0.03), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
      #Re-code NA's in the taxa.id data.frame as "Other" so that taxa with a phylum of 'NA' is included in the phylum-specific plot:
        taxa.id.temp <- taxa.id
        taxa.id.temp[is.na(taxa.id.temp)] <- "Other"
      phyla <- levels(data.melted$phylum)
      phyla.cols <- cols.for.phyla$col
      ranks <- order(all.rates$ape.boot.median, decreasing=TRUE)
      tax.order <- data.frame(axis.loc=seq(1, dim(all.rates)[1]), ranks)
      y.mins.r <- c(all.rates$ape.boot.median, all.rates$ape.boot.CI.L)
      y.maxs.r <- c(all.rates$ape.boot.median, all.rates$ape.boot.CI.U)
      y.min.r <- min(y.mins.r[y.mins.r != -Inf], na.rm=TRUE)
      y.max.r <- max(y.maxs.r[y.maxs.r != Inf], na.rm=TRUE)

      eaf.add.panel <- function(DATA){
        plot(y=DATA$ape.boot.median[tax.order$ranks], x=1:dim(DATA)[1], type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", ylim=c(y.min.r, y.max.r), main="")
        counter <- 1
        for (p in 1:length(phyla)){
          curr.comp.ranked <- DATA[tax.order$ranks,]
          mids <- curr.comp.ranked$ape.boot.median[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
          lowers <- curr.comp.ranked$ape.boot.CI.L[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
          uppers <- curr.comp.ranked$ape.boot.CI.U[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
          tax.nums <- counter:(counter+length(mids)-1)
          counter <- counter+length(mids)
          arrows(x0=tax.nums, y0=lowers, x1=tax.nums, y1=uppers, length=0, angle=90, code=3, col=as.character(phyla.cols[p]))
          # points(x=tax.nums, y=mids, pch=21, cex=0.6, col=as.character(phyla.cols[p]), bg=as.character(phyla.cols[p]))
          points(x=tax.nums, y=mids, pch=21, cex=0.6, col="gray30", bg=as.character(phyla.cols[p]))
          par(xpd=NA)
          if(!is.element(p, c(3,4,5,6,8,9,10,11,13,14,16,17,18,21,22))){
            text(x=mean(tax.nums), y=y.max.r, labels=phyla[p], adj=c(0,0.5), pos=3, offset=1.5, cex=0.6, col=as.character(phyla.cols)[p])
          }
          else if (is.element(p, c(5,8,13))){
            text(x=mean(tax.nums), y=y.max.r, labels=phyla[p], adj=c(0,0.5), pos=3, offset=2.0, cex=0.6, col=as.character(phyla.cols)[p])
          }
          else if (is.element(p, c(21))){     #shift taxon 21 to the left so as to fit on the page:
            text(x=mean(tax.nums)-2, y=y.max.r, labels=phyla[p], adj=c(0,0.5), pos=3, offset=2.0, cex=0.6, col=as.character(phyla.cols)[p])
          }
          else if (is.element(p, c(4,9,14,22))){
            text(x=mean(tax.nums), y=y.max.r, labels=phyla[p], adj=c(0,0.5), pos=3, offset=1.0, cex=0.6, col=as.character(phyla.cols)[p])
          }
          else if (is.element(p, c(3,10,16))){
            text(x=mean(tax.nums), y=y.max.r, labels=phyla[p], adj=c(0,0.5), pos=3, offset=2.5, cex=0.6, col=as.character(phyla.cols)[p])
          }
          else if (is.element(p, c(6,11,17))){
            text(x=mean(tax.nums), y=y.max.r, labels=phyla[p], adj=c(0,0.5), pos=3, offset=0.5, cex=0.6, col=as.character(phyla.cols)[p])
          }
          else if (is.element(p, c(18))){
            text(x=mean(tax.nums), y=y.max.r, labels=phyla[p], adj=c(0,0.5), pos=3, offset=3.0, cex=0.6, col=as.character(phyla.cols)[p])
          }
          par(xpd=FALSE)
        }
        abline(h=0, col="black", lty=2)
        par(mgp=c(3,0,0))			        #set spacing for the x-axis labels (the second value)
        axis(side=1, at=1:dim(DATA)[1], labels=rep(NA, 103), tck=-0.015, las=1, cex.axis=0.6)
        mtext("Genus", side=1, line=0.4, cex=0.75)
        par(mgp=c(3,0.45,0))				  #set spacing for the y-axis labels (the second value)
        axis(side=2, tck=-0.015, las=1, cex.axis=0.6)
        mtext(expression(paste("Excess atom fraction "^18, "O", sep="")), side=2, line=1.7, at=NA, cex=0.75)
      }

      eaf.add.panel(DATA=all.rates)
      rm(eaf.add.panel, taxa.id.temp)
      par(mfrow=c(1,1))

      dev.off()
      
      #Phyla represented in the data ('NA' is included in 'Other' in the graph):
      phyla
      length(phyla)
  #}


#Graph the bootstrapped excess atom fraction of 18O in .eps format (horizontal layout, publication width):_____________________________________________________
  #{
    #Growth rates - in units of r:
      graphics.off()
      # dev.new(width=6, height=4)
      setEPS(width=6, height=4)
      postscript(file="qSIP_output/Figures/TM_Taxa_PhylumGroups_eaf18O_horizontal_pub_width.eps")
      par(mai=c(0.29,0.47,0.60,0.03), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
      #Re-code NA's in the taxa.id data.frame as "Other" so that taxa with a phylum of 'NA' is included in the phylum-specific plot:
        taxa.id.temp <- taxa.id
        taxa.id.temp[is.na(taxa.id.temp)] <- "Other"
      phyla <- levels(data.melted$phylum)
      phyla.cols <- cols.for.phyla$col
      ranks <- order(all.rates$ape.boot.median, decreasing=TRUE)
      tax.order <- data.frame(axis.loc=seq(1, dim(all.rates)[1]), ranks)
      y.mins.r <- c(all.rates$ape.boot.median, all.rates$ape.boot.CI.L)
      y.maxs.r <- c(all.rates$ape.boot.median, all.rates$ape.boot.CI.U)
      y.min.r <- min(y.mins.r[y.mins.r != -Inf], na.rm=TRUE)
      y.max.r <- max(y.maxs.r[y.maxs.r != Inf], na.rm=TRUE)

      eaf.add.panel <- function(DATA){
        plot(y=DATA$ape.boot.median[tax.order$ranks], x=1:dim(DATA)[1], type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", ylim=c(y.min.r, y.max.r), main="")
        counter <- 1
        for (p in 1:length(phyla)){
          curr.comp.ranked <- DATA[tax.order$ranks,]
          mids <- curr.comp.ranked$ape.boot.median[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
          lowers <- curr.comp.ranked$ape.boot.CI.L[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
          uppers <- curr.comp.ranked$ape.boot.CI.U[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
          tax.nums <- counter:(counter+length(mids)-1)
          counter <- counter+length(mids)
          arrows(x0=tax.nums, y0=lowers, x1=tax.nums, y1=uppers, length=0, angle=90, code=3, col=as.character(phyla.cols[p]))
          # points(x=tax.nums, y=mids, pch=21, cex=0.6, col=as.character(phyla.cols[p]), bg=as.character(phyla.cols[p]))
          points(x=tax.nums, y=mids, pch=21, cex=0.6, col="gray30", bg=as.character(phyla.cols[p]))
          par(xpd=NA)
          if(!is.element(p, c(3,4,5,6,8,9,10,11,13,14,16,17,18,21,22))){
            text(x=mean(tax.nums), y=y.max.r, labels=phyla[p], adj=c(0,0.5), pos=3, offset=1.5, cex=0.6, col=as.character(phyla.cols)[p])
          }
          else if (is.element(p, c(5,8,13))){
            text(x=mean(tax.nums), y=y.max.r, labels=phyla[p], adj=c(0,0.5), pos=3, offset=2.0, cex=0.6, col=as.character(phyla.cols)[p])
          }
          else if (is.element(p, c(21))){     #shift taxon 21 to the left so as to fit on the page:
            text(x=mean(tax.nums)-2, y=y.max.r, labels=phyla[p], adj=c(0,0.5), pos=3, offset=2.0, cex=0.6, col=as.character(phyla.cols)[p])
          }
          else if (is.element(p, c(4,9,14,22))){
            text(x=mean(tax.nums), y=y.max.r, labels=phyla[p], adj=c(0,0.5), pos=3, offset=1.0, cex=0.6, col=as.character(phyla.cols)[p])
          }
          else if (is.element(p, c(3,10,16))){
            text(x=mean(tax.nums), y=y.max.r, labels=phyla[p], adj=c(0,0.5), pos=3, offset=2.5, cex=0.6, col=as.character(phyla.cols)[p])
          }
          else if (is.element(p, c(6,11,17))){
            text(x=mean(tax.nums), y=y.max.r, labels=phyla[p], adj=c(0,0.5), pos=3, offset=0.5, cex=0.6, col=as.character(phyla.cols)[p])
          }
          else if (is.element(p, c(18))){
            text(x=mean(tax.nums), y=y.max.r, labels=phyla[p], adj=c(0,0.5), pos=3, offset=3.0, cex=0.6, col=as.character(phyla.cols)[p])
          }
          par(xpd=FALSE)
        }
        abline(h=0, col="black", lty=2)
        par(mgp=c(3,0,0))			        #set spacing for the x-axis labels (the second value)
        axis(side=1, at=1:dim(DATA)[1], labels=rep(NA, 103), tck=-0.015, las=1, cex.axis=0.6)
        mtext("Genus", side=1, line=0.4, cex=0.75)
        par(mgp=c(3,0.45,0))				  #set spacing for the y-axis labels (the second value)
        axis(side=2, tck=-0.015, las=1, cex.axis=0.6)
        mtext(expression(paste("Excess atom fraction "^18, "O", sep="")), side=2, line=1.7, at=NA, cex=0.75)
      }

      eaf.add.panel(DATA=all.rates)
      rm(eaf.add.panel, taxa.id.temp)
      par(mfrow=c(1,1))

      dev.off()
    
      #Phyla represented in the data ('NA' is included in 'Other' in the graph):
      phyla
      length(phyla)
  #}


#Graph the bootstrapped birth, death, and net growth and flux rates:___________________________________________________________________________________________
  #{
    #Growth & C flux rates - in units of r and C flux:
      graphics.off()
      # dev.new(width=7.5, height=7.5)
      pdf(file="qSIP_output/Figures/TM_Taxa_same_order_PhylumGroups_growth&flux.pdf", width=7.5, height=7.5)
      par(mfcol=c(1,2))
      par(mai=c(0.44,0.42,0.2,0.05), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
      #Re-code NA's in the taxa.id data.frame as "Other" so that taxa with a phylum of 'NA' is included in the phylum-specific plot:
        taxa.id.temp <- taxa.id
        taxa.id.temp[is.na(taxa.id.temp)] <- "Other"
      phyla <- levels(data.melted$phylum)
      phyla.cols <- cols.for.phyla$col
      ranks <- order(all.rates$r.net.boot.median)
      tax.order <- data.frame(axis.loc=seq(1, dim(all.rates)[1]), ranks)
      x.mins.r <- c(all.rates$r.net.boot.median, all.rates$r.net.boot.CI.L, all.rates$b.boot.CI.L, all.rates$d.boot.CI.L)
      x.maxs.r <- c(all.rates$r.net.boot.median, all.rates$r.net.boot.CI.U, all.rates$b.boot.CI.U, all.rates$d.boot.CI.U)
      x.mins.f <- c(all.rates$f.net.boot.median, all.rates$f.net.boot.CI.L, all.rates$f.gross.boot.CI.L, all.rates$f.death.boot.CI.L)
      x.maxs.f <- c(all.rates$f.net.boot.median, all.rates$f.net.boot.CI.U, all.rates$f.gross.boot.CI.U, all.rates$f.death.boot.CI.U)
      x.min.r <- min(x.mins.r[x.mins.r != -Inf], na.rm=TRUE)
      x.max.r <- max(x.maxs.r[x.maxs.r != Inf], na.rm=TRUE)
      x.min.f <- min(x.mins.f[x.mins.f != -Inf], na.rm=TRUE)/1000
      x.max.f <- max(x.maxs.f[x.maxs.f != Inf], na.rm=TRUE)/1000
  
      r.add.panel <- function(DATA){
        plot(y=1:dim(DATA)[1], x=DATA$r.net.boot.median[tax.order$ranks], type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min.r, x.max.r), main="")
        counter <- 1
        for (p in 1:length(phyla)){
          curr.comp.ranked <- DATA[tax.order$ranks,]
          mids.net <- curr.comp.ranked$r.net.boot.median[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
          lowers.net <- curr.comp.ranked$r.net.boot.CI.L[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
          uppers.net <- curr.comp.ranked$r.net.boot.CI.U[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
          tax.nums <- counter:(counter+length(mids.net)-1)
          counter <- counter+length(mids.net)
          points(x=mids.net, y=tax.nums, pch=21, cex=0.6, col=as.character(phyla.cols[p]), bg=as.character(phyla.cols[p]))
          arrows(x0=lowers.net, y0=tax.nums, x1=uppers.net, y1=tax.nums, length=0, angle=90, code=3, col=as.character(phyla.cols[p]))
          par(xpd=NA)
          if(!is.element(p, c(5,6,9,10,16,17,20))){
            text(x=x.max.r, y=mean(tax.nums), labels=phyla[p], adj=c(0,0.5), pos=4, offset=-1.5, cex=0.6, col=as.character(phyla.cols)[p])
          }
          else if (is.element(p, c(5,10,20))){
            text(x=x.max.r, y=mean(tax.nums), labels=phyla[p], adj=c(0,0.5), pos=4, offset=0, cex=0.6, col=as.character(phyla.cols)[p])
          }
          else if (is.element(p, c(16))){
            text(x=x.max.r, y=mean(tax.nums), labels=phyla[p], adj=c(0,0.5), pos=4, offset=1.5, cex=0.6, col=as.character(phyla.cols)[p])
          }
          else if (is.element(p, c(6,9,17))){
            text(x=x.max.r, y=mean(tax.nums), labels=phyla[p], adj=c(0,0.5), pos=4, offset=3, cex=0.6, col=as.character(phyla.cols)[p])
          }
          par(xpd=FALSE)
          mids.b <- curr.comp.ranked$b.boot.median[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
          lowers.b <- curr.comp.ranked$b.boot.CI.L[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
          uppers.b <- curr.comp.ranked$b.boot.CI.U[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
          arrows(x0=lowers.b, y0=tax.nums, x1=uppers.b, y1=tax.nums, length=0, angle=90, code=3, col=as.character(phyla.cols[p]), lwd=0.5)
          points(x=mids.b, y=tax.nums, pch=21, cex=0.6, col=as.character(phyla.cols[p]), bg="white")
          mids.d <- curr.comp.ranked$d.boot.median[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
          lowers.d <- curr.comp.ranked$d.boot.CI.L[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
          uppers.d <- curr.comp.ranked$d.boot.CI.U[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
          arrows(x0=lowers.d, y0=tax.nums, x1=uppers.d, y1=tax.nums, length=0, angle=90, code=3, col=as.character(phyla.cols[p]), lwd=0.5)
          points(x=mids.d, y=tax.nums, pch=21, cex=0.6, col="gray30", bg=as.character(phyla.cols[p]))
        }
        abline(v=0, col="black")
        par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
        axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
        mtext(expression(paste(italic("r")~" (day"^-1, ")", sep="")), side=1, line=1.4, cex=0.75)
        par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
        axis(side=2, at=1:dim(DATA)[1], labels=rep(NA, 103), tck=-0.015, las=1, cex.axis=0.6)
        mtext("Taxon", side=2, line=1.35, cex=0.75)
      }
  
      f.add.panel <- function(DATA, log=FALSE){
        if (log == FALSE){
          plot(y=1:dim(DATA)[1], x=DATA$f.net.boot.median[tax.order$ranks]/1000, type="n", bty="n", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min.f, x.max.f), main="")
          counter <- 1
          for (p in 1:length(phyla)){
            curr.comp.ranked <- DATA[tax.order$ranks,]
            mids <- curr.comp.ranked$f.net.boot.median[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
            lowers <- curr.comp.ranked$f.net.boot.CI.L[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
            uppers <- curr.comp.ranked$f.net.boot.CI.U[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
            tax.nums <- counter:(counter+length(mids)-1)
            counter <- counter+length(mids)
            points(x=mids/1000, y=tax.nums, pch=21, cex=0.6, col=as.character(phyla.cols[p]), bg=as.character(phyla.cols[p]))
            arrows(x0=lowers/1000, y0=tax.nums, x1=uppers/1000, y1=tax.nums, length=0, angle=90, code=3, col=as.character(phyla.cols[p]))
            mids.gross <- curr.comp.ranked$f.gross.boot.median[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
            lowers.gross <- curr.comp.ranked$f.gross.boot.CI.L[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
            uppers.gross <- curr.comp.ranked$f.gross.boot.CI.U[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
            arrows(x0=lowers.gross/1000, y0=tax.nums, x1=uppers.gross/1000, y1=tax.nums, length=0, angle=90, code=3, col=as.character(phyla.cols[p]), lwd=0.5)
            points(x=mids.gross/1000, y=tax.nums, pch=21, cex=0.6, col=as.character(phyla.cols[p]), bg="white")
            mids.death <- curr.comp.ranked$f.death.boot.median[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
            lowers.death <- curr.comp.ranked$f.death.boot.CI.L[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
            uppers.death <- curr.comp.ranked$f.death.boot.CI.U[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
            arrows(x0=lowers.death/1000, y0=tax.nums, x1=uppers.death/1000, y1=tax.nums, length=0, angle=90, code=3, col=as.character(phyla.cols[p]), lwd=0.5)
            points(x=mids.death/1000, y=tax.nums, pch=21, cex=0.6, col="gray30", bg=as.character(phyla.cols[p]))
          }
          legend(x=x.max.f, y=dim(DATA)[1], legend=c("net", "gross", "death"), pch=21, col=c(as.character(phyla.cols[12]), as.character(phyla.cols[12]), "gray30"), pt.bg=c(as.character(phyla.cols[12]), "white", as.character(phyla.cols[12])), bty="o", cex=0.6, xjust=1, yjust=1)
          abline(v=0, col="black")
          par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
          axis(side=1, at=c(x.min.f-((x.max.f-x.min.f)*0.04), x.max.f+((x.max.f-x.min.f)*0.04)), labels=FALSE, tck=0, las=1, cex.axis=0.6)
          axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
          mtext(expression(paste("C flux into biomass (ng C g soil"^-1, " day"^-1, ")", sep="")), side=1, line=1.4, cex=0.75)
          par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
          # axis(side=2, at=1:dim(DATA)[1], labels=rep(NA, 103), tck=-0.015, las=1, cex.axis=0.6)
          # mtext("Taxon", side=2, line=1.35, cex=0.75)
        }
        if (log == TRUE){
          x.min.f <- min(x.mins.f[x.mins.f != -Inf & x.mins.f > x.min.f*1000], na.rm=TRUE)/1000   #reset the minimum to the second-smallest number          
          plot(y=1:dim(DATA)[1], x=log10((DATA$f.net.boot.median[tax.order$ranks]/1000)+100000), type="n", bty="n", xlab="", ylab="", xaxt="n", yaxt="n", xlim=log10(c(x.min.f+100000, x.max.f+100000)), main="")
          counter <- 1
          for (p in 1:length(phyla)){
            curr.comp.ranked <- DATA[tax.order$ranks,]
            mids <- curr.comp.ranked$f.net.boot.median[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
            lowers <- curr.comp.ranked$f.net.boot.CI.L[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
            uppers <- curr.comp.ranked$f.net.boot.CI.U[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
            tax.nums <- counter:(counter+length(mids)-1)
            counter <- counter+length(mids)
            points(x=log10((mids/1000)+100000), y=tax.nums, pch=21, cex=0.6, col=as.character(phyla.cols[p]), bg=as.character(phyla.cols[p]))
            arrows(x0=log10((lowers/1000)+100000), y0=tax.nums, x1=log10((uppers/1000)+100000), y1=tax.nums, length=0, angle=90, code=3, col=as.character(phyla.cols[p]))
            mids.gross <- curr.comp.ranked$f.gross.boot.median[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
            lowers.gross <- curr.comp.ranked$f.gross.boot.CI.L[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
            uppers.gross <- curr.comp.ranked$f.gross.boot.CI.U[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
            arrows(x0=log10((lowers.gross/1000)+100000), y0=tax.nums, x1=log10((uppers.gross/1000)+100000), y1=tax.nums, length=0, angle=90, code=3, col=as.character(phyla.cols[p]), lwd=0.5)
            points(x=log10((mids.gross/1000)+100000), y=tax.nums, pch=21, cex=0.6, col=as.character(phyla.cols[p]), bg="white")
            mids.death <- curr.comp.ranked$f.death.boot.median[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
            lowers.death <- curr.comp.ranked$f.death.boot.CI.L[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
            uppers.death <- curr.comp.ranked$f.death.boot.CI.U[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
            arrows(x0=log10((lowers.death/1000)+100000), y0=tax.nums, x1=log10((uppers.death/1000)+100000), y1=tax.nums, length=0, angle=90, code=3, col=as.character(phyla.cols[p]), lwd=0.5)
            points(x=log10((mids.death/1000)+100000), y=tax.nums, pch=21, cex=0.6, col="gray30", bg=as.character(phyla.cols[p]))
          }
          legend(x=log10(x.max.f+100000), y=dim(DATA)[1], legend=c("net", "gross", "death"), pch=21, col=c(as.character(phyla.cols[12]), as.character(phyla.cols[12]), "gray30"), pt.bg=c(as.character(phyla.cols[12]), "white", as.character(phyla.cols[12])), bty="o", cex=0.6, xjust=1, yjust=1)
          # abline(v=0, col="black")
          par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
          axis(side=1, at=c(log10(x.min.f+100000)-((log10(x.max.f+100000)-log10(x.min.f+100000))*0.04), log10(x.max.f+100000)+((log10(x.max.f+100000)-log10(x.min.f+100000))*0.04)), labels=FALSE, tck=0, las=1, cex.axis=0.6)
          AT <- c(10,20,30,40,50,60,70,80,90,100,200,300,400,500,600,700,800,900,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,20000,30000,40000,50000,60000,70000,80000,90000,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000)
          axis(side=1, at=log10(AT), labels=AT, tck=-0.015, las=1, cex.axis=0.6)
          mtext(expression(paste("C flux into biomass (ng C g soil"^-1, " day"^-1, ")  [log"[10], " + 100000]", sep="")), side=1, line=1.4, cex=0.75)
          par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
          # axis(side=2, at=1:dim(DATA)[1], labels=rep(NA, 103), tck=-0.015, las=1, cex.axis=0.6)
          # mtext("Taxon", side=2, line=1.35, cex=0.75)
        }
      }

      r.add.panel(DATA=all.rates)
      # legend(x=0.5*x.max.r, y=0.6*dim(all.rates)[1], legend=phyla, bty="n", text.col=as.character(phyla.cols), cex=0.6)
      f.add.panel(DATA=all.rates, log=FALSE)
      rm(r.add.panel, f.add.panel, taxa.id.temp)
      par(mfrow=c(1,1))

      dev.off()
  #}


#Graph the observed birth, death, and (bootstrapped) net growth and flux rates:________________________________________________________________________________
  #{
    #Growth & C flux rates - in units of r and C flux:
      graphics.off()
      # dev.new(width=7.5, height=7.5)
      pdf(file="qSIP_output/Figures/TM_Taxa_same_order_PhylumGroups_growth&flux_observed.pdf", width=7.5, height=7.5)
      par(mfcol=c(1,2))
      par(mai=c(0.44,0.42,0.2,0.05), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
      #Re-code NA's in the taxa.id data.frame as "Other" so that taxa with a phylum of 'NA' is included in the phylum-specific plot:
        taxa.id.temp <- taxa.id
        taxa.id.temp[is.na(taxa.id.temp)] <- "Other"
      phyla <- levels(data.melted$phylum)
      phyla.cols <- cols.for.phyla$col
      ranks <- order(all.rates$r.net.boot.median)
      tax.order <- data.frame(axis.loc=seq(1, dim(all.rates)[1]), ranks)
      x.mins.r <- c(all.rates$r.net.boot.median, all.rates$r.net.boot.CI.L, all.rates$b.obs, all.rates$d.obs)
      x.maxs.r <- c(all.rates$r.net.boot.median, all.rates$r.net.boot.CI.U, all.rates$b.obs, all.rates$d.obs)
      x.mins.f <- c(all.rates$f.net.boot.median, all.rates$f.net.boot.CI.L, all.rates$f.gross.obs, all.rates$f.death.obs)
      x.maxs.f <- c(all.rates$f.net.boot.median, all.rates$f.net.boot.CI.U, all.rates$f.gross.obs, all.rates$f.death.obs)
      x.min.r <- min(x.mins.r[x.mins.r != -Inf], na.rm=TRUE)
      x.max.r <- max(x.maxs.r[x.maxs.r != Inf], na.rm=TRUE)
      x.min.f <- min(x.mins.f[x.mins.f != -Inf], na.rm=TRUE)/1000
      x.max.f <- max(x.maxs.f[x.maxs.f != Inf], na.rm=TRUE)/1000

      r.add.panel <- function(DATA){
        plot(y=1:dim(DATA)[1], x=DATA$r.net.boot.median[tax.order$ranks], type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min.r, x.max.r), main="")
        counter <- 1
        for (p in 1:length(phyla)){
          curr.comp.ranked <- DATA[tax.order$ranks,]
          mids.net <- curr.comp.ranked$r.net.boot.median[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
          lowers.net <- curr.comp.ranked$r.net.boot.CI.L[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
          uppers.net <- curr.comp.ranked$r.net.boot.CI.U[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
          tax.nums <- counter:(counter+length(mids.net)-1)
          counter <- counter+length(mids.net)
          points(x=mids.net, y=tax.nums, pch=21, cex=0.6, col=as.character(phyla.cols[p]), bg=as.character(phyla.cols[p]))
          arrows(x0=lowers.net, y0=tax.nums, x1=uppers.net, y1=tax.nums, length=0, angle=90, code=3, col=as.character(phyla.cols[p]))
          par(xpd=NA)
          if(!is.element(p, c(5,6,9,10,16,17,20))){
            text(x=x.max.r, y=mean(tax.nums), labels=phyla[p], adj=c(0,0.5), pos=4, offset=-1.5, cex=0.6, col=as.character(phyla.cols)[p])
          }
          else if (is.element(p, c(5,10,20))){
            text(x=x.max.r, y=mean(tax.nums), labels=phyla[p], adj=c(0,0.5), pos=4, offset=0, cex=0.6, col=as.character(phyla.cols)[p])
          }
          else if (is.element(p, c(16))){
            text(x=x.max.r, y=mean(tax.nums), labels=phyla[p], adj=c(0,0.5), pos=4, offset=1.5, cex=0.6, col=as.character(phyla.cols)[p])
          }
          else if (is.element(p, c(6,9,17))){
            text(x=x.max.r, y=mean(tax.nums), labels=phyla[p], adj=c(0,0.5), pos=4, offset=3, cex=0.6, col=as.character(phyla.cols)[p])
          }
          par(xpd=FALSE)
          mids.b <- curr.comp.ranked$b.obs[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
          points(x=mids.b, y=tax.nums, pch=21, cex=0.6, col=as.character(phyla.cols[p]), bg="white")
          mids.d <- curr.comp.ranked$d.obs[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
          points(x=mids.d, y=tax.nums, pch=21, cex=0.6, col="gray30", bg=as.character(phyla.cols[p]))
        }
        abline(v=0, col="black")
        par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
        axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
        mtext(expression(paste(italic("r")~" (day"^-1, ")", sep="")), side=1, line=1.4, cex=0.75)
        par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
        axis(side=2, at=1:dim(DATA)[1], labels=rep(NA, 103), tck=-0.015, las=1, cex.axis=0.6)
        mtext("Taxon", side=2, line=1.35, cex=0.75)
      }

      f.add.panel <- function(DATA, log=FALSE){
        if (log == FALSE){
          plot(y=1:dim(DATA)[1], x=DATA$f.net.boot.median[tax.order$ranks]/1000, type="n", bty="n", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min.f, x.max.f), main="")
          counter <- 1
          for (p in 1:length(phyla)){
            curr.comp.ranked <- DATA[tax.order$ranks,]
            mids <- curr.comp.ranked$f.net.boot.median[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
            lowers <- curr.comp.ranked$f.net.boot.CI.L[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
            uppers <- curr.comp.ranked$f.net.boot.CI.U[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
            tax.nums <- counter:(counter+length(mids)-1)
            counter <- counter+length(mids)
            points(x=mids/1000, y=tax.nums, pch=21, cex=0.6, col=as.character(phyla.cols[p]), bg=as.character(phyla.cols[p]))
            arrows(x0=lowers/1000, y0=tax.nums, x1=uppers/1000, y1=tax.nums, length=0, angle=90, code=3, col=as.character(phyla.cols[p]))
            mids.gross <- curr.comp.ranked$f.gross.obs[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
            points(x=mids.gross/1000, y=tax.nums, pch=21, cex=0.6, col=as.character(phyla.cols[p]), bg="white")
            mids.death <- curr.comp.ranked$f.death.obs[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
            points(x=mids.death/1000, y=tax.nums, pch=21, cex=0.6, col="gray30", bg=as.character(phyla.cols[p]))
          }
          legend(x=x.max.f, y=dim(DATA)[1], legend=c("net", "gross", "death"), pch=21, col=c(as.character(phyla.cols[12]), as.character(phyla.cols[12]), "gray30"), pt.bg=c(as.character(phyla.cols[12]), "white", as.character(phyla.cols[12])), bty="o", cex=0.6, xjust=1, yjust=1)
          abline(v=0, col="black")
          par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
          axis(side=1, at=c(x.min.f-((x.max.f-x.min.f)*0.04), x.max.f+((x.max.f-x.min.f)*0.04)), labels=FALSE, tck=0, las=1, cex.axis=0.6)
          axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
          mtext(expression(paste("C flux into biomass (ng C g soil"^-1, " day"^-1, ")", sep="")), side=1, line=1.4, cex=0.75)
          par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
          # axis(side=2, at=1:dim(DATA)[1], labels=rep(NA, 103), tck=-0.015, las=1, cex.axis=0.6)
          # mtext("Taxon", side=2, line=1.35, cex=0.75)
        }
        if (log == TRUE){
          x.min.f <- min(x.mins.f[x.mins.f != -Inf & x.mins.f > x.min.f*1000], na.rm=TRUE)/1000   #reset the minimum to the second-smallest number          
          plot(y=1:dim(DATA)[1], x=log10((DATA$f.net.boot.median[tax.order$ranks]/1000)+100000), type="n", bty="n", xlab="", ylab="", xaxt="n", yaxt="n", xlim=log10(c(x.min.f+100000, x.max.f+100000)), main="")
          counter <- 1
          for (p in 1:length(phyla)){
            curr.comp.ranked <- DATA[tax.order$ranks,]
            mids <- curr.comp.ranked$f.net.boot.median[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
            lowers <- curr.comp.ranked$f.net.boot.CI.L[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
            uppers <- curr.comp.ranked$f.net.boot.CI.U[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
            tax.nums <- counter:(counter+length(mids)-1)
            counter <- counter+length(mids)
            points(x=log10((mids/1000)+100000), y=tax.nums, pch=21, cex=0.6, col=as.character(phyla.cols[p]), bg=as.character(phyla.cols[p]))
            arrows(x0=log10((lowers/1000)+100000), y0=tax.nums, x1=log10((uppers/1000)+100000), y1=tax.nums, length=0, angle=90, code=3, col=as.character(phyla.cols[p]))
            mids.gross <- curr.comp.ranked$f.gross.obs[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
            points(x=log10((mids.gross/1000)+100000), y=tax.nums, pch=21, cex=0.6, col=as.character(phyla.cols[p]), bg="white")
            mids.death <- curr.comp.ranked$f.death.obs[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
            points(x=log10((mids.death/1000)+100000), y=tax.nums, pch=21, cex=0.6, col="gray30", bg=as.character(phyla.cols[p]))
          }
          legend(x=log10(x.max.f+100000), y=dim(DATA)[1], legend=c("net", "gross", "death"), pch=21, col=c(as.character(phyla.cols[12]), as.character(phyla.cols[12]), "gray30"), pt.bg=c(as.character(phyla.cols[12]), "white", as.character(phyla.cols[12])), bty="o", cex=0.6, xjust=1, yjust=1)
          # abline(v=0, col="black")
          par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
          axis(side=1, at=c(log10(x.min.f+100000)-((log10(x.max.f+100000)-log10(x.min.f+100000))*0.04), log10(x.max.f+100000)+((log10(x.max.f+100000)-log10(x.min.f+100000))*0.04)), labels=FALSE, tck=0, las=1, cex.axis=0.6)
          AT <- c(10,20,30,40,50,60,70,80,90,100,200,300,400,500,600,700,800,900,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,20000,30000,40000,50000,60000,70000,80000,90000,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000)
          axis(side=1, at=log10(AT), labels=AT, tck=-0.015, las=1, cex.axis=0.6)
          mtext(expression(paste("C flux into biomass (ng C g soil"^-1, " day"^-1, ")  [log"[10], " + 100000]", sep="")), side=1, line=1.4, cex=0.75)
          par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
          # axis(side=2, at=1:dim(DATA)[1], labels=rep(NA, 103), tck=-0.015, las=1, cex.axis=0.6)
          # mtext("Taxon", side=2, line=1.35, cex=0.75)
        }
      }

      r.add.panel(DATA=all.rates)
      # legend(x=0.5*x.max.r, y=0.6*dim(all.rates)[1], legend=phyla, bty="n", text.col=as.character(phyla.cols), cex=0.6)
      f.add.panel(DATA=all.rates, log=FALSE)
      rm(r.add.panel, f.add.panel, taxa.id.temp)
      par(mfrow=c(1,1))

      dev.off()
  #}


#Graph the bootstrapped birth, death, and net growth rates (no flux panel):____________________________________________________________________________________
  #{
    #Growth rates - in units of r:
      graphics.off()
      # dev.new(width=7.5, height=7.5)
      pdf(file="qSIP_output/Figures/TM_Taxa_PhylumGroups_growth.pdf", width=7.5, height=7.5)
      par(mfcol=c(1,2))
      par(mai=c(0.44,0.42,0.2,0.05), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
      #Re-code NA's in the taxa.id data.frame as "Other" so that taxa with a phylum of 'NA' is included in the phylum-specific plot:
        taxa.id.temp <- taxa.id
        taxa.id.temp[is.na(taxa.id.temp)] <- "Other"
      phyla <- levels(data.melted$phylum)
      phyla.cols <- cols.for.phyla$col
      ranks <- order(all.rates$r.net.boot.median)
      tax.order <- data.frame(axis.loc=seq(1, dim(all.rates)[1]), ranks)
      x.mins.r <- c(all.rates$r.net.boot.median, all.rates$r.net.boot.CI.L, all.rates$b.boot.CI.L, all.rates$d.boot.CI.L)
      x.maxs.r <- c(all.rates$r.net.boot.median, all.rates$r.net.boot.CI.U, all.rates$b.boot.CI.U, all.rates$d.boot.CI.U)
      x.min.r <- min(x.mins.r[x.mins.r != -Inf], na.rm=TRUE)
      x.max.r <- max(x.maxs.r[x.maxs.r != Inf], na.rm=TRUE)

      r.add.panel <- function(DATA){
        plot(y=1:dim(DATA)[1], x=DATA$r.net.boot.median[tax.order$ranks], type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min.r, x.max.r), main="")
        counter <- 1
        for (p in 1:length(phyla)){
          curr.comp.ranked <- DATA[tax.order$ranks,]
          mids.net <- curr.comp.ranked$r.net.boot.median[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
          lowers.net <- curr.comp.ranked$r.net.boot.CI.L[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
          uppers.net <- curr.comp.ranked$r.net.boot.CI.U[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
          tax.nums <- counter:(counter+length(mids.net)-1)
          counter <- counter+length(mids.net)
          points(x=mids.net, y=tax.nums, pch=21, cex=0.6, col=as.character(phyla.cols[p]), bg=as.character(phyla.cols[p]))
          arrows(x0=lowers.net, y0=tax.nums, x1=uppers.net, y1=tax.nums, length=0, angle=90, code=3, col=as.character(phyla.cols[p]))
          par(xpd=NA)
          if(!is.element(p, c(5,6,9,10,16,17,20))){
            text(x=x.max.r, y=mean(tax.nums), labels=phyla[p], adj=c(0,0.5), pos=4, offset=-1.5, cex=0.6, col=as.character(phyla.cols)[p])
          }
          else if (is.element(p, c(5,10,20))){
            text(x=x.max.r, y=mean(tax.nums), labels=phyla[p], adj=c(0,0.5), pos=4, offset=0, cex=0.6, col=as.character(phyla.cols)[p])
          }
          else if (is.element(p, c(16))){
            text(x=x.max.r, y=mean(tax.nums), labels=phyla[p], adj=c(0,0.5), pos=4, offset=1.5, cex=0.6, col=as.character(phyla.cols)[p])
          }
          else if (is.element(p, c(6,9,17))){
            text(x=x.max.r, y=mean(tax.nums), labels=phyla[p], adj=c(0,0.5), pos=4, offset=3, cex=0.6, col=as.character(phyla.cols)[p])
          }
          par(xpd=FALSE)
          mids.b <- curr.comp.ranked$b.boot.median[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
          lowers.b <- curr.comp.ranked$b.boot.CI.L[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
          uppers.b <- curr.comp.ranked$b.boot.CI.U[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
          arrows(x0=lowers.b, y0=tax.nums, x1=uppers.b, y1=tax.nums, length=0, angle=90, code=3, col=as.character(phyla.cols[p]), lwd=0.5)
          points(x=mids.b, y=tax.nums, pch=21, cex=0.6, col=as.character(phyla.cols[p]), bg="white")
          mids.d <- curr.comp.ranked$d.boot.median[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
          lowers.d <- curr.comp.ranked$d.boot.CI.L[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
          uppers.d <- curr.comp.ranked$d.boot.CI.U[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
          arrows(x0=lowers.d, y0=tax.nums, x1=uppers.d, y1=tax.nums, length=0, angle=90, code=3, col=as.character(phyla.cols[p]), lwd=0.5)
          points(x=mids.d, y=tax.nums, pch=21, cex=0.6, col="gray30", bg=as.character(phyla.cols[p]))
        }
        legend(x=1.04*x.min.r, y=1.04*(dim(DATA)[1]-1), legend=c("net", "gross", "death"), pch=21, col=c(as.character(phyla.cols[12]), as.character(phyla.cols[12]), "gray30"), pt.bg=c(as.character(phyla.cols[12]), "white", as.character(phyla.cols[12])), bty="o", cex=0.6, xjust=0, yjust=1, xpd=NA)
        abline(v=0, col="black")
        par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
        axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
        mtext(expression(paste(italic("r")~" (day"^-1, ")", sep="")), side=1, line=1.4, cex=0.75)
        par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
        axis(side=2, at=1:dim(DATA)[1], labels=rep(NA, 103), tck=-0.015, las=1, cex.axis=0.6)
        mtext("Taxon", side=2, line=1.35, cex=0.75)
      }

      r.add.panel(DATA=all.rates)
      rm(r.add.panel, taxa.id.temp)
      par(mfrow=c(1,1))

      dev.off()
  #}


#Graph the observed birth, death, and (bootstrapped) net growth rates (no flux panel):_________________________________________________________________________
  #{
    #Growth rates - in units of r:
      graphics.off()
      # dev.new(width=7.5, height=7.5)
      pdf(file="qSIP_output/Figures/TM_Taxa_PhylumGroups_growth_observed.pdf", width=7.5, height=7.5)
      par(mfcol=c(1,2))
      par(mai=c(0.44,0.42,0.2,0.05), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
      #Re-code NA's in the taxa.id data.frame as "Other" so that taxa with a phylum of 'NA' is included in the phylum-specific plot:
        taxa.id.temp <- taxa.id
        taxa.id.temp[is.na(taxa.id.temp)] <- "Other"
      phyla <- levels(data.melted$phylum)
      phyla.cols <- cols.for.phyla$col
      ranks <- order(all.rates$r.net.boot.median)
      tax.order <- data.frame(axis.loc=seq(1, dim(all.rates)[1]), ranks)
      x.mins.r <- c(all.rates$r.net.boot.median, all.rates$r.net.boot.CI.L, all.rates$b.obs, all.rates$d.obs)
      x.maxs.r <- c(all.rates$r.net.boot.median, all.rates$r.net.boot.CI.U, all.rates$b.obs, all.rates$d.obs)
      x.min.r <- min(x.mins.r[x.mins.r != -Inf], na.rm=TRUE)
      x.max.r <- max(x.maxs.r[x.maxs.r != Inf], na.rm=TRUE)

      r.add.panel <- function(DATA){
        plot(y=1:dim(DATA)[1], x=DATA$r.net.boot.median[tax.order$ranks], type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min.r, x.max.r), main="")
        counter <- 1
        for (p in 1:length(phyla)){
          curr.comp.ranked <- DATA[tax.order$ranks,]
          mids.net <- curr.comp.ranked$r.net.boot.median[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
          lowers.net <- curr.comp.ranked$r.net.boot.CI.L[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
          uppers.net <- curr.comp.ranked$r.net.boot.CI.U[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
          tax.nums <- counter:(counter+length(mids.net)-1)
          counter <- counter+length(mids.net)
          points(x=mids.net, y=tax.nums, pch=21, cex=0.6, col=as.character(phyla.cols[p]), bg=as.character(phyla.cols[p]))
          arrows(x0=lowers.net, y0=tax.nums, x1=uppers.net, y1=tax.nums, length=0, angle=90, code=3, col=as.character(phyla.cols[p]))
          par(xpd=NA)
          if(!is.element(p, c(5,6,9,10,16,17,20))){
            text(x=x.max.r, y=mean(tax.nums), labels=phyla[p], adj=c(0,0.5), pos=4, offset=-1.5, cex=0.6, col=as.character(phyla.cols)[p])
          }
          else if (is.element(p, c(5,10,20))){
            text(x=x.max.r, y=mean(tax.nums), labels=phyla[p], adj=c(0,0.5), pos=4, offset=0, cex=0.6, col=as.character(phyla.cols)[p])
          }
          else if (is.element(p, c(16))){
            text(x=x.max.r, y=mean(tax.nums), labels=phyla[p], adj=c(0,0.5), pos=4, offset=1.5, cex=0.6, col=as.character(phyla.cols)[p])
          }
          else if (is.element(p, c(6,9,17))){
            text(x=x.max.r, y=mean(tax.nums), labels=phyla[p], adj=c(0,0.5), pos=4, offset=3, cex=0.6, col=as.character(phyla.cols)[p])
          }
          par(xpd=FALSE)
          mids.b <- curr.comp.ranked$b.obs[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
          points(x=mids.b, y=tax.nums, pch=21, cex=0.6, col=as.character(phyla.cols[p]), bg="white")
          mids.d <- curr.comp.ranked$d.obs[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
          points(x=mids.d, y=tax.nums, pch=21, cex=0.6, col="gray30", bg=as.character(phyla.cols[p]))
        }
        legend(x=1.04*x.min.r, y=1.04*(dim(DATA)[1]-1), legend=c("net", "gross", "death"), pch=21, col=c(as.character(phyla.cols[12]), as.character(phyla.cols[12]), "gray30"), pt.bg=c(as.character(phyla.cols[12]), "white", as.character(phyla.cols[12])), bty="o", cex=0.6, xjust=0, yjust=1, xpd=NA)
        abline(v=0, col="black")
        par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
        axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
        mtext(expression(paste(italic("r")~" (day"^-1, ")", sep="")), side=1, line=1.4, cex=0.75)
        par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
        axis(side=2, at=1:dim(DATA)[1], labels=rep(NA, 103), tck=-0.015, las=1, cex.axis=0.6)
        mtext("Taxon", side=2, line=1.35, cex=0.75)
      }

      r.add.panel(DATA=all.rates)
      rm(r.add.panel, taxa.id.temp)
      par(mfrow=c(1,1))

      dev.off()
  #}


#Graph the bootstrapped birth, death, and net growth rates only (horizontal layout):___________________________________________________________________________
  #{
    #Growth rates - in units of r:
      graphics.off()
      # dev.new(width=7.5, height=4)
      pdf(file="qSIP_output/Figures/TM_Taxa_PhylumGroups_growth_horizontal.pdf", width=7.5, height=4)
      par(mai=c(0.29,0.47,0.60,0.03), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
      #Re-code NA's in the taxa.id data.frame as "Other" so that taxa with a phylum of 'NA' is included in the phylum-specific plot:
        taxa.id.temp <- taxa.id
        taxa.id.temp[is.na(taxa.id.temp)] <- "Other"
      phyla <- levels(data.melted$phylum)
      phyla.cols <- cols.for.phyla$col
      ranks <- order(all.rates$r.net.boot.median, decreasing=TRUE)
      tax.order <- data.frame(axis.loc=seq(1, dim(all.rates)[1]), ranks)
      y.mins.r <- c(all.rates$r.net.boot.median, all.rates$r.net.boot.CI.L, all.rates$b.boot.CI.L, all.rates$d.boot.CI.L)
      y.maxs.r <- c(all.rates$r.net.boot.median, all.rates$r.net.boot.CI.U, all.rates$b.boot.CI.U, all.rates$d.boot.CI.U)
      y.min.r <- min(y.mins.r[y.mins.r != -Inf], na.rm=TRUE)
      y.max.r <- max(y.maxs.r[y.maxs.r != Inf], na.rm=TRUE)

      r.add.panel <- function(DATA){
        plot(y=DATA$r.net.boot.median[tax.order$ranks], x=1:dim(DATA)[1], type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", ylim=c(y.min.r, y.max.r), main="")
        counter <- 1
        for (p in 1:length(phyla)){
          curr.comp.ranked <- DATA[tax.order$ranks,]
          mids.net <- curr.comp.ranked$r.net.boot.median[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
          lowers.net <- curr.comp.ranked$r.net.boot.CI.L[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
          uppers.net <- curr.comp.ranked$r.net.boot.CI.U[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
          tax.nums <- counter:(counter+length(mids.net)-1)
          counter <- counter+length(mids.net)
          points(x=tax.nums, y=mids.net, pch=21, cex=0.6, col=as.character(phyla.cols[p]), bg=as.character(phyla.cols[p]))
          arrows(x0=tax.nums, y0=lowers.net, x1=tax.nums, y1=uppers.net, length=0, angle=90, code=3, col=as.character(phyla.cols[p]))
          par(xpd=NA)
          if(!is.element(p, c(3,5,6,8,9,10,11,13,14,16,17,18,21,22))){
            text(x=mean(tax.nums), y=y.max.r, labels=phyla[p], adj=c(0,0.5), pos=3, offset=1.5, cex=0.6, col=as.character(phyla.cols)[p])
          }
          else if (is.element(p, c(3,5,8,13,21))){
            text(x=mean(tax.nums), y=y.max.r, labels=phyla[p], adj=c(0,0.5), pos=3, offset=2.0, cex=0.6, col=as.character(phyla.cols)[p])
          }
          else if (is.element(p, c(6,9,14,22))){
            text(x=mean(tax.nums), y=y.max.r, labels=phyla[p], adj=c(0,0.5), pos=3, offset=1.0, cex=0.6, col=as.character(phyla.cols)[p])
          }
          else if (is.element(p, c(10,16))){
            text(x=mean(tax.nums), y=y.max.r, labels=phyla[p], adj=c(0,0.5), pos=3, offset=2.5, cex=0.6, col=as.character(phyla.cols)[p])
          }
          else if (is.element(p, c(11,17))){
            text(x=mean(tax.nums), y=y.max.r, labels=phyla[p], adj=c(0,0.5), pos=3, offset=0.5, cex=0.6, col=as.character(phyla.cols)[p])
          }
          else if (is.element(p, c(18))){
            text(x=mean(tax.nums), y=y.max.r, labels=phyla[p], adj=c(0,0.5), pos=3, offset=3.0, cex=0.6, col=as.character(phyla.cols)[p])
          }
          par(xpd=FALSE)
          mids.b <- curr.comp.ranked$b.boot.median[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
          lowers.b <- curr.comp.ranked$b.boot.CI.L[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
          uppers.b <- curr.comp.ranked$b.boot.CI.U[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
          arrows(x0=tax.nums, y0=lowers.b, x1=tax.nums, y1=uppers.b, length=0, angle=90, code=3, col=as.character(phyla.cols[p]), lwd=0.5)
          points(x=tax.nums, y=mids.b, pch=21, cex=0.6, col=as.character(phyla.cols[p]), bg="white")
          mids.d <- curr.comp.ranked$d.boot.median[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
          lowers.d <- curr.comp.ranked$d.boot.CI.L[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
          uppers.d <- curr.comp.ranked$d.boot.CI.U[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
          arrows(x0=tax.nums, y0=lowers.d, x1=tax.nums, y1=uppers.d, length=0, angle=90, code=3, col=as.character(phyla.cols[p]), lwd=0.5)
          points(x=tax.nums, y=mids.d, pch=21, cex=0.6, col="gray30", bg=as.character(phyla.cols[p]))
        }
        legend(x=1-(0.02*(dim(DATA)[1]-1)), y=y.min.r, legend=c("gross", "net", "death"), pch=21, col=c(as.character(phyla.cols[2]), as.character(phyla.cols[2]), "gray30"), pt.bg=c("white", as.character(phyla.cols[2]), as.character(phyla.cols[2])), bty="o", cex=0.6, xjust=0, yjust=0, xpd=NA)
        abline(h=0, col="black")
        par(mgp=c(3,0,0))			#set spacing for the x-axis labels (the second value)
        axis(side=1, at=1:dim(DATA)[1], labels=rep(NA, 103), tck=-0.015, las=1, cex.axis=0.6)
        mtext("Taxon", side=1, line=0.4, cex=0.75)
        par(mgp=c(3,0.45,0))				  #set spacing for the y-axis labels (the second value)
        axis(side=2, tck=-0.015, las=1, cex.axis=0.6)
        mtext(expression(paste(italic("r")~" (day"^-1, ")", sep="")), side=2, line=1.7, at=0, cex=0.75)
      }

      r.add.panel(DATA=all.rates)
      rm(r.add.panel, taxa.id.temp)
      par(mfrow=c(1,1))

      dev.off()
  #}


#Graph the bootstrapped birth, death, and net growth rates only (horizontal layout - formatted for powerpoint):________________________________________________
  #{
    #Growth rates - in units of r:
      graphics.off()
      # dev.new(width=10, height=6)
      pdf(file="qSIP_output/Figures/TM_Taxa_PhylumGroups_growth_horizontal_slide.pdf", width=10, height=6)
      par(mai=c(0.49,0.80,1.03,0.03), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
      #Re-code NA's in the taxa.id data.frame as "Other" so that taxa with a phylum of 'NA' is included in the phylum-specific plot:
        taxa.id.temp <- taxa.id
        taxa.id.temp[is.na(taxa.id.temp)] <- "Other"
      phyla <- levels(data.melted$phylum)
      phyla.cols <- cols.for.phyla$col
      ranks <- order(all.rates$r.net.boot.median, decreasing=TRUE)
      tax.order <- data.frame(axis.loc=seq(1, dim(all.rates)[1]), ranks)
      y.mins.r <- c(all.rates$r.net.boot.median, all.rates$r.net.boot.CI.L, all.rates$b.boot.CI.L, all.rates$d.boot.CI.L)
      y.maxs.r <- c(all.rates$r.net.boot.median, all.rates$r.net.boot.CI.U, all.rates$b.boot.CI.U, all.rates$d.boot.CI.U)
      y.min.r <- min(y.mins.r[y.mins.r != -Inf], na.rm=TRUE)
      y.max.r <- max(y.maxs.r[y.maxs.r != Inf], na.rm=TRUE)

      r.add.panel <- function(DATA){
        plot(y=DATA$r.net.boot.median[tax.order$ranks], x=1:dim(DATA)[1], type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", ylim=c(y.min.r, y.max.r), main="")
        counter <- 1
        for (p in 1:length(phyla)){
          curr.comp.ranked <- DATA[tax.order$ranks,]
          mids.net <- curr.comp.ranked$r.net.boot.median[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
          lowers.net <- curr.comp.ranked$r.net.boot.CI.L[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
          uppers.net <- curr.comp.ranked$r.net.boot.CI.U[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
          tax.nums <- counter:(counter+length(mids.net)-1)
          counter <- counter+length(mids.net)
          par(xpd=NA)
          if(!is.element(p, c(2,3,5,6,7,8,10,11,12,13,14,16,17,18,21,22))){
            text(x=mean(tax.nums), y=y.max.r, labels=phyla[p], adj=c(0,0.5), pos=3, offset=2.10, cex=1.1, col=as.character(phyla.cols)[p])
          }
          else if (is.element(p, c(3,7,13,22))){
            text(x=mean(tax.nums), y=y.max.r, labels=phyla[p], adj=c(0,0.5), pos=3, offset=2.90, cex=1.1, col=as.character(phyla.cols)[p])
          }
          else if (is.element(p, c(2,5,11,14))){
            text(x=mean(tax.nums), y=y.max.r, labels=phyla[p], adj=c(0,0.5), pos=3, offset=1.30, cex=1.1, col=as.character(phyla.cols)[p])
          }
          else if (is.element(p, c(21))){
            text(x=mean(tax.nums)-5, y=y.max.r, labels=phyla[p], adj=c(0,0.5), pos=3, offset=1.30, cex=1.1, col=as.character(phyla.cols)[p])
          }
          else if (is.element(p, c(8,16))){
            text(x=mean(tax.nums), y=y.max.r, labels=phyla[p], adj=c(0,0.5), pos=3, offset=3.70, cex=1.1, col=as.character(phyla.cols)[p])
          }
          else if (is.element(p, c(6,12,17))){
            text(x=mean(tax.nums), y=y.max.r, labels=phyla[p], adj=c(0,0.5), pos=3, offset=0.50, cex=1.1, col=as.character(phyla.cols)[p])
          }
          else if (is.element(p, c(10,18))){
            text(x=mean(tax.nums), y=y.max.r, labels=phyla[p], adj=c(0,0.5), pos=3, offset=4.50, cex=1.1, col=as.character(phyla.cols)[p])
          }
          par(xpd=FALSE)
          mids.b <- curr.comp.ranked$b.boot.median[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
          lowers.b <- curr.comp.ranked$b.boot.CI.L[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
          uppers.b <- curr.comp.ranked$b.boot.CI.U[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
          mids.d <- curr.comp.ranked$d.boot.median[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
          lowers.d <- curr.comp.ranked$d.boot.CI.L[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]
          uppers.d <- curr.comp.ranked$d.boot.CI.U[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id.temp$taxon[taxa.id.temp$phylum == phyla[p]])]

          arrows(x0=tax.nums, y0=lowers.b, x1=tax.nums, y1=uppers.b, length=0, angle=90, code=3, col=as.character(phyla.cols[p]), lwd=1.5)
          arrows(x0=tax.nums, y0=lowers.d, x1=tax.nums, y1=uppers.d, length=0, angle=90, code=3, col=as.character(phyla.cols[p]), lwd=1.5)
          arrows(x0=tax.nums, y0=lowers.net, x1=tax.nums, y1=uppers.net, length=0, angle=90, code=3, col=as.character(phyla.cols[p]), lwd=1.5)
          points(x=tax.nums, y=mids.b, pch=21, cex=1.5, col=as.character(phyla.cols[p]), bg="white")
          points(x=tax.nums, y=mids.net, pch=21, cex=1.5, col=as.character(phyla.cols[p]), bg=as.character(phyla.cols[p]))
          points(x=tax.nums, y=mids.d, pch=21, cex=1.5, col="gray30", bg=as.character(phyla.cols[p]))
        }
        legend(x=dim(DATA)[1]+(0.02*(dim(DATA)[1]-1)), y=y.max.r+(0.02*(y.max.r-y.min.r)), legend=c("reproduction", "net growth", "death"), pch=21, col=c(as.character(phyla.cols[19]), as.character(phyla.cols[19]), "gray30"), pt.bg=c("white", as.character(phyla.cols[19]), as.character(phyla.cols[19])), bty="o", pt.cex=1.5, cex=1.0, xjust=1, yjust=1, xpd=NA)
        abline(h=0, col="black", lty=1, lwd=1.5)
        par(mgp=c(3,0,0))			#set spacing for the x-axis labels (the second value)
        axis(side=1, at=1:dim(DATA)[1], labels=rep(NA, 103), tck=-0.010, las=1, cex.axis=1.0)
        mtext("Taxon", side=1, line=1.4, cex=1.6)
        par(mgp=c(3,0.45,0))				  #set spacing for the y-axis labels (the second value)
        axis(side=2, tck=-0.010, las=1, cex.axis=1.0)
        mtext(expression(paste(italic("r")~" (day"^-1, ")", sep="")), side=2, line=2.6, at=0, cex=1.6)
      }

      r.add.panel(DATA=all.rates)
      rm(r.add.panel, taxa.id.temp)
      par(mfrow=c(1,1))

      dev.off()
  #}


#Graph the bootstrapped birth, death, and net growth rates only, without phyla (horizontal layout):____________________________________________________________
  #{
    #Growth rates - in units of r:
      graphics.off()
      # dev.new(width=7.5, height=3.5)
      pdf(file="qSIP_output/Figures/TM_Taxa_growth_horizontal.pdf", width=7.5, height=3.5)
      par(mai=c(0.29,0.47,0.10,0.03), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
      phyla.cols <- cols.for.phyla$col
      ranks <- order(all.rates$r.net.boot.median, decreasing=TRUE)
      tax.order <- data.frame(axis.loc=seq(1, dim(all.rates)[1]), ranks)
      y.mins.r <- c(all.rates$r.net.boot.median, all.rates$r.net.boot.CI.L, all.rates$b.boot.CI.L, all.rates$d.boot.CI.L)
      y.maxs.r <- c(all.rates$r.net.boot.median, all.rates$r.net.boot.CI.U, all.rates$b.boot.CI.U, all.rates$d.boot.CI.U)
      y.min.r <- min(y.mins.r[y.mins.r != -Inf], na.rm=TRUE)
      y.max.r <- max(y.maxs.r[y.maxs.r != Inf], na.rm=TRUE)

      r.add.panel <- function(DATA){
        plot(y=DATA$r.net.boot.median[tax.order$ranks], x=1:dim(DATA)[1], type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", ylim=c(y.min.r, y.max.r), main="")
        curr.comp.ranked <- DATA[tax.order$ranks,]
        mids.net <- curr.comp.ranked$r.net.boot.median
        lowers.net <- curr.comp.ranked$r.net.boot.CI.L
        uppers.net <- curr.comp.ranked$r.net.boot.CI.U
        tax.nums <- 1:length(mids.net)
        points(x=tax.nums, y=mids.net, pch=21, cex=0.6, col=phyla.cols[7], bg=phyla.cols[7])
        arrows(x0=tax.nums, y0=lowers.net, x1=tax.nums, y1=uppers.net, length=0, angle=90, code=3, col=phyla.cols[7])
        mids.b <- curr.comp.ranked$b.boot.median
        lowers.b <- curr.comp.ranked$b.boot.CI.L
        uppers.b <- curr.comp.ranked$b.boot.CI.U
        arrows(x0=tax.nums, y0=lowers.b, x1=tax.nums, y1=uppers.b, length=0, angle=90, code=3, col="black", lwd=0.5)
        points(x=tax.nums, y=mids.b, pch=21, cex=0.6, col="black", bg="white")
        mids.d <- curr.comp.ranked$d.boot.median
        lowers.d <- curr.comp.ranked$d.boot.CI.L
        uppers.d <- curr.comp.ranked$d.boot.CI.U
        arrows(x0=tax.nums, y0=lowers.d, x1=tax.nums, y1=uppers.d, length=0, angle=90, code=3, col="black", lwd=0.5)
        points(x=tax.nums, y=mids.d, pch=21, cex=0.6, col="black", bg="black")
        legend(x=dim(DATA)[1]+(0.02*(dim(DATA)[1]-1)), y=y.max.r+(0.02*(y.max.r-y.min.r)), legend=c("gross", "net", "death"), pch=21, col=c("black", phyla.cols[7], "black"), pt.bg=c("white", phyla.cols[7], "black"), bty="o", cex=0.6, xjust=1, yjust=1, xpd=NA)
        abline(h=0, col="red", lty=2)
        par(mgp=c(3,0,0))			#set spacing for the x-axis labels (the second value)
        axis(side=1, at=1:dim(DATA)[1], labels=rep(NA, 103), tck=-0.015, las=1, cex.axis=0.6)
        mtext("Ranked genus", side=1, line=0.4, cex=0.75)
        par(mgp=c(3,0.45,0))				  #set spacing for the y-axis labels (the second value)
        axis(side=2, tck=-0.015, las=1, cex.axis=0.6)
        mtext(expression(paste(italic("r")~" (day"^-1, ")", sep="")), side=2, line=1.7, at=0, cex=0.75)
      }

      r.add.panel(DATA=all.rates)
      rm(r.add.panel)
      par(mfrow=c(1,1))

      dev.off()
  #}


#Graph the bootstrapped birth, death, and net growth rates only, without phyla (horizontal layout, publication width):_________________________________________
  #{
    #Growth rates - in units of r:
      graphics.off()
      # dev.new(width=6.0, height=3.5)
      pdf(file="qSIP_output/Figures/TM_Taxa_growth_horizontal_pub_width.pdf", width=6.0, height=3.5)
      par(mai=c(0.29,0.47,0.10,0.03), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
      phyla.cols <- cols.for.phyla$col
      ranks <- order(all.rates$r.net.boot.median, decreasing=TRUE)
      tax.order <- data.frame(axis.loc=seq(1, dim(all.rates)[1]), ranks)
      y.mins.r <- c(all.rates$r.net.boot.median, all.rates$r.net.boot.CI.L, all.rates$b.boot.CI.L, all.rates$d.boot.CI.L)
      y.maxs.r <- c(all.rates$r.net.boot.median, all.rates$r.net.boot.CI.U, all.rates$b.boot.CI.U, all.rates$d.boot.CI.U)
      y.min.r <- min(y.mins.r[y.mins.r != -Inf], na.rm=TRUE)
      y.max.r <- max(y.maxs.r[y.maxs.r != Inf], na.rm=TRUE)

      r.add.panel <- function(DATA){
        plot(y=DATA$r.net.boot.median[tax.order$ranks], x=1:dim(DATA)[1], type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", ylim=c(y.min.r, y.max.r), main="")
        curr.comp.ranked <- DATA[tax.order$ranks,]
        mids.net <- curr.comp.ranked$r.net.boot.median
        lowers.net <- curr.comp.ranked$r.net.boot.CI.L
        uppers.net <- curr.comp.ranked$r.net.boot.CI.U
        tax.nums <- 1:length(mids.net)
        points(x=tax.nums, y=mids.net, pch=21, cex=0.6, col=phyla.cols[7], bg=phyla.cols[7])
        arrows(x0=tax.nums, y0=lowers.net, x1=tax.nums, y1=uppers.net, length=0, angle=90, code=3, col=phyla.cols[7])
        mids.b <- curr.comp.ranked$b.boot.median
        lowers.b <- curr.comp.ranked$b.boot.CI.L
        uppers.b <- curr.comp.ranked$b.boot.CI.U
        arrows(x0=tax.nums, y0=lowers.b, x1=tax.nums, y1=uppers.b, length=0, angle=90, code=3, col="black", lwd=0.5)
        points(x=tax.nums, y=mids.b, pch=21, cex=0.6, col="black", bg="white")
        mids.d <- curr.comp.ranked$d.boot.median
        lowers.d <- curr.comp.ranked$d.boot.CI.L
        uppers.d <- curr.comp.ranked$d.boot.CI.U
        arrows(x0=tax.nums, y0=lowers.d, x1=tax.nums, y1=uppers.d, length=0, angle=90, code=3, col="black", lwd=0.5)
        points(x=tax.nums, y=mids.d, pch=21, cex=0.6, col="black", bg="black")
        legend(x=dim(DATA)[1]+(0.02*(dim(DATA)[1]-1)), y=y.max.r+(0.02*(y.max.r-y.min.r)), legend=c("gross", "net", "death"), pch=21, col=c("black", phyla.cols[7], "black"), pt.bg=c("white", phyla.cols[7], "black"), bty="o", cex=0.6, xjust=1, yjust=1, xpd=NA)
        abline(h=0, col="red", lty=2)
        par(mgp=c(3,0,0))			#set spacing for the x-axis labels (the second value)
        axis(side=1, at=1:dim(DATA)[1], labels=rep(NA, 103), tck=-0.015, las=1, cex.axis=0.6)
        mtext("Ranked genus", side=1, line=0.4, cex=0.75)
        par(mgp=c(3,0.45,0))				  #set spacing for the y-axis labels (the second value)
        axis(side=2, tck=-0.015, las=1, cex.axis=0.6)
        mtext(expression(paste(italic("r")~" (day"^-1, ")", sep="")), side=2, line=1.7, at=0, cex=0.75)
      }

      r.add.panel(DATA=all.rates)
      rm(r.add.panel)
      par(mfrow=c(1,1))

      dev.off()
  #}


#Graph the bootstrapped abundance at Time t, without phyla (horizontal layout):________________________________________________________________________________
  #{
    #Growth rates - in units of r:
      graphics.off()
      # dev.new(width=7.5, height=3.5)
      pdf(file="qSIP_output/Figures/TM_Taxa_abundance_horizontal.pdf", width=7.5, height=3.5)
      par(mai=c(0.29,0.47,0.10,0.03), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
      phyla.cols <- cols.for.phyla$col
      ranks <- order(all.rates$N.tot.Tt.boot.median, decreasing=TRUE)
      tax.order <- data.frame(axis.loc=seq(1, dim(all.rates)[1]), ranks)
      y.mins.r <- c(all.rates$N.tot.Tt.boot.median, all.rates$N.tot.Tt.boot.CI.L)/(10^6)
      y.maxs.r <- c(all.rates$N.tot.Tt.boot.median, all.rates$N.tot.Tt.boot.CI.U)/(10^6)
      y.min.r <- min(y.mins.r[y.mins.r != -Inf], na.rm=TRUE)
      y.max.r <- max(y.maxs.r[y.maxs.r != Inf], na.rm=TRUE)

      N.add.panel <- function(DATA){
        plot(y=DATA$N.tot.Tt.boot.median[tax.order$ranks]/(10^6), x=1:dim(DATA)[1], type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", ylim=c(y.min.r, y.max.r), main="")
        curr.comp.ranked <- DATA[tax.order$ranks,]
        mids <- curr.comp.ranked$N.tot.Tt.boot.median/(10^6)
        lowers <- curr.comp.ranked$N.tot.Tt.boot.CI.L/(10^6)
        uppers <- curr.comp.ranked$N.tot.Tt.boot.CI.U/(10^6)
        tax.nums <- 1:length(mids)
        arrows(x0=tax.nums, y0=lowers, x1=tax.nums, y1=uppers, length=0, angle=90, code=3, col=as.character(phyla.cols[2]))
        points(x=tax.nums, y=mids, pch=21, cex=0.6, col=as.character(phyla.cols[2]), bg=as.character(phyla.cols[2]))
        par(mgp=c(3,0,0))			        #set spacing for the x-axis labels (the second value)
        axis(side=1, at=1:dim(DATA)[1], labels=rep(NA, 103), tck=-0.015, las=1, cex.axis=0.6)
        mtext("Ranked genus", side=1, line=0.4, cex=0.75)
        par(mgp=c(3,0.45,0))				  #set spacing for the y-axis labels (the second value)
        axis(side=2, tck=-0.015, las=1, cex.axis=0.6)
        mtext(expression(paste("Abundance (10"^6, " 16S copies g soil"^-1, ")", sep="")), side=2, line=1.7, at=NA, cex=0.75)
      }

      N.add.panel(DATA=all.rates)
      rm(N.add.panel)
      par(mfrow=c(1,1))

      dev.off()
  #}


#Graph the bootstrapped abundance at Time t, without phyla (horizontal layout - formatted for powerpoint):_____________________________________________________
  #{
    #Growth rates - in units of r:
      graphics.off()
      # dev.new(width=10, height=5)
      pdf(file="qSIP_output/Figures/TM_Taxa_abundance_horizontal_slide.pdf", width=10, height=5)
      par(mai=c(0.49,0.80,0.06,0.03), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
      phyla.cols <- cols.for.phyla$col
      ranks <- order(all.rates$N.tot.Tt.boot.median, decreasing=TRUE)
      tax.order <- data.frame(axis.loc=seq(1, dim(all.rates)[1]), ranks)
      y.mins.r <- c(all.rates$N.tot.Tt.boot.median, all.rates$N.tot.Tt.boot.CI.L)/(10^6)
      y.maxs.r <- c(all.rates$N.tot.Tt.boot.median, all.rates$N.tot.Tt.boot.CI.U)/(10^6)
      y.min.r <- min(y.mins.r[y.mins.r != -Inf], na.rm=TRUE)
      y.max.r <- max(y.maxs.r[y.maxs.r != Inf], na.rm=TRUE)

      N.add.panel <- function(DATA){
        plot(y=DATA$N.tot.Tt.boot.median[tax.order$ranks]/(10^6), x=1:dim(DATA)[1], type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", ylim=c(y.min.r, y.max.r), main="")
        curr.comp.ranked <- DATA[tax.order$ranks,]
        mids <- curr.comp.ranked$N.tot.Tt.boot.median/(10^6)
        lowers <- curr.comp.ranked$N.tot.Tt.boot.CI.L/(10^6)
        uppers <- curr.comp.ranked$N.tot.Tt.boot.CI.U/(10^6)
        tax.nums <- 1:length(mids)
        arrows(x0=tax.nums, y0=lowers, x1=tax.nums, y1=uppers, length=0, angle=90, code=3, col=as.character(phyla.cols[2]), lwd=1.5)
        points(x=tax.nums, y=mids, pch=21, cex=1.5, col=as.character(phyla.cols[2]), bg=as.character(phyla.cols[2]))
        par(mgp=c(3,0,0))			        #set spacing for the x-axis labels (the second value)
        axis(side=1, at=1:dim(DATA)[1], labels=rep(NA, 103), tck=-0.010, las=1, cex.axis=1.0)
        mtext("Ranked genus", side=1, line=1.4, cex=1.6)
        par(mgp=c(3,0.45,0))				  #set spacing for the y-axis labels (the second value)
        axis(side=2, tck=-0.010, las=1, cex.axis=1.0)
        mtext(expression(paste("Abundance (10"^6, " 16S copies g soil"^-1, ")", sep="")), side=2, line=2.6, at=NA, cex=1.6)
      }

      N.add.panel(DATA=all.rates)
      rm(N.add.panel)
      par(mfrow=c(1,1))

      dev.off()
  #}


#Graph the bootstrapped abundance at Time t, without phyla (horizontal layout, publication width):_____________________________________________________________
  #{
    #Growth rates - in units of r:
      graphics.off()
      # dev.new(width=6.0, height=3.5)
      pdf(file="qSIP_output/Figures/TM_Taxa_abundance_horizontal_pub_width.pdf", width=6.0, height=3.5)
      par(mai=c(0.29,0.47,0.10,0.03), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
      phyla.cols <- cols.for.phyla$col
      ranks <- order(all.rates$N.tot.Tt.boot.median, decreasing=TRUE)
      tax.order <- data.frame(axis.loc=seq(1, dim(all.rates)[1]), ranks)
      y.mins.r <- c(all.rates$N.tot.Tt.boot.median, all.rates$N.tot.Tt.boot.CI.L)/(10^6)
      y.maxs.r <- c(all.rates$N.tot.Tt.boot.median, all.rates$N.tot.Tt.boot.CI.U)/(10^6)
      y.min.r <- min(y.mins.r[y.mins.r != -Inf], na.rm=TRUE)
      y.max.r <- max(y.maxs.r[y.maxs.r != Inf], na.rm=TRUE)

      N.add.panel <- function(DATA){
        plot(y=DATA$N.tot.Tt.boot.median[tax.order$ranks]/(10^6), x=1:dim(DATA)[1], type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", ylim=c(y.min.r, y.max.r), main="")
        curr.comp.ranked <- DATA[tax.order$ranks,]
        mids <- curr.comp.ranked$N.tot.Tt.boot.median/(10^6)
        lowers <- curr.comp.ranked$N.tot.Tt.boot.CI.L/(10^6)
        uppers <- curr.comp.ranked$N.tot.Tt.boot.CI.U/(10^6)
        tax.nums <- 1:length(mids)
        arrows(x0=tax.nums, y0=lowers, x1=tax.nums, y1=uppers, length=0, angle=90, code=3, col=as.character(phyla.cols[2]))
        points(x=tax.nums, y=mids, pch=21, cex=0.6, col=as.character(phyla.cols[2]), bg=as.character(phyla.cols[2]))
        par(mgp=c(3,0,0))			        #set spacing for the x-axis labels (the second value)
        axis(side=1, at=1:dim(DATA)[1], labels=rep(NA, 103), tck=-0.015, las=1, cex.axis=0.6)
        mtext("Ranked genus", side=1, line=0.4, cex=0.75)
        par(mgp=c(3,0.45,0))				  #set spacing for the y-axis labels (the second value)
        axis(side=2, tck=-0.015, las=1, cex.axis=0.6)
        mtext(expression(paste("Abundance (10"^6, " 16S copies g soil"^-1, ")", sep="")), side=2, line=1.7, at=NA, cex=0.75)
      }

      N.add.panel(DATA=all.rates)
      rm(N.add.panel)
      par(mfrow=c(1,1))

      dev.off()
  #}


#Create a two-paneled graph of the bootstrapped abundances (Tt) & growth rates, without phyla (horizontal layout):_____________________________________________
  #{
    #Abundances & growth rates - in units of copies/g soil and r (day-1):
      graphics.off()
      # dev.new(width=7.5, height=7.0)
      pdf(file="qSIP_output/Figures/TM_Taxa_abundance&growth_horizontal.pdf", width=7.5, height=7.0)
      par(mfrow=c(2,1))
      #Abundance panel (Time t):
        par(mai=c(0.29,0.47,0.10,0.03), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
        phyla.cols <- cols.for.phyla$col
        ranks <- order(all.rates$N.tot.Tt.boot.median, decreasing=TRUE)
        tax.order <- data.frame(axis.loc=seq(1, dim(all.rates)[1]), ranks)
        y.mins.N <- c(all.rates$N.tot.Tt.boot.median, all.rates$N.tot.Tt.boot.CI.L)/(10^6)
        y.maxs.N <- c(all.rates$N.tot.Tt.boot.median, all.rates$N.tot.Tt.boot.CI.U)/(10^6)
        y.min.N <- min(y.mins.N[y.mins.N != -Inf], na.rm=TRUE)
        y.max.N <- max(y.maxs.N[y.maxs.N != Inf], na.rm=TRUE)

        N.add.panel <- function(DATA){
          plot(y=DATA$N.tot.Tt.boot.median[tax.order$ranks]/(10^6), x=1:dim(DATA)[1], type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", ylim=c(y.min.N, y.max.N), main="")
          curr.comp.ranked <- DATA[tax.order$ranks,]
          mids <- curr.comp.ranked$N.tot.Tt.boot.median/(10^6)
          lowers <- curr.comp.ranked$N.tot.Tt.boot.CI.L/(10^6)
          uppers <- curr.comp.ranked$N.tot.Tt.boot.CI.U/(10^6)
          tax.nums <- 1:length(mids)
          arrows(x0=tax.nums, y0=lowers, x1=tax.nums, y1=uppers, length=0, angle=90, code=3, col=as.character(phyla.cols[2]))
          points(x=tax.nums, y=mids, pch=21, cex=0.6, col=as.character(phyla.cols[2]), bg=as.character(phyla.cols[2]))
          par(mgp=c(3,0,0))			        #set spacing for the x-axis labels (the second value)
          axis(side=1, at=1:dim(DATA)[1], labels=rep(NA, 103), tck=-0.015, las=1, cex.axis=0.6)
          mtext("Ranked genus", side=1, line=0.4, cex=0.75)
          par(mgp=c(3,0.45,0))				  #set spacing for the y-axis labels (the second value)
          axis(side=2, tck=-0.015, las=1, cex.axis=0.6)
          mtext(expression(paste("Abundance (10"^6, " 16S copies g soil"^-1, ")", sep="")), side=2, line=1.7, at=NA, cex=0.75)
          text(x=(min(1:dim(DATA)[1])+(((max(1:dim(DATA)[1])-min(1:dim(DATA)[1]))))*1.03), y=y.max.N*1.03, "A", adj=c(1,1), cex=0.75)
        }

        N.add.panel(DATA=all.rates)

      #Growth panel:
        par(mai=c(0.29,0.47,0.10,0.03), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
        phyla.cols <- cols.for.phyla$col
        ranks <- order(all.rates$r.net.boot.median, decreasing=TRUE)
        tax.order <- data.frame(axis.loc=seq(1, dim(all.rates)[1]), ranks)
        y.mins.r <- c(all.rates$r.net.boot.median, all.rates$r.net.boot.CI.L, all.rates$b.boot.CI.L, all.rates$d.boot.CI.L)
        y.maxs.r <- c(all.rates$r.net.boot.median, all.rates$r.net.boot.CI.U, all.rates$b.boot.CI.U, all.rates$d.boot.CI.U)
        y.min.r <- min(y.mins.r[y.mins.r != -Inf], na.rm=TRUE)
        y.max.r <- max(y.maxs.r[y.maxs.r != Inf], na.rm=TRUE)

        r.add.panel <- function(DATA){
          plot(y=DATA$r.net.boot.median[tax.order$ranks], x=1:dim(DATA)[1], type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", ylim=c(y.min.r, y.max.r), main="")
          curr.comp.ranked <- DATA[tax.order$ranks,]
          mids.net <- curr.comp.ranked$r.net.boot.median
          lowers.net <- curr.comp.ranked$r.net.boot.CI.L
          uppers.net <- curr.comp.ranked$r.net.boot.CI.U
          tax.nums <- 1:length(mids.net)
          points(x=tax.nums, y=mids.net, pch=21, cex=0.6, col=phyla.cols[7], bg=phyla.cols[7])
          arrows(x0=tax.nums, y0=lowers.net, x1=tax.nums, y1=uppers.net, length=0, angle=90, code=3, col=phyla.cols[7])
          mids.b <- curr.comp.ranked$b.boot.median
          lowers.b <- curr.comp.ranked$b.boot.CI.L
          uppers.b <- curr.comp.ranked$b.boot.CI.U
          arrows(x0=tax.nums, y0=lowers.b, x1=tax.nums, y1=uppers.b, length=0, angle=90, code=3, col="black", lwd=0.5)
          points(x=tax.nums, y=mids.b, pch=21, cex=0.6, col="black", bg="white")
          mids.d <- curr.comp.ranked$d.boot.median
          lowers.d <- curr.comp.ranked$d.boot.CI.L
          uppers.d <- curr.comp.ranked$d.boot.CI.U
          arrows(x0=tax.nums, y0=lowers.d, x1=tax.nums, y1=uppers.d, length=0, angle=90, code=3, col="black", lwd=0.5)
          points(x=tax.nums, y=mids.d, pch=21, cex=0.6, col="black", bg="black")
          legend(x=1-(0.02*(dim(DATA)[1]-1)), y=y.min.r, legend=c("gross", "net", "death"), pch=21, col=c("black", phyla.cols[7], "black"), pt.bg=c("white", phyla.cols[7], "black"), bty="o", cex=0.6, xjust=0, yjust=0, xpd=NA)
          abline(h=0, col="red", lty=2)
          par(mgp=c(3,0,0))			#set spacing for the x-axis labels (the second value)
          axis(side=1, at=1:dim(DATA)[1], labels=rep(NA, 103), tck=-0.015, las=1, cex.axis=0.6)
          mtext("Ranked genus", side=1, line=0.4, cex=0.75)
          par(mgp=c(3,0.45,0))				  #set spacing for the y-axis labels (the second value)
          axis(side=2, tck=-0.015, las=1, cex.axis=0.6)
          mtext(expression(paste(italic("r")~" (day"^-1, ")", sep="")), side=2, line=1.7, at=0, cex=0.75)
          text(x=(min(1:dim(DATA)[1])+(((max(1:dim(DATA)[1])-min(1:dim(DATA)[1]))))*1.03), y=y.max.r*1.03, "B", adj=c(1,1), cex=0.75)
        }

      r.add.panel(DATA=all.rates)
      rm(N.add.panel, r.add.panel)
      par(mfrow=c(1,1))

      dev.off()
  #}


#Create a two-paneled graph of the bootstrapped abundances (Tt) & growth rates, without phyla (horizontal layout, publication width):__________________________
  #{
    #Abundances & growth rates - in units of copies/g soil and r (day-1):
      graphics.off()
      # dev.new(width=6.0, height=7.0)
      pdf(file="qSIP_output/Figures/TM_Taxa_abundance&growth_horizontal_pub_width.pdf", width=6.0, height=7.0)
      par(mfrow=c(2,1))
      #Abundance panel (Time t):
        par(mai=c(0.29,0.47,0.10,0.03), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
        phyla.cols <- cols.for.phyla$col
        ranks <- order(all.rates$N.tot.Tt.boot.median, decreasing=TRUE)
        tax.order <- data.frame(axis.loc=seq(1, dim(all.rates)[1]), ranks)
        y.mins.N <- c(all.rates$N.tot.Tt.boot.median, all.rates$N.tot.Tt.boot.CI.L)/(10^6)
        y.maxs.N <- c(all.rates$N.tot.Tt.boot.median, all.rates$N.tot.Tt.boot.CI.U)/(10^6)
        y.min.N <- min(y.mins.N[y.mins.N != -Inf], na.rm=TRUE)
        y.max.N <- max(y.maxs.N[y.maxs.N != Inf], na.rm=TRUE)

        N.add.panel <- function(DATA){
          plot(y=DATA$N.tot.Tt.boot.median[tax.order$ranks]/(10^6), x=1:dim(DATA)[1], type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", ylim=c(y.min.N, y.max.N), main="")
          curr.comp.ranked <- DATA[tax.order$ranks,]
          mids <- curr.comp.ranked$N.tot.Tt.boot.median/(10^6)
          lowers <- curr.comp.ranked$N.tot.Tt.boot.CI.L/(10^6)
          uppers <- curr.comp.ranked$N.tot.Tt.boot.CI.U/(10^6)
          tax.nums <- 1:length(mids)
          arrows(x0=tax.nums, y0=lowers, x1=tax.nums, y1=uppers, length=0, angle=90, code=3, col=as.character(phyla.cols[2]))
          points(x=tax.nums, y=mids, pch=21, cex=0.6, col=as.character(phyla.cols[2]), bg=as.character(phyla.cols[2]))
          par(mgp=c(3,0,0))			        #set spacing for the x-axis labels (the second value)
          axis(side=1, at=1:dim(DATA)[1], labels=rep(NA, 103), tck=-0.015, las=1, cex.axis=0.6)
          mtext("Ranked genus", side=1, line=0.4, cex=0.75)
          par(mgp=c(3,0.45,0))				  #set spacing for the y-axis labels (the second value)
          axis(side=2, tck=-0.015, las=1, cex.axis=0.6)
          mtext(expression(paste("Abundance (10"^6, " 16S copies g soil"^-1, ")", sep="")), side=2, line=1.7, at=NA, cex=0.75)
          text(x=(min(1:dim(DATA)[1])+(((max(1:dim(DATA)[1])-min(1:dim(DATA)[1]))))*1.03), y=y.max.N*1.03, "A", adj=c(1,1), cex=0.75)
        }

        N.add.panel(DATA=all.rates)

      #Growth panel:
        par(mai=c(0.29,0.47,0.10,0.03), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
        phyla.cols <- cols.for.phyla$col
        ranks <- order(all.rates$r.net.boot.median, decreasing=TRUE)
        tax.order <- data.frame(axis.loc=seq(1, dim(all.rates)[1]), ranks)
        y.mins.r <- c(all.rates$r.net.boot.median, all.rates$r.net.boot.CI.L, all.rates$b.boot.CI.L, all.rates$d.boot.CI.L)
        y.maxs.r <- c(all.rates$r.net.boot.median, all.rates$r.net.boot.CI.U, all.rates$b.boot.CI.U, all.rates$d.boot.CI.U)
        y.min.r <- min(y.mins.r[y.mins.r != -Inf], na.rm=TRUE)
        y.max.r <- max(y.maxs.r[y.maxs.r != Inf], na.rm=TRUE)

        r.add.panel <- function(DATA){
          plot(y=DATA$r.net.boot.median[tax.order$ranks], x=1:dim(DATA)[1], type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", ylim=c(y.min.r, y.max.r), main="")
          curr.comp.ranked <- DATA[tax.order$ranks,]
          mids.net <- curr.comp.ranked$r.net.boot.median
          lowers.net <- curr.comp.ranked$r.net.boot.CI.L
          uppers.net <- curr.comp.ranked$r.net.boot.CI.U
          tax.nums <- 1:length(mids.net)
          points(x=tax.nums, y=mids.net, pch=21, cex=0.6, col=phyla.cols[7], bg=phyla.cols[7])
          arrows(x0=tax.nums, y0=lowers.net, x1=tax.nums, y1=uppers.net, length=0, angle=90, code=3, col=phyla.cols[7])
          mids.b <- curr.comp.ranked$b.boot.median
          lowers.b <- curr.comp.ranked$b.boot.CI.L
          uppers.b <- curr.comp.ranked$b.boot.CI.U
          arrows(x0=tax.nums, y0=lowers.b, x1=tax.nums, y1=uppers.b, length=0, angle=90, code=3, col="black", lwd=0.5)
          points(x=tax.nums, y=mids.b, pch=21, cex=0.6, col="black", bg="white")
          mids.d <- curr.comp.ranked$d.boot.median
          lowers.d <- curr.comp.ranked$d.boot.CI.L
          uppers.d <- curr.comp.ranked$d.boot.CI.U
          arrows(x0=tax.nums, y0=lowers.d, x1=tax.nums, y1=uppers.d, length=0, angle=90, code=3, col="black", lwd=0.5)
          points(x=tax.nums, y=mids.d, pch=21, cex=0.6, col="black", bg="black")
          legend(x=1-(0.02*(dim(DATA)[1]-1)), y=y.min.r, legend=c("gross", "net", "death"), pch=21, col=c("black", phyla.cols[7], "black"), pt.bg=c("white", phyla.cols[7], "black"), bty="o", cex=0.6, xjust=0, yjust=0, xpd=NA)
          abline(h=0, col="red", lty=2)
          par(mgp=c(3,0,0))			#set spacing for the x-axis labels (the second value)
          axis(side=1, at=1:dim(DATA)[1], labels=rep(NA, 103), tck=-0.015, las=1, cex.axis=0.6)
          mtext("Ranked genus", side=1, line=0.4, cex=0.75)
          par(mgp=c(3,0.45,0))				  #set spacing for the y-axis labels (the second value)
          axis(side=2, tck=-0.015, las=1, cex.axis=0.6)
          mtext(expression(paste(italic("r")~" (day"^-1, ")", sep="")), side=2, line=1.7, at=0, cex=0.75)
          text(x=(min(1:dim(DATA)[1])+(((max(1:dim(DATA)[1])-min(1:dim(DATA)[1]))))*1.03), y=y.max.r*1.03, "B", adj=c(1,1), cex=0.75)
        }

      r.add.panel(DATA=all.rates)
      rm(N.add.panel, r.add.panel)
      par(mfrow=c(1,1))

      dev.off()
  #}


#Create a two-paneled graph of the bootstrapped abundances (Tt) & growth rates, without phyla in .eps format (horizontal layout, publication width):__________________________
  #{
    #Abundances & growth rates - in units of copies/g soil and r (day-1):
      graphics.off()
      # dev.new(width=6.0, height=7.0)
      setEPS(width=6.0, height=7.0)
      postscript(file="qSIP_output/Figures/TM_Taxa_abundance&growth_horizontal_pub_width.eps")
      par(mfrow=c(2,1))
      #Abundance panel (Time t):
        par(mai=c(0.29,0.47,0.10,0.03), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
        phyla.cols <- cols.for.phyla$col
        ranks <- order(all.rates$N.tot.Tt.boot.median, decreasing=TRUE)
        tax.order <- data.frame(axis.loc=seq(1, dim(all.rates)[1]), ranks)
        y.mins.N <- c(all.rates$N.tot.Tt.boot.median, all.rates$N.tot.Tt.boot.CI.L)/(10^6)
        y.maxs.N <- c(all.rates$N.tot.Tt.boot.median, all.rates$N.tot.Tt.boot.CI.U)/(10^6)
        y.min.N <- min(y.mins.N[y.mins.N != -Inf], na.rm=TRUE)
        y.max.N <- max(y.maxs.N[y.maxs.N != Inf], na.rm=TRUE)

        N.add.panel <- function(DATA){
          plot(y=DATA$N.tot.Tt.boot.median[tax.order$ranks]/(10^6), x=1:dim(DATA)[1], type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", ylim=c(y.min.N, y.max.N), main="")
          curr.comp.ranked <- DATA[tax.order$ranks,]
          mids <- curr.comp.ranked$N.tot.Tt.boot.median/(10^6)
          lowers <- curr.comp.ranked$N.tot.Tt.boot.CI.L/(10^6)
          uppers <- curr.comp.ranked$N.tot.Tt.boot.CI.U/(10^6)
          tax.nums <- 1:length(mids)
          arrows(x0=tax.nums, y0=lowers, x1=tax.nums, y1=uppers, length=0, angle=90, code=3, col=as.character(phyla.cols[2]))
          points(x=tax.nums, y=mids, pch=21, cex=0.6, col=as.character(phyla.cols[2]), bg=as.character(phyla.cols[2]))
          par(mgp=c(3,0,0))			        #set spacing for the x-axis labels (the second value)
          axis(side=1, at=1:dim(DATA)[1], labels=rep(NA, 103), tck=-0.015, las=1, cex.axis=0.6)
          mtext("Ranked genus", side=1, line=0.4, cex=0.75)
          par(mgp=c(3,0.45,0))				  #set spacing for the y-axis labels (the second value)
          axis(side=2, tck=-0.015, las=1, cex.axis=0.6)
          # mtext(expression(paste("Abundance (10"^6, " 16S copies g soil"^-1, ")", sep="")), side=2, line=1.7, at=NA, cex=0.75)     #alternative version of axis-label
          mtext(expression(paste("Abundance (10"^6, " 16S copies / g soil)", sep="")), side=2, line=1.7, at=NA, cex=0.75)
          text(x=(min(1:dim(DATA)[1])+(((max(1:dim(DATA)[1])-min(1:dim(DATA)[1]))))*1.03), y=y.max.N*1.03, "A", adj=c(1,1), cex=0.75)
        }

        N.add.panel(DATA=all.rates)

      #Growth panel:
        par(mai=c(0.29,0.47,0.10,0.03), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
        phyla.cols <- cols.for.phyla$col
        ranks <- order(all.rates$r.net.boot.median, decreasing=TRUE)
        tax.order <- data.frame(axis.loc=seq(1, dim(all.rates)[1]), ranks)
        y.mins.r <- c(all.rates$r.net.boot.median, all.rates$r.net.boot.CI.L, all.rates$b.boot.CI.L, all.rates$d.boot.CI.L)
        y.maxs.r <- c(all.rates$r.net.boot.median, all.rates$r.net.boot.CI.U, all.rates$b.boot.CI.U, all.rates$d.boot.CI.U)
        y.min.r <- min(y.mins.r[y.mins.r != -Inf], na.rm=TRUE)
        y.max.r <- max(y.maxs.r[y.maxs.r != Inf], na.rm=TRUE)

        r.add.panel <- function(DATA){
          plot(y=DATA$r.net.boot.median[tax.order$ranks], x=1:dim(DATA)[1], type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", ylim=c(y.min.r, y.max.r), main="")
          curr.comp.ranked <- DATA[tax.order$ranks,]
          mids.net <- curr.comp.ranked$r.net.boot.median
          lowers.net <- curr.comp.ranked$r.net.boot.CI.L
          uppers.net <- curr.comp.ranked$r.net.boot.CI.U
          tax.nums <- 1:length(mids.net)
          points(x=tax.nums, y=mids.net, pch=21, cex=0.6, col=phyla.cols[7], bg=phyla.cols[7])
          arrows(x0=tax.nums, y0=lowers.net, x1=tax.nums, y1=uppers.net, length=0, angle=90, code=3, col=phyla.cols[7])
          mids.b <- curr.comp.ranked$b.boot.median
          lowers.b <- curr.comp.ranked$b.boot.CI.L
          uppers.b <- curr.comp.ranked$b.boot.CI.U
          arrows(x0=tax.nums, y0=lowers.b, x1=tax.nums, y1=uppers.b, length=0, angle=90, code=3, col="black", lwd=0.5)
          points(x=tax.nums, y=mids.b, pch=21, cex=0.6, col="black", bg="white")
          mids.d <- curr.comp.ranked$d.boot.median
          lowers.d <- curr.comp.ranked$d.boot.CI.L
          uppers.d <- curr.comp.ranked$d.boot.CI.U
          arrows(x0=tax.nums, y0=lowers.d, x1=tax.nums, y1=uppers.d, length=0, angle=90, code=3, col="black", lwd=0.5)
          points(x=tax.nums, y=mids.d, pch=21, cex=0.6, col="black", bg="black")
          legend(x=1-(0.02*(dim(DATA)[1]-1)), y=y.min.r, legend=c("gross", "net", "death"), pch=21, col=c("black", phyla.cols[7], "black"), pt.bg=c("white", phyla.cols[7], "black"), bty="o", cex=0.6, xjust=0, yjust=0, xpd=NA)
          abline(h=0, col="red", lty=2)
          par(mgp=c(3,0,0))			#set spacing for the x-axis labels (the second value)
          axis(side=1, at=1:dim(DATA)[1], labels=rep(NA, 103), tck=-0.015, las=1, cex.axis=0.6)
          mtext("Ranked genus", side=1, line=0.4, cex=0.75)
          par(mgp=c(3,0.45,0))				  #set spacing for the y-axis labels (the second value)
          axis(side=2, tck=-0.015, las=1, cex.axis=0.6)
          mtext(expression(paste(italic("r")~" (day"^-1, ")", sep="")), side=2, line=1.7, at=0, cex=0.75)
          text(x=(min(1:dim(DATA)[1])+(((max(1:dim(DATA)[1])-min(1:dim(DATA)[1]))))*1.03), y=y.max.r*1.03, "B", adj=c(1,1), cex=0.75)
        }

      r.add.panel(DATA=all.rates)
      rm(N.add.panel, r.add.panel)
      par(mfrow=c(1,1))

      dev.off()
  #}


#Graph the bootstrapped birth, death, and net growth rates only, without phyla (horizontal layout - formatted for powerpoint):_________________________________
  #{
    #Growth rates - in units of r:
      graphics.off()
      # dev.new(width=10, height=5)
      pdf(file="qSIP_output/Figures/TM_Taxa_growth_horizontal_slide.pdf", width=10, height=5)
      par(mai=c(0.49,0.80,0.06,0.03), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
      phyla.cols <- cols.for.phyla$col
      ranks <- order(all.rates$r.net.boot.median, decreasing=TRUE)
      tax.order <- data.frame(axis.loc=seq(1, dim(all.rates)[1]), ranks)
      y.mins.r <- c(all.rates$r.net.boot.median, all.rates$r.net.boot.CI.L, all.rates$b.boot.CI.L, all.rates$d.boot.CI.L)
      y.maxs.r <- c(all.rates$r.net.boot.median, all.rates$r.net.boot.CI.U, all.rates$b.boot.CI.U, all.rates$d.boot.CI.U)
      y.min.r <- min(y.mins.r[y.mins.r != -Inf], na.rm=TRUE)
      y.max.r <- max(y.maxs.r[y.maxs.r != Inf], na.rm=TRUE)

      r.add.panel <- function(DATA){
        plot(y=DATA$r.net.boot.median[tax.order$ranks], x=1:dim(DATA)[1], type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", ylim=c(y.min.r, y.max.r), main="")
        curr.comp.ranked <- DATA[tax.order$ranks,]
        mids.net <- curr.comp.ranked$r.net.boot.median
        lowers.net <- curr.comp.ranked$r.net.boot.CI.L
        uppers.net <- curr.comp.ranked$r.net.boot.CI.U
        tax.nums <- 1:length(mids.net)
        mids.b <- curr.comp.ranked$b.boot.median
        lowers.b <- curr.comp.ranked$b.boot.CI.L
        uppers.b <- curr.comp.ranked$b.boot.CI.U
        mids.d <- curr.comp.ranked$d.boot.median
        lowers.d <- curr.comp.ranked$d.boot.CI.L
        uppers.d <- curr.comp.ranked$d.boot.CI.U

        arrows(x0=tax.nums, y0=lowers.b, x1=tax.nums, y1=uppers.b, length=0, angle=90, code=3, col="black", lwd=1.5)
        arrows(x0=tax.nums, y0=lowers.d, x1=tax.nums, y1=uppers.d, length=0, angle=90, code=3, col="red", lwd=1.5)
        arrows(x0=tax.nums, y0=lowers.net, x1=tax.nums, y1=uppers.net, length=0, angle=90, code=3, col=phyla.cols[7], lwd=1.5)
        points(x=tax.nums, y=mids.b, pch=21, cex=1.5, col="black", bg="white")
        points(x=tax.nums, y=mids.net, pch=21, cex=1.5, col="black", bg=phyla.cols[7])
        points(x=tax.nums, y=mids.d, pch=21, cex=1.5, col="black", bg="red")

        legend(x=dim(DATA)[1]+(0.02*(dim(DATA)[1]-1)), y=y.max.r+(0.02*(y.max.r-y.min.r)), legend=c("reproduction", "net growth", "death"), pch=21, col=c("black", "black", "black"), pt.bg=c("white", phyla.cols[7], "red"), bty="o", pt.cex=1.5, cex=1.0, xjust=1, yjust=1, xpd=NA)
        abline(h=0, col=as.character(phyla.cols[2]), lty=1, lwd=1.5)
        par(mgp=c(3,0,0))			#set spacing for the x-axis labels (the second value)
        axis(side=1, at=1:dim(DATA)[1], labels=rep(NA, 103), tck=-0.010, las=1, cex.axis=1.0)
        mtext("Ranked genus", side=1, line=1.4, cex=1.6)
        par(mgp=c(3,0.45,0))				  #set spacing for the y-axis labels (the second value)
        axis(side=2, tck=-0.010, las=1, cex.axis=1.0)
        mtext(expression(paste(italic("r")~" (day"^-1, ")", sep="")), side=2, line=2.6, at=0, cex=1.6)
      }

      r.add.panel(DATA=all.rates)
      rm(r.add.panel)
      par(mfrow=c(1,1))

      dev.off()
  #}


#Graph the bootstrapped birth rates only, without phyla (horizontal layout - formatted for powerpoint):________________________________________________________
  #{
    #Growth rates - in units of r:
      graphics.off()
      # dev.new(width=10, height=5)
      pdf(file="qSIP_output/Figures/TM_Taxa_birth_horizontal_slide.pdf", width=10, height=5)
      par(mai=c(0.49,0.80,0.06,0.03), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
      phyla.cols <- cols.for.phyla$col
      ranks <- order(all.rates$b.boot.median, decreasing=TRUE)
      tax.order <- data.frame(axis.loc=seq(1, dim(all.rates)[1]), ranks)
      y.mins.r <- c(all.rates$r.net.boot.median, all.rates$r.net.boot.CI.L, all.rates$b.boot.CI.L, all.rates$d.boot.CI.L)
      y.maxs.r <- c(all.rates$r.net.boot.median, all.rates$r.net.boot.CI.U, all.rates$b.boot.CI.U, all.rates$d.boot.CI.U)
      y.min.r <- min(y.mins.r[y.mins.r != -Inf], na.rm=TRUE)
      y.max.r <- max(y.maxs.r[y.maxs.r != Inf], na.rm=TRUE)

      r.add.panel <- function(DATA){
        plot(y=DATA$r.net.boot.median[tax.order$ranks], x=1:dim(DATA)[1], type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", ylim=c(y.min.r, y.max.r), main="")
        curr.comp.ranked <- DATA[tax.order$ranks,]
        mids.net <- curr.comp.ranked$r.net.boot.median
        lowers.net <- curr.comp.ranked$r.net.boot.CI.L
        uppers.net <- curr.comp.ranked$r.net.boot.CI.U
        tax.nums <- 1:length(mids.net)
        mids.b <- curr.comp.ranked$b.boot.median
        lowers.b <- curr.comp.ranked$b.boot.CI.L
        uppers.b <- curr.comp.ranked$b.boot.CI.U
        mids.d <- curr.comp.ranked$d.boot.median
        lowers.d <- curr.comp.ranked$d.boot.CI.L
        uppers.d <- curr.comp.ranked$d.boot.CI.U

        arrows(x0=tax.nums, y0=lowers.b, x1=tax.nums, y1=uppers.b, length=0, angle=90, code=3, col="black", lwd=1.5)
        points(x=tax.nums, y=mids.b, pch=21, cex=1.5, col="black", bg="white")

        abline(h=0, col=as.character(phyla.cols[2]), lty=1, lwd=1.5)
        par(mgp=c(3,0,0))			#set spacing for the x-axis labels (the second value)
        axis(side=1, at=1:dim(DATA)[1], labels=rep(NA, 103), tck=-0.010, las=1, cex.axis=1.0)
        mtext("Ranked genus", side=1, line=1.4, cex=1.6)
        par(mgp=c(3,0.45,0))				  #set spacing for the y-axis labels (the second value)
        axis(side=2, tck=-0.010, las=1, cex.axis=1.0)
        mtext(expression(paste("reproduction (day"^-1, ")", sep="")), side=2, line=2.6, at=0, cex=1.6)
      }

      r.add.panel(DATA=all.rates)
      rm(r.add.panel)
      par(mfrow=c(1,1))

      dev.off()
  #}


#Graph the bootstrapped death rates only, without phyla (horizontal layout - formatted for powerpoint):________________________________________________________
  #{
    #Growth rates - in units of r:
      graphics.off()
      # dev.new(width=10, height=5)
      pdf(file="qSIP_output/Figures/TM_Taxa_death_horizontal_slide.pdf", width=10, height=5)
      par(mai=c(0.49,0.80,0.06,0.03), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
      phyla.cols <- cols.for.phyla$col
      ranks <- order(all.rates$d.boot.median, decreasing=TRUE)
      tax.order <- data.frame(axis.loc=seq(1, dim(all.rates)[1]), ranks)
      y.mins.r <- c(all.rates$r.net.boot.median, all.rates$r.net.boot.CI.L, all.rates$b.boot.CI.L, all.rates$d.boot.CI.L)
      y.maxs.r <- c(all.rates$r.net.boot.median, all.rates$r.net.boot.CI.U, all.rates$b.boot.CI.U, all.rates$d.boot.CI.U)
      y.min.r <- min(y.mins.r[y.mins.r != -Inf], na.rm=TRUE)
      y.max.r <- max(y.maxs.r[y.maxs.r != Inf], na.rm=TRUE)

      r.add.panel <- function(DATA){
        plot(y=DATA$r.net.boot.median[tax.order$ranks], x=1:dim(DATA)[1], type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", ylim=c(y.min.r, y.max.r), main="")
        curr.comp.ranked <- DATA[tax.order$ranks,]
        mids.net <- curr.comp.ranked$r.net.boot.median
        lowers.net <- curr.comp.ranked$r.net.boot.CI.L
        uppers.net <- curr.comp.ranked$r.net.boot.CI.U
        tax.nums <- 1:length(mids.net)
        mids.b <- curr.comp.ranked$b.boot.median
        lowers.b <- curr.comp.ranked$b.boot.CI.L
        uppers.b <- curr.comp.ranked$b.boot.CI.U
        mids.d <- curr.comp.ranked$d.boot.median
        lowers.d <- curr.comp.ranked$d.boot.CI.L
        uppers.d <- curr.comp.ranked$d.boot.CI.U

        arrows(x0=tax.nums, y0=lowers.d, x1=tax.nums, y1=uppers.d, length=0, angle=90, code=3, col="red", lwd=1.5)
        points(x=tax.nums, y=mids.d, pch=21, cex=1.5, col="black", bg="red")

        abline(h=0, col=as.character(phyla.cols[2]), lty=1, lwd=1.5)
        par(mgp=c(3,0,0))			#set spacing for the x-axis labels (the second value)
        axis(side=1, at=1:dim(DATA)[1], labels=rep(NA, 103), tck=-0.010, las=1, cex.axis=1.0)
        mtext("Ranked genus", side=1, line=1.4, cex=1.6)
        par(mgp=c(3,0.45,0))				  #set spacing for the y-axis labels (the second value)
        axis(side=2, tck=-0.010, las=1, cex.axis=1.0)
        mtext(expression(paste("mortality (day"^-1, ")", sep="")), side=2, line=2.6, at=0, cex=1.6)
      }

      r.add.panel(DATA=all.rates)
      rm(r.add.panel)
      par(mfrow=c(1,1))

      dev.off()
  #}


#Graph the bootstrapped net growth rates only, without phyla (horizontal layout - formatted for powerpoint):___________________________________________________
  #{
    #Growth rates - in units of r:
      graphics.off()
      # dev.new(width=10, height=5)
      pdf(file="qSIP_output/Figures/TM_Taxa_net_growth_horizontal_slide.pdf", width=10, height=5)
      par(mai=c(0.49,0.80,0.06,0.03), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
      phyla.cols <- cols.for.phyla$col
      ranks <- order(all.rates$r.net.boot.median, decreasing=TRUE)
      tax.order <- data.frame(axis.loc=seq(1, dim(all.rates)[1]), ranks)
      y.mins.r <- c(all.rates$r.net.boot.median, all.rates$r.net.boot.CI.L, all.rates$b.boot.CI.L, all.rates$d.boot.CI.L)
      y.maxs.r <- c(all.rates$r.net.boot.median, all.rates$r.net.boot.CI.U, all.rates$b.boot.CI.U, all.rates$d.boot.CI.U)
      y.min.r <- min(y.mins.r[y.mins.r != -Inf], na.rm=TRUE)
      y.max.r <- max(y.maxs.r[y.maxs.r != Inf], na.rm=TRUE)

      r.add.panel <- function(DATA){
        plot(y=DATA$r.net.boot.median[tax.order$ranks], x=1:dim(DATA)[1], type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", ylim=c(y.min.r, y.max.r), main="")
        curr.comp.ranked <- DATA[tax.order$ranks,]
        mids.net <- curr.comp.ranked$r.net.boot.median
        lowers.net <- curr.comp.ranked$r.net.boot.CI.L
        uppers.net <- curr.comp.ranked$r.net.boot.CI.U
        tax.nums <- 1:length(mids.net)
        mids.b <- curr.comp.ranked$b.boot.median
        lowers.b <- curr.comp.ranked$b.boot.CI.L
        uppers.b <- curr.comp.ranked$b.boot.CI.U
        mids.d <- curr.comp.ranked$d.boot.median
        lowers.d <- curr.comp.ranked$d.boot.CI.L
        uppers.d <- curr.comp.ranked$d.boot.CI.U

        arrows(x0=tax.nums, y0=lowers.net, x1=tax.nums, y1=uppers.net, length=0, angle=90, code=3, col=phyla.cols[7], lwd=1.5)
        points(x=tax.nums, y=mids.net, pch=21, cex=1.5, col="black", bg=phyla.cols[7])

        abline(h=0, col=as.character(phyla.cols[2]), lty=1, lwd=1.5)
        par(mgp=c(3,0,0))			#set spacing for the x-axis labels (the second value)
        axis(side=1, at=1:dim(DATA)[1], labels=rep(NA, 103), tck=-0.010, las=1, cex.axis=1.0)
        mtext("Ranked genus", side=1, line=1.4, cex=1.6)
        par(mgp=c(3,0.45,0))				  #set spacing for the y-axis labels (the second value)
        axis(side=2, tck=-0.010, las=1, cex.axis=1.0)
        mtext(expression(paste(italic("r")~" (day"^-1, ")", sep="")), side=2, line=2.6, at=0, cex=1.6)
      }

      r.add.panel(DATA=all.rates)
      rm(r.add.panel)
      par(mfrow=c(1,1))

      dev.off()
  #}


#Graph the bootstrapped C fluxes only, without phyla (horizontal layout - formatted for powerpoint):___________________________________________________
  #{
    #Fluxes - in units of micrograms C per gram soil per day (originally was in picograms C per gram soil per day):
      graphics.off()
      # dev.new(width=10, height=5)
      pdf(file="qSIP_output/Figures/TM_Taxa_flux_horizontal_slide.pdf", width=10, height=5)
      par(mai=c(0.49,0.88,0.06,0.03), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
      phyla.cols <- cols.for.phyla$col
      ranks <- order(all.rates$f.gross.boot.median, decreasing=TRUE)
      tax.order <- data.frame(axis.loc=seq(1, dim(all.rates)[1]), ranks)
      y.mins.f <- c(all.rates$f.gross.boot.median, all.rates$f.gross.boot.CI.L)
      y.maxs.f <- c(all.rates$f.gross.boot.median, all.rates$f.gross.boot.CI.U)
      y.min.f <- min(y.mins.f[y.mins.f != -Inf], na.rm=TRUE)
      y.max.f <- max(y.maxs.f[y.maxs.f != Inf], na.rm=TRUE)

      f.add.panel <- function(DATA){
        plot(y=DATA$f.gross.boot.median[tax.order$ranks]/1000000, x=1:dim(DATA)[1], type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", ylim=c(y.min.f/1000000, y.max.f/1000000), main="")
        curr.comp.ranked <- DATA[tax.order$ranks,]
        mids.f <- curr.comp.ranked$f.gross.boot.median/1000000
        lowers.f <- curr.comp.ranked$f.gross.boot.CI.L/1000000
        uppers.f <- curr.comp.ranked$f.gross.boot.CI.U/1000000
        tax.nums <- 1:length(mids.f)

        arrows(x0=tax.nums, y0=lowers.f, x1=tax.nums, y1=uppers.f, length=0, angle=90, code=3, col="black", lwd=1.5)
        points(x=tax.nums, y=mids.f, pch=21, cex=1.5, col="black", bg="orange")

        abline(h=0, col=as.character(phyla.cols[2]), lty=1, lwd=1.5)
        par(mgp=c(3,0,0))			#set spacing for the x-axis labels (the second value)
        axis(side=1, at=1:dim(DATA)[1], labels=rep(NA, 103), tck=-0.010, las=1, cex.axis=1.0)
        mtext("Ranked genus", side=1, line=1.4, cex=1.6)
        par(mgp=c(3,0.45,0))				  #set spacing for the y-axis labels (the second value)
        axis(side=2, tck=-0.010, las=1, cex.axis=1.0)
        mtext(expression(paste("Production (", mu, "g C g soil"^-1, " day"^-1, ")", sep="")), side=2, line=3.25, cex=1.6)
      }

      f.add.panel(DATA=all.rates)
      rm(f.add.panel)
      par(mfrow=c(1,1))

      dev.off()
  #}


#Plot rough graph of rates vs. number of copies/g.soil for all 3 rates (r.net, b, d) at both levels of abundance (T0 & Tt):____________________________________
  #{
    #Check normality of variables in regressions:
      #Histograms of rates:
        graphics.off()
        par(mfrow=c(3,3))
        hist(all.rates$d.obs)
        hist(all.rates$b.obs)
        hist(all.rates$r.net.obs)
        hist(log10(all.rates$d.obs+1))
        hist(log10(all.rates$b.obs+1))
        hist(log10(all.rates$r.net.obs+1))
        hist(sqrt(all.rates$d.obs+1))
        hist(sqrt(all.rates$b.obs+1))
        hist(sqrt(all.rates$r.net.obs+1))
        par(mfrow=c(1,1))

      #QQ plots of rates:
        dev.new()
        par(mfrow=c(3,3))
        qqnorm(all.rates$d.obs)
        mtext(text="d.obs", side=3, line=0, cex=0.6)
        qqnorm(all.rates$b.obs)
        mtext(text="b.obs", side=3, line=0, cex=0.6)
        qqnorm(all.rates$r.net.obs)
        mtext(text="r.net.obs", side=3, line=0, cex=0.6)
        qqnorm(log10(all.rates$d.obs+1))
        mtext(text="log10(d.obs+1)", side=3, line=0, cex=0.6)
        qqnorm(log10(all.rates$b.obs+1))
        mtext(text="log10(b.obs+1)", side=3, line=0, cex=0.6)
        qqnorm(log10(all.rates$r.net.obs+1))
        mtext(text="log10(r.net.obs+1)", side=3, line=0, cex=0.6)
        qqnorm(sqrt(all.rates$d.obs+1))
        mtext(text="sqrt(d.obs+1)", side=3, line=0, cex=0.6)
        qqnorm(sqrt(all.rates$b.obs+1))
        mtext(text="sqrt(b.obs+1)", side=3, line=0, cex=0.6)
        qqnorm(sqrt(all.rates$r.net.obs+1))
        mtext(text="sqrt(r.net.obs+1)", side=3, line=0, cex=0.6)
        par(mfrow=c(1,1))

      #Histograms of abundances:
        dev.new()
        par(mfrow=c(2,2))
        hist(all.rates$N.tot.T0.obs)
        hist(all.rates$N.tot.Tt.obs)
        hist(log10(all.rates$N.tot.T0.obs))
        hist(log10(all.rates$N.tot.Tt.obs))
        par(mfrow=c(1,1))
    
      #QQ plots of abundances:
        dev.new()
        par(mfrow=c(2,2))
        qqnorm(all.rates$N.tot.T0.obs)
        mtext(text="N.tot.T0.obs", side=3, line=0, cex=0.6)
        qqnorm(all.rates$N.tot.Tt.obs)
        mtext(text="N.tot.Tt.obs", side=3, line=0, cex=0.6)
        qqnorm(log10(all.rates$N.tot.T0.obs))
        mtext(text="log10(N.tot.T0.obs)", side=3, line=0, cex=0.6)
        qqnorm(log10(all.rates$N.tot.Tt.obs))
        mtext(text="log10(N.tot.Tt.obs)", side=3, line=0, cex=0.6)
        par(mfrow=c(1,1))

        graphics.off()

    #Create the plots (only transform the abundances (log10)):
      #Note that abundances are imperfect indices of abundance because the number of 16S copies per cell may not be equivalent for all taxa
      #Based on above examinations of normality, log10 transform the abundances; do not transform the rates
        graphics.off()
        min(c(all.rates$d.obs, all.rates$b.obs, all.rates$r.net.obs))
        max(c(all.rates$d.obs, all.rates$b.obs, all.rates$r.net.obs))
        min(c(log10(all.rates$N.tot.T0.obs), log10(all.rates$N.tot.Tt.obs)))
        max(c(log10(all.rates$N.tot.T0.obs), log10(all.rates$N.tot.Tt.obs)))      
        par(mfrow=c(2,3))
        plot(y=all.rates$d.obs, x=log10(all.rates$N.tot.T0.obs), bty="l", type="p", pch=21, col="black", bg="black", ylab="d (day-1)", xlab="log10(abundance) (copies g soil-1 at T0)", xlim=c(4,10), ylim=c(-0.5, 0.5))
        regLine(lm(all.rates$d.obs~log10(all.rates$N.tot.T0.obs)), col="red", lwd=3)
        plot(y=all.rates$b.obs, x=log10(all.rates$N.tot.T0.obs), bty="l", type="p", pch=21, col="black", bg="black", ylab="b (day-1)", xlab="log10(abundance) (copies g soil-1 at T0)", xlim=c(4,10), ylim=c(-0.5, 0.5))
        plot(y=all.rates$r.net.obs, x=log10(all.rates$N.tot.T0.obs), bty="l", type="p", pch=21, col="black", bg="black", ylab="r.net (day-1)", xlab="log10(abundance) (copies g soil-1 at T0)", xlim=c(4,10), ylim=c(-0.5, 0.5))
        regLine(lm(all.rates$r.net.obs~log10(all.rates$N.tot.T0.obs)), col="red", lwd=3)
        plot(y=all.rates$d.obs, x=log10(all.rates$N.tot.Tt.obs), bty="l", type="p", pch=21, col="black", bg="black", ylab="d (day-1)", xlab="log10(abundance) (copies g soil-1 at Tt)", xlim=c(4,10), ylim=c(-0.5, 0.5))
        plot(y=all.rates$b.obs, x=log10(all.rates$N.tot.Tt.obs), bty="l", type="p", pch=21, col="black", bg="black", ylab="b (day-1)", xlab="log10(abundance) (copies g soil-1 at Tt)", xlim=c(4,10), ylim=c(-0.5, 0.5))
        regLine(lm(all.rates$b.obs~log10(all.rates$N.tot.Tt.obs)), col="gray60", lwd=3)
        plot(y=all.rates$r.net.obs, x=log10(all.rates$N.tot.Tt.obs), bty="l", type="p", pch=21, col="black", bg="black", ylab="r (day-1)", xlab="log10(abundance) (copies g soil-1 at Tt)", xlim=c(4,10), ylim=c(-0.5, 0.5))
        par(mfrow=c(1,1))
      #Look at regressions:
        # d vs copies.T0:
          X.temp <- log10(all.rates$N.tot.T0.obs)
          Y.temp <- all.rates$d.obs
          summary(lm(Y.temp~X.temp))  #significant: taxa with the highest death rates during rewetting tended to be those with the most 16S copies present at the start of the incubation (i.e., dominant taxa died rapidly)
        # b vs copies.T0:
          X.temp <- log10(all.rates$N.tot.T0.obs)
          Y.temp <- all.rates$b.obs
          summary(lm(Y.temp~X.temp))
        # r.net vs copies.T0:
          X.temp <- log10(all.rates$N.tot.T0.obs)
          Y.temp <- all.rates$r.net.obs
          summary(lm(Y.temp~X.temp))  #significant: taxa with the greatest net rates of decline during rewetting tended to be those with the most 16S copies present at the start of the incubation (because these taxa also had the highest death rates)
        # d vs copies.Tt:
          X.temp <- log10(all.rates$N.tot.Tt.obs)
          Y.temp <- all.rates$d.obs
          summary(lm(Y.temp~X.temp))
        # b vs copies.Tt:
          X.temp <- log10(all.rates$N.tot.Tt.obs)
          Y.temp <- all.rates$b.obs
          summary(lm(Y.temp~X.temp))  #marginally significant: taxa with the highest 'birth' rates during rewetting tended to have the most 16S copies present by the end of the incubation (i.e. rapid growers became more dominant)
        # r.net vs copies.Tt:
          X.temp <- log10(all.rates$N.tot.Tt.obs)
          Y.temp <- all.rates$r.net.obs
          summary(lm(Y.temp~X.temp))

    #Look at the relationship between death rate and net growth rate:
      graphics.off()
      plot(y=all.rates$r.net.obs, x=all.rates$d.obs, bty="l", type="p", pch=21, col="black", bg="black", ylab="r.net (day-1)", xlab="d (day-1)")
      regLine(lm(all.rates$r.net.obs~all.rates$d.obs), col="red", lwd=3)
      summary(lm(all.rates$r.net.obs~all.rates$d.obs))
      #death rates strongly drove net population growth

    #Look at the relationship between birth rate and net growth rate:
      graphics.off()
      plot(y=all.rates$r.net.obs, x=all.rates$b.obs, bty="l", type="p", pch=21, col="black", bg="black", ylab="r.net (day-1)", xlab="b (day-1)")
      regLine(lm(all.rates$r.net.obs~all.rates$b.obs), col="red", lwd=3)
      summary(lm(all.rates$r.net.obs~all.rates$b.obs))
      #birth rates had only a weak influence on net population growth

    #Look at the relationship between birth rate and death rate:
      graphics.off()
      plot(y=all.rates$b.obs, x=all.rates$d.obs, bty="l", type="p", pch=21, col="black", bg="black", ylab="b (day-1)", xlab="d (day-1)")
      summary(lm(all.rates$b.obs~all.rates$d.obs))
      length(all.rates$b.obs[!is.na(all.rates$b.obs)])      #number of y observations
      length(all.rates$d.obs[!is.na(all.rates$d.obs)])      #number of x observations
      #there was no relationship between 'birth' rate and death rate

      graphics.off()
  #}


#Plot 2-paneled graph of rates vs. number of copies/g.soil (b vs N.Tt & d vs N.T0):____________________________________________________________________________
  #{
    graphics.off()
    # dev.new(width=3.0, height=6.0)
    pdf(file="qSIP_output/Figures/TM_Taxa_birth&death_vs_abundance.pdf", width=3.0, height=6.0)
    par(mfrow=c(2,1))
    #Birth vs. abundance panel (Time t):
      par(mai=c(0.47,0.47,0.02,0.02), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
      x.mins.N <- log10(c(all.rates$N.tot.T0.obs, all.rates$N.tot.Tt.obs))
      x.maxs.N <- log10(c(all.rates$N.tot.T0.obs, all.rates$N.tot.Tt.obs))
      x.min.N <- min(x.mins.N[x.mins.N != -Inf], na.rm=TRUE)
      x.max.N <- max(x.maxs.N[x.maxs.N != Inf], na.rm=TRUE)
      y.mins.b <- all.rates$b.obs
      y.maxs.b <- all.rates$b.obs
      y.min.b <- min(y.mins.b[y.mins.b != -Inf], na.rm=TRUE)
      y.max.b <- max(y.maxs.b[y.maxs.b != Inf], na.rm=TRUE)
      AT <- seq(5, 9, length.out=5)
      AT.minor <- log10(c(seq(0.01,0.1,0.01), seq(0.1,1,0.1), seq(1,10,1), seq(10,100,10), seq(100,1000,100), seq(1000,10000,1000))*10^6)

      b.add.panel <- function(DATA){
        plot(y=DATA$b.obs, x=log10(DATA$N.tot.Tt.obs), type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min.N, x.max.N), ylim=c(y.min.b, y.max.b), main="")
        # regLine(lm(DATA$b.obs~log10(DATA$N.tot.Tt.obs)), col="gray60", lwd=3)
        points(x=log10(DATA$N.tot.Tt.obs), y=DATA$b.obs, pch=21, cex=0.6, col="black", bg="black")
        par(mgp=c(3,0,0))			        #set spacing for the x-axis labels (the second value)
        axis(side=1, at=AT, labels=(10^AT)/(10^6), tck=-0.015, las=1, cex.axis=0.6)
        axis(side=1, at=AT.minor, labels=NA, tck=-0.0075, las=1, cex.axis=0.6)
        mtext(expression(paste("Day 10 abundance (10"^6, " 16S copies g soil"^-1, ")", sep="")), side=1, line=1.4, cex=0.75)
        par(mgp=c(3,0.45,0))				  #set spacing for the y-axis labels (the second value)
        axis(side=2, tck=-0.015, las=1, cex.axis=0.6)
        mtext(expression(paste(italic("b")~" (day"^-1, ")", sep="")), side=2, line=1.7, at=NA, cex=0.75)
        text(x=(x.min.N+((x.max.N-x.min.N)*1.03)), y=(y.min.b+((y.max.b-y.min.b)*1.03)), "A", adj=c(1,1), cex=0.75)
      }

      b.add.panel(DATA=all.rates)

    #Death vs. abundance panel (Time 0):
      par(mai=c(0.42,0.47,0.07,0.02), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
      x.mins.N <- log10(c(all.rates$N.tot.T0.obs, all.rates$N.tot.Tt.obs))
      x.maxs.N <- log10(c(all.rates$N.tot.T0.obs, all.rates$N.tot.Tt.obs))
      x.min.N <- min(x.mins.N[x.mins.N != -Inf], na.rm=TRUE)
      x.max.N <- max(x.maxs.N[x.maxs.N != Inf], na.rm=TRUE)
      y.mins.d <- all.rates$d.obs
      y.maxs.d <- all.rates$d.obs
      y.min.d <- min(y.mins.d[y.mins.d != -Inf], na.rm=TRUE)
      y.max.d <- max(y.maxs.d[y.maxs.d != Inf], na.rm=TRUE)
      AT <- seq(5, 9, length.out=5)
      AT.minor <- log10(c(seq(0.01,0.1,0.01), seq(0.1,1,0.1), seq(1,10,1), seq(10,100,10), seq(100,1000,100), seq(1000,10000,1000))*10^6)

      d.add.panel <- function(DATA){
        plot(y=DATA$d.obs, x=log10(DATA$N.tot.T0.obs), type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min.N, x.max.N), ylim=c(y.min.d, y.max.d), main="")
        regLine(lm(DATA$d.obs~log10(DATA$N.tot.T0.obs)), col="gray60", lwd=3)
        points(x=log10(DATA$N.tot.T0.obs), y=DATA$d.obs, pch=21, cex=0.6, col="black", bg="black")
        par(mgp=c(3,0,0))			        #set spacing for the x-axis labels (the second value)
        axis(side=1, at=AT, labels=(10^AT)/(10^6), tck=-0.015, las=1, cex.axis=0.6)
        axis(side=1, at=AT.minor, labels=NA, tck=-0.0075, las=1, cex.axis=0.6)
        mtext(expression(paste("Day 0 abundance (10"^6, " 16S copies g soil"^-1, ")", sep="")), side=1, line=1.4, cex=0.75)
        par(mgp=c(3,0.45,0))				  #set spacing for the y-axis labels (the second value)
        axis(side=2, tck=-0.015, las=1, cex.axis=0.6)
        mtext(expression(paste(italic("d")~" (day"^-1, ")", sep="")), side=2, line=1.7, at=NA, cex=0.75)
        text(x=(x.min.N+((x.max.N-x.min.N)*1.03)), y=(y.min.d+((y.max.d-y.min.d)*1.03)), "B", adj=c(1,1), cex=0.75)
      }

    d.add.panel(DATA=all.rates)
    rm(b.add.panel, d.add.panel)
    par(mfrow=c(1,1))

    dev.off()
  
    #Statistics for the regressions:
    summary(lm(all.rates$b.obs~log10(all.rates$N.tot.Tt.obs)))
      length(all.rates$b.obs[!is.na(all.rates$b.obs)])                                #number of y observations
      length(log10(all.rates$N.tot.Tt.obs)[!is.na(log10(all.rates$N.tot.Tt.obs))])    #number of x observations
    summary(lm(all.rates$d.obs~log10(all.rates$N.tot.T0.obs)))
      length(all.rates$d.obs[!is.na(all.rates$d.obs)])                                #number of y observations
      length(log10(all.rates$N.tot.T0.obs)[!is.na(log10(all.rates$N.tot.T0.obs))])    #number of x observations
  #}


#Plot 2-paneled graph of rates vs. number of copies/g.soil in .eps format (b vs N.Tt & d vs N.T0):____________________________________________________________________________
  #{
    graphics.off()
    # dev.new(width=3.0, height=6.0)
    setEPS(width=3.0, height=6.0)
    postscript(file="qSIP_output/Figures/TM_Taxa_birth&death_vs_abundance.eps")
    par(mfrow=c(2,1))
    #Birth vs. abundance panel (Time t):
      par(mai=c(0.47,0.47,0.02,0.02), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
      x.mins.N <- log10(c(all.rates$N.tot.T0.obs, all.rates$N.tot.Tt.obs))
      x.maxs.N <- log10(c(all.rates$N.tot.T0.obs, all.rates$N.tot.Tt.obs))
      x.min.N <- min(x.mins.N[x.mins.N != -Inf], na.rm=TRUE)
      x.max.N <- max(x.maxs.N[x.maxs.N != Inf], na.rm=TRUE)
      y.mins.b <- all.rates$b.obs
      y.maxs.b <- all.rates$b.obs
      y.min.b <- min(y.mins.b[y.mins.b != -Inf], na.rm=TRUE)
      y.max.b <- max(y.maxs.b[y.maxs.b != Inf], na.rm=TRUE)
      AT <- seq(5, 9, length.out=5)
      AT.minor <- log10(c(seq(0.01,0.1,0.01), seq(0.1,1,0.1), seq(1,10,1), seq(10,100,10), seq(100,1000,100), seq(1000,10000,1000))*10^6)

      b.add.panel <- function(DATA){
        plot(y=DATA$b.obs, x=log10(DATA$N.tot.Tt.obs), type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min.N, x.max.N), ylim=c(y.min.b, y.max.b), main="")
        # regLine(lm(DATA$b.obs~log10(DATA$N.tot.Tt.obs)), col="gray60", lwd=3)
        points(x=log10(DATA$N.tot.Tt.obs), y=DATA$b.obs, pch=21, cex=0.6, col="black", bg="black")
        par(mgp=c(3,0,0))			        #set spacing for the x-axis labels (the second value)
        axis(side=1, at=AT, labels=(10^AT)/(10^6), tck=-0.015, las=1, cex.axis=0.6)
        axis(side=1, at=AT.minor, labels=NA, tck=-0.0075, las=1, cex.axis=0.6)
        # mtext(expression(paste("Day 10 abundance (10"^6, " 16S copies g soil"^-1, ")", sep="")), side=1, line=1.4, cex=0.75)     #alternative version of axis-label
        mtext(expression(paste("Day 10 abundance (10"^6, " 16S copies / g soil)", sep="")), side=1, line=1.4, cex=0.75)
        par(mgp=c(3,0.45,0))				  #set spacing for the y-axis labels (the second value)
        axis(side=2, tck=-0.015, las=1, cex.axis=0.6)
        mtext(expression(paste(italic("b")~" (day"^-1, ")", sep="")), side=2, line=1.7, at=NA, cex=0.75)
        text(x=(x.min.N+((x.max.N-x.min.N)*1.03)), y=(y.min.b+((y.max.b-y.min.b)*1.03)), "A", adj=c(1,1), cex=0.75)
      }

      b.add.panel(DATA=all.rates)

    #Death vs. abundance panel (Time 0):
      par(mai=c(0.42,0.47,0.07,0.02), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
      x.mins.N <- log10(c(all.rates$N.tot.T0.obs, all.rates$N.tot.Tt.obs))
      x.maxs.N <- log10(c(all.rates$N.tot.T0.obs, all.rates$N.tot.Tt.obs))
      x.min.N <- min(x.mins.N[x.mins.N != -Inf], na.rm=TRUE)
      x.max.N <- max(x.maxs.N[x.maxs.N != Inf], na.rm=TRUE)
      y.mins.d <- all.rates$d.obs
      y.maxs.d <- all.rates$d.obs
      y.min.d <- min(y.mins.d[y.mins.d != -Inf], na.rm=TRUE)
      y.max.d <- max(y.maxs.d[y.maxs.d != Inf], na.rm=TRUE)
      AT <- seq(5, 9, length.out=5)
      AT.minor <- log10(c(seq(0.01,0.1,0.01), seq(0.1,1,0.1), seq(1,10,1), seq(10,100,10), seq(100,1000,100), seq(1000,10000,1000))*10^6)

      d.add.panel <- function(DATA){
        plot(y=DATA$d.obs, x=log10(DATA$N.tot.T0.obs), type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min.N, x.max.N), ylim=c(y.min.d, y.max.d), main="")
        regLine(lm(DATA$d.obs~log10(DATA$N.tot.T0.obs)), col="gray60", lwd=3)
        points(x=log10(DATA$N.tot.T0.obs), y=DATA$d.obs, pch=21, cex=0.6, col="black", bg="black")
        par(mgp=c(3,0,0))			        #set spacing for the x-axis labels (the second value)
        axis(side=1, at=AT, labels=(10^AT)/(10^6), tck=-0.015, las=1, cex.axis=0.6)
        axis(side=1, at=AT.minor, labels=NA, tck=-0.0075, las=1, cex.axis=0.6)
        # mtext(expression(paste("Day 0 abundance (10"^6, " 16S copies g soil"^-1, ")", sep="")), side=1, line=1.4, cex=0.75)     #alternative version of axis-label
        mtext(expression(paste("Day 0 abundance (10"^6, " 16S copies / g soil)", sep="")), side=1, line=1.4, cex=0.75)
        par(mgp=c(3,0.45,0))				  #set spacing for the y-axis labels (the second value)
        axis(side=2, tck=-0.015, las=1, cex.axis=0.6)
        mtext(expression(paste(italic("d")~" (day"^-1, ")", sep="")), side=2, line=1.7, at=NA, cex=0.75)
        text(x=(x.min.N+((x.max.N-x.min.N)*1.03)), y=(y.min.d+((y.max.d-y.min.d)*1.03)), "B", adj=c(1,1), cex=0.75)
      }

    d.add.panel(DATA=all.rates)
    rm(b.add.panel, d.add.panel)
    par(mfrow=c(1,1))

    dev.off()

    #Statistics for the regressions:
    summary(lm(all.rates$b.obs~log10(all.rates$N.tot.Tt.obs)))
      length(all.rates$b.obs[!is.na(all.rates$b.obs)])                                #number of y observations
      length(log10(all.rates$N.tot.Tt.obs)[!is.na(log10(all.rates$N.tot.Tt.obs))])    #number of x observations
    summary(lm(all.rates$d.obs~log10(all.rates$N.tot.T0.obs)))
      length(all.rates$d.obs[!is.na(all.rates$d.obs)])                                #number of y observations
      length(log10(all.rates$N.tot.T0.obs)[!is.na(log10(all.rates$N.tot.T0.obs))])    #number of x observations
  #}


#Plot 2-paneled graph of rates vs. number of copies/g.soil (b vs N.Tt & d vs N.T0) - formatted for powerpoint:_________________________________________________
  #{
    graphics.off()
    # dev.new(width=10.0, height=5.0)
    pdf(file="qSIP_output/Figures/TM_Taxa_birth&death_vs_abundance_slide.pdf", width=10.0, height=5.0)
    par(mfrow=c(1,2))
    #Birth vs. abundance panel (Time t):
      par(mai=c(0.58,0.80,0.03,0.03), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
      x.mins.N <- log10(c(all.rates$N.tot.T0.obs, all.rates$N.tot.Tt.obs))
      x.maxs.N <- log10(c(all.rates$N.tot.T0.obs, all.rates$N.tot.Tt.obs))
      x.min.N <- min(x.mins.N[x.mins.N != -Inf], na.rm=TRUE)
      x.max.N <- max(x.maxs.N[x.maxs.N != Inf], na.rm=TRUE)
      y.mins.b <- all.rates$b.obs
      y.maxs.b <- all.rates$b.obs
      y.min.b <- min(y.mins.b[y.mins.b != -Inf], na.rm=TRUE)
      y.max.b <- max(y.maxs.b[y.maxs.b != Inf], na.rm=TRUE)
      AT <- seq(5, 9, length.out=5)
      AT.minor <- log10(c(seq(0.01,0.1,0.01), seq(0.1,1,0.1), seq(1,10,1), seq(10,100,10), seq(100,1000,100), seq(1000,10000,1000))*10^6)

      b.add.panel <- function(DATA){
        plot(y=DATA$b.obs, x=log10(DATA$N.tot.Tt.obs), type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min.N, x.max.N), ylim=c(y.min.b, y.max.b), main="")
        # regLine(lm(DATA$b.obs~log10(DATA$N.tot.Tt.obs)), col="gray60", lwd=5)
        points(x=log10(DATA$N.tot.Tt.obs), y=DATA$b.obs, pch=21, cex=1.5, col="black", bg="black")
        par(mgp=c(3,0.2,0))			        #set spacing for the x-axis labels (the second value)
        axis(side=1, at=AT, labels=(10^AT)/(10^6), tck=-0.010, las=1, cex.axis=1.0)
        axis(side=1, at=AT.minor, labels=NA, tck=-0.005, las=1, cex.axis=1.0)
        mtext(expression(paste("10"^6, " 16S copies g soil"^-1, " (day 10)", sep="")), side=1, line=2.6, cex=1.6)
        par(mgp=c(3,0.45,0))				  #set spacing for the y-axis labels (the second value)
        axis(side=2, tck=-0.010, las=1, cex.axis=1.0)
        mtext(expression(paste(italic("b")~" (day"^-1, ")", sep="")), side=2, line=2.8, at=NA, cex=1.6)
      }

      b.add.panel(DATA=all.rates)

    #Death vs. abundance panel (Time 0):
      par(mai=c(0.58,0.80,0.03,0.03), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
      x.mins.N <- log10(c(all.rates$N.tot.T0.obs, all.rates$N.tot.Tt.obs))
      x.maxs.N <- log10(c(all.rates$N.tot.T0.obs, all.rates$N.tot.Tt.obs))
      x.min.N <- min(x.mins.N[x.mins.N != -Inf], na.rm=TRUE)
      x.max.N <- max(x.maxs.N[x.maxs.N != Inf], na.rm=TRUE)
      y.mins.d <- all.rates$d.obs
      y.maxs.d <- all.rates$d.obs
      y.min.d <- min(y.mins.d[y.mins.d != -Inf], na.rm=TRUE)
      y.max.d <- max(y.maxs.d[y.maxs.d != Inf], na.rm=TRUE)
      AT <- seq(5, 9, length.out=5)
      AT.minor <- log10(c(seq(0.01,0.1,0.01), seq(0.1,1,0.1), seq(1,10,1), seq(10,100,10), seq(100,1000,100), seq(1000,10000,1000))*10^6)

      d.add.panel <- function(DATA){
        plot(y=DATA$d.obs, x=log10(DATA$N.tot.T0.obs), type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min.N, x.max.N), ylim=c(y.min.d, y.max.d), main="")
        regLine(lm(DATA$d.obs~log10(DATA$N.tot.T0.obs)), col="gray60", lwd=5)
        points(x=log10(DATA$N.tot.T0.obs), y=DATA$d.obs, pch=21, cex=1.5, col="black", bg="black")
        par(mgp=c(3,0.2,0))			        #set spacing for the x-axis labels (the second value)
        axis(side=1, at=AT, labels=(10^AT)/(10^6), tck=-0.010, las=1, cex.axis=1.0)
        axis(side=1, at=AT.minor, labels=NA, tck=-0.005, las=1, cex.axis=1.0)
        mtext(expression(paste("10"^6, " 16S copies g soil"^-1, " (day 0)", sep="")), side=1, line=2.6, cex=1.6)
        par(mgp=c(3,0.45,0))				  #set spacing for the y-axis labels (the second value)
        axis(side=2, tck=-0.010, las=1, cex.axis=1.0)
        mtext(expression(paste(italic("d")~" (day"^-1, ")", sep="")), side=2, line=2.8, at=NA, cex=1.6)
      }

    d.add.panel(DATA=all.rates)
    rm(b.add.panel, d.add.panel)
    par(mfrow=c(1,1))

    dev.off()

    #Statistics for the regressions:
    summary(lm(all.rates$b.obs~log10(all.rates$N.tot.Tt.obs)))
      length(all.rates$b.obs[!is.na(all.rates$b.obs)])                                #number of y observations
      length(log10(all.rates$N.tot.Tt.obs)[!is.na(log10(all.rates$N.tot.Tt.obs))])    #number of x observations
    summary(lm(all.rates$d.obs~log10(all.rates$N.tot.T0.obs)))
      length(all.rates$d.obs[!is.na(all.rates$d.obs)])                                #number of y observations
      length(log10(all.rates$N.tot.T0.obs)[!is.na(log10(all.rates$N.tot.T0.obs))])    #number of x observations
  #}


#Plot graph of birth rate vs. death rate (b vs d) - formatted for powerpoint:__________________________________________________________________________________
  #{
    graphics.off()
    # dev.new(width=5.02, height=5.0)
    pdf(file="qSIP_output/Figures/TM_Taxa_birth_vs_death_slide.pdf", width=5.02, height=5.0)
    #Birth vs. death panel:
      par(mai=c(0.58,0.80,0.03,0.05), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
      x.min <- min(all.rates$d.obs[all.rates$d.obs != -Inf], na.rm=TRUE)
      x.max <- max(all.rates$d.obs[all.rates$d.obs != Inf], na.rm=TRUE)
      y.min <- min(all.rates$b.obs[all.rates$b.obs != -Inf], na.rm=TRUE)
      y.max <- max(all.rates$b.obs[all.rates$b.obs != Inf], na.rm=TRUE)

      bd.add.panel <- function(DATA){
        plot(y=DATA$b.obs, x=DATA$d.obs, type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min, x.max), ylim=c(y.min, y.max), main="")
        # regLine(lm(DATA$b.obs~DATA$d.obs), col="gray60", lwd=5)
        points(x=DATA$d.obs, y=DATA$b.obs, pch=21, cex=1.5, col="black", bg="black")
        par(mgp=c(3,0.2,0))			        #set spacing for the x-axis labels (the second value)
        axis(side=1, tck=-0.010, las=1, cex.axis=1.0)
        mtext(expression(paste(italic("d")~" (day"^-1, ")", sep="")), side=1, line=2.6, cex=1.6)
        par(mgp=c(3,0.45,0))				  #set spacing for the y-axis labels (the second value)
        axis(side=2, tck=-0.010, las=1, cex.axis=1.0)
        mtext(expression(paste(italic("b")~" (day"^-1, ")", sep="")), side=2, line=2.8, at=NA, cex=1.6)
      }

    bd.add.panel(DATA=all.rates)
    rm(bd.add.panel)

    dev.off()

    #Statistics for the regressions:
    summary(lm(all.rates$b.obs~all.rates$d.obs))
      length(all.rates$b.obs[!is.na(all.rates$b.obs)])    #number of y observations
      length(all.rates$d.obs[!is.na(all.rates$d.obs)])    #number of x observations
  #}


#Estimate assemblage-level turnover rates and compare to quick & dirty literature values using radiolabeled nucleotides/amino acids:___________________________
  #{
    #Calculate assemblage-level turnover using death rates and initial abundances (convert death rates to positive values):
      graphics.off()
      assemblage.turnover.d.obs <- sum(-all.rates$d.obs * all.rates$N.tot.T0.obs) / sum(all.rates$N.tot.T0.obs)
      assemblage.turnover.d.boots <- apply(-d.boots.only * N.T0.boots.only, 2, sum) / apply(N.T0.boots.only, 2, sum)
      assemblage.turnover.d.median <- median(assemblage.turnover.d.boots, na.rm=TRUE)
      assemblage.turnover.d.CI.L <- quantile(assemblage.turnover.d.boots, probs=0.05, na.rm=TRUE)
      assemblage.turnover.d.CI.U <- quantile(assemblage.turnover.d.boots, probs=0.95, na.rm=TRUE)
    #Calculate assemblage-level turnover using birth rates and initial abundances (may be better to use time t abundances; but death rate is more appropriate anyway):
      assemblage.turnover.b.obs <- sum(all.rates$b.obs * all.rates$N.tot.T0.obs) / sum(all.rates$N.tot.T0.obs)
      assemblage.turnover.b.boots <- apply(b.boots.only * N.T0.boots.only, 2, sum) / apply(N.T0.boots.only, 2, sum)
      assemblage.turnover.b.median <- median(assemblage.turnover.b.boots, na.rm=TRUE)
      assemblage.turnover.b.CI.L <- quantile(assemblage.turnover.b.boots, probs=0.05, na.rm=TRUE)
      assemblage.turnover.b.CI.U <- quantile(assemblage.turnover.b.boots, probs=0.95, na.rm=TRUE)

    #Quick & dirty literature values for bacterial assemblage-level turnover from Rousk & Baath 2011 (Table 1) (this literature data was imported earlier):
      summary(Lit.turnover)
    #Graph quick & dirty literature turnover values to compare with those from these data (rough graph):
      x.values.lit <- (1:dim(Lit.turnover)[1])-1
      plot(y=Lit.turnover$mean.turnover.rate.days, x=x.values.lit, bty="l", type="p", pch=21, col="black", bg="black", xlim=c(0,100), ylim=c(0,max(Lit.turnover[,4:6], na.rm=TRUE)))
      arrows(x0=x.values.lit, y0=Lit.turnover$min.turnover.rate.days, x1=x.values.lit, y1=Lit.turnover$max.turnover.rate.days, length=0, angle=90, code=3)
      arrows(x0=15, y0=assemblage.turnover.d.CI.L, x1=15, y1=assemblage.turnover.d.CI.U, length=0, angle=90, code=3, col="red")
      points(x=15, y=assemblage.turnover.d.median, pch=21, col="black", bg="red")    #death rate is more reflective of true turnover than is birth rate
      arrows(x0=16, y0=assemblage.turnover.b.CI.L, x1=16, y1=assemblage.turnover.b.CI.U, length=0, angle=90, code=3, col="green")
      points(x=16, y=assemblage.turnover.b.median, pch=21, col="black", bg="green")
      legend(x=100, y=1.04*max(Lit.turnover[,4:6], na.rm=TRUE), legend=c("assemblage turnover (d)", "assemblage turnover (b)"), pch=c(21,21), col="black", pt.bg=c("red", "green"), bty="o", cex=1.0, xjust=1, yjust=1)

    #Turnover values to report in the text:
      #This study (qSIP):
        assemblage.turnover.d.median
        assemblage.turnover.d.CI.L
        assemblage.turnover.d.CI.U
      #Quick & dirty literature values (bulk radiolabeling of the bacterial assemblage):
        mean(Lit.turnover$mean.turnover.rate.days, na.rm=TRUE)
        median(Lit.turnover$mean.turnover.rate.days, na.rm=TRUE)
        min(Lit.turnover$min.turnover.rate.days, na.rm=TRUE)
        max(Lit.turnover$max.turnover.rate.days, na.rm=TRUE)

    graphics.off()
  #}


#Graph the quick & dirty literature bacterial assemblage turnover rates in comparison to those estimated here):________________________________________________
  #{
    #Turnover rates - in units of day-1:
      graphics.off()
      # dev.new(width=3.0, height=3.0)
      pdf(file="qSIP_output/Figures/TM_assemblage_turnover_comparison.pdf", width=3.0, height=3.0)
      par(mai=c(0.29,0.47,0.02,0.02), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
      y.mins <- c(0, as.numeric(unlist(Lit.turnover[,4:6])), assemblage.turnover.d.CI.L)
      y.maxs <- c(0, as.numeric(unlist(Lit.turnover[,4:6])), assemblage.turnover.d.CI.U)
      y.min <- min(y.mins[y.mins != -Inf], na.rm=TRUE)
      y.max <- max(y.maxs[y.maxs != Inf], na.rm=TRUE)
      x.values.lit <- (1:dim(Lit.turnover)[1])-1

      plot(y=c(1,1), x=c(0,15), type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", ylim=c(y.min, y.max), main="")
      arrows(x0=x.values.lit, y0=Lit.turnover$min.turnover.rate.days, x1=x.values.lit, y1=Lit.turnover$max.turnover.rate.days, length=0, angle=90, code=3, col="black", lwd=1.5)
      points(x=x.values.lit, y=Lit.turnover$mean.turnover.rate.days, pch=21, cex=1.0, col="black", bg="black", lwd=1.5)
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
  #}


#Graph the quick & dirty literature bacterial assemblage turnover rates in comparison to those estimated here - formatted for powerpoint:______________________
  #{
    #Turnover rates - in units of day-1:
      graphics.off()
      # dev.new(width=3.5, height=5)
      pdf(file="qSIP_output/Figures/TM_assemblage_turnover_comparison_slide.pdf", width=3.5, height=5)
      par(mai=c(0.49,0.80,0.03,0.03), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
      y.mins <- c(0, as.numeric(unlist(Lit.turnover[,4:6])), assemblage.turnover.d.CI.L)
      y.maxs <- c(0, as.numeric(unlist(Lit.turnover[,4:6])), assemblage.turnover.d.CI.U)
      y.min <- min(y.mins[y.mins != -Inf], na.rm=TRUE)
      y.max <- max(y.maxs[y.maxs != Inf], na.rm=TRUE)
      x.values.lit <- (1:dim(Lit.turnover)[1])-1

      plot(y=c(1,1), x=c(0,15), type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", ylim=c(y.min, y.max), main="")
      arrows(x0=x.values.lit, y0=Lit.turnover$min.turnover.rate.days, x1=x.values.lit, y1=Lit.turnover$max.turnover.rate.days, length=0, angle=90, code=3, col="black", lwd=2)
      points(x=x.values.lit, y=Lit.turnover$mean.turnover.rate.days, pch=21, cex=1.5, col="black", bg="black", lwd=1.5)
      arrows(x0=13.3, y0=assemblage.turnover.d.CI.L, x1=13.3, y1=assemblage.turnover.d.CI.U, length=0, angle=90, code=3, col="red", lwd=2)
      points(x=13.3, y=assemblage.turnover.d.median, pch=21, cex=1.5, col="black", bg="red", lwd=1.5)    #death rate is more reflective of true turnover than is birth rate

      par(mgp=c(3,0,0))			#set spacing for the x-axis labels (the second value)
      axis(side=1, at=0:15, labels=NA, tck=0, las=1, cex.axis=1.0)
      mtext("Literature values", side=1, line=1.4, at=5.5, cex=1.6)
      par(mgp=c(3,0.45,0))				  #set spacing for the y-axis labels (the second value)
      axis(side=2, tck=-0.010, las=1, cex.axis=1.0)
      mtext(expression(paste("Assemblage turnover rate (day"^-1, ")", sep="")), side=2, line=2.6, cex=1.6)

      dev.off()
  #}


#Graph the quick & dirty literature bacterial assemblage turnover rates only - formatted for powerpoint:_______________________________________________________
  #{
    #Turnover rates - in units of day-1:
      graphics.off()
      # dev.new(width=3.5, height=5)
      pdf(file="qSIP_output/Figures/TM_assemblage_turnover_literature_only_slide.pdf", width=3.5, height=5)
      par(mai=c(0.49,0.80,0.03,0.03), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
      y.mins <- c(0, as.numeric(unlist(Lit.turnover[,4:6])), assemblage.turnover.d.CI.L)
      y.maxs <- c(0, as.numeric(unlist(Lit.turnover[,4:6])), assemblage.turnover.d.CI.U)
      y.min <- min(y.mins[y.mins != -Inf], na.rm=TRUE)
      y.max <- max(y.maxs[y.maxs != Inf], na.rm=TRUE)
      x.values.lit <- (1:dim(Lit.turnover)[1])-1

      plot(y=c(1,1), x=c(0,15), type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", ylim=c(y.min, y.max), main="")
      arrows(x0=x.values.lit, y0=Lit.turnover$min.turnover.rate.days, x1=x.values.lit, y1=Lit.turnover$max.turnover.rate.days, length=0, angle=90, code=3, col="black", lwd=2)
      points(x=x.values.lit, y=Lit.turnover$mean.turnover.rate.days, pch=21, cex=1.5, col="black", bg="black", lwd=1.5)

      par(mgp=c(3,0,0))			#set spacing for the x-axis labels (the second value)
      axis(side=1, at=0:15, labels=NA, tck=0, las=1, cex.axis=1.0)
      mtext("Literature values", side=1, line=1.4, at=5.5, cex=1.6)
      par(mgp=c(3,0.45,0))				  #set spacing for the y-axis labels (the second value)
      axis(side=2, tck=-0.010, las=1, cex.axis=1.0)
      mtext(expression(paste("Assemblage turnover rate (day"^-1, ")", sep="")), side=2, line=2.6, cex=1.6)

      dev.off()
  #}


#Estimate assemblage-level turnover rates and compare to a more comprehensive selected collection of literature values using various methods:__________________
  #{
    #Import a more comprehensive set of selected literature values for bacterial assemblage-level turnover:
      Lit.turnover2 <- read.table("qSIP_data/Lit_turnover_for_R_selected.txt", header=TRUE, sep="\t")
      summary(Lit.turnover2)
    #Calculate median and 90%CIs of turnover rate for all literature values:
      Lit.turnover2$median.turnover.rate.d <- NA
      Lit.turnover2$CI.L.turnover.rate.d <- NA
      Lit.turnover2$CI.U.turnover.rate.d <- NA
      set.seed(100)
      for (i in 1:dim(Lit.turnover2)[1]){
        if (!is.na(Lit.turnover2$mean.turnover.time.d[i]) & Lit.turnover2$author[i] != "Baath" & Lit.turnover2$year[i] != 1998){   #Baath 1998 has turnover times and rates already provided
          temp.times.d <- rnorm(1000, mean=Lit.turnover2$mean.turnover.time.d[i], sd=Lit.turnover2$sd.turnover.time.d[i])
          temp.rates.d <- 1/temp.times.d
          Lit.turnover2$mean.turnover.rate.d[i] <- mean(temp.rates.d, na.rm=TRUE)
          Lit.turnover2$sd.turnover.rate.d[i] <- sd(temp.rates.d, na.rm=TRUE)
          Lit.turnover2$median.turnover.rate.d[i] <- median(temp.rates.d, na.rm=TRUE)
          Lit.turnover2$CI.L.turnover.rate.d[i] <- quantile(temp.rates.d, probs=0.05, na.rm=TRUE)
          Lit.turnover2$CI.U.turnover.rate.d[i] <- quantile(temp.rates.d, probs=0.95, na.rm=TRUE)
        } else if (!is.na(Lit.turnover2$mean.turnover.rate.d[i])){
          if (Lit.turnover2$author[i] == "Hunt et al."){   #Hunt et al. only have one estimate -- set this as the median
            Lit.turnover2$mean.turnover.rate.d[i] <- 1.20/365   #ALSO CHANGE THIS VALUE BACK TO THE ORIGINAL VALUE REPORTED IN THEIR PAPER (EVEN THOUGH THAT IS WRONG) TO AVOID HAVING TO INCLUDE ALL THE TEXT EXPLAINING WHY IT IS WRONG AND HOW IT SHOULD BE CORRECTED; SAVE THAT FOR THE TURNOVER META-ANALYSIS THAT USES ALL DATA
            temp.rates.d <- rnorm(1000, mean=Lit.turnover2$mean.turnover.rate.d[i], sd=Lit.turnover2$sd.turnover.rate.d[i])
            Lit.turnover2$median.turnover.rate.d[i] <- Lit.turnover2$mean.turnover.rate.d[i]
          } else if (Lit.turnover2$author[i] != "Hunt et al."){
            temp.rates.d <- rnorm(1000, mean=Lit.turnover2$mean.turnover.rate.d[i], sd=Lit.turnover2$sd.turnover.rate.d[i])
            Lit.turnover2$median.turnover.rate.d[i] <- median(temp.rates.d, na.rm=TRUE)
            Lit.turnover2$CI.L.turnover.rate.d[i] <- quantile(temp.rates.d, probs=0.05, na.rm=TRUE)
            Lit.turnover2$CI.U.turnover.rate.d[i] <- quantile(temp.rates.d, probs=0.95, na.rm=TRUE)
          }
        }
      }

    #Graph the more comprehensive set of selected literature bacterial assemblage turnover rates in comparison to those estimated here:
        #Turnover rates - in units of day-1:
          graphics.off()
          # dev.new(width=3.0, height=3.0)
          pdf(file="qSIP_output/Figures/TM_assemblage_turnover_comparison_selected.pdf", width=3.0, height=3.0)
          par(mai=c(0.29,0.47,0.02,0.02), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
          y.mins <- c(0, as.numeric(unlist(Lit.turnover2[,26:28])), assemblage.turnover.d.CI.L)
          y.maxs <- c(0, as.numeric(unlist(Lit.turnover2[,26:28])), assemblage.turnover.d.CI.U)
          y.min <- min(y.mins[y.mins != -Inf], na.rm=TRUE)
          y.max <- max(y.maxs[y.maxs != Inf], na.rm=TRUE)
          x.values.lit <- (1:dim(Lit.turnover2)[1])-1
          y.medians.lit <- Lit.turnover2$median.turnover.rate.d[order(Lit.turnover2$median.turnover.rate.d, decreasing=FALSE)]
          y.CI.L.lit <- Lit.turnover2$CI.L.turnover.rate.d[order(Lit.turnover2$median.turnover.rate.d, decreasing=FALSE)]
          y.CI.U.lit <- Lit.turnover2$CI.U.turnover.rate.d[order(Lit.turnover2$median.turnover.rate.d, decreasing=FALSE)]

          plot(y=c(1,1), x=c(0, max(x.values.lit)+7), type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", ylim=c(y.min, y.max), main="")
          arrows(x0=x.values.lit, y0=y.CI.L.lit, x1=x.values.lit, y1=y.CI.U.lit, length=0, angle=90, code=3, col="black", lwd=1.5)
          points(x=x.values.lit, y=y.medians.lit, pch=21, cex=1.0, col="black", bg="black", lwd=1.5)
          arrows(x0=33.22, y0=assemblage.turnover.d.CI.L, x1=33.22, y1=assemblage.turnover.d.CI.U, length=0, angle=90, code=3, col="black", lwd=1.5)
          points(x=33.22, y=assemblage.turnover.d.median, pch=21, cex=1.0, col="black", bg="white", lwd=1.5)    #death rate is more reflective of true turnover than is birth rate

          par(mgp=c(3,0,0))			#set spacing for the x-axis labels (the second value)
          axis(side=1, at=0:15, labels=NA, tck=0, las=1, cex.axis=1.0)
          mtext("Literature values", side=1, line=0.4, at=29/2, cex=0.75)
          mtext("This study", side=1, line=0.4, at=33.22, cex=0.75)
          par(mgp=c(3,0.45,0))				  #set spacing for the y-axis labels (the second value)
          axis(side=2, tck=-0.015, las=1, cex.axis=0.6)
          mtext(expression(paste("Assemblage turnover rate (day"^-1, ")", sep="")), side=2, line=1.7, cex=0.75)

          dev.off()


    #Graph the more comprehensive set of selected literature bacterial assemblage turnover rates in comparison to those estimated here -- in .eps format:
        #Turnover rates - in units of day-1:
          graphics.off()
          # dev.new(width=3.0, height=3.0)
          setEPS(width=3.0, height=3.0)
          postscript(file="qSIP_output/Figures/TM_assemblage_turnover_comparison_selected.eps", width=3.0, height=3.0)
          par(mai=c(0.29,0.47,0.02,0.02), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
          y.mins <- c(0, as.numeric(unlist(Lit.turnover2[,26:28])), assemblage.turnover.d.CI.L)
          y.maxs <- c(0, as.numeric(unlist(Lit.turnover2[,26:28])), assemblage.turnover.d.CI.U)
          y.min <- min(y.mins[y.mins != -Inf], na.rm=TRUE)
          y.max <- max(y.maxs[y.maxs != Inf], na.rm=TRUE)
          x.values.lit <- (1:dim(Lit.turnover2)[1])-1
          y.medians.lit <- Lit.turnover2$median.turnover.rate.d[order(Lit.turnover2$median.turnover.rate.d, decreasing=FALSE)]
          y.CI.L.lit <- Lit.turnover2$CI.L.turnover.rate.d[order(Lit.turnover2$median.turnover.rate.d, decreasing=FALSE)]
          y.CI.U.lit <- Lit.turnover2$CI.U.turnover.rate.d[order(Lit.turnover2$median.turnover.rate.d, decreasing=FALSE)]

          plot(y=c(1,1), x=c(0, max(x.values.lit)+7), type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", ylim=c(y.min, y.max), main="")
          arrows(x0=x.values.lit, y0=y.CI.L.lit, x1=x.values.lit, y1=y.CI.U.lit, length=0, angle=90, code=3, col="black", lwd=1.5)
          points(x=x.values.lit, y=y.medians.lit, pch=21, cex=1.0, col="black", bg="black", lwd=1.5)
          arrows(x0=33.22, y0=assemblage.turnover.d.CI.L, x1=33.22, y1=assemblage.turnover.d.CI.U, length=0, angle=90, code=3, col="black", lwd=1.5)
          points(x=33.22, y=assemblage.turnover.d.median, pch=21, cex=1.0, col="black", bg="white", lwd=1.5)    #death rate is more reflective of true turnover than is birth rate

          par(mgp=c(3,0,0))			#set spacing for the x-axis labels (the second value)
          axis(side=1, at=0:15, labels=NA, tck=0, las=1, cex.axis=1.0)
          mtext("Literature values", side=1, line=0.4, at=29/2, cex=0.75)
          mtext("This study", side=1, line=0.4, at=33.22, cex=0.75)
          par(mgp=c(3,0.45,0))				  #set spacing for the y-axis labels (the second value)
          axis(side=2, tck=-0.015, las=1, cex.axis=0.6)
          mtext(expression(paste("Assemblage turnover rate (day"^-1, ")", sep="")), side=2, line=1.7, cex=0.75)

          dev.off()


    #Graph the more comprehensive set of selected literature bacterial assemblage turnover rates in comparison to those estimated here - formatted for powerpoint:__________________
      #{
        #Turnover rates - in units of day-1:
          graphics.off()
          # dev.new(width=3.5, height=5)
          pdf(file="qSIP_output/Figures/TM_assemblage_turnover_comparison_selected_slide.pdf", width=3.5, height=5)
          par(mai=c(0.49,0.80,0.03,0.03), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
          y.mins <- c(0, as.numeric(unlist(Lit.turnover2[,26:28])), assemblage.turnover.d.CI.L)
          y.maxs <- c(0, as.numeric(unlist(Lit.turnover2[,26:28])), assemblage.turnover.d.CI.U)
          y.min <- min(y.mins[y.mins != -Inf], na.rm=TRUE)
          y.max <- max(y.maxs[y.maxs != Inf], na.rm=TRUE)
          x.values.lit <- (1:dim(Lit.turnover2)[1])-1

          plot(y=c(1,1), x=c(0, max(x.values.lit)+7), type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", ylim=c(y.min, y.max), main="")
          arrows(x0=x.values.lit, y0=y.CI.L.lit, x1=x.values.lit, y1=y.CI.U.lit, length=0, angle=90, code=3, col="black", lwd=2)
          points(x=x.values.lit, y=y.medians.lit, pch=21, cex=1.5, col="black", bg="black", lwd=1.5)
          arrows(x0=33.22, y0=assemblage.turnover.d.CI.L, x1=33.22, y1=assemblage.turnover.d.CI.U, length=0, angle=90, code=3, col="red", lwd=2)
          points(x=33.22, y=assemblage.turnover.d.median, pch=21, cex=1.5, col="black", bg="red", lwd=1.5)    #death rate is more reflective of true turnover than is birth rate

          par(mgp=c(3,0,0))			#set spacing for the x-axis labels (the second value)
          axis(side=1, at=0:15, labels=NA, tck=0, las=1, cex.axis=1.0)
          mtext("Literature values", side=1, line=1.4, at=29/2, cex=1.6)
          par(mgp=c(3,0.45,0))				  #set spacing for the y-axis labels (the second value)
          axis(side=2, tck=-0.010, las=1, cex.axis=1.0)
          mtext(expression(paste("Assemblage turnover rate (day"^-1, ")", sep="")), side=2, line=2.6, cex=1.6)

          dev.off()
      #}


    #Graph the more comprehensive set of selected literature bacterial assemblage turnover rates only - formatted for powerpoint:___________________________________________________
      #{
        #Turnover rates - in units of day-1:
          graphics.off()
          # dev.new(width=3.5, height=5)
          pdf(file="qSIP_output/Figures/TM_assemblage_turnover_literature_only_selected_slide.pdf", width=3.5, height=5)
          par(mai=c(0.49,0.80,0.03,0.03), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
          y.mins <- c(0, as.numeric(unlist(Lit.turnover2[,26:28])), assemblage.turnover.d.CI.L)
          y.maxs <- c(0, as.numeric(unlist(Lit.turnover2[,26:28])), assemblage.turnover.d.CI.U)
          y.min <- min(y.mins[y.mins != -Inf], na.rm=TRUE)
          y.max <- max(y.maxs[y.maxs != Inf], na.rm=TRUE)
          x.values.lit <- (1:dim(Lit.turnover2)[1])-1

          plot(y=c(1,1), x=c(0, max(x.values.lit)+7), type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", ylim=c(y.min, y.max), main="")
          arrows(x0=x.values.lit, y0=y.CI.L.lit, x1=x.values.lit, y1=y.CI.U.lit, length=0, angle=90, code=3, col="black", lwd=2)
          points(x=x.values.lit, y=y.medians.lit, pch=21, cex=1.5, col="black", bg="black", lwd=1.5)

          # plot(y=c(1,1), x=c(0,15), type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", ylim=c(y.min, y.max), main="")
          # arrows(x0=x.values.lit, y0=Lit.turnover$min.turnover.rate.days, x1=x.values.lit, y1=Lit.turnover$max.turnover.rate.days, length=0, angle=90, code=3, col="black", lwd=2)
          # points(x=x.values.lit, y=Lit.turnover$mean.turnover.rate.days, pch=21, cex=1.5, col="black", bg="black", lwd=1.5)

          par(mgp=c(3,0,0))			#set spacing for the x-axis labels (the second value)
          axis(side=1, at=0:15, labels=NA, tck=0, las=1, cex.axis=1.0)
          mtext("Literature values", side=1, line=1.4, at=29/2, cex=1.6)
          par(mgp=c(3,0.45,0))				  #set spacing for the y-axis labels (the second value)
          axis(side=2, tck=-0.010, las=1, cex.axis=1.0)
          mtext(expression(paste("Assemblage turnover rate (day"^-1, ")", sep="")), side=2, line=2.6, cex=1.6)

          dev.off()
      #}


    #Calculate the mean, median, minimum, and maximum values of the more comprehensive selected collection of literature values:
      mean(Lit.turnover2$median.turnover.rate.d)
      median(Lit.turnover2$median.turnover.rate.d)
      min(Lit.turnover2$median.turnover.rate.d)
      max(Lit.turnover2$median.turnover.rate.d)

    #Calculate the minimum and maximum values of just the radiolabeled-derived values of the more comprehensive selected collection of literature values:
      Lit.turnover2[c(4:13,20:22),c(1:2,8,9,12,13,26:28)]
      Lit.turnover2$median.turnover.rate.d[c(4:13,20:22)]
      min(Lit.turnover2$median.turnover.rate.d[c(4:13,20:22)])
      max(Lit.turnover2$median.turnover.rate.d[c(4:13,20:22)])

    #Calculate the minimum and maximum values of just the Baath/Bloem radiolabeled-derived values of the more comprehensive selected collection of literature values:
      Lit.turnover2[c(4:13),c(1:2,8,9,12,13,26:28)]
      Lit.turnover2$median.turnover.rate.d[c(4:13)]
      min(Lit.turnover2$median.turnover.rate.d[c(4:13)])
      max(Lit.turnover2$median.turnover.rate.d[c(4:13)])

  #}




#Save workspace:
  save(list=ls(all.names=TRUE), file="qSIP_workspaces/TM_01-02-03-04/.RData", envir=.GlobalEnv)



