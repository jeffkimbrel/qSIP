# This code performs various analyses and produces various figures from the growth, flux, and excess atom fraction results


#_0a_Calculate the number of taxa that grew, grew with priming, grew on glucose with priming:
  #{
  #Week 1:
    #Total number of taxa:
    N.taxa <- length(levels(factor(all.comparisons$taxonID)))
    N.taxa
    #Number of taxa that grew with no added glucose:
    inds <- all.comparisons$comparisonID == 1
    sum(all.comparisons$r.boot.CI.L[inds] > 0, na.rm=TRUE)  
    #Number of taxa that grew with added glucose:
    inds <- all.comparisons$comparisonID == 2
    sum(all.comparisons$r.boot.CI.L[inds] > 0, na.rm=TRUE)  
    #Number of taxa that grew on added glucose:
    inds <- all.comparisons$comparisonID == 3
    sum(all.comparisons$r.boot.CI.L[inds] > 0, na.rm=TRUE)  

  #Week 6:
    #Total number of taxa:
    N.taxa <- length(levels(factor(all.comparisons$taxonID)))
    N.taxa
    #Number of taxa that grew with no added glucose:
    inds <- all.comparisons$comparisonID == 5
    sum(all.comparisons$r.boot.CI.L[inds] > 0, na.rm=TRUE)  
    #Number of taxa that grew with added glucose:
    inds <- all.comparisons$comparisonID == 6
    sum(all.comparisons$r.boot.CI.L[inds] > 0, na.rm=TRUE)  
    #Number of taxa that grew on added glucose:
    inds <- all.comparisons$comparisonID == 7
  sum(all.comparisons$r.boot.CI.L[inds] > 0, na.rm=TRUE)  
  #}




#_0b_Graph output [wad, eaf, r, f] (taxa in same order for all 6 panels within a graph, not colored phylogenetically):
  #{
  #difference in WADS:
    dev.off()
    # dev.new(width=10, height=7.5)
    pdf(file="qSIP_output/Figures/BH_Taxa_same_order_NoPhyla_wad.diff.pdf", width=10, height=7.5)
    par(mfrow=c(2,3))
    par(mai=c(0.44,0.42,0.2,0.05), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
    inds1 <- all.comparisons$comparisonID == 1
    ranks <- order(all.comparisons$wad.diff.boot.median[inds1])
    x.min <- min(all.comparisons$wad.diff.boot.CI.L, na.rm=TRUE)
    x.max <- max(all.comparisons$wad.diff.boot.CI.U, na.rm=TRUE)

    wad.diff.add.panel <- function(comparisonID){
      inds <- all.comparisons$comparisonID == comparisonID
      plot(1:N.taxa~all.comparisons$wad.diff.boot.median[inds][ranks], type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min, x.max), main="")
      mtext(paste(all.comparisons$trt.code.1[inds][1], "v", all.comparisons$trt.code.2[inds][1], sep=" "), side=3, line=-0.1, cex=0.75)
      points(x=all.comparisons$wad.diff.boot.median[inds][ranks], y=1:N.taxa, pch=21, cex=0.6, bg="black")
      arrows(x0=all.comparisons$wad.diff.boot.CI.L[inds][ranks], y0=1:N.taxa, x1=all.comparisons$wad.diff.boot.CI.U[inds][ranks], y1=1:N.taxa, length=0, angle=90, code=3)
      pts.with.message <- c(1:N.taxa)[all.comparisons$message[inds][ranks] != "none"]
      arrows(x0=all.comparisons$wad.diff.boot.CI.U[inds][ranks][pts.with.message]+((x.max-x.min)*0.01), y0=pts.with.message, x1=all.comparisons$wad.diff.boot.CI.U[inds][ranks][pts.with.message]+((x.max-x.min)*0.1), y1=pts.with.message, length=0.06, angle=30, code=1, col="red")
      abline(v=0, col="red")
      par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
      axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
      mtext(expression(paste("Difference in weighted average density (g mL"^-1, ")", sep="")), side=1, line=1.4, cex=0.75)
      par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
      axis(side=2, at=1:N.taxa, labels=rep(NA, 103), tck=-0.015, las=1, cex.axis=0.6)
      mtext("Taxon", side=2, line=1.35, cex=0.75)
    }

    wad.diff.add.panel(comparisonID=1)
    wad.diff.add.panel(comparisonID=2)
    wad.diff.add.panel(comparisonID=3)
    wad.diff.add.panel(comparisonID=5)
    wad.diff.add.panel(comparisonID=6)
    wad.diff.add.panel(comparisonID=7)
    rm(wad.diff.add.panel)
    par(mfrow=c(1,1))

    dev.off()


  #excess atom fraction (atom percent excess, expressed as a fraction rather than as a percentage):
    dev.off()
    # dev.new(width=10, height=7.5)
    pdf(file="qSIP_output/Figures/BH_Taxa_same_order_NoPhyla_ape.pdf", width=10, height=7.5)
    par(mfrow=c(2,3))
    par(mai=c(0.44,0.42,0.2,0.05), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
    inds1 <- all.comparisons$comparisonID == 1
    ranks <- order(all.comparisons$ape.boot.median[inds1])
    x.min <- min(all.comparisons$ape.boot.CI.L, na.rm=TRUE)
    x.max <- max(all.comparisons$ape.boot.CI.U, na.rm=TRUE)

    ape.add.panel <- function(comparisonID){
      inds <- all.comparisons$comparisonID == comparisonID
      plot(1:N.taxa~all.comparisons$ape.boot.median[inds][ranks], type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min, x.max), main="")
      mtext(paste(all.comparisons$trt.code.1[inds][1], "v", all.comparisons$trt.code.2[inds][1], sep=" "), side=3, line=-0.1, cex=0.75)
      points(x=all.comparisons$ape.boot.median[inds][ranks], y=1:N.taxa, pch=21, cex=0.6, bg="black")
      arrows(x0=all.comparisons$ape.boot.CI.L[inds][ranks], y0=1:N.taxa, x1=all.comparisons$ape.boot.CI.U[inds][ranks], y1=1:N.taxa, length=0, angle=90, code=3)
      pts.with.message <- c(1:N.taxa)[all.comparisons$message[inds][ranks] != "none"]
      arrows(x0=all.comparisons$ape.boot.CI.U[inds][ranks][pts.with.message]+((x.max-x.min)*0.01), y0=pts.with.message, x1=all.comparisons$ape.boot.CI.U[inds][ranks][pts.with.message]+((x.max-x.min)*0.1), y1=pts.with.message, length=0.06, angle=30, code=1, col="red")
      abline(v=0, col="red")
      par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
      axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
      mtext("Excess atom fraction", side=1, line=1.4, cex=0.75)
      par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
      axis(side=2, at=1:N.taxa, labels=rep(NA, 103), tck=-0.015, las=1, cex.axis=0.6)
      mtext("Taxon", side=2, line=1.35, cex=0.75)
    }

    ape.add.panel(comparisonID=1)
    ape.add.panel(comparisonID=2)
    ape.add.panel(comparisonID=3)
    ape.add.panel(comparisonID=5)
    ape.add.panel(comparisonID=6)
    ape.add.panel(comparisonID=7)
    rm(ape.add.panel)
    par(mfrow=c(1,1))

    dev.off()


  #intrinsic rate of increase (r):
    dev.off()
    # dev.new(width=10, height=7.5)
    pdf(file="qSIP_output/Figures/BH_Taxa_same_order_NoPhyla_r.pdf", width=10, height=7.5)
    par(mfrow=c(2,3))
    par(mai=c(0.44,0.42,0.2,0.05), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
    inds1 <- all.comparisons$comparisonID == 1
    ranks <- order(all.comparisons$r.boot.median[inds1])
    x.min <- min(all.comparisons$r.boot.CI.L, na.rm=TRUE)
    x.max <- max(all.comparisons$r.boot.CI.U, na.rm=TRUE)

    r.add.panel <- function(comparisonID){
      inds <- all.comparisons$comparisonID == comparisonID
      plot(1:N.taxa~all.comparisons$r.boot.median[inds][ranks], type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min, x.max), main="")
      mtext(paste(all.comparisons$trt.code.1[inds][1], "v", all.comparisons$trt.code.2[inds][1], sep=" "), side=3, line=-0.1, cex=0.75)
      points(x=all.comparisons$r.boot.median[inds][ranks], y=1:N.taxa, pch=21, cex=0.6, bg="black")
      arrows(x0=all.comparisons$r.boot.CI.L[inds][ranks], y0=1:N.taxa, x1=all.comparisons$r.boot.CI.U[inds][ranks], y1=1:N.taxa, length=0, angle=90, code=3)
      pts.with.message <- c(1:N.taxa)[all.comparisons$message[inds][ranks] != "none"]
      arrows(x0=all.comparisons$r.boot.CI.U[inds][ranks][pts.with.message]+((x.max-x.min)*0.01), y0=pts.with.message, x1=all.comparisons$r.boot.CI.U[inds][ranks][pts.with.message]+((x.max-x.min)*0.1), y1=pts.with.message, length=0.06, angle=30, code=1, col="red")
      abline(v=0, col="red")
      par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
      axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
      mtext(expression(paste("r (day"^-1, ")", sep="")), side=1, line=1.4, cex=0.75)
      par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
      axis(side=2, at=1:N.taxa, labels=rep(NA, 103), tck=-0.015, las=1, cex.axis=0.6)
      mtext("Taxon", side=2, line=1.35, cex=0.75)
    }

    r.add.panel(comparisonID=1)
    r.add.panel(comparisonID=2)
    r.add.panel(comparisonID=3)
    r.add.panel(comparisonID=5)
    r.add.panel(comparisonID=6)
    r.add.panel(comparisonID=7)
    rm(r.add.panel)
    par(mfrow=c(1,1))
    
    dev.off()


  #flux of carbon into biomass (convert picograms C per g soil per day to nanograms C per g soil per day):
    dev.off()
    # dev.new(width=10, height=7.5)
    pdf(file="qSIP_output/Figures/BH_Taxa_same_order_NoPhyla_f.pdf", width=10, height=7.5)
    par(mfrow=c(2,3))
    par(mai=c(0.44,0.42,0.2,0.05), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
    inds1 <- all.comparisons$comparisonID == 1
    ranks <- order(all.comparisons$f.boot.median[inds1]/1000)
    x.min <- min(all.comparisons$f.boot.CI.L/1000, na.rm=TRUE)
    x.max <- max(all.comparisons$f.boot.CI.U/1000, na.rm=TRUE)

    f.add.panel <- function(comparisonID){
      inds <- all.comparisons$comparisonID == comparisonID
      plot(y=1:N.taxa, x=all.comparisons$f.boot.median[inds][ranks]/1000, type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min, x.max), main="")
      mtext(paste(all.comparisons$trt.code.1[inds][1], "v", all.comparisons$trt.code.2[inds][1], sep=" "), side=3, line=-0.1, cex=0.75)
      points(x=all.comparisons$f.boot.median[inds][ranks]/1000, y=1:N.taxa, pch=21, cex=0.6, bg="black")
      arrows(x0=all.comparisons$f.boot.CI.L[inds][ranks]/1000, y0=1:N.taxa, x1=all.comparisons$f.boot.CI.U[inds][ranks]/1000, y1=1:N.taxa, length=0, angle=90, code=3)
      pts.with.message <- c(1:N.taxa)[all.comparisons$message[inds][ranks] != "none"]
      arrows(x0=(all.comparisons$f.boot.CI.U[inds][ranks][pts.with.message]/1000)+((x.max-x.min)*0.01), y0=pts.with.message, x1=(all.comparisons$f.boot.CI.U[inds][ranks][pts.with.message]/1000)+((x.max-x.min)*0.1), y1=pts.with.message, length=0.06, angle=30, code=1, col="red")
      abline(v=0, col="red")
      par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
      axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
      mtext(expression(paste("C flux into biomass (ng C g soil"^-1, " day"^-1, ")", sep="")), side=1, line=1.4, cex=0.75)
      par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
      axis(side=2, at=1:N.taxa, labels=rep(NA, 103), tck=-0.015, las=1, cex.axis=0.6)
      mtext("Taxon", side=2, line=1.35, cex=0.75)
    }

    f.add.panel(comparisonID=1)
    f.add.panel(comparisonID=2)
    f.add.panel(comparisonID=3)
    f.add.panel(comparisonID=5)
    f.add.panel(comparisonID=6)
    f.add.panel(comparisonID=7)
    rm(f.add.panel)
    par(mfrow=c(1,1))
    
    dev.off()
  #}




#_0c_Graph output [wad, eaf, r, f] (taxa in same order for all 6 panels within a graph, colored by phylum):
  #{
  #difference in WADS:
    dev.off()
    # dev.new(width=10, height=7.5)
    pdf(file="qSIP_output/Figures/BH_Taxa_same_order_Phyla_wad.diff.pdf", width=10, height=7.5)
    par(mfrow=c(2,3))
    par(mai=c(0.44,0.42,0.2,0.05), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
    phyla <- levels(data.melted$phylum)
    phyla.cols <- rainbow(length(phyla))
    inds1 <- all.comparisons$comparisonID == 1
    ranks <- order(all.comparisons$wad.diff.boot.median[inds1])
    tax.order <- data.frame(axis.loc=seq(1, N.taxa), ranks)
    x.min <- min(all.comparisons$wad.diff.boot.CI.L, na.rm=TRUE)
    x.max <- max(all.comparisons$wad.diff.boot.CI.U, na.rm=TRUE)

    wad.diff.add.panel <- function(comparisonID){
      inds <- all.comparisons$comparisonID == comparisonID
      plot(1:N.taxa~all.comparisons$wad.diff.boot.median[inds][tax.order$ranks], type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min, x.max), main="")
      mtext(paste(all.comparisons$trt.code.1[inds][1], "v", all.comparisons$trt.code.2[inds][1], sep=" "), side=3, line=-0.1, cex=0.75)
      for (p in 1:length(phyla)){
        current.comparison <- all.comparisons[inds,]
        curr.comp.ranked <- current.comparison[tax.order$ranks,]
        tax.nums <- tax.order$axis.loc[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        mids <- curr.comp.ranked$wad.diff.boot.median[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        lowers <- curr.comp.ranked$wad.diff.boot.CI.L[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        uppers <- curr.comp.ranked$wad.diff.boot.CI.U[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        points(x=mids, y=tax.nums, pch=21, cex=0.6, col=phyla.cols[p], bg=phyla.cols[p])
        arrows(x0=lowers, y0=tax.nums, x1=uppers, y1=tax.nums, length=0, angle=90, code=3, col=phyla.cols[p])
      }
      abline(v=0, col="black")
      par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
      axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
      mtext(expression(paste("Difference in weighted average density (g mL"^-1, ")", sep="")), side=1, line=1.4, cex=0.75)
      par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
      axis(side=2, at=1:N.taxa, labels=rep(NA, 103), tck=-0.015, las=1, cex.axis=0.6)
      mtext("Taxon", side=2, line=1.35, cex=0.75)
    }

    wad.diff.add.panel(comparisonID=1)
    legend(x=1.2*x.min, y=1.05*N.taxa, legend=phyla, bty="n", text.col=phyla.cols, cex=0.4)
    wad.diff.add.panel(comparisonID=2)
    wad.diff.add.panel(comparisonID=3)
    wad.diff.add.panel(comparisonID=5)
    wad.diff.add.panel(comparisonID=6)
    wad.diff.add.panel(comparisonID=7)
    rm(wad.diff.add.panel)
    par(mfrow=c(1,1))

    dev.off()


  #excess atom fraction (atom percent excess, expressed as a fraction rather than as a percentage):
    dev.off()
    # dev.new(width=10, height=7.5)
    pdf(file="qSIP_output/Figures/BH_Taxa_same_order_Phyla_ape.pdf", width=10, height=7.5)
    par(mfrow=c(2,3))
    par(mai=c(0.44,0.42,0.2,0.05), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
    phyla <- levels(data.melted$phylum)
    phyla.cols <- rainbow(length(phyla))
    inds1 <- all.comparisons$comparisonID == 1
    ranks <- order(all.comparisons$ape.boot.median[inds1])
    tax.order <- data.frame(axis.loc=seq(1, N.taxa), ranks)
    x.min <- min(all.comparisons$ape.boot.CI.L, na.rm=TRUE)
    x.max <- max(all.comparisons$ape.boot.CI.U, na.rm=TRUE)

    ape.add.panel <- function(comparisonID){
      inds <- all.comparisons$comparisonID == comparisonID
      plot(1:N.taxa~all.comparisons$ape.boot.median[inds][tax.order$ranks], type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min, x.max), main="")
      mtext(paste(all.comparisons$trt.code.1[inds][1], "v", all.comparisons$trt.code.2[inds][1], sep=" "), side=3, line=-0.1, cex=0.75)
      for (p in 1:length(phyla)){
        current.comparison <- all.comparisons[inds,]
        curr.comp.ranked <- current.comparison[tax.order$ranks,]
        tax.nums <- tax.order$axis.loc[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        mids <- curr.comp.ranked$ape.boot.median[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        lowers <- curr.comp.ranked$ape.boot.CI.L[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        uppers <- curr.comp.ranked$ape.boot.CI.U[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        points(x=mids, y=tax.nums, pch=21, cex=0.6, col=phyla.cols[p], bg=phyla.cols[p])
        arrows(x0=lowers, y0=tax.nums, x1=uppers, y1=tax.nums, length=0, angle=90, code=3, col=phyla.cols[p])
      }
      abline(v=0, col="black")
      par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
      axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
      mtext("Excess atom fraction", side=1, line=1.4, cex=0.75)
      par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
      axis(side=2, at=1:N.taxa, labels=rep(NA, 103), tck=-0.015, las=1, cex.axis=0.6)
      mtext("Taxon", side=2, line=1.35, cex=0.75)
    }

    ape.add.panel(comparisonID=1)
    legend(x=1.2*x.min, y=1.05*N.taxa, legend=phyla, bty="n", text.col=phyla.cols, cex=0.4)
    ape.add.panel(comparisonID=2)
    ape.add.panel(comparisonID=3)
    ape.add.panel(comparisonID=5)
    ape.add.panel(comparisonID=6)
    ape.add.panel(comparisonID=7)
    rm(ape.add.panel)
    par(mfrow=c(1,1))

    dev.off()


  #intrinsic rate of increase (r):
    dev.off()
    # dev.new(width=10, height=7.5)
    pdf(file="qSIP_output/Figures/BH_Taxa_same_order_Phyla_r.pdf", width=10, height=7.5)
    par(mfrow=c(2,3))
    par(mai=c(0.44,0.42,0.2,0.05), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
    phyla <- levels(data.melted$phylum)
    phyla.cols <- rainbow(length(phyla))
    inds1 <- all.comparisons$comparisonID == 1
    ranks <- order(all.comparisons$r.boot.median[inds1])
    tax.order <- data.frame(axis.loc=seq(1, N.taxa), ranks)
    x.min <- min(all.comparisons$r.boot.CI.L, na.rm=TRUE)
    x.max <- max(all.comparisons$r.boot.CI.U, na.rm=TRUE)

    r.add.panel <- function(comparisonID){
      inds <- all.comparisons$comparisonID == comparisonID
      plot(1:N.taxa~all.comparisons$r.boot.median[inds][tax.order$ranks], type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min, x.max), main="")
      mtext(paste(all.comparisons$trt.code.1[inds][1], "v", all.comparisons$trt.code.2[inds][1], sep=" "), side=3, line=-0.1, cex=0.75)
      for (p in 1:length(phyla)){
        current.comparison <- all.comparisons[inds,]
        curr.comp.ranked <- current.comparison[tax.order$ranks,]
        tax.nums <- tax.order$axis.loc[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        mids <- curr.comp.ranked$r.boot.median[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        lowers <- curr.comp.ranked$r.boot.CI.L[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        uppers <- curr.comp.ranked$r.boot.CI.U[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        points(x=mids, y=tax.nums, pch=21, cex=0.6, col=phyla.cols[p], bg=phyla.cols[p])
        arrows(x0=lowers, y0=tax.nums, x1=uppers, y1=tax.nums, length=0, angle=90, code=3, col=phyla.cols[p])
      }
      abline(v=0, col="black")
      par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
      axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
      mtext(expression(paste("r (day"^-1, ")", sep="")), side=1, line=1.4, cex=0.75)
      par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
      axis(side=2, at=1:N.taxa, labels=rep(NA, 103), tck=-0.015, las=1, cex.axis=0.6)
      mtext("Taxon", side=2, line=1.35, cex=0.75)
    }

    r.add.panel(comparisonID=1)
    legend(x=0.5*x.max, y=0.8*N.taxa, legend=phyla, bty="n", text.col=phyla.cols, cex=0.4)
    r.add.panel(comparisonID=2)
    r.add.panel(comparisonID=3)
    r.add.panel(comparisonID=5)
    r.add.panel(comparisonID=6)
    r.add.panel(comparisonID=7)
    rm(r.add.panel)
    par(mfrow=c(1,1))

    dev.off()


  #flux of carbon into biomass (convert picograms C per g soil per day to nanograms C per g soil per day):
    dev.off()
    # dev.new(width=10, height=7.5)
    pdf(file="qSIP_output/Figures/BH_Taxa_same_order_Phyla_f.pdf", width=10, height=7.5)
    par(mfrow=c(2,3))
    par(mai=c(0.44,0.42,0.2,0.05), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
    phyla <- levels(data.melted$phylum)
    phyla.cols <- rainbow(length(phyla))
    inds1 <- all.comparisons$comparisonID == 1
    ranks <- order(all.comparisons$f.boot.median[inds1]/1000)
    tax.order <- data.frame(axis.loc=seq(1, N.taxa), ranks)
    x.min <- min(all.comparisons$f.boot.CI.L/1000, na.rm=TRUE)
    x.max <- max(all.comparisons$f.boot.CI.U/1000, na.rm=TRUE)

    f.add.panel <- function(comparisonID){
      inds <- all.comparisons$comparisonID == comparisonID
      plot(y=1:N.taxa, x=all.comparisons$f.boot.median[inds][tax.order$ranks]/1000, type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min, x.max), main="")
      mtext(paste(all.comparisons$trt.code.1[inds][1], "v", all.comparisons$trt.code.2[inds][1], sep=" "), side=3, line=-0.1, cex=0.75)
      for (p in 1:length(phyla)){
        current.comparison <- all.comparisons[inds,]
        curr.comp.ranked <- current.comparison[tax.order$ranks,]
        tax.nums <- tax.order$axis.loc[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        mids <- curr.comp.ranked$f.boot.median[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        lowers <- curr.comp.ranked$f.boot.CI.L[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        uppers <- curr.comp.ranked$f.boot.CI.U[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        points(x=mids/1000, y=tax.nums, pch=21, cex=0.6, col=phyla.cols[p], bg=phyla.cols[p])
        arrows(x0=lowers/1000, y0=tax.nums, x1=uppers/1000, y1=tax.nums, length=0, angle=90, code=3, col=phyla.cols[p])
      }
      abline(v=0, col="black")
      par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
      axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
      mtext(expression(paste("C flux into biomass (ng C g soil"^-1, " day"^-1, ")", sep="")), side=1, line=1.4, cex=0.75)
      par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
      axis(side=2, at=1:N.taxa, labels=rep(NA, 103), tck=-0.015, las=1, cex.axis=0.6)
      mtext("Taxon", side=2, line=1.35, cex=0.75)
    }

    f.add.panel(comparisonID=1)
    legend(x=0.5*x.max, y=0.8*N.taxa, legend=phyla, bty="n", text.col=phyla.cols, cex=0.4)
    f.add.panel(comparisonID=2)
    f.add.panel(comparisonID=3)
    f.add.panel(comparisonID=5)
    f.add.panel(comparisonID=6)
    f.add.panel(comparisonID=7)
    rm(f.add.panel)
    par(mfrow=c(1,1))

    dev.off()
  #}




#_0d_Graph output [wad, eaf, r, f] (taxa re-sorted highest to lowest in each of the 6 panels within a graph, colored by phylum):
  #{
  #difference in WADS:
    dev.off()
    # dev.new(width=10, height=7.5)
    pdf(file="qSIP_output/Figures/BH_Taxa_resorted_Phyla_wad.diff.pdf", width=10, height=7.5)
    par(mfrow=c(2,3))
    par(mai=c(0.44,0.42,0.2,0.05), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
    phyla <- levels(data.melted$phylum)
    phyla.cols <- rainbow(length(phyla))
    x.min <- min(all.comparisons$wad.diff.boot.CI.L, na.rm=TRUE)
    x.max <- max(all.comparisons$wad.diff.boot.CI.U, na.rm=TRUE)

    wad.diff.add.panel <- function(comparisonID){
      inds <- all.comparisons$comparisonID == comparisonID
      ranks <- order(all.comparisons$wad.diff.boot.median[inds])
      tax.order <- data.frame(axis.loc=seq(1, N.taxa), ranks)
      plot(1:N.taxa~all.comparisons$wad.diff.boot.median[inds][tax.order$ranks], type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min, x.max), main="")
      mtext(paste(all.comparisons$trt.code.1[inds][1], "v", all.comparisons$trt.code.2[inds][1], sep=" "), side=3, line=-0.1, cex=0.75)
      for (p in 1:length(phyla)){
        current.comparison <- all.comparisons[inds,]
        curr.comp.ranked <- current.comparison[tax.order$ranks,]
        tax.nums <- tax.order$axis.loc[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        mids <- curr.comp.ranked$wad.diff.boot.median[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        lowers <- curr.comp.ranked$wad.diff.boot.CI.L[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        uppers <- curr.comp.ranked$wad.diff.boot.CI.U[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        points(x=mids, y=tax.nums, pch=21, cex=0.6, col=phyla.cols[p], bg=phyla.cols[p])
        arrows(x0=lowers, y0=tax.nums, x1=uppers, y1=tax.nums, length=0, angle=90, code=3, col=phyla.cols[p])
      }
      abline(v=0, col="black")
      par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
      axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
      mtext(expression(paste("Difference in weighted average density (g mL"^-1, ")", sep="")), side=1, line=1.4, cex=0.75)
      par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
      axis(side=2, at=1:N.taxa, labels=rep(NA, 103), tck=-0.015, las=1, cex.axis=0.6)
      mtext("Taxon", side=2, line=1.35, cex=0.75)
    }

    wad.diff.add.panel(comparisonID=1)
    legend(x=1.2*x.min, y=1.05*N.taxa, legend=phyla, bty="n", text.col=phyla.cols, cex=0.4)
    wad.diff.add.panel(comparisonID=2)
    wad.diff.add.panel(comparisonID=3)
    wad.diff.add.panel(comparisonID=5)
    wad.diff.add.panel(comparisonID=6)
    wad.diff.add.panel(comparisonID=7)
    rm(wad.diff.add.panel)
    par(mfrow=c(1,1))

    dev.off()


  #excess atom fraction (atom percent excess, expressed as a fraction rather than as a percentage):
    dev.off()
    # dev.new(width=10, height=7.5)
    pdf(file="qSIP_output/Figures/BH_Taxa_resorted_Phyla_ape.pdf", width=10, height=7.5)
    par(mfrow=c(2,3))
    par(mai=c(0.44,0.42,0.2,0.05), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
    phyla <- levels(data.melted$phylum)
    phyla.cols <- rainbow(length(phyla))
    x.min <- min(all.comparisons$ape.boot.CI.L, na.rm=TRUE)
    x.max <- max(all.comparisons$ape.boot.CI.U, na.rm=TRUE)

    ape.add.panel <- function(comparisonID){
      inds <- all.comparisons$comparisonID == comparisonID
      ranks <- order(all.comparisons$ape.boot.median[inds])
      tax.order <- data.frame(axis.loc=seq(1, N.taxa), ranks)
      plot(1:N.taxa~all.comparisons$ape.boot.median[inds][tax.order$ranks], type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min, x.max), main="")
      mtext(paste(all.comparisons$trt.code.1[inds][1], "v", all.comparisons$trt.code.2[inds][1], sep=" "), side=3, line=-0.1, cex=0.75)
      for (p in 1:length(phyla)){
        current.comparison <- all.comparisons[inds,]
        curr.comp.ranked <- current.comparison[tax.order$ranks,]
        tax.nums <- tax.order$axis.loc[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        mids <- curr.comp.ranked$ape.boot.median[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        lowers <- curr.comp.ranked$ape.boot.CI.L[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        uppers <- curr.comp.ranked$ape.boot.CI.U[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        points(x=mids, y=tax.nums, pch=21, cex=0.6, col=phyla.cols[p], bg=phyla.cols[p])
        arrows(x0=lowers, y0=tax.nums, x1=uppers, y1=tax.nums, length=0, angle=90, code=3, col=phyla.cols[p])
      }
      abline(v=0, col="black")
      par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
      axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
      mtext("Excess atom fraction", side=1, line=1.4, cex=0.75)
      par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
      axis(side=2, at=1:N.taxa, labels=rep(NA, 103), tck=-0.015, las=1, cex.axis=0.6)
      mtext("Taxon", side=2, line=1.35, cex=0.75)
    }

    ape.add.panel(comparisonID=1)
    legend(x=1.2*x.min, y=1.05*N.taxa, legend=phyla, bty="n", text.col=phyla.cols, cex=0.4)
    ape.add.panel(comparisonID=2)
    ape.add.panel(comparisonID=3)
    ape.add.panel(comparisonID=5)
    ape.add.panel(comparisonID=6)
    ape.add.panel(comparisonID=7)
    rm(ape.add.panel)
    par(mfrow=c(1,1))

    dev.off()


  #intrinsic rate of increase (r):
    dev.off()
    # dev.new(width=10, height=7.5)
    pdf(file="qSIP_output/Figures/BH_Taxa_resorted_Phyla_r.pdf", width=10, height=7.5)
    par(mfrow=c(2,3))
    par(mai=c(0.44,0.42,0.2,0.05), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
    phyla <- levels(data.melted$phylum)
    phyla.cols <- rainbow(length(phyla))
    x.min <- min(all.comparisons$r.boot.CI.L, na.rm=TRUE)
    x.max <- max(all.comparisons$r.boot.CI.U, na.rm=TRUE)

    r.add.panel <- function(comparisonID){
      inds <- all.comparisons$comparisonID == comparisonID
      ranks <- order(all.comparisons$r.boot.median[inds])
      tax.order <- data.frame(axis.loc=seq(1, N.taxa), ranks)
      plot(1:N.taxa~all.comparisons$r.boot.median[inds][tax.order$ranks], type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min, x.max), main="")
      mtext(paste(all.comparisons$trt.code.1[inds][1], "v", all.comparisons$trt.code.2[inds][1], sep=" "), side=3, line=-0.1, cex=0.75)
      for (p in 1:length(phyla)){
        current.comparison <- all.comparisons[inds,]
        curr.comp.ranked <- current.comparison[tax.order$ranks,]
        tax.nums <- tax.order$axis.loc[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        mids <- curr.comp.ranked$r.boot.median[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        lowers <- curr.comp.ranked$r.boot.CI.L[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        uppers <- curr.comp.ranked$r.boot.CI.U[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        points(x=mids, y=tax.nums, pch=21, cex=0.6, col=phyla.cols[p], bg=phyla.cols[p])
        arrows(x0=lowers, y0=tax.nums, x1=uppers, y1=tax.nums, length=0, angle=90, code=3, col=phyla.cols[p])
      }
      abline(v=0, col="black")
      par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
      axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
      mtext(expression(paste("r (day"^-1, ")", sep="")), side=1, line=1.4, cex=0.75)
      par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
      axis(side=2, at=1:N.taxa, labels=rep(NA, 103), tck=-0.015, las=1, cex.axis=0.6)
      mtext("Taxon", side=2, line=1.35, cex=0.75)
    }

    r.add.panel(comparisonID=1)
    legend(x=0.5*x.max, y=0.8*N.taxa, legend=phyla, bty="n", text.col=phyla.cols, cex=0.4)
    r.add.panel(comparisonID=2)
    r.add.panel(comparisonID=3)
    r.add.panel(comparisonID=5)
    r.add.panel(comparisonID=6)
    r.add.panel(comparisonID=7)
    rm(r.add.panel)
    par(mfrow=c(1,1))

    dev.off()


  #flux of carbon into biomass (convert picograms C per g soil per day to nanograms C per g soil per day):
    dev.off()
    # dev.new(width=10, height=7.5)
    pdf(file="qSIP_output/Figures/BH_Taxa_resorted_Phyla_f.pdf", width=10, height=7.5)
    par(mfrow=c(2,3))
    par(mai=c(0.44,0.42,0.2,0.05), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
    phyla <- levels(data.melted$phylum)
    phyla.cols <- rainbow(length(phyla))
    x.min <- min(all.comparisons$f.boot.CI.L/1000, na.rm=TRUE)
    x.max <- max(all.comparisons$f.boot.CI.U/1000, na.rm=TRUE)

    f.add.panel <- function(comparisonID){
      inds <- all.comparisons$comparisonID == comparisonID
      ranks <- order(all.comparisons$f.boot.median[inds]/1000)
      tax.order <- data.frame(axis.loc=seq(1, N.taxa), ranks)
      plot(y=1:N.taxa, x=all.comparisons$f.boot.median[inds][tax.order$ranks]/1000, type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min, x.max), main="")
      mtext(paste(all.comparisons$trt.code.1[inds][1], "v", all.comparisons$trt.code.2[inds][1], sep=" "), side=3, line=-0.1, cex=0.75)
      for (p in 1:length(phyla)){
        current.comparison <- all.comparisons[inds,]
        curr.comp.ranked <- current.comparison[tax.order$ranks,]
        tax.nums <- tax.order$axis.loc[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        mids <- curr.comp.ranked$f.boot.median[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        lowers <- curr.comp.ranked$f.boot.CI.L[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        uppers <- curr.comp.ranked$f.boot.CI.U[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        points(x=mids/1000, y=tax.nums, pch=21, cex=0.6, col=phyla.cols[p], bg=phyla.cols[p])
        arrows(x0=lowers/1000, y0=tax.nums, x1=uppers/1000, y1=tax.nums, length=0, angle=90, code=3, col=phyla.cols[p])
      }
      abline(v=0, col="black")
      par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
      axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
      mtext(expression(paste("C flux into biomass (ng C g soil"^-1, " day"^-1, ")", sep="")), side=1, line=1.4, cex=0.75)
      par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
      axis(side=2, at=1:N.taxa, labels=rep(NA, 103), tck=-0.015, las=1, cex.axis=0.6)
      mtext("Taxon", side=2, line=1.35, cex=0.75)
    }

    f.add.panel(comparisonID=1)
    legend(x=0.5*x.max, y=0.8*N.taxa, legend=phyla, bty="n", text.col=phyla.cols, cex=0.4)
    f.add.panel(comparisonID=2)
    f.add.panel(comparisonID=3)
    f.add.panel(comparisonID=5)
    f.add.panel(comparisonID=6)
    f.add.panel(comparisonID=7)
    rm(f.add.panel)
    par(mfrow=c(1,1))

    dev.off()
  #}
    



#_0e_Graph output [wad, eaf, r, f] (taxa sorted the same way (and grouped by phylum) in all 6 panels within a graph, colored by phylum):
  #{
  #difference in WADS:
    dev.off()
    # dev.new(width=10, height=7.5)
    pdf(file="qSIP_output/Figures/BH_Taxa_same_order_PhylumGroups_wad.diff.pdf", width=10, height=7.5)
    par(mfrow=c(2,3))
    par(mai=c(0.44,0.42,0.2,0.05), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
    phyla <- levels(data.melted$phylum)
    phyla.cols <- rainbow(length(phyla))
    inds1 <- all.comparisons$comparisonID == 1
    ranks <- order(all.comparisons$wad.diff.boot.median[inds1])
    tax.order <- data.frame(axis.loc=seq(1, N.taxa), ranks)
    x.min <- min(all.comparisons$wad.diff.boot.CI.L, na.rm=TRUE)
    x.max <- max(all.comparisons$wad.diff.boot.CI.U, na.rm=TRUE)

    wad.diff.add.panel <- function(comparisonID){
      inds <- all.comparisons$comparisonID == comparisonID
      plot(1:N.taxa~all.comparisons$wad.diff.boot.median[inds][tax.order$ranks], type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min, x.max), main="")
      mtext(paste(all.comparisons$trt.code.1[inds][1], "v", all.comparisons$trt.code.2[inds][1], sep=" "), side=3, line=-0.1, cex=0.75)
      counter <- 1
      for (p in 1:length(phyla)){
        current.comparison <- all.comparisons[inds,]
        curr.comp.ranked <- current.comparison[tax.order$ranks,]
        mids <- curr.comp.ranked$wad.diff.boot.median[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        lowers <- curr.comp.ranked$wad.diff.boot.CI.L[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        uppers <- curr.comp.ranked$wad.diff.boot.CI.U[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        tax.nums <- counter:(counter+length(mids)-1)
        counter <- counter+length(mids)
        points(x=mids, y=tax.nums, pch=21, cex=0.6, col=phyla.cols[p], bg=phyla.cols[p])
        arrows(x0=lowers, y0=tax.nums, x1=uppers, y1=tax.nums, length=0, angle=90, code=3, col=phyla.cols[p])
      }
      abline(v=0, col="black")
      par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
      axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
      mtext(expression(paste("Difference in weighted average density (g mL"^-1, ")", sep="")), side=1, line=1.4, cex=0.75)
      par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
      axis(side=2, at=1:N.taxa, labels=rep(NA, 103), tck=-0.015, las=1, cex.axis=0.6)
      mtext("Taxon", side=2, line=1.35, cex=0.75)
    }

    wad.diff.add.panel(comparisonID=1)
    legend(x=1.2*x.min, y=1.05*N.taxa, legend=phyla, bty="n", text.col=phyla.cols, cex=0.4)
    wad.diff.add.panel(comparisonID=2)
    wad.diff.add.panel(comparisonID=3)
    wad.diff.add.panel(comparisonID=5)
    wad.diff.add.panel(comparisonID=6)
    wad.diff.add.panel(comparisonID=7)
    rm(wad.diff.add.panel)
    par(mfrow=c(1,1))

    dev.off()


  #excess atom fraction (atom percent excess, expressed as a fraction rather than as a percentage):
    dev.off()
    # dev.new(width=10, height=7.5)
    pdf(file="qSIP_output/Figures/BH_Taxa_same_order_PhylumGroups_ape.pdf", width=10, height=7.5)
    par(mfrow=c(2,3))
    par(mai=c(0.44,0.42,0.2,0.05), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
    phyla <- levels(data.melted$phylum)
    phyla.cols <- rainbow(length(phyla))
    inds1 <- all.comparisons$comparisonID == 1
    ranks <- order(all.comparisons$ape.boot.median[inds1])
    tax.order <- data.frame(axis.loc=seq(1, N.taxa), ranks)
    x.min <- min(all.comparisons$ape.boot.CI.L, na.rm=TRUE)
    x.max <- max(all.comparisons$ape.boot.CI.U, na.rm=TRUE)

    ape.add.panel <- function(comparisonID){
      inds <- all.comparisons$comparisonID == comparisonID
      plot(1:N.taxa~all.comparisons$ape.boot.median[inds][tax.order$ranks], type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min, x.max), main="")
      mtext(paste(all.comparisons$trt.code.1[inds][1], "v", all.comparisons$trt.code.2[inds][1], sep=" "), side=3, line=-0.1, cex=0.75)
      counter <- 1
      for (p in 1:length(phyla)){
        current.comparison <- all.comparisons[inds,]
        curr.comp.ranked <- current.comparison[tax.order$ranks,]
        mids <- curr.comp.ranked$ape.boot.median[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        lowers <- curr.comp.ranked$ape.boot.CI.L[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        uppers <- curr.comp.ranked$ape.boot.CI.U[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        tax.nums <- counter:(counter+length(mids)-1)
        counter <- counter+length(mids)
        points(x=mids, y=tax.nums, pch=21, cex=0.6, col=phyla.cols[p], bg=phyla.cols[p])
        arrows(x0=lowers, y0=tax.nums, x1=uppers, y1=tax.nums, length=0, angle=90, code=3, col=phyla.cols[p])
      }
      abline(v=0, col="black")
      par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
      axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
      mtext("Excess atom fraction", side=1, line=1.4, cex=0.75)
      par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
      axis(side=2, at=1:N.taxa, labels=rep(NA, 103), tck=-0.015, las=1, cex.axis=0.6)
      mtext("Taxon", side=2, line=1.35, cex=0.75)
    }

    ape.add.panel(comparisonID=1)
    legend(x=1.2*x.min, y=1.05*N.taxa, legend=phyla, bty="n", text.col=phyla.cols, cex=0.4)
    ape.add.panel(comparisonID=2)
    ape.add.panel(comparisonID=3)
    ape.add.panel(comparisonID=5)
    ape.add.panel(comparisonID=6)
    ape.add.panel(comparisonID=7)
    rm(ape.add.panel)
    par(mfrow=c(1,1))

    dev.off()


  #intrinsic rate of increase (r):
    dev.off()
    # dev.new(width=10, height=7.5)
    pdf(file="qSIP_output/Figures/BH_Taxa_same_order_PhylumGroups_r.pdf", width=10, height=7.5)
    par(mfrow=c(2,3))
    par(mai=c(0.44,0.42,0.2,0.05), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
    phyla <- levels(data.melted$phylum)
    phyla.cols <- rainbow(length(phyla))
    inds1 <- all.comparisons$comparisonID == 1
    ranks <- order(all.comparisons$r.boot.median[inds1])
    tax.order <- data.frame(axis.loc=seq(1, N.taxa), ranks)
    x.min <- min(all.comparisons$r.boot.CI.L, na.rm=TRUE)
    x.max <- max(all.comparisons$r.boot.CI.U, na.rm=TRUE)

    r.add.panel <- function(comparisonID){
      inds <- all.comparisons$comparisonID == comparisonID
      plot(1:N.taxa~all.comparisons$r.boot.median[inds][tax.order$ranks], type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min, x.max), main="")
      mtext(paste(all.comparisons$trt.code.1[inds][1], "v", all.comparisons$trt.code.2[inds][1], sep=" "), side=3, line=-0.1, cex=0.75)
      counter <- 1
      for (p in 1:length(phyla)){
        current.comparison <- all.comparisons[inds,]
        curr.comp.ranked <- current.comparison[tax.order$ranks,]
        mids <- curr.comp.ranked$r.boot.median[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        lowers <- curr.comp.ranked$r.boot.CI.L[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        uppers <- curr.comp.ranked$r.boot.CI.U[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        tax.nums <- counter:(counter+length(mids)-1)
        counter <- counter+length(mids)
        points(x=mids, y=tax.nums, pch=21, cex=0.6, col=phyla.cols[p], bg=phyla.cols[p])
        arrows(x0=lowers, y0=tax.nums, x1=uppers, y1=tax.nums, length=0, angle=90, code=3, col=phyla.cols[p])
      }
      abline(v=0, col="black")
      par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
      axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
      mtext(expression(paste("r (day"^-1, ")", sep="")), side=1, line=1.4, cex=0.75)
      par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
      axis(side=2, at=1:N.taxa, labels=rep(NA, 103), tck=-0.015, las=1, cex.axis=0.6)
      mtext("Taxon", side=2, line=1.35, cex=0.75)
    }

    r.add.panel(comparisonID=1)
    legend(x=0.5*x.max, y=0.8*N.taxa, legend=phyla, bty="n", text.col=phyla.cols, cex=0.4)
    r.add.panel(comparisonID=2)
    r.add.panel(comparisonID=3)
    r.add.panel(comparisonID=5)
    r.add.panel(comparisonID=6)
    r.add.panel(comparisonID=7)
    rm(r.add.panel)
    par(mfrow=c(1,1))

    dev.off()


  #flux of carbon into biomass (convert picograms C per g soil per day to nanograms C per g soil per day):
    dev.off()
    # dev.new(width=10, height=7.5)
    pdf(file="qSIP_output/Figures/BH_Taxa_same_order_PhylumGroups_f.pdf", width=10, height=7.5)
    par(mfrow=c(2,3))
    par(mai=c(0.44,0.42,0.2,0.05), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
    phyla <- levels(data.melted$phylum)
    phyla.cols <- rainbow(length(phyla))
    inds1 <- all.comparisons$comparisonID == 1
    ranks <- order(all.comparisons$f.boot.median[inds1]/1000)
    tax.order <- data.frame(axis.loc=seq(1, N.taxa), ranks)
    x.min <- min(all.comparisons$f.boot.CI.L/1000, na.rm=TRUE)
    x.max <- max(all.comparisons$f.boot.CI.U/1000, na.rm=TRUE)

    f.add.panel <- function(comparisonID){
      inds <- all.comparisons$comparisonID == comparisonID
      plot(y=1:N.taxa, x=all.comparisons$f.boot.median[inds][tax.order$ranks]/1000, type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min, x.max), main="")
      mtext(paste(all.comparisons$trt.code.1[inds][1], "v", all.comparisons$trt.code.2[inds][1], sep=" "), side=3, line=-0.1, cex=0.75)
      counter <- 1
      for (p in 1:length(phyla)){
        current.comparison <- all.comparisons[inds,]
        curr.comp.ranked <- current.comparison[tax.order$ranks,]
        mids <- curr.comp.ranked$f.boot.median[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        lowers <- curr.comp.ranked$f.boot.CI.L[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        uppers <- curr.comp.ranked$f.boot.CI.U[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        tax.nums <- counter:(counter+length(mids)-1)
        counter <- counter+length(mids)
        points(x=mids/1000, y=tax.nums, pch=21, cex=0.6, col=phyla.cols[p], bg=phyla.cols[p])
        arrows(x0=lowers/1000, y0=tax.nums, x1=uppers/1000, y1=tax.nums, length=0, angle=90, code=3, col=phyla.cols[p])
      }
      abline(v=0, col="black")
      par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
      axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
      mtext(expression(paste("C flux into biomass (ng C g soil"^-1, " day"^-1, ")", sep="")), side=1, line=1.4, cex=0.75)
      par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
      axis(side=2, at=1:N.taxa, labels=rep(NA, 103), tck=-0.015, las=1, cex.axis=0.6)
      mtext("Taxon", side=2, line=1.35, cex=0.75)
    }

    f.add.panel(comparisonID=1)
    legend(x=0.5*x.max, y=0.8*N.taxa, legend=phyla, bty="n", text.col=phyla.cols, cex=0.4)
    f.add.panel(comparisonID=2)
    f.add.panel(comparisonID=3)
    f.add.panel(comparisonID=5)
    f.add.panel(comparisonID=6)
    f.add.panel(comparisonID=7)
    rm(f.add.panel)
    par(mfrow=c(1,1))

    dev.off()
  #}




#_0f_Plot primed-glucose growth rates (r determined by 13C; units=day-1) vs. primed growth rates (r determined by 18O; units=day-1):
#(taxa not colored phylogenetically)
  #{
  #Week 1:
    dev.off()
    # dev.new(width=6, height=6)
    pdf(file="qSIP_output/Figures/BH_GlucGrowth_v_Growth_NoPhyla_Week1.pdf", width=6, height=6)
    par(mai=c(0.9, 0.9, 0.4, 0.4))
    inds2 <- all.comparisons$comparisonID == 2            #growth with priming
    ranks <- order(all.comparisons$r.boot.median[inds1])
    inds3 <- all.comparisons$comparisonID == 3            #growth on glucose with priming
    x.min <- min(c(all.comparisons$r.boot.CI.L[inds2], all.comparisons$r.boot.CI.L[inds3]), na.rm=TRUE)
    x.max <- max(c(all.comparisons$r.boot.CI.U[inds2], all.comparisons$r.boot.CI.U[inds3]), na.rm=TRUE)
    plot(all.comparisons$r.boot.median[inds3][ranks]~all.comparisons$r.boot.median[inds2][ranks], type="n", bty="l", xlab=expression(paste("Growth rate (r) with priming (day"^-1, ")", sep="")), ylab=expression(paste("Growth rate (r) on glucose with priming (day"^-1, ")", sep="")), xlim=c(x.min, x.max), ylim=c(x.min, x.max), main="Week 1")
    #plot(y=c(0,x.max), x=c(0,x.max), type="n", bty="l", xlab=expression(paste("Growth rate (r) with priming (day"^-1, ")", sep="")), ylab=expression(paste("Growth rate (r) on glucose with priming (day"^-1, ")", sep="")), xaxs="i", yaxs="i", main="Week 1")
    abline(h=0, col="red")
    abline(v=0, col="red")
    abline(a=0, b=1, lty=2, col="blue")
    points(x=all.comparisons$r.boot.median[inds2][ranks], y=all.comparisons$r.boot.median[inds3][ranks], pch=21, cex=0.6, bg="black")
    arrows(x0=all.comparisons$r.boot.CI.L[inds2][ranks], y0=all.comparisons$r.boot.median[inds3][ranks], x1=all.comparisons$r.boot.CI.U[inds2][ranks], y1=all.comparisons$r.boot.median[inds3][ranks], length=0, angle=90, code=3)
    arrows(x0=all.comparisons$r.boot.median[inds2][ranks], y0=all.comparisons$r.boot.CI.L[inds3][ranks], x1=all.comparisons$r.boot.median[inds2][ranks], y1=all.comparisons$r.boot.CI.U[inds3][ranks], length=0, angle=90, code=3)

    dev.off()


  #Week 6:
    dev.off()
    # dev.new(width=6, height=6)
    pdf(file="qSIP_output/Figures/BH_GlucGrowth_v_Growth_NoPhyla_Week6.pdf", width=6, height=6)
    par(mai=c(0.9, 0.9, 0.4, 0.4))
    inds6 <- all.comparisons$comparisonID == 6            #growth with priming
    ranks <- order(all.comparisons$r.boot.median[inds1])
    inds7 <- all.comparisons$comparisonID == 7            #growth on glucose with priming
    x.min <- min(c(all.comparisons$r.boot.CI.L[inds6], all.comparisons$r.boot.CI.L[inds7]), na.rm=TRUE)
    x.max <- max(c(all.comparisons$r.boot.CI.U[inds6], all.comparisons$r.boot.CI.U[inds7]), na.rm=TRUE)
    plot(all.comparisons$r.boot.median[inds7][ranks]~all.comparisons$r.boot.median[inds6][ranks], type="n", bty="l", xlab=expression(paste("Growth rate (r) with priming (day"^-1, ")", sep="")), ylab=expression(paste("Growth rate (r) on glucose with priming (day"^-1, ")", sep="")), xlim=c(x.min, x.max), ylim=c(x.min, x.max), main="Week 6")
    #plot(y=c(0,x.max), x=c(0,x.max), type="n", bty="l", xlab=expression(paste("Growth rate (r) with priming (day"^-1, ")", sep="")), ylab=expression(paste("Growth rate (r) on glucose with priming (day"^-1, ")", sep="")), xaxs="i", yaxs="i", main="Week 1")
    abline(h=0, col="red")
    abline(v=0, col="red")
    abline(a=0, b=1, lty=2, col="blue")
    points(x=all.comparisons$r.boot.median[inds6][ranks], y=all.comparisons$r.boot.median[inds7][ranks], pch=21, cex=0.6, bg="black")
    arrows(x0=all.comparisons$r.boot.CI.L[inds6][ranks], y0=all.comparisons$r.boot.median[inds7][ranks], x1=all.comparisons$r.boot.CI.U[inds6][ranks], y1=all.comparisons$r.boot.median[inds7][ranks], length=0, angle=90, code=3)
    arrows(x0=all.comparisons$r.boot.median[inds6][ranks], y0=all.comparisons$r.boot.CI.L[inds7][ranks], x1=all.comparisons$r.boot.median[inds6][ranks], y1=all.comparisons$r.boot.CI.U[inds7][ranks], length=0, angle=90, code=3)

    dev.off()
  #}




#_0g_Plot primed-glucose growth rates (r determined by 13C; units=day-1) vs. primed growth rates (r determined by 18O; units=day-1):
#(taxa colored by phylum)
  #{
  #Week 1:
    dev.off()
    # dev.new(width=6, height=6)
    pdf(file="qSIP_output/Figures/BH_GlucGrowth_v_Growth_Phyla_Week1.pdf", width=6, height=6)
    par(mai=c(0.9, 0.9, 0.4, 0.4))
    phyla <- levels(data.melted$phylum)
    phyla.cols <- rainbow(length(phyla))
    inds2 <- all.comparisons$comparisonID == 2            #growth with priming
    ranks <- order(all.comparisons$r.boot.median[inds1])
    tax.order <- data.frame(axis.loc=seq(1, N.taxa), ranks)
    inds3 <- all.comparisons$comparisonID == 3            #growth on glucose with priming
    x.min <- min(c(all.comparisons$r.boot.CI.L[inds2], all.comparisons$r.boot.CI.L[inds3]), na.rm=TRUE)
    x.max <- max(c(all.comparisons$r.boot.CI.U[inds2], all.comparisons$r.boot.CI.U[inds3]), na.rm=TRUE)
    plot(all.comparisons$r.boot.median[inds3][ranks]~all.comparisons$r.boot.median[inds2][ranks], type="n", bty="l", xlab=expression(paste("Growth rate (r) with priming (day"^-1, ")", sep="")), ylab=expression(paste("Growth rate (r) on glucose with priming (day"^-1, ")", sep="")), xlim=c(x.min, x.max), ylim=c(x.min, x.max), main="Week 1")
    #plot(y=c(0,x.max), x=c(0,x.max), type="n", bty="l", xlab=expression(paste("Growth rate (r) with priming (day"^-1, ")", sep="")), ylab=expression(paste("Growth rate (r) on glucose with priming (day"^-1, ")", sep="")), xaxs="i", yaxs="i", main="Week 1")
    abline(h=0, col="red")
    abline(v=0, col="red")
    abline(a=0, b=1, lty=2, col="blue")
    legend(x=0, y=1.05*x.max, legend=phyla, bty="n", text.col=phyla.cols, cex=0.6)
    for (p in 1:length(phyla)){
      comparison.x <- all.comparisons[inds2,]
      comp.ranked.x <- comparison.x[tax.order$ranks,]
      mids.x <- comp.ranked.x$r.boot.median[as.character(comp.ranked.x$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
      lowers.x <- comp.ranked.x$r.boot.CI.L[as.character(comp.ranked.x$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
      uppers.x <- comp.ranked.x$r.boot.CI.U[as.character(comp.ranked.x$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
      comparison.y <- all.comparisons[inds3,]
      comp.ranked.y <- comparison.y[tax.order$ranks,]
      mids.y <- comp.ranked.y$r.boot.median[as.character(comp.ranked.y$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
      lowers.y <- comp.ranked.y$r.boot.CI.L[as.character(comp.ranked.y$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
      uppers.y <- comp.ranked.y$r.boot.CI.U[as.character(comp.ranked.y$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
      points(x=mids.x, y=mids.y, pch=21, cex=0.6, col=phyla.cols[p], bg=phyla.cols[p])
      arrows(x0=lowers.x, y0=mids.y, x1=uppers.x, y1=mids.y, length=0, angle=90, code=3, col=phyla.cols[p])
      arrows(x0=mids.x, y0=lowers.y, x1=mids.x, y1=uppers.y, length=0, angle=90, code=3, col=phyla.cols[p])
    }

    dev.off()


  #Week 6:
    dev.off()
    # dev.new(width=6, height=6)
    pdf(file="qSIP_output/Figures/BH_GlucGrowth_v_Growth_Phyla_Week6.pdf", width=6, height=6)
    par(mai=c(0.9, 0.9, 0.4, 0.4))
    phyla <- levels(data.melted$phylum)
    phyla.cols <- rainbow(length(phyla))
    inds6 <- all.comparisons$comparisonID == 6            #growth with priming
    ranks <- order(all.comparisons$r.boot.median[inds1])
    tax.order <- data.frame(axis.loc=seq(1, N.taxa), ranks)
    inds7 <- all.comparisons$comparisonID == 7            #growth on glucose with priming
    x.min <- min(c(all.comparisons$r.boot.CI.L[inds6], all.comparisons$r.boot.CI.L[inds7]), na.rm=TRUE)
    x.max <- max(c(all.comparisons$r.boot.CI.U[inds6], all.comparisons$r.boot.CI.U[inds7]), na.rm=TRUE)
    plot(all.comparisons$r.boot.median[inds7][ranks]~all.comparisons$r.boot.median[inds6][ranks], type="n", bty="l", xlab=expression(paste("Growth rate (r) with priming (day"^-1, ")", sep="")), ylab=expression(paste("Growth rate (r) on glucose with priming (day"^-1, ")", sep="")), xlim=c(x.min, x.max), ylim=c(x.min, x.max), main="Week 6")
    #plot(y=c(0,x.max), x=c(0,x.max), type="n", bty="l", xlab=expression(paste("Growth rate (r) with priming (day"^-1, ")", sep="")), ylab=expression(paste("Growth rate (r) on glucose with priming (day"^-1, ")", sep="")), xaxs="i", yaxs="i", main="Week 6")
    abline(h=0, col="red")
    abline(v=0, col="red")
    abline(a=0, b=1, lty=2, col="blue")
    legend(x=0.18, y=1.05*x.max, legend=phyla, bty="n", text.col=phyla.cols, cex=0.6)
    for (p in 1:length(phyla)){
      comparison.x <- all.comparisons[inds6,]
      comp.ranked.x <- comparison.x[tax.order$ranks,]
      mids.x <- comp.ranked.x$r.boot.median[as.character(comp.ranked.x$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
      lowers.x <- comp.ranked.x$r.boot.CI.L[as.character(comp.ranked.x$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
      uppers.x <- comp.ranked.x$r.boot.CI.U[as.character(comp.ranked.x$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
      comparison.y <- all.comparisons[inds7,]
      comp.ranked.y <- comparison.y[tax.order$ranks,]
      mids.y <- comp.ranked.y$r.boot.median[as.character(comp.ranked.y$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
      lowers.y <- comp.ranked.y$r.boot.CI.L[as.character(comp.ranked.y$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
      uppers.y <- comp.ranked.y$r.boot.CI.U[as.character(comp.ranked.y$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
      points(x=mids.x, y=mids.y, pch=21, cex=0.6, col=phyla.cols[p], bg=phyla.cols[p])
      arrows(x0=lowers.x, y0=mids.y, x1=uppers.x, y1=mids.y, length=0, angle=90, code=3, col=phyla.cols[p])
      arrows(x0=mids.x, y0=lowers.y, x1=mids.x, y1=uppers.y, length=0, angle=90, code=3, col=phyla.cols[p])
    }

    dev.off()
  #}
    



#_0h_Calculate & plot the fraction of the growth rate with priming (r determined by 18O) attributable to glucose (r determined by 13C):
  #{
    #Week 1:
      dev.off()
      pdf(file="qSIP_output/Figures/BH_FractionGrowthFromGlucose_Week1.pdf", width=6, height=6)
      inds2 <- all.comparisons$comparisonID == 2  #with priming
      inds3 <- all.comparisons$comparisonID == 3  #on glucose with priming
      plot(y=(all.comparisons$r.boot.median[inds3] / all.comparisons$r.boot.median[inds2]), x=1:N.taxa, bty="l", type="p", pch=21, bg="black", xlab="Taxon", ylab="r on glucose with priming / r with priming", main="Week 1")
      dev.off()
    #Week 6:
      dev.off()
      pdf(file="qSIP_output/Figures/BH_FractionGrowthFromGlucose_Week6.pdf", width=6, height=6)
      inds6 <- all.comparisons$comparisonID == 6  #with priming
      inds7 <- all.comparisons$comparisonID == 7  #on glucose with priming
      plot(y=(all.comparisons$r.boot.median[inds7] / all.comparisons$r.boot.median[inds6]), x=1:N.taxa, bty="l", type="p", pch=21, bg="black", xlab="Taxon", ylab="r on glucose with priming / r with priming", main="Week 6")
      dev.off()

    #Weeks 1 & 6 plotted on same scale:
      dev.off()
      # dev.new(width=12, height=6)
      pdf(file="qSIP_output/Figures/BH_FractionGrowthFromGlucose_Week1+6.pdf", width=12, height=6)
      par(mfrow=c(1,2))
      inds2 <- all.comparisons$comparisonID == 2  #with priming; Week 1
      inds3 <- all.comparisons$comparisonID == 3  #on glucose with priming; Week 1
      inds6 <- all.comparisons$comparisonID == 6  #with priming; Week 6
      inds7 <- all.comparisons$comparisonID == 7  #on glucose with priming; Week 6
      y.min <- min(c((all.comparisons$r.boot.median[inds3] / all.comparisons$r.boot.median[inds2]), (all.comparisons$r.boot.median[inds7] / all.comparisons$r.boot.median[inds6])), na.rm=TRUE)
      y.max <- max(c((all.comparisons$r.boot.median[inds3] / all.comparisons$r.boot.median[inds2]), (all.comparisons$r.boot.median[inds7] / all.comparisons$r.boot.median[inds6])), na.rm=TRUE)
      plot(y=(all.comparisons$r.boot.median[inds3] / all.comparisons$r.boot.median[inds2]), x=1:N.taxa, bty="l", type="p", pch=21, bg="black", ylim=c(y.min, y.max), xlab="Taxon", ylab="r on glucose with priming / r with priming", main="Week 1")
      plot(y=(all.comparisons$r.boot.median[inds7] / all.comparisons$r.boot.median[inds6]), x=1:N.taxa, bty="l", type="p", pch=21, bg="black", ylim=c(y.min, y.max), xlab="Taxon", ylab="r on glucose with priming / r with priming", main="Week 6")

      dev.off()
  #}




#_1a_Create a refined plot containing both wad.diff and ape on same graph using two axes (Weeks 1 & 6):
  #{
  #first, determine the empirical relationship between wad.diff and ape for 13C and 18O comparisons:
    plot(y=all.comparisons$wad.diff.obs, x=all.comparisons$ape.obs) #includes both 18O and 13C comparisons
      #18O-only data points:
      plot(y=all.comparisons$wad.diff.obs[all.comparisons$comparisonID %in% c(1,2,5,6)], x=all.comparisons$ape.obs[all.comparisons$comparisonID %in% c(1,2,5,6)]) #only 18O comparisons
      ape.obs.18O <- all.comparisons$ape.obs[all.comparisons$comparisonID %in% c(1,2,5,6)]
      wad.diff.obs.18O <- all.comparisons$wad.diff.obs[all.comparisons$comparisonID %in% c(1,2,5,6)]
      summary(lm(wad.diff.obs.18O~ape.obs.18O))
      coef(lm(wad.diff.obs.18O~ape.obs.18O))
      #13C-only data points:
      plot(y=all.comparisons$wad.diff.obs[all.comparisons$comparisonID %in% c(3,7)], x=all.comparisons$ape.obs[all.comparisons$comparisonID %in% c(3,7)]) #only 13C comparisons
      ape.obs.13C <- all.comparisons$ape.obs[all.comparisons$comparisonID %in% c(3,7)]
      wad.diff.obs.13C <- all.comparisons$wad.diff.obs[all.comparisons$comparisonID %in% c(3,7)]
      summary(lm(wad.diff.obs.13C~ape.obs.13C))
      coef(lm(wad.diff.obs.13C~ape.obs.13C))

  #plot difference in WADS & APE on same graph with two axes:
  #(NOTE THAT THERE IS ONE TAXON WHICH HAS PHYLUM == NA AND IS THUS NOT PLOTTED ON THIS FIGURE)
    dev.off()
    dimensions <- c((18.3/2.54), (13.725/2.54))
    # dev.new(width=dimensions[1], height=dimensions[2])
    pdf(file="qSIP_output/Figures/BH_Pretty_waddiff&ape.pdf", width=dimensions[1], height=dimensions[2])
    layout(mat=matrix(c(1,2,3,4,5,6), 2, 3, byrow=TRUE), widths=c(1/3,1/3,1/3), heights=c(1/2,1/2))
    phyla <- levels(data.melted$phylum)
    phyla.cols <- rainbow(length(phyla))
    # x.min <- min(all.comparisons$wad.diff.boot.CI.L, na.rm=TRUE)
    # x.max <- max(all.comparisons$wad.diff.boot.CI.U, na.rm=TRUE)
    x.min <- -0.02
    x.max <- 0.0425
    AT.wad.diff <- seq(-0.02, 0.04, 0.02)
    AT.ape.18O <- format(seq(-0.4, 0.8, 0.2), digits=1)
    AT.18O <- (as.numeric(AT.ape.18O)*coef(lm(wad.diff.obs.18O~ape.obs.18O))[2])+coef(lm(wad.diff.obs.18O~ape.obs.18O))[1]
    AT.ape.13C <- format(seq(-0.4, 0.8, 0.2), digits=1)
    AT.13C <- (as.numeric(AT.ape.13C)*coef(lm(wad.diff.obs.13C~ape.obs.13C))[2])+coef(lm(wad.diff.obs.13C~ape.obs.13C))[1]

    wad.diff.add.panel <- function(comparisonID, AT.iso, y.axis){
      inds <- all.comparisons$comparisonID == comparisonID
      ranks <- order(all.comparisons$wad.diff.boot.median[inds])
      tax.order <- data.frame(axis.loc=seq(1, N.taxa), ranks)
      plot(1:N.taxa~all.comparisons$wad.diff.boot.median[inds][tax.order$ranks], type="n", bty="n", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min, x.max), main="")
      for (p in 1:length(phyla)){
        current.comparison <- all.comparisons[inds,]
        curr.comp.ranked <- current.comparison[tax.order$ranks,]
        tax.nums <- tax.order$axis.loc[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        mids <- curr.comp.ranked$wad.diff.boot.median[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        lowers <- curr.comp.ranked$wad.diff.boot.CI.L[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        uppers <- curr.comp.ranked$wad.diff.boot.CI.U[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        points(x=mids, y=tax.nums, pch=21, cex=0.2, col=phyla.cols[p], bg=phyla.cols[p])
        arrows(x0=lowers, y0=tax.nums, x1=uppers, y1=tax.nums, length=0, angle=90, code=3, col=phyla.cols[p], lwd=0.5)
      }
      axis(side=1, at=c(x.min, x.max)*1000, labels=FALSE, tck=0)
      axis(side=3, at=c(x.min, x.max)*1000, labels=FALSE, tck=0)
      abline(v=0, col="black", lty=2)
      par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
      axis(side=1, at=AT.wad.diff, labels=FALSE, tck=-0.015, las=1, cex.axis=0.6)
      par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
      axis(side=3, at=AT.iso, labels=FALSE, tck=-0.015, las=1, cex.axis=0.6)
      if (y.axis == TRUE){
        par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
        axis(side=2, at=c(-1000, 1000), labels=FALSE, tck=0)
      }
    }

    par(mai=c(0.1,0.25,0.44,0.05), cex=1, mex=0.75)
    wad.diff.add.panel(comparisonID=1, AT.iso=AT.18O, y.axis=TRUE)
    par(mgp=c(3,0.25,0))				  #set spacing for the x-axis labels (the second value)
    axis(side=3, at=AT.18O, labels=AT.ape.18O, tck=-0.015, las=1, cex.axis=0.6)
    par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
    axis(side=2, at=1:N.taxa, labels=rep(NA, N.taxa), tck=-0.015, las=1, cex.axis=0.6)
    mtext("Ranked taxon", side=2, line=0.45, cex=0.75)
    text(x=x.min*1.03, y=N.taxa*1.01, labels="Week 1", adj=c(0,1), cex=0.6)
    text(x=x.min*1.03, y=N.taxa*0.95, labels="control", adj=c(0,1), cex=0.6)
    text(x=x.min*1.03, y=N.taxa*0.89, labels=expression(paste(""^18, "O", sep="")), adj=c(0,1), cex=0.6)
    par(xpd=NA)
    legend(x=0.038, y=0.90*N.taxa, xjust=0.5, yjust=1, legend=phyla, bty="n", text.col=phyla.cols, cex=0.6)
    par(xpd=FALSE)
    par(mai=c(0.1,0.25,0.44,0.05), cex=1, mex=0.75)
    wad.diff.add.panel(comparisonID=2, AT.iso=AT.18O, y.axis=FALSE)
    par(mgp=c(3,0.25,0))				  #set spacing for the x-axis labels (the second value)
    axis(side=3, at=AT.18O, labels=AT.ape.18O, tck=-0.015, las=1, cex.axis=0.6)
    text(x=x.min*1.03, y=N.taxa*1.01, labels="Week 1", adj=c(0,1), cex=0.6)
    text(x=x.min*1.03, y=N.taxa*0.95, labels="added carbon", adj=c(0,1), cex=0.6)
    text(x=x.min*1.03, y=N.taxa*0.89, labels=expression(paste(""^18, "O", sep="")), adj=c(0,1), cex=0.6)
    par(mai=c(0.1,0.25,0.44,0.05), cex=1, mex=0.75)
    mtext(expression(paste("Excess atom fraction ("^18, "O or "^13, "C)", sep="")), side=3, line=1.3, cex=0.75)
    wad.diff.add.panel(comparisonID=3, AT.iso=AT.13C, y.axis=FALSE)
    par(mgp=c(3,0.25,0))				  #set spacing for the x-axis labels (the second value)
    axis(side=3, at=AT.13C, labels=AT.ape.13C, tck=-0.015, las=1, cex.axis=0.6)
    text(x=x.min*1.03, y=N.taxa*1.01, labels="Week 1", adj=c(0,1), cex=0.6)
    text(x=x.min*1.03, y=N.taxa*0.95, labels="added carbon", adj=c(0,1), cex=0.6)
    text(x=x.min*1.03, y=N.taxa*0.89, labels=expression(paste(""^13, "C", sep="")), adj=c(0,1), cex=0.6)
    par(mai=c(0.44,0.25,0.1,0.05), cex=1, mex=0.75)
    wad.diff.add.panel(comparisonID=5, AT.iso=AT.18O, y.axis=TRUE)
    par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
    axis(side=2, at=1:N.taxa, labels=rep(NA, N.taxa), tck=-0.015, las=1, cex.axis=0.6)
    mtext("Ranked taxon", side=2, line=0.45, cex=0.75)
    text(x=x.min*1.03, y=N.taxa*1.01, labels="Week 6", adj=c(0,1), cex=0.6)
    text(x=x.min*1.03, y=N.taxa*0.95, labels="control", adj=c(0,1), cex=0.6)
    text(x=x.min*1.03, y=N.taxa*0.89, labels=expression(paste(""^18, "O", sep="")), adj=c(0,1), cex=0.6)
    par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
    axis(side=1, at=AT.wad.diff, labels=AT.wad.diff, tck=-0.015, las=1, cex.axis=0.6)
    par(mai=c(0.44,0.25,0.1,0.05), cex=1, mex=0.75)
    wad.diff.add.panel(comparisonID=6, AT.iso=AT.18O, y.axis=FALSE)
    text(x=x.min*1.03, y=N.taxa*1.01, labels="Week 6", adj=c(0,1), cex=0.6)
    text(x=x.min*1.03, y=N.taxa*0.95, labels="added carbon", adj=c(0,1), cex=0.6)
    text(x=x.min*1.03, y=N.taxa*0.89, labels=expression(paste(""^18, "O", sep="")), adj=c(0,1), cex=0.6)
    mtext(expression(paste("Difference in weighted average density (g mL"^-1, ")", sep="")), side=1, line=1.4, cex=0.75)
    par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
    axis(side=1, at=AT.wad.diff, labels=AT.wad.diff, tck=-0.015, las=1, cex.axis=0.6)
    par(mai=c(0.44,0.25,0.1,0.05), cex=1, mex=0.75)
    wad.diff.add.panel(comparisonID=7, AT.iso=AT.13C, y.axis=FALSE)
    text(x=x.min*1.03, y=N.taxa*1.01, labels="Week 6", adj=c(0,1), cex=0.6)
    text(x=x.min*1.03, y=N.taxa*0.95, labels="added carbon", adj=c(0,1), cex=0.6)
    text(x=x.min*1.03, y=N.taxa*0.89, labels=expression(paste(""^13, "C", sep="")), adj=c(0,1), cex=0.6)
    par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
    axis(side=1, at=AT.wad.diff, labels=AT.wad.diff, tck=-0.015, las=1, cex.axis=0.6)
    rm(wad.diff.add.panel)
    par(mfrow=c(1,1))

    dev.off()
  #}
    

#Import colors for phyla for use below in 1b & 5b:
  cols.for.phyla <- read.table("qSIP_data/BH_PhylaColors.txt", header=T, sep="\t", comment.char="")


#_1b_Create a refined plot containing both wad.diff and ape on same graph using two axes (Week 1 only):
  #{
  #first, determine the empirical relationship between wad.diff and ape for 13C and 18O comparisons:
    plot(y=all.comparisons$wad.diff.obs, x=all.comparisons$ape.obs) #includes both 18O and 13C comparisons
      #18O-only data points:
      plot(y=all.comparisons$wad.diff.obs[all.comparisons$comparisonID %in% c(1,2,5,6)], x=all.comparisons$ape.obs[all.comparisons$comparisonID %in% c(1,2,5,6)]) #only 18O comparisons
      ape.obs.18O <- all.comparisons$ape.obs[all.comparisons$comparisonID %in% c(1,2,5,6)]
      wad.diff.obs.18O <- all.comparisons$wad.diff.obs[all.comparisons$comparisonID %in% c(1,2,5,6)]
      summary(lm(wad.diff.obs.18O~ape.obs.18O))
      coef(lm(wad.diff.obs.18O~ape.obs.18O))
      #13C-only data points:
      plot(y=all.comparisons$wad.diff.obs[all.comparisons$comparisonID %in% c(3,7)], x=all.comparisons$ape.obs[all.comparisons$comparisonID %in% c(3,7)]) #only 13C comparisons
      ape.obs.13C <- all.comparisons$ape.obs[all.comparisons$comparisonID %in% c(3,7)]
      wad.diff.obs.13C <- all.comparisons$wad.diff.obs[all.comparisons$comparisonID %in% c(3,7)]
      summary(lm(wad.diff.obs.13C~ape.obs.13C))
      coef(lm(wad.diff.obs.13C~ape.obs.13C))

  #plot difference in WADS & APE on same graph with two axes:
  #(NOTE THAT THERE IS ONE TAXON WHICH HAS PHYLUM == NA AND IS THUS NOT PLOTTED ON THIS FIGURE)
    dev.off()
    dimensions <- c((18.3/2.54), (7.726101/2.54))
    # dev.new(width=dimensions[1], height=dimensions[2])
    pdf(file="qSIP_output/Figures/BH_Pretty_waddiff&ape_Week1_colorbrewer.pdf", width=dimensions[1], height=dimensions[2])
    layout(mat=matrix(c(1,2,3), 1, 3, byrow=TRUE), widths=c(1/3,1/3,1/3))
    phyla <- levels(data.melted$phylum)
    # phyla.cols <- rainbow(length(phyla))
    phyla.cols <- cols.for.phyla$col
    # x.min <- min(all.comparisons$wad.diff.boot.CI.L[all.comparisons$comparisonID %in% c(1,2,3)], na.rm=TRUE)
    # x.max <- max(all.comparisons$wad.diff.boot.CI.U[all.comparisons$comparisonID %in% c(1,2,3)], na.rm=TRUE)
    x.min <- -0.02
    x.max <- 0.0425
    AT.wad.diff <- seq(-0.02, 0.04, 0.02)
    AT.ape.18O <- format(seq(-0.4, 0.8, 0.2), digits=1)
    AT.18O <- (as.numeric(AT.ape.18O)*coef(lm(wad.diff.obs.18O~ape.obs.18O))[2])+coef(lm(wad.diff.obs.18O~ape.obs.18O))[1]
    AT.ape.13C <- format(seq(-0.4, 0.8, 0.2), digits=1)
    AT.13C <- (as.numeric(AT.ape.13C)*coef(lm(wad.diff.obs.13C~ape.obs.13C))[2])+coef(lm(wad.diff.obs.13C~ape.obs.13C))[1]

    wad.diff.add.panel <- function(comparisonID, AT.iso, y.axis){
      inds <- all.comparisons$comparisonID == comparisonID
      ranks <- order(all.comparisons$wad.diff.boot.median[inds])
      tax.order <- data.frame(axis.loc=seq(1, N.taxa), ranks)
      plot(1:N.taxa~all.comparisons$wad.diff.boot.median[inds][tax.order$ranks], type="n", bty="n", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min, x.max), main="")
      for (p in 1:length(phyla)){
        current.comparison <- all.comparisons[inds,]
        curr.comp.ranked <- current.comparison[tax.order$ranks,]
        tax.nums <- tax.order$axis.loc[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        mids <- curr.comp.ranked$wad.diff.boot.median[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        lowers <- curr.comp.ranked$wad.diff.boot.CI.L[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        uppers <- curr.comp.ranked$wad.diff.boot.CI.U[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        points(x=mids, y=tax.nums, pch=21, cex=0.15, col=as.character(phyla.cols)[p], bg=as.character(phyla.cols)[p])
        arrows(x0=lowers, y0=tax.nums, x1=uppers, y1=tax.nums, length=0, angle=90, code=3, col=as.character(phyla.cols)[p], lwd=0.5)
      }
      axis(side=1, at=c(x.min, x.max)*1000, labels=FALSE, tck=0)
      axis(side=3, at=c(x.min, x.max)*1000, labels=FALSE, tck=0)
      abline(v=0, col="black", lty=2)
      if (y.axis == TRUE){
        par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
        axis(side=2, at=c(-1000, 1000), labels=FALSE, tck=0)
      }
    }

    par(mai=c(0.44,0.25,0.44,0.05), cex=1, mex=0.75)
    wad.diff.add.panel(comparisonID=1, AT.iso=AT.18O, y.axis=TRUE)
    par(mgp=c(3,0.25,0))				  #set spacing for the x-axis labels (the second value)
    axis(side=3, at=AT.18O, labels=AT.ape.18O, tck=-0.015, las=1, cex.axis=0.6)
    par(mgp=c(3,-0.08,0))				  #set spacing for the x-axis labels (the second value)
    axis(side=1, at=AT.wad.diff, labels=AT.wad.diff, tck=-0.015, las=1, cex.axis=0.6)
    par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
    axis(side=2, at=1:N.taxa, labels=rep(NA, N.taxa), tck=-0.015, las=1, cex.axis=0.6)
    mtext("Ranked taxon", side=2, line=0.45, cex=0.75)
    text(x=x.min*1.03, y=N.taxa*1.01, labels="control", adj=c(0,1), cex=0.6)
    text(x=x.min*1.03, y=N.taxa*0.95, labels=expression(paste(""^18, "O", sep="")), adj=c(0,1), cex=0.6)
    par(xpd=NA)
    legend(x=0.035, y=0.82*N.taxa, xjust=0.5, yjust=1, legend=phyla[1:ceiling(0.5*length(phyla))], bty="n", text.col=as.character(phyla.cols)[1:ceiling(0.5*length(phyla))], cex=0.4)
    legend(x=0.056, y=0.82*N.taxa, xjust=0.5, yjust=1, legend=phyla[(ceiling(0.5*length(phyla))+1):length(phyla)], bty="n", text.col=as.character(phyla.cols)[(ceiling(0.5*length(phyla))+1):length(phyla)], cex=0.4)
    par(xpd=FALSE)
    par(mai=c(0.44,0.25,0.44,0.05), cex=1, mex=0.75)
    wad.diff.add.panel(comparisonID=2, AT.iso=AT.18O, y.axis=FALSE)
    par(mgp=c(3,0.25,0))				  #set spacing for the x-axis labels (the second value)
    axis(side=3, at=AT.18O, labels=AT.ape.18O, tck=-0.015, las=1, cex.axis=0.6)
    par(mgp=c(3,-0.08,0))				  #set spacing for the x-axis labels (the second value)
    axis(side=1, at=AT.wad.diff, labels=AT.wad.diff, tck=-0.015, las=1, cex.axis=0.6)
    text(x=x.min*1.03, y=N.taxa*1.01, labels="added carbon", adj=c(0,1), cex=0.6)
    text(x=x.min*1.03, y=N.taxa*0.95, labels=expression(paste(""^18, "O", sep="")), adj=c(0,1), cex=0.6)
    par(mai=c(0.44,0.25,0.44,0.05), cex=1, mex=0.75)
    mtext(expression(paste("Atom fraction excess ("^18, "O or "^13, "C)", sep="")), side=3, line=1.3, cex=0.75)
    mtext(expression(paste("Difference in density (g cm"^-3, ")", sep="")), side=1, line=1.4, cex=0.75)
    wad.diff.add.panel(comparisonID=3, AT.iso=AT.13C, y.axis=FALSE)
    par(mgp=c(3,0.25,0))				  #set spacing for the x-axis labels (the second value)
    axis(side=3, at=AT.13C, labels=AT.ape.13C, tck=-0.015, las=1, cex.axis=0.6)
    par(mgp=c(3,-0.08,0))				  #set spacing for the x-axis labels (the second value)
    axis(side=1, at=AT.wad.diff, labels=AT.wad.diff, tck=-0.015, las=1, cex.axis=0.6)
    text(x=x.min*1.03, y=N.taxa*1.01, labels="added carbon", adj=c(0,1), cex=0.6)
    text(x=x.min*1.03, y=N.taxa*0.95, labels=expression(paste(""^13, "C", sep="")), adj=c(0,1), cex=0.6)
    rm(wad.diff.add.panel)
    par(mfrow=c(1,1))

    dev.off()
  #}


#_1c_Create a refined plot containing both wad.diff and ape on same graph using two axes (Week 1 only) with phyla sorted:
  #{
  #first, determine the empirical relationship between wad.diff and ape for 13C and 18O comparisons:
    plot(y=all.comparisons$wad.diff.obs, x=all.comparisons$ape.obs) #includes both 18O and 13C comparisons
      #18O-only data points:
      plot(y=all.comparisons$wad.diff.obs[all.comparisons$comparisonID %in% c(1,2,5,6)], x=all.comparisons$ape.obs[all.comparisons$comparisonID %in% c(1,2,5,6)]) #only 18O comparisons
      ape.obs.18O <- all.comparisons$ape.obs[all.comparisons$comparisonID %in% c(1,2,5,6)]
      wad.diff.obs.18O <- all.comparisons$wad.diff.obs[all.comparisons$comparisonID %in% c(1,2,5,6)]
      summary(lm(wad.diff.obs.18O~ape.obs.18O))
      coef(lm(wad.diff.obs.18O~ape.obs.18O))
      #13C-only data points:
      plot(y=all.comparisons$wad.diff.obs[all.comparisons$comparisonID %in% c(3,7)], x=all.comparisons$ape.obs[all.comparisons$comparisonID %in% c(3,7)]) #only 13C comparisons
      ape.obs.13C <- all.comparisons$ape.obs[all.comparisons$comparisonID %in% c(3,7)]
      wad.diff.obs.13C <- all.comparisons$wad.diff.obs[all.comparisons$comparisonID %in% c(3,7)]
      summary(lm(wad.diff.obs.13C~ape.obs.13C))
      coef(lm(wad.diff.obs.13C~ape.obs.13C))

  #plot difference in WADS & APE on same graph with two axes:
  #(NOTE THAT THERE IS ONE TAXON WHICH HAS PHYLUM == NA AND IS THUS NOT PLOTTED ON THIS FIGURE)
    dev.off()
    dimensions <- c((18.3/2.54), (7.726101/2.54))
    # dev.new(width=dimensions[1], height=dimensions[2])
    pdf(file="qSIP_output/Figures/BH_Pretty_waddiff&ape_Week1_phyla_sorted_colorbrewer.pdf", width=dimensions[1], height=dimensions[2])
    layout(mat=matrix(c(1,2,3), 1, 3, byrow=TRUE), widths=c(1/3,1/3,1/3))
    phyla <- levels(data.melted$phylum)
    phyla <- sort(phyla, decreasing=TRUE)
    # phyla.cols <- rainbow(length(phyla))
    phyla.cols <- cols.for.phyla$col
    phyla.cols <- phyla.cols[length(phyla.cols):1]
    # x.min <- min(all.comparisons$wad.diff.boot.CI.L[all.comparisons$comparisonID %in% c(1,2,3)], na.rm=TRUE)
    # x.max <- max(all.comparisons$wad.diff.boot.CI.U[all.comparisons$comparisonID %in% c(1,2,3)], na.rm=TRUE)
    x.min <- -0.02
    x.max <- 0.0425
    AT.wad.diff <- seq(-0.02, 0.04, 0.02)
    AT.ape.18O <- format(seq(-0.4, 0.8, 0.2), digits=1)
    AT.18O <- (as.numeric(AT.ape.18O)*coef(lm(wad.diff.obs.18O~ape.obs.18O))[2])+coef(lm(wad.diff.obs.18O~ape.obs.18O))[1]
    AT.ape.13C <- format(seq(-0.4, 0.8, 0.2), digits=1)
    AT.13C <- (as.numeric(AT.ape.13C)*coef(lm(wad.diff.obs.13C~ape.obs.13C))[2])+coef(lm(wad.diff.obs.13C~ape.obs.13C))[1]

    wad.diff.add.panel <- function(comparisonID, AT.iso, y.axis, label=FALSE){
      inds <- all.comparisons$comparisonID == comparisonID
      ranks <- order(all.comparisons$wad.diff.boot.median[inds])
      tax.order <- data.frame(axis.loc=seq(1, N.taxa), ranks)
      plot(1:N.taxa~all.comparisons$wad.diff.boot.median[inds][tax.order$ranks], type="n", bty="n", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min, x.max), main="")
      counter <- 1
      for (p in 1:length(phyla)){
        current.comparison <- all.comparisons[inds,]
        curr.comp.ranked <- current.comparison[tax.order$ranks,]
        # tax.nums <- tax.order$axis.loc[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        mids <- curr.comp.ranked$wad.diff.boot.median[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        lowers <- curr.comp.ranked$wad.diff.boot.CI.L[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        uppers <- curr.comp.ranked$wad.diff.boot.CI.U[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        tax.nums <- counter:(counter+length(mids)-1)
        counter <- counter+length(mids)
        arrows(x0=lowers, y0=tax.nums, x1=uppers, y1=tax.nums, length=0, angle=90, code=3, col=as.character(phyla.cols)[p], lwd=0.5)
        points(x=mids, y=tax.nums, pch=21, cex=0.15, lwd=0.125, col="gray30", bg=as.character(phyla.cols)[p])
        if (label == TRUE){
          par(xpd=NA)
          if(is.element(p, c(1,17))){
            if (p == 17){   #shift Elusimicrobia down one unit
              text(x=x.max, y=mean(tax.nums)-1, labels=phyla[p], adj=c(0,0.5), pos=4, offset=-3, cex=0.35, col=as.character(phyla.cols)[p])              
            }
            else{
              text(x=x.max, y=mean(tax.nums), labels=phyla[p], adj=c(0,0.5), pos=4, offset=-3, cex=0.35, col=as.character(phyla.cols)[p])
            }
          }
          else if (is.element(p, c(2,5,6,7,11,14,19,20,24,26,27,28,29))){
            text(x=x.max, y=mean(tax.nums), labels=phyla[p], adj=c(0,0.5), pos=4, offset=-2.25, cex=0.35, col=as.character(phyla.cols)[p])
          }
          else if (is.element(p, c(3,25))){
            text(x=x.max, y=mean(tax.nums), labels=phyla[p], adj=c(0,0.5), pos=4, offset=-1.5, cex=0.35, col=as.character(phyla.cols)[p])
          }
          else if (is.element(p, c(12,16,21))){
            text(x=x.max, y=mean(tax.nums), labels=phyla[p], adj=c(0,0.5), pos=4, offset=-0.75, cex=0.35, col=as.character(phyla.cols)[p])
          }
          else if (is.element(p, c(4,9,13,18))){
            if (p == 18){   #shift Cyanobacteria up one unit
              text(x=x.max, y=mean(tax.nums)+1, labels=phyla[p], adj=c(0,0.5), pos=4, offset=0, cex=0.35, col=as.character(phyla.cols)[p])              
            }
            else{
              text(x=x.max, y=mean(tax.nums), labels=phyla[p], adj=c(0,0.5), pos=4, offset=0, cex=0.35, col=as.character(phyla.cols)[p])
            }
          }
          else if (is.element(p, c(8,15,22))){
            text(x=x.max, y=mean(tax.nums), labels=phyla[p], adj=c(0,0.5), pos=4, offset=0.75, cex=0.35, col=as.character(phyla.cols)[p])
          }
          else if (is.element(p, c(10,23))){
            text(x=x.max, y=mean(tax.nums), labels=phyla[p], adj=c(0,0.5), pos=4, offset=1.5, cex=0.35, col=as.character(phyla.cols)[p])
          }
          else if (is.element(p, c(99))){
            text(x=x.max, y=mean(tax.nums), labels=phyla[p], adj=c(0,0.5), pos=4, offset=2.25, cex=0.35, col=as.character(phyla.cols)[p])
          }
          else if (is.element(p, c(99))){
            text(x=x.max, y=mean(tax.nums), labels=phyla[p], adj=c(0,0.5), pos=4, offset=3, cex=0.35, col=as.character(phyla.cols)[p])
          }
          par(xpd=FALSE)
        }
      }
      axis(side=1, at=c(x.min, x.max)*1000, labels=FALSE, tck=0)
      axis(side=3, at=c(x.min, x.max)*1000, labels=FALSE, tck=0)
      abline(v=0, col="black", lty=2)
      if (y.axis == TRUE){
        par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
        axis(side=2, at=c(-1000, 1000), labels=FALSE, tck=0)
      }
    }

    par(mai=c(0.44,0.25,0.44,0.05), cex=1, mex=0.75)
    wad.diff.add.panel(comparisonID=1, AT.iso=AT.18O, y.axis=TRUE, label=TRUE)
    par(mgp=c(3,0.25,0))				  #set spacing for the x-axis labels (the second value)
    axis(side=3, at=AT.18O, labels=AT.ape.18O, tck=-0.015, las=1, cex.axis=0.6)
    par(mgp=c(3,-0.08,0))				  #set spacing for the x-axis labels (the second value)
    axis(side=1, at=AT.wad.diff, labels=AT.wad.diff, tck=-0.015, las=1, cex.axis=0.6)
    par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
    axis(side=2, at=1:N.taxa, labels=rep(NA, N.taxa), tck=-0.015, las=1, cex.axis=0.6)
    mtext("Ranked taxon", side=2, line=0.45, cex=0.75)
    text(x=x.min*1.03, y=N.taxa*1.01, labels="control", adj=c(0,1), cex=0.6)
    text(x=x.min*1.03, y=N.taxa*0.95, labels=expression(paste(""^18, "O", sep="")), adj=c(0,1), cex=0.6)
    # par(xpd=NA)
    # legend(x=0.035, y=0.82*N.taxa, xjust=0.5, yjust=1, legend=phyla[1:ceiling(0.5*length(phyla))], bty="n", text.col=as.character(phyla.cols)[1:ceiling(0.5*length(phyla))], cex=0.4)
    # legend(x=0.056, y=0.82*N.taxa, xjust=0.5, yjust=1, legend=phyla[(ceiling(0.5*length(phyla))+1):length(phyla)], bty="n", text.col=as.character(phyla.cols)[(ceiling(0.5*length(phyla))+1):length(phyla)], cex=0.4)
    # par(xpd=FALSE)
    par(mai=c(0.44,0.25,0.44,0.05), cex=1, mex=0.75)
    wad.diff.add.panel(comparisonID=2, AT.iso=AT.18O, y.axis=FALSE, label=FALSE)
    par(mgp=c(3,0.25,0))				  #set spacing for the x-axis labels (the second value)
    axis(side=3, at=AT.18O, labels=AT.ape.18O, tck=-0.015, las=1, cex.axis=0.6)
    par(mgp=c(3,-0.08,0))				  #set spacing for the x-axis labels (the second value)
    axis(side=1, at=AT.wad.diff, labels=AT.wad.diff, tck=-0.015, las=1, cex.axis=0.6)
    text(x=x.min*1.03, y=N.taxa*1.01, labels="added carbon", adj=c(0,1), cex=0.6)
    text(x=x.min*1.03, y=N.taxa*0.95, labels=expression(paste(""^18, "O", sep="")), adj=c(0,1), cex=0.6)
    par(mai=c(0.44,0.25,0.44,0.05), cex=1, mex=0.75)
    mtext(expression(paste("Atom fraction excess ("^18, "O or "^13, "C)", sep="")), side=3, line=1.3, cex=0.75)
    mtext(expression(paste("Difference in density (g cm"^-3, ")", sep="")), side=1, line=1.4, cex=0.75)
    wad.diff.add.panel(comparisonID=3, AT.iso=AT.13C, y.axis=FALSE, label=FALSE)
    par(mgp=c(3,0.25,0))				  #set spacing for the x-axis labels (the second value)
    axis(side=3, at=AT.13C, labels=AT.ape.13C, tck=-0.015, las=1, cex.axis=0.6)
    par(mgp=c(3,-0.08,0))				  #set spacing for the x-axis labels (the second value)
    axis(side=1, at=AT.wad.diff, labels=AT.wad.diff, tck=-0.015, las=1, cex.axis=0.6)
    text(x=x.min*1.03, y=N.taxa*1.01, labels="added carbon", adj=c(0,1), cex=0.6)
    text(x=x.min*1.03, y=N.taxa*0.95, labels=expression(paste(""^13, "C", sep="")), adj=c(0,1), cex=0.6)
    rm(wad.diff.add.panel)
    par(mfrow=c(1,1))

    dev.off()
  #}


#_1d_Create a refined plot containing both wad.diff and ape on same graph using two axes (Week 1 only) with all taxa colored gray:
  #{
  #first, determine the empirical relationship between wad.diff and ape for 13C and 18O comparisons:
    plot(y=all.comparisons$wad.diff.obs, x=all.comparisons$ape.obs) #includes both 18O and 13C comparisons
      #18O-only data points:
      plot(y=all.comparisons$wad.diff.obs[all.comparisons$comparisonID %in% c(1,2,5,6)], x=all.comparisons$ape.obs[all.comparisons$comparisonID %in% c(1,2,5,6)]) #only 18O comparisons
      ape.obs.18O <- all.comparisons$ape.obs[all.comparisons$comparisonID %in% c(1,2,5,6)]
      wad.diff.obs.18O <- all.comparisons$wad.diff.obs[all.comparisons$comparisonID %in% c(1,2,5,6)]
      summary(lm(wad.diff.obs.18O~ape.obs.18O))
      coef(lm(wad.diff.obs.18O~ape.obs.18O))
      #13C-only data points:
      plot(y=all.comparisons$wad.diff.obs[all.comparisons$comparisonID %in% c(3,7)], x=all.comparisons$ape.obs[all.comparisons$comparisonID %in% c(3,7)]) #only 13C comparisons
      ape.obs.13C <- all.comparisons$ape.obs[all.comparisons$comparisonID %in% c(3,7)]
      wad.diff.obs.13C <- all.comparisons$wad.diff.obs[all.comparisons$comparisonID %in% c(3,7)]
      summary(lm(wad.diff.obs.13C~ape.obs.13C))
      coef(lm(wad.diff.obs.13C~ape.obs.13C))

  #plot difference in WADS & APE on same graph with two axes:
  #(NOTE THAT THERE IS ONE TAXON WHICH HAS PHYLUM == NA AND IS THUS NOT PLOTTED ON THIS FIGURE)
    dev.off()
    dimensions <- c((18.3/2.54), (7.726101/2.54))
    # dev.new(width=dimensions[1], height=dimensions[2])
    pdf(file="qSIP_output/Figures/BH_Pretty_waddiff&ape_Week1_b&w.pdf", width=dimensions[1], height=dimensions[2])
    layout(mat=matrix(c(1,2,3), 1, 3, byrow=TRUE), widths=c(1/3,1/3,1/3))
    phyla <- levels(data.melted$phylum)
    # phyla.cols <- rainbow(length(phyla))
    # phyla.cols <- cols.for.phyla$col
    # x.min <- min(all.comparisons$wad.diff.boot.CI.L[all.comparisons$comparisonID %in% c(1,2,3)], na.rm=TRUE)
    # x.max <- max(all.comparisons$wad.diff.boot.CI.U[all.comparisons$comparisonID %in% c(1,2,3)], na.rm=TRUE)
    x.min <- -0.02
    x.max <- 0.0425
    AT.wad.diff <- seq(-0.02, 0.04, 0.02)
    AT.ape.18O <- format(seq(-0.4, 0.8, 0.2), digits=1)
    AT.18O <- (as.numeric(AT.ape.18O)*coef(lm(wad.diff.obs.18O~ape.obs.18O))[2])+coef(lm(wad.diff.obs.18O~ape.obs.18O))[1]
    AT.ape.13C <- format(seq(-0.4, 0.8, 0.2), digits=1)
    AT.13C <- (as.numeric(AT.ape.13C)*coef(lm(wad.diff.obs.13C~ape.obs.13C))[2])+coef(lm(wad.diff.obs.13C~ape.obs.13C))[1]

    wad.diff.add.panel <- function(comparisonID, AT.iso, y.axis){
      inds <- all.comparisons$comparisonID == comparisonID
      ranks <- order(all.comparisons$wad.diff.boot.median[inds])
      tax.order <- data.frame(axis.loc=seq(1, N.taxa), ranks)
      plot(1:N.taxa~all.comparisons$wad.diff.boot.median[inds][tax.order$ranks], type="n", bty="n", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min, x.max), main="")
      for (p in 1:length(phyla)){
        current.comparison <- all.comparisons[inds,]
        curr.comp.ranked <- current.comparison[tax.order$ranks,]
        tax.nums <- tax.order$axis.loc[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        mids <- curr.comp.ranked$wad.diff.boot.median[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        lowers <- curr.comp.ranked$wad.diff.boot.CI.L[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        uppers <- curr.comp.ranked$wad.diff.boot.CI.U[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        arrows(x0=lowers, y0=tax.nums, x1=uppers, y1=tax.nums, length=0, angle=90, code=3, col="gray20", lwd=0.5)
        points(x=mids, y=tax.nums, pch=21, cex=0.15, lwd=0.125, col="gray20", bg="white")
      }
      axis(side=1, at=c(x.min, x.max)*1000, labels=FALSE, tck=0)
      axis(side=3, at=c(x.min, x.max)*1000, labels=FALSE, tck=0)
      abline(v=0, col="red", lty=2)
      if (y.axis == TRUE){
        par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
        axis(side=2, at=c(-1000, 1000), labels=FALSE, tck=0)
      }
    }

    par(mai=c(0.44,0.25,0.44,0.05), cex=1, mex=0.75)
    wad.diff.add.panel(comparisonID=1, AT.iso=AT.18O, y.axis=TRUE)
    par(mgp=c(3,0.25,0))				  #set spacing for the x-axis labels (the second value)
    axis(side=3, at=AT.18O, labels=AT.ape.18O, tck=-0.015, las=1, cex.axis=0.6)
    par(mgp=c(3,-0.08,0))				  #set spacing for the x-axis labels (the second value)
    axis(side=1, at=AT.wad.diff, labels=AT.wad.diff, tck=-0.015, las=1, cex.axis=0.6)
    par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
    axis(side=2, at=1:N.taxa, labels=rep(NA, N.taxa), tck=-0.015, las=1, cex.axis=0.6)
    mtext("Ranked taxon", side=2, line=0.45, cex=0.75)
    text(x=x.min*1.03, y=N.taxa*1.01, labels="control", adj=c(0,1), cex=0.6)
    text(x=x.min*1.03, y=N.taxa*0.95, labels=expression(paste(""^18, "O", sep="")), adj=c(0,1), cex=0.6)
    # par(xpd=NA)
    # legend(x=0.035, y=0.82*N.taxa, xjust=0.5, yjust=1, legend=phyla[1:ceiling(0.5*length(phyla))], bty="n", text.col=as.character(phyla.cols)[1:ceiling(0.5*length(phyla))], cex=0.4)
    # legend(x=0.056, y=0.82*N.taxa, xjust=0.5, yjust=1, legend=phyla[(ceiling(0.5*length(phyla))+1):length(phyla)], bty="n", text.col=as.character(phyla.cols)[(ceiling(0.5*length(phyla))+1):length(phyla)], cex=0.4)
    # par(xpd=FALSE)
    par(mai=c(0.44,0.25,0.44,0.05), cex=1, mex=0.75)
    wad.diff.add.panel(comparisonID=2, AT.iso=AT.18O, y.axis=FALSE)
    par(mgp=c(3,0.25,0))				  #set spacing for the x-axis labels (the second value)
    axis(side=3, at=AT.18O, labels=AT.ape.18O, tck=-0.015, las=1, cex.axis=0.6)
    par(mgp=c(3,-0.08,0))				  #set spacing for the x-axis labels (the second value)
    axis(side=1, at=AT.wad.diff, labels=AT.wad.diff, tck=-0.015, las=1, cex.axis=0.6)
    text(x=x.min*1.03, y=N.taxa*1.01, labels="added carbon", adj=c(0,1), cex=0.6)
    text(x=x.min*1.03, y=N.taxa*0.95, labels=expression(paste(""^18, "O", sep="")), adj=c(0,1), cex=0.6)
    par(mai=c(0.44,0.25,0.44,0.05), cex=1, mex=0.75)
    mtext(expression(paste("Atom fraction excess ("^18, "O or "^13, "C)", sep="")), side=3, line=1.3, cex=0.75)
    mtext(expression(paste("Difference in density (g cm"^-3, ")", sep="")), side=1, line=1.4, cex=0.75)
    wad.diff.add.panel(comparisonID=3, AT.iso=AT.13C, y.axis=FALSE)
    par(mgp=c(3,0.25,0))				  #set spacing for the x-axis labels (the second value)
    axis(side=3, at=AT.13C, labels=AT.ape.13C, tck=-0.015, las=1, cex.axis=0.6)
    par(mgp=c(3,-0.08,0))				  #set spacing for the x-axis labels (the second value)
    axis(side=1, at=AT.wad.diff, labels=AT.wad.diff, tck=-0.015, las=1, cex.axis=0.6)
    text(x=x.min*1.03, y=N.taxa*1.01, labels="added carbon", adj=c(0,1), cex=0.6)
    text(x=x.min*1.03, y=N.taxa*0.95, labels=expression(paste(""^13, "C", sep="")), adj=c(0,1), cex=0.6)
    rm(wad.diff.add.panel)
    par(mfrow=c(1,1))

    dev.off()
  #}


#_2a_Create a refined plot of primed-glucose growth rates (r determined by 13C; units=day-1) vs. primed growth rates (r determined by 18O; units=day-1)
#(graph versions with both square and squashed axes)
  #{
  #Plot 2 panels: week 1 & week 6 - SQUARE AXES
    dev.off()
    dimensions <- c((8.9/2.54), (17.8/2.54))
    # dev.new(width=dimensions[1], height=dimensions[2])
    pdf(file="qSIP_output/Figures/BH_Pretty_GlucGrowth_v_Growth_Square.pdf", width=dimensions[1], height=dimensions[2])
    layout(mat=matrix(c(1,2), 2, 1), heights=c(1/2,1/2))
    phyla <- levels(data.melted$phylum)
    phyla.cols <- rainbow(length(phyla))
    inds1 <- all.comparisons$comparisonID == 1
    inds2 <- all.comparisons$comparisonID == 2            #growth with priming (week 1)
    inds3 <- all.comparisons$comparisonID == 3            #growth on glucose with priming (week 1)
    inds6 <- all.comparisons$comparisonID == 6            #growth with priming (week 6)
    inds7 <- all.comparisons$comparisonID == 7            #growth on glucose with priming (week 6)
    ranks <- order(all.comparisons$r.boot.median[inds1])
    tax.order <- data.frame(axis.loc=seq(1, N.taxa), ranks)
    min(c(all.comparisons$r.boot.CI.L[inds2], all.comparisons$r.boot.CI.L[inds3], all.comparisons$r.boot.CI.L[inds6], all.comparisons$r.boot.CI.L[inds7]), na.rm=TRUE)
    max(c(all.comparisons$r.boot.CI.U[inds2], all.comparisons$r.boot.CI.U[inds3], all.comparisons$r.boot.CI.U[inds6], all.comparisons$r.boot.CI.U[inds7]), na.rm=TRUE)
    max(c(all.comparisons$r.boot.median[inds2], all.comparisons$r.boot.median[inds6]), na.rm=TRUE)
    max(c(all.comparisons$r.boot.median[inds3], all.comparisons$r.boot.median[inds7]), na.rm=TRUE)
    x.min <- 0
    x.max <- 0.33
    y.min <- 0
    y.max <- 0.33
  
    #Week 1:
    par(mai=c(0.44,0.54,0.2,0.05), cex=1, mex=0.75)
    plot(all.comparisons$r.boot.median[inds3][ranks]~all.comparisons$r.boot.median[inds2][ranks], type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min, x.max), ylim=c(y.min, y.max), main="")
    par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
    axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
    par(mgp=c(3,0.5,0))			#set spacing for the y-axis labels (the second value)
    axis(side=2, tck=-0.015, las=1, cex.axis=0.6)
    mtext(expression(paste("Growth on added carbon (day"^-1, ")", sep="")), side=2, line=2.1, cex=0.75)
    mtext("Week 1", side=3, line=-0.1, cex=0.75)
    abline(a=0, b=1, lty=2)
    legend(x=-0.02, y=1.05*y.max, xjust=0, yjust=1, legend=phyla, bty="n", text.col=phyla.cols, cex=0.35)
    for (p in 1:length(phyla)){
      comparison.x <- all.comparisons[inds2,]
      comp.ranked.x <- comparison.x[tax.order$ranks,]
      mids.x <- comp.ranked.x$r.boot.median[as.character(comp.ranked.x$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
      lowers.x <- comp.ranked.x$r.boot.CI.L[as.character(comp.ranked.x$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
      uppers.x <- comp.ranked.x$r.boot.CI.U[as.character(comp.ranked.x$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
      comparison.y <- all.comparisons[inds3,]
      comp.ranked.y <- comparison.y[tax.order$ranks,]
      mids.y <- comp.ranked.y$r.boot.median[as.character(comp.ranked.y$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
      lowers.y <- comp.ranked.y$r.boot.CI.L[as.character(comp.ranked.y$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
      uppers.y <- comp.ranked.y$r.boot.CI.U[as.character(comp.ranked.y$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
      points(x=mids.x, y=mids.y, pch=21, cex=0.4, col=phyla.cols[p], bg=phyla.cols[p])
      arrows(x0=lowers.x, y0=mids.y, x1=uppers.x, y1=mids.y, length=0, angle=90, code=3, col=phyla.cols[p])
      arrows(x0=mids.x, y0=lowers.y, x1=mids.x, y1=uppers.y, length=0, angle=90, code=3, col=phyla.cols[p])
    }

    #Week 6:
    par(mai=c(0.44,0.54,0.2,0.05), cex=1, mex=0.75)
    plot(all.comparisons$r.boot.median[inds7][ranks]~all.comparisons$r.boot.median[inds6][ranks], type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min, x.max), ylim=c(y.min, y.max), main="")
    par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
    axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
    mtext(expression(paste("Growth in the presence of added carbon (day"^-1, ")", sep="")), side=1, line=1.4, cex=0.75)
    par(mgp=c(3,0.5,0))			#set spacing for the y-axis labels (the second value)
    axis(side=2, tck=-0.015, las=1, cex.axis=0.6)
    mtext(expression(paste("Growth on added carbon (day"^-1, ")", sep="")), side=2, line=2.1, cex=0.75)
    mtext("Week 6", side=3, line=-0.1, cex=0.75)
    abline(a=0, b=1, lty=2)
    for (p in 1:length(phyla)){
      comparison.x <- all.comparisons[inds6,]
      comp.ranked.x <- comparison.x[tax.order$ranks,]
      mids.x <- comp.ranked.x$r.boot.median[as.character(comp.ranked.x$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
      lowers.x <- comp.ranked.x$r.boot.CI.L[as.character(comp.ranked.x$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
      uppers.x <- comp.ranked.x$r.boot.CI.U[as.character(comp.ranked.x$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
      comparison.y <- all.comparisons[inds7,]
      comp.ranked.y <- comparison.y[tax.order$ranks,]
      mids.y <- comp.ranked.y$r.boot.median[as.character(comp.ranked.y$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
      lowers.y <- comp.ranked.y$r.boot.CI.L[as.character(comp.ranked.y$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
      uppers.y <- comp.ranked.y$r.boot.CI.U[as.character(comp.ranked.y$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
      points(x=mids.x, y=mids.y, pch=21, cex=0.4, col=phyla.cols[p], bg=phyla.cols[p])
      arrows(x0=lowers.x, y0=mids.y, x1=uppers.x, y1=mids.y, length=0, angle=90, code=3, col=phyla.cols[p])
      arrows(x0=mids.x, y0=lowers.y, x1=mids.x, y1=uppers.y, length=0, angle=90, code=3, col=phyla.cols[p])
    }

    dev.off()


  #Plot 2 panels: week 1 & week 6 - SQUASHED AXES
    dev.off()
    dimensions <- c((8.9/2.54), (17.8/2.54))
    # dev.new(width=dimensions[1], height=dimensions[2])
    pdf(file="qSIP_output/Figures/BH_Pretty_GlucGrowth_v_Growth_Squashed.pdf", width=dimensions[1], height=dimensions[2])
    layout(mat=matrix(c(1,2), 2, 1), heights=c(1/2,1/2))
    phyla <- levels(data.melted$phylum)
    phyla.cols <- rainbow(length(phyla))
    inds1 <- all.comparisons$comparisonID == 1
    inds2 <- all.comparisons$comparisonID == 2            #growth with priming (week 1)
    inds3 <- all.comparisons$comparisonID == 3            #growth on glucose with priming (week 1)
    inds6 <- all.comparisons$comparisonID == 6            #growth with priming (week 6)
    inds7 <- all.comparisons$comparisonID == 7            #growth on glucose with priming (week 6)
    ranks <- order(all.comparisons$r.boot.median[inds1])
    tax.order <- data.frame(axis.loc=seq(1, N.taxa), ranks)
    min(c(all.comparisons$r.boot.CI.L[inds2], all.comparisons$r.boot.CI.L[inds3], all.comparisons$r.boot.CI.L[inds6], all.comparisons$r.boot.CI.L[inds7]), na.rm=TRUE)
    max(c(all.comparisons$r.boot.CI.U[inds2], all.comparisons$r.boot.CI.U[inds3], all.comparisons$r.boot.CI.U[inds6], all.comparisons$r.boot.CI.U[inds7]), na.rm=TRUE)
    max(c(all.comparisons$r.boot.median[inds2], all.comparisons$r.boot.median[inds6]), na.rm=TRUE)
    max(c(all.comparisons$r.boot.median[inds3], all.comparisons$r.boot.median[inds7]), na.rm=TRUE)
    x.min <- 0
    x.max <- 0.33
    y.min <- 0
    y.max <- 0.11

    #Week 1:
    par(mai=c(0.44,0.54,0.2,0.05), cex=1, mex=0.75)
    plot(all.comparisons$r.boot.median[inds3][ranks]~all.comparisons$r.boot.median[inds2][ranks], type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min, x.max), ylim=c(y.min, y.max), main="")
    par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
    axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
    par(mgp=c(3,0.5,0))			#set spacing for the y-axis labels (the second value)
    axis(side=2, tck=-0.015, las=1, cex.axis=0.6)
    mtext(expression(paste("Growth on added carbon (day"^-1, ")", sep="")), side=2, line=2.1, cex=0.75)
    mtext("Week 1", side=3, line=-0.1, cex=0.75)
    abline(a=0, b=1, lty=2)
    for (p in 1:length(phyla)){
      comparison.x <- all.comparisons[inds2,]
      comp.ranked.x <- comparison.x[tax.order$ranks,]
      mids.x <- comp.ranked.x$r.boot.median[as.character(comp.ranked.x$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
      lowers.x <- comp.ranked.x$r.boot.CI.L[as.character(comp.ranked.x$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
      uppers.x <- comp.ranked.x$r.boot.CI.U[as.character(comp.ranked.x$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
      comparison.y <- all.comparisons[inds3,]
      comp.ranked.y <- comparison.y[tax.order$ranks,]
      mids.y <- comp.ranked.y$r.boot.median[as.character(comp.ranked.y$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
      lowers.y <- comp.ranked.y$r.boot.CI.L[as.character(comp.ranked.y$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
      uppers.y <- comp.ranked.y$r.boot.CI.U[as.character(comp.ranked.y$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
      points(x=mids.x, y=mids.y, pch=21, cex=0.4, col=phyla.cols[p], bg=phyla.cols[p])
      arrows(x0=lowers.x, y0=mids.y, x1=uppers.x, y1=mids.y, length=0, angle=90, code=3, col=phyla.cols[p])
      arrows(x0=mids.x, y0=lowers.y, x1=mids.x, y1=uppers.y, length=0, angle=90, code=3, col=phyla.cols[p])
    }

    #Week 6:
    par(mai=c(0.44,0.54,0.2,0.05), cex=1, mex=0.75)
    plot(all.comparisons$r.boot.median[inds7][ranks]~all.comparisons$r.boot.median[inds6][ranks], type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min, x.max), ylim=c(y.min, y.max), main="")
    par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
    axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
    mtext(expression(paste("Growth in the presence of added carbon (day"^-1, ")", sep="")), side=1, line=1.4, cex=0.75)
    par(mgp=c(3,0.5,0))			#set spacing for the y-axis labels (the second value)
    axis(side=2, tck=-0.015, las=1, cex.axis=0.6)
    mtext(expression(paste("Growth on added carbon (day"^-1, ")", sep="")), side=2, line=2.1, cex=0.75)
    mtext("Week 6", side=3, line=-0.1, cex=0.75)
    abline(a=0, b=1, lty=2)
    legend(x=0.143, y=1.05*y.max, xjust=0, yjust=1, legend=phyla, bty="n", text.col=phyla.cols, cex=0.35)
    for (p in 1:length(phyla)){
      comparison.x <- all.comparisons[inds6,]
      comp.ranked.x <- comparison.x[tax.order$ranks,]
      mids.x <- comp.ranked.x$r.boot.median[as.character(comp.ranked.x$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
      lowers.x <- comp.ranked.x$r.boot.CI.L[as.character(comp.ranked.x$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
      uppers.x <- comp.ranked.x$r.boot.CI.U[as.character(comp.ranked.x$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
      comparison.y <- all.comparisons[inds7,]
      comp.ranked.y <- comparison.y[tax.order$ranks,]
      mids.y <- comp.ranked.y$r.boot.median[as.character(comp.ranked.y$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
      lowers.y <- comp.ranked.y$r.boot.CI.L[as.character(comp.ranked.y$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
      uppers.y <- comp.ranked.y$r.boot.CI.U[as.character(comp.ranked.y$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
      points(x=mids.x, y=mids.y, pch=21, cex=0.4, col=phyla.cols[p], bg=phyla.cols[p])
      arrows(x0=lowers.x, y0=mids.y, x1=uppers.x, y1=mids.y, length=0, angle=90, code=3, col=phyla.cols[p])
      arrows(x0=mids.x, y0=lowers.y, x1=mids.x, y1=uppers.y, length=0, angle=90, code=3, col=phyla.cols[p])
    }

    dev.off()
  #}


#_2b_Create a refined plot of primed-glucose flux rates (f determined by 13C; units = pgC g soil-1 day-1) vs. primed flux rates (f determined by 18O; units = pgC g soil-1 day-1)
#(graph versions with both square and squashed axes)
  #{
  #NOTE: flux units are converted to ngC g soil-1 day-1 in plots below:
  #Plot 2 panels: week 1 & week 6 - SQUARE AXES
    dev.off()
    dimensions <- c((8.9/2.54), (17.8/2.54))
    # dev.new(width=dimensions[1], height=dimensions[2])
    pdf(file="qSIP_output/Figures/BH_Pretty_GlucFlux_v_Flux_Square.pdf", width=dimensions[1], height=dimensions[2])
    layout(mat=matrix(c(1,2), 2, 1), heights=c(1/2,1/2))
    phyla <- levels(data.melted$phylum)
    phyla.cols <- rainbow(length(phyla))
    inds1 <- all.comparisons$comparisonID == 1
    inds2 <- all.comparisons$comparisonID == 2            #flux with priming (week 1)
    inds3 <- all.comparisons$comparisonID == 3            #flux on glucose with priming (week 1)
    inds6 <- all.comparisons$comparisonID == 6            #flux with priming (week 6)
    inds7 <- all.comparisons$comparisonID == 7            #flux on glucose with priming (week 6)
    ranks <- order(all.comparisons$f.boot.median[inds1])
    tax.order <- data.frame(axis.loc=seq(1, N.taxa), ranks)
    min(c(all.comparisons$f.boot.CI.L[inds2], all.comparisons$f.boot.CI.L[inds3], all.comparisons$f.boot.CI.L[inds6], all.comparisons$f.boot.CI.L[inds7]), na.rm=TRUE)*(1/1000)
    max(c(all.comparisons$f.boot.CI.U[inds2], all.comparisons$f.boot.CI.U[inds3], all.comparisons$f.boot.CI.U[inds6], all.comparisons$f.boot.CI.U[inds7]), na.rm=TRUE)*(1/1000)
    max(c(all.comparisons$f.boot.median[inds2], all.comparisons$f.boot.median[inds6]), na.rm=TRUE)*(1/1000)
    max(c(all.comparisons$f.boot.median[inds3], all.comparisons$f.boot.median[inds7]), na.rm=TRUE)*(1/1000)
    x.min <- -0.5
    x.max <- 17
    y.min <- -0.5
    y.max <- 17

    #Week 1:
    par(mai=c(0.44,0.48,0.2,0.05), cex=1, mex=0.75)
    plot(y=(all.comparisons$f.boot.median[inds3][ranks]*(1/1000)), x=(all.comparisons$f.boot.median[inds2][ranks]*(1/1000)), type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min, x.max), ylim=c(y.min, y.max), main="")
    par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
    axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
    par(mgp=c(3,0.5,0))			#set spacing for the y-axis labels (the second value)
    axis(side=2, tck=-0.015, las=1, cex.axis=0.6)
    mtext(expression(paste("Flux on added carbon (ng C g soil"^-1, " day"^-1, ")", sep="")), side=2, line=1.5, cex=0.75)
    mtext("Week 1", side=3, line=-0.1, cex=0.75)
    abline(a=0, b=1, lty=2)
    legend(x=-1.5, y=1.05*y.max, xjust=0, yjust=1, legend=phyla, bty="n", text.col=phyla.cols, cex=0.35)
    for (p in 1:length(phyla)){
      comparison.x <- all.comparisons[inds2,]
      comp.ranked.x <- comparison.x[tax.order$ranks,]
      mids.x <- (1/1000)*comp.ranked.x$f.boot.median[as.character(comp.ranked.x$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
      lowers.x <- (1/1000)*comp.ranked.x$f.boot.CI.L[as.character(comp.ranked.x$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
      uppers.x <- (1/1000)*comp.ranked.x$f.boot.CI.U[as.character(comp.ranked.x$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
      comparison.y <- all.comparisons[inds3,]
      comp.ranked.y <- comparison.y[tax.order$ranks,]
      mids.y <- (1/1000)*comp.ranked.y$f.boot.median[as.character(comp.ranked.y$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
      lowers.y <- (1/1000)*comp.ranked.y$f.boot.CI.L[as.character(comp.ranked.y$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
      uppers.y <- (1/1000)*comp.ranked.y$f.boot.CI.U[as.character(comp.ranked.y$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
      points(x=mids.x, y=mids.y, pch=21, cex=0.4, col=phyla.cols[p], bg=phyla.cols[p])
      arrows(x0=lowers.x, y0=mids.y, x1=uppers.x, y1=mids.y, length=0, angle=90, code=3, col=phyla.cols[p])
      arrows(x0=mids.x, y0=lowers.y, x1=mids.x, y1=uppers.y, length=0, angle=90, code=3, col=phyla.cols[p])
    }

    #Week 6:
    par(mai=c(0.44,0.48,0.2,0.05), cex=1, mex=0.75)
    plot(y=(all.comparisons$f.boot.median[inds7][ranks]*(1/1000)), x=(all.comparisons$f.boot.median[inds6][ranks]*(1/1000)), type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min, x.max), ylim=c(y.min, y.max), main="")
    par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
    axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
    mtext(expression(paste("Flux in presence of added carbon (ng C g soil"^-1, " day"^-1, ")", sep="")), side=1, line=1.4, cex=0.75)
    par(mgp=c(3,0.5,0))			#set spacing for the y-axis labels (the second value)
    axis(side=2, tck=-0.015, las=1, cex.axis=0.6)
    mtext(expression(paste("Flux on added carbon (ng C g soil"^-1, " day"^-1, ")", sep="")), side=2, line=1.5, cex=0.75)
    mtext("Week 6", side=3, line=-0.1, cex=0.75)
    abline(a=0, b=1, lty=2)
    for (p in 1:length(phyla)){
      comparison.x <- all.comparisons[inds6,]
      comp.ranked.x <- comparison.x[tax.order$ranks,]
      mids.x <- (1/1000)*comp.ranked.x$f.boot.median[as.character(comp.ranked.x$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
      lowers.x <- (1/1000)*comp.ranked.x$f.boot.CI.L[as.character(comp.ranked.x$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
      uppers.x <- (1/1000)*comp.ranked.x$f.boot.CI.U[as.character(comp.ranked.x$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
      comparison.y <- all.comparisons[inds7,]
      comp.ranked.y <- comparison.y[tax.order$ranks,]
      mids.y <- (1/1000)*comp.ranked.y$f.boot.median[as.character(comp.ranked.y$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
      lowers.y <- (1/1000)*comp.ranked.y$f.boot.CI.L[as.character(comp.ranked.y$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
      uppers.y <- (1/1000)*comp.ranked.y$f.boot.CI.U[as.character(comp.ranked.y$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
      points(x=mids.x, y=mids.y, pch=21, cex=0.4, col=phyla.cols[p], bg=phyla.cols[p])
      arrows(x0=lowers.x, y0=mids.y, x1=uppers.x, y1=mids.y, length=0, angle=90, code=3, col=phyla.cols[p])
      arrows(x0=mids.x, y0=lowers.y, x1=mids.x, y1=uppers.y, length=0, angle=90, code=3, col=phyla.cols[p])
    }

    dev.off()


  #Plot 2 panels: week 1 & week 6 - SQUASHED AXES
    dev.off()
    dimensions <- c((8.9/2.54), (17.8/2.54))
    # dev.new(width=dimensions[1], height=dimensions[2])
    pdf(file="qSIP_output/Figures/BH_Pretty_GlucFlux_v_Flux_Squashed.pdf", width=dimensions[1], height=dimensions[2])
    layout(mat=matrix(c(1,2), 2, 1), heights=c(1/2,1/2))
    phyla <- levels(data.melted$phylum)
    phyla.cols <- rainbow(length(phyla))
    inds1 <- all.comparisons$comparisonID == 1
    inds2 <- all.comparisons$comparisonID == 2            #flux with priming (week 1)
    inds3 <- all.comparisons$comparisonID == 3            #flux on glucose with priming (week 1)
    inds6 <- all.comparisons$comparisonID == 6            #flux with priming (week 6)
    inds7 <- all.comparisons$comparisonID == 7            #flux on glucose with priming (week 6)
    ranks <- order(all.comparisons$f.boot.median[inds1])
    tax.order <- data.frame(axis.loc=seq(1, N.taxa), ranks)
    min(c(all.comparisons$f.boot.CI.L[inds2], all.comparisons$f.boot.CI.L[inds3], all.comparisons$f.boot.CI.L[inds6], all.comparisons$f.boot.CI.L[inds7]), na.rm=TRUE)*(1/1000)
    max(c(all.comparisons$f.boot.CI.U[inds2], all.comparisons$f.boot.CI.U[inds3], all.comparisons$f.boot.CI.U[inds6], all.comparisons$f.boot.CI.U[inds7]), na.rm=TRUE)*(1/1000)
    max(c(all.comparisons$f.boot.median[inds2], all.comparisons$f.boot.median[inds6]), na.rm=TRUE)*(1/1000)
    max(c(all.comparisons$f.boot.median[inds3], all.comparisons$f.boot.median[inds7]), na.rm=TRUE)*(1/1000)
    x.min <- -0.5
    x.max <- 17
    y.min <- -0.5
    y.max <- 9

    #Week 1:
    par(mai=c(0.44,0.48,0.2,0.05), cex=1, mex=0.75)
    plot(y=(all.comparisons$f.boot.median[inds3][ranks]*(1/1000)), x=(all.comparisons$f.boot.median[inds2][ranks]*(1/1000)), type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min, x.max), ylim=c(y.min, y.max), main="")
    par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
    axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
    par(mgp=c(3,0.5,0))			#set spacing for the y-axis labels (the second value)
    axis(side=2, tck=-0.015, las=1, cex.axis=0.6)
    mtext(expression(paste("Flux on added carbon (ng C g soil"^-1, " day"^-1, ")", sep="")), side=2, line=1.5, cex=0.75)
    mtext("Week 1", side=3, line=-0.1, cex=0.75)
    abline(a=0, b=1, lty=2)
    for (p in 1:length(phyla)){
      comparison.x <- all.comparisons[inds2,]
      comp.ranked.x <- comparison.x[tax.order$ranks,]
      mids.x <- (1/1000)*comp.ranked.x$f.boot.median[as.character(comp.ranked.x$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
      lowers.x <- (1/1000)*comp.ranked.x$f.boot.CI.L[as.character(comp.ranked.x$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
      uppers.x <- (1/1000)*comp.ranked.x$f.boot.CI.U[as.character(comp.ranked.x$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
      comparison.y <- all.comparisons[inds3,]
      comp.ranked.y <- comparison.y[tax.order$ranks,]
      mids.y <- (1/1000)*comp.ranked.y$f.boot.median[as.character(comp.ranked.y$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
      lowers.y <- (1/1000)*comp.ranked.y$f.boot.CI.L[as.character(comp.ranked.y$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
      uppers.y <- (1/1000)*comp.ranked.y$f.boot.CI.U[as.character(comp.ranked.y$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
      points(x=mids.x, y=mids.y, pch=21, cex=0.4, col=phyla.cols[p], bg=phyla.cols[p])
      arrows(x0=lowers.x, y0=mids.y, x1=uppers.x, y1=mids.y, length=0, angle=90, code=3, col=phyla.cols[p])
      arrows(x0=mids.x, y0=lowers.y, x1=mids.x, y1=uppers.y, length=0, angle=90, code=3, col=phyla.cols[p])
    }

    #Week 6:
    par(mai=c(0.44,0.48,0.2,0.05), cex=1, mex=0.75)
    plot(y=(all.comparisons$f.boot.median[inds7][ranks]*(1/1000)), x=(all.comparisons$f.boot.median[inds6][ranks]*(1/1000)), type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min, x.max), ylim=c(y.min, y.max), main="")
    par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
    axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
    mtext(expression(paste("Flux in presence of added carbon (ng C g soil"^-1, " day"^-1, ")", sep="")), side=1, line=1.4, cex=0.75)
    par(mgp=c(3,0.5,0))			#set spacing for the y-axis labels (the second value)
    axis(side=2, tck=-0.015, las=1, cex.axis=0.6)
    mtext(expression(paste("Flux on added carbon (ng C g soil"^-1, " day"^-1, ")", sep="")), side=2, line=1.5, cex=0.75)
    mtext("Week 6", side=3, line=-0.1, cex=0.75)
    abline(a=0, b=1, lty=2)
    legend(x=13, y=1.05*y.max, xjust=0, yjust=1, legend=phyla, bty="n", text.col=phyla.cols, cex=0.35)
    for (p in 1:length(phyla)){
      comparison.x <- all.comparisons[inds6,]
      comp.ranked.x <- comparison.x[tax.order$ranks,]
      mids.x <- (1/1000)*comp.ranked.x$f.boot.median[as.character(comp.ranked.x$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
      lowers.x <- (1/1000)*comp.ranked.x$f.boot.CI.L[as.character(comp.ranked.x$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
      uppers.x <- (1/1000)*comp.ranked.x$f.boot.CI.U[as.character(comp.ranked.x$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
      comparison.y <- all.comparisons[inds7,]
      comp.ranked.y <- comparison.y[tax.order$ranks,]
      mids.y <- (1/1000)*comp.ranked.y$f.boot.median[as.character(comp.ranked.y$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
      lowers.y <- (1/1000)*comp.ranked.y$f.boot.CI.L[as.character(comp.ranked.y$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
      uppers.y <- (1/1000)*comp.ranked.y$f.boot.CI.U[as.character(comp.ranked.y$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
      points(x=mids.x, y=mids.y, pch=21, cex=0.4, col=phyla.cols[p], bg=phyla.cols[p])
      arrows(x0=lowers.x, y0=mids.y, x1=uppers.x, y1=mids.y, length=0, angle=90, code=3, col=phyla.cols[p])
      arrows(x0=mids.x, y0=lowers.y, x1=mids.x, y1=uppers.y, length=0, angle=90, code=3, col=phyla.cols[p])
    }

    dev.off()
  #}


#_2c_For those taxa that exhibit growth on glucose, calculate the range in the ratio of (growth on glucose/growth in presence of glucose):
  #{
    #Week 1:
      inds2 <- all.comparisons$comparisonID == 2            #growth with priming (week 1)
      inds3 <- all.comparisons$comparisonID == 3            #growth on glucose with priming (week 1)
      ranks <- order(all.comparisons$r.boot.CI.L[inds3])
      w1.growth.on.glucose <- all.comparisons$r.boot.median[inds3][ranks]
      w1.growth.with.glucose <- all.comparisons$r.boot.median[inds2][ranks]
      w1.ratio.for.gluc.growers <- w1.growth.on.glucose[all.comparisons$r.boot.CI.L[inds3][ranks] > 0] / w1.growth.with.glucose[all.comparisons$r.boot.CI.L[inds3][ranks] > 0]
      range(w1.ratio.for.gluc.growers)
      mean(w1.ratio.for.gluc.growers)
      median(w1.ratio.for.gluc.growers)

    #Week 6:
      inds6 <- all.comparisons$comparisonID == 6            #growth with priming (week 6)
      inds7 <- all.comparisons$comparisonID == 7            #growth on glucose with priming (week 6)
      ranks <- order(all.comparisons$r.boot.CI.L[inds7])
      w6.growth.on.glucose <- all.comparisons$r.boot.median[inds7][ranks]
      w6.growth.with.glucose <- all.comparisons$r.boot.median[inds6][ranks]
      w6.ratio.for.gluc.growers <- w6.growth.on.glucose[all.comparisons$r.boot.CI.L[inds7][ranks] > 0] / w6.growth.with.glucose[all.comparisons$r.boot.CI.L[inds7][ranks] > 0]
      range(w6.ratio.for.gluc.growers)
      mean(w6.ratio.for.gluc.growers)
      median(w6.ratio.for.gluc.growers)

      dev.off()
      # dev.new(width=6, height=6)
      pdf(file="qSIP_output/Figures/BH_Ugly_GlucGrowth_to_Growth_for_GlucGrowers.pdf", width=6, height=6)
      boxplot(list(w1.ratio.for.gluc.growers, w6.ratio.for.gluc.growers), names=c("Week 1", "Week 6"), ylim=c(0,1), ylab="growth on glucose / growth with glucose", main="for only those taxa with positive growth on glucose")

      dev.off()
  #}


#_3a_Create a plot of average growth rates for each phylum:
  #{
  #Plot 6 panels: week 1 & week 6
  #(NOTE THAT THERE IS ONE TAXON WHICH HAS PHYLUM == NA AND IS THUS NOT INCLUDED IN THIS FIGURE)
    dev.off()
    dimensions <- c((18.3/2.54), (13.725/2.54))
    # dev.new(width=dimensions[1], height=dimensions[2])
    pdf(file="qSIP_output/Figures/BH_Pretty_Growth_by_Phylum_6panel.pdf", width=dimensions[1], height=dimensions[2])
    layout(mat=matrix(c(1,2,3,4,5,6), 2, 3, byrow=TRUE), widths=lcm(c(((dimensions[1]-0.27)/3)+0.27,(dimensions[1]-0.27)/3,(dimensions[1]-0.27)/3)*2.54), heights=lcm(c((dimensions[2]-0.08)/2,((dimensions[2]-0.08)/2)+0.08)*2.54))
    phyla <- levels(data.melted$phylum)
    phyla.cols <- rainbow(length(phyla))
    phykey <- data.frame(phy.num=1:length(phyla), abbr=substr(phyla, 1, 3), phyla=phyla)  #re-assign the 'phy.num' variable to different phyla to change the order for which phyla are plotted
    phykey$abbr <- as.character(phykey$abbr)
    phykey$abbr[phykey$phyla == "[Thermi]"] <- "[The]"  #recode some ambiguous abbreviations
    phykey$abbr[phykey$phyla == "Chlamydiae"] <- "Chlam"
    phykey$abbr[phykey$phyla == "Chlorobi"] <- "Chlbi"
    phykey$abbr[phykey$phyla == "Chloroflexi"] <- "Chlfx"    
    phykey$abbr <- factor(phykey$abbr)
    taxa.id.present <- taxa.id[as.character(taxa.id$taxon) %in% levels(data.melted$taxon),]
      taxa.id.present$code <- factor(taxa.id.present$code)
      taxa.id.present$kingdom <- factor(taxa.id.present$kingdom)
      taxa.id.present$phylum <- factor(taxa.id.present$phylum)
      taxa.id.present$class <- factor(taxa.id.present$class)
      taxa.id.present$order <- factor(taxa.id.present$order)
      taxa.id.present$family <- factor(taxa.id.present$family)
      taxa.id.present$genus <- factor(taxa.id.present$genus)
    PHY <- data.frame(taxonID=taxa.id.present$taxon, phyla=taxa.id.present$phylum)
    PHY <- PHY[!is.na(PHY$phyla),]  #EXCLUDE THE TAXON THAT IS "NA" for PHYLUM
    for (i in 1:dim(PHY)[1]){
      PHY$abbr[i] <- as.character(phykey$abbr[phykey$phyla == PHY$phyla[i]])
      PHY$phy.num[i] <- phykey$phy.num[phykey$phyla == PHY$phyla[i]]
    }
    PHY$abbr <- factor(PHY$abbr)
    y.min <- min(all.comparisons$r.boot.median, na.rm=TRUE)
    y.max <- max(all.comparisons$r.boot.median, na.rm=TRUE)

    r.by.phylum.add.panel <- function(comparisonID, y.axis, x.axis.lab=FALSE, overlay=FALSE){
      inds <- all.comparisons$comparisonID == comparisonID
      if (overlay != TRUE){
        plot(y=rep(x.max, length(phykey$abbr)), x=1:length(phykey$abbr), type="n", bty="n", xlab="", ylab="", xaxt="n", yaxt="n", ylim=c(y.min, y.max), main="")
      }
      for (p in 1:length(phykey$abbr)){
        AT <- unique(PHY$phy.num[PHY$abbr == phykey$abbr[p]])
        LAB <- phykey$abbr[p]
        y <- all.comparisons$r.boot.median[inds & all.comparisons$taxonID %in% PHY$taxonID[PHY$abbr == phykey$abbr[p]]]
        x <- rep(AT, length(y))
        N <- length(y)
        boxplot(y~x, add=TRUE, at=AT, bty="l", xlab="", ylab="", xaxt="n", yaxt="n", col="gray96", outpch=21, medlwd=2, boxcol=phyla.cols[p], medcol=phyla.cols[p], whiskcol=phyla.cols[p], staplecol=phyla.cols[p], outcol=phyla.cols[p], outbg="gray96", boxwex=0.85, cex=0.4, outline=TRUE, frame=FALSE)
        par(xpd=NA)
        text(x=x, y=max(y, na.rm=TRUE)+0.01, labels=N, adj=c(0.5,0), cex=0.35, col="gray60")
        par(xpd=FALSE)
      }
      axis(side=1, at=c(x.min, x.max)*1000, labels=FALSE, tck=0)
       if (x.axis.lab == TRUE){
        par(mgp=c(3,0.5,0))				  #set spacing for the x-axis labels (the second value)
        axis(side=1, at=1:length(phykey$abbr), labels=phykey$abbr, tck=-0.015, las=3, cex.axis=0.4)
      }
      else {
        par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
        axis(side=1, at=1:length(phykey$abbr), labels=FALSE, tck=-0.015, las=1, cex.axis=0.6)
      }
      if (y.axis == TRUE){
        par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
        axis(side=2, at=c(-1000, 1000), labels=FALSE, tck=0)
      }
    }

    par(mai=c(0.27,0.52,0.1,0.05), cex=1, mex=0.75)
    r.by.phylum.add.panel(comparisonID=1, y.axis=TRUE, x.axis.lab=FALSE)
    par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
    axis(side=2, tck=-0.015, las=1, cex.axis=0.6)
    mtext(expression(paste(italic(r), " (day"^-1, ")", sep="")), side=2, line=1.95, cex=0.75)
    text(x=1*1.03, y=y.max*1.01, labels="Week 1", adj=c(0,1), cex=0.6)
    text(x=1*1.03, y=y.max*0.95, labels="control", adj=c(0,1), cex=0.6)
    text(x=1*1.03, y=y.max*0.89, labels=expression(paste(""^18, "O", sep="")), adj=c(0,1), cex=0.6)
    par(xpd=NA)
    segments(x0=1-(29*0.04), y0=0, x1=99.5, y1=0, col="black", lty=2)
    par(xpd=FALSE)
    r.by.phylum.add.panel(comparisonID=1, y.axis=TRUE, x.axis.lab=FALSE, overlay=TRUE)

    par(mai=c(0.27,0.25,0.1,0.05), cex=1, mex=0.75)
    r.by.phylum.add.panel(comparisonID=2, y.axis=FALSE, x.axis.lab=FALSE)
    text(x=1*1.03, y=y.max*1.01, labels="Week 1", adj=c(0,1), cex=0.6)
    text(x=1*1.03, y=y.max*0.95, labels="added carbon", adj=c(0,1), cex=0.6)
    text(x=1*1.03, y=y.max*0.89, labels=expression(paste(""^18, "O", sep="")), adj=c(0,1), cex=0.6)

    par(mai=c(0.27,0.25,0.1,0.05), cex=1, mex=0.75)
    r.by.phylum.add.panel(comparisonID=3, y.axis=FALSE, x.axis.lab=FALSE)
    text(x=1*1.03, y=y.max*1.01, labels="Week 1", adj=c(0,1), cex=0.6)
    text(x=1*1.03, y=y.max*0.95, labels="added carbon", adj=c(0,1), cex=0.6)
    text(x=1*1.03, y=y.max*0.89, labels=expression(paste(""^13, "C", sep="")), adj=c(0,1), cex=0.6)

    par(mai=c(0.44,0.52,0.1,0.05), cex=1, mex=0.75)
    r.by.phylum.add.panel(comparisonID=5, y.axis=TRUE, x.axis.lab=TRUE)
    par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
    axis(side=2, tck=-0.015, las=1, cex.axis=0.6)
    mtext(expression(paste(italic(r), " (day"^-1, ")", sep="")), side=2, line=1.95, cex=0.75)
    text(x=1*1.03, y=y.max*1.01, labels="Week 6", adj=c(0,1), cex=0.6)
    text(x=1*1.03, y=y.max*0.95, labels="control", adj=c(0,1), cex=0.6)
    text(x=1*1.03, y=y.max*0.89, labels=expression(paste(""^18, "O", sep="")), adj=c(0,1), cex=0.6)
    par(xpd=NA)
    segments(x0=1-(29*0.04), y0=0, x1=99.5, y1=0, col="black", lty=2)
    par(xpd=FALSE)
    r.by.phylum.add.panel(comparisonID=5, y.axis=TRUE, x.axis.lab=TRUE, overlay=TRUE)

    par(mai=c(0.44,0.25,0.1,0.05), cex=1, mex=0.75)
    r.by.phylum.add.panel(comparisonID=6, y.axis=FALSE, x.axis.lab=TRUE)
    text(x=1*1.03, y=y.max*1.01, labels="Week 6", adj=c(0,1), cex=0.6)
    text(x=1*1.03, y=y.max*0.95, labels="added carbon", adj=c(0,1), cex=0.6)
    text(x=1*1.03, y=y.max*0.89, labels=expression(paste(""^18, "O", sep="")), adj=c(0,1), cex=0.6)

    par(mai=c(0.44,0.25,0.1,0.05), cex=1, mex=0.75)
    r.by.phylum.add.panel(comparisonID=7, y.axis=FALSE, x.axis.lab=TRUE)
    text(x=1*1.03, y=y.max*1.01, labels="Week 6", adj=c(0,1), cex=0.6)
    text(x=1*1.03, y=y.max*0.95, labels="added carbon", adj=c(0,1), cex=0.6)
    text(x=1*1.03, y=y.max*0.89, labels=expression(paste(""^13, "C", sep="")), adj=c(0,1), cex=0.6)
    par(xpd=NA)
    legend(x=19, y=1.22*y.max, xjust=0, yjust=1, legend=phyla, bty="n", text.col=phyla.cols, cex=0.35)
    par(xpd=FALSE)

    dev.off()
  #}


#_3b_Create a plot of average growth rates for each phylum (exclude 13C growth):
  #{
  #Plot 4 panels: week 1 & week 6
  #(NOTE THAT THERE IS ONE TAXON WHICH HAS PHYLUM == NA AND IS THUS NOT INCLUDED IN THIS FIGURE)
    dev.off()
    dimensions <- c(((18.3-5.8714)/2.54), (13.725/2.54))
    # dev.new(width=dimensions[1], height=dimensions[2])
    pdf(file="qSIP_output/Figures/BH_Pretty_Growth_by_Phylum_4panel.pdf", width=dimensions[1], height=dimensions[2])
    layout(mat=matrix(c(1,2,3,4), 2, 2, byrow=TRUE), widths=lcm(c(((dimensions[1]-0.27)/2)+0.27,(dimensions[1]-0.27)/2)*2.54), heights=lcm(c((dimensions[2]-0.08)/2,((dimensions[2]-0.08)/2)+0.08)*2.54))
    phyla <- levels(data.melted$phylum)
    phyla.cols <- rainbow(length(phyla))
    phykey <- data.frame(phy.num=1:length(phyla), abbr=substr(phyla, 1, 3), phyla=phyla)  #re-assign the 'phy.num' variable to different phyla to change the order for which phyla are plotted
    phykey$abbr <- as.character(phykey$abbr)
    phykey$abbr[phykey$phyla == "[Thermi]"] <- "[The]"  #recode some ambiguous abbreviations
    phykey$abbr[phykey$phyla == "Chlamydiae"] <- "Chlam"
    phykey$abbr[phykey$phyla == "Chlorobi"] <- "Chlbi"
    phykey$abbr[phykey$phyla == "Chloroflexi"] <- "Chlfx"    
    phykey$abbr <- factor(phykey$abbr)
    taxa.id.present <- taxa.id[as.character(taxa.id$taxon) %in% levels(data.melted$taxon),]
      taxa.id.present$code <- factor(taxa.id.present$code)
      taxa.id.present$kingdom <- factor(taxa.id.present$kingdom)
      taxa.id.present$phylum <- factor(taxa.id.present$phylum)
      taxa.id.present$class <- factor(taxa.id.present$class)
      taxa.id.present$order <- factor(taxa.id.present$order)
      taxa.id.present$family <- factor(taxa.id.present$family)
      taxa.id.present$genus <- factor(taxa.id.present$genus)
    PHY <- data.frame(taxonID=taxa.id.present$taxon, phyla=taxa.id.present$phylum)
    PHY <- PHY[!is.na(PHY$phyla),]  #EXCLUDE THE TAXON THAT IS "NA" for PHYLUM
    for (i in 1:dim(PHY)[1]){
      PHY$abbr[i] <- as.character(phykey$abbr[phykey$phyla == PHY$phyla[i]])
      PHY$phy.num[i] <- phykey$phy.num[phykey$phyla == PHY$phyla[i]]
    }
    PHY$abbr <- factor(PHY$abbr)
    y.min <- min(all.comparisons$r.boot.median, na.rm=TRUE)
    y.max <- max(all.comparisons$r.boot.median, na.rm=TRUE)

    r.by.phylum.add.panel <- function(comparisonID, y.axis, x.axis.lab=FALSE, overlay=FALSE){
      inds <- all.comparisons$comparisonID == comparisonID
      if (overlay != TRUE){
        plot(y=rep(x.max, length(phykey$abbr)), x=1:length(phykey$abbr), type="n", bty="n", xlab="", ylab="", xaxt="n", yaxt="n", ylim=c(y.min, y.max), main="")
      }
      for (p in 1:length(phykey$abbr)){
        AT <- unique(PHY$phy.num[PHY$abbr == phykey$abbr[p]])
        LAB <- phykey$abbr[p]
        y <- all.comparisons$r.boot.median[inds & all.comparisons$taxonID %in% PHY$taxonID[PHY$abbr == phykey$abbr[p]]]
        x <- rep(AT, length(y))
        N <- length(y)
        boxplot(y~x, add=TRUE, at=AT, bty="l", xlab="", ylab="", xaxt="n", yaxt="n", col="gray96", outpch=21, medlwd=2, boxcol=phyla.cols[p], medcol=phyla.cols[p], whiskcol=phyla.cols[p], staplecol=phyla.cols[p], outcol=phyla.cols[p], outbg="gray96", boxwex=0.85, cex=0.4, outline=TRUE, frame=FALSE)
        par(xpd=NA)
        text(x=x, y=max(y, na.rm=TRUE)+0.01, labels=N, adj=c(0.5,0), cex=0.35, col="gray60")
        par(xpd=FALSE)
      }
      axis(side=1, at=c(x.min, x.max)*1000, labels=FALSE, tck=0)
       if (x.axis.lab == TRUE){
        par(mgp=c(3,0.5,0))				  #set spacing for the x-axis labels (the second value)
        axis(side=1, at=1:length(phykey$abbr), labels=phykey$abbr, tck=-0.015, las=3, cex.axis=0.4)
      }
      else {
        par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
        axis(side=1, at=1:length(phykey$abbr), labels=FALSE, tck=-0.015, las=1, cex.axis=0.6)
      }
      if (y.axis == TRUE){
        par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
        axis(side=2, at=c(-1000, 1000), labels=FALSE, tck=0)
      }
    }

    par(mai=c(0.27,0.52,0.1,0.05), cex=1, mex=0.75)
    r.by.phylum.add.panel(comparisonID=1, y.axis=TRUE, x.axis.lab=FALSE)
    par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
    axis(side=2, tck=-0.015, las=1, cex.axis=0.6)
    mtext(expression(paste(italic(r), " (day"^-1, ")", sep="")), side=2, line=1.95, cex=0.75)
    text(x=1*1.03, y=y.max*1.01, labels="Week 1", adj=c(0,1), cex=0.6)
    text(x=1*1.03, y=y.max*0.95, labels="control", adj=c(0,1), cex=0.6)
    text(x=1*1.03, y=y.max*0.89, labels=expression(paste(""^18, "O", sep="")), adj=c(0,1), cex=0.6)
    par(xpd=NA)
    segments(x0=1-(29*0.04), y0=0, x1=64.5, y1=0, col="black", lty=2)
    par(xpd=FALSE)
    r.by.phylum.add.panel(comparisonID=1, y.axis=TRUE, x.axis.lab=FALSE, overlay=TRUE)

    par(mai=c(0.27,0.25,0.1,0.05), cex=1, mex=0.75)
    r.by.phylum.add.panel(comparisonID=2, y.axis=FALSE, x.axis.lab=FALSE)
    text(x=1*1.03, y=y.max*1.01, labels="Week 1", adj=c(0,1), cex=0.6)
    text(x=1*1.03, y=y.max*0.95, labels="added carbon", adj=c(0,1), cex=0.6)
    text(x=1*1.03, y=y.max*0.89, labels=expression(paste(""^18, "O", sep="")), adj=c(0,1), cex=0.6)

    par(mai=c(0.44,0.52,0.1,0.05), cex=1, mex=0.75)
    r.by.phylum.add.panel(comparisonID=5, y.axis=TRUE, x.axis.lab=TRUE)
    par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
    axis(side=2, tck=-0.015, las=1, cex.axis=0.6)
    mtext(expression(paste(italic(r), " (day"^-1, ")", sep="")), side=2, line=1.95, cex=0.75)
    text(x=1*1.03, y=y.max*1.01, labels="Week 6", adj=c(0,1), cex=0.6)
    text(x=1*1.03, y=y.max*0.95, labels="control", adj=c(0,1), cex=0.6)
    text(x=1*1.03, y=y.max*0.89, labels=expression(paste(""^18, "O", sep="")), adj=c(0,1), cex=0.6)
    par(xpd=NA)
    segments(x0=1-(29*0.04), y0=0, x1=64.5, y1=0, col="black", lty=2)
    par(xpd=FALSE)
    r.by.phylum.add.panel(comparisonID=5, y.axis=TRUE, x.axis.lab=TRUE, overlay=TRUE)

    par(mai=c(0.44,0.25,0.1,0.05), cex=1, mex=0.75)
    r.by.phylum.add.panel(comparisonID=6, y.axis=FALSE, x.axis.lab=TRUE)
    text(x=1*1.03, y=y.max*1.01, labels="Week 6", adj=c(0,1), cex=0.6)
    text(x=1*1.03, y=y.max*0.95, labels="added carbon", adj=c(0,1), cex=0.6)
    text(x=1*1.03, y=y.max*0.89, labels=expression(paste(""^18, "O", sep="")), adj=c(0,1), cex=0.6)
    par(xpd=NA)
    legend(x=13, y=1.24*y.max, xjust=0, yjust=1, legend=phyla, bty="n", text.col=phyla.cols, cex=0.35)
    par(xpd=FALSE)

    dev.off()
  #}


#_3c_Determine whether any taxa that grew on glucose did not grow in the presence of glucose (ideally, there shouldn't be any):
  #{
  #Week 1: calculate the number of taxa that grew, grew with priming, grew on glucose with priming:
    #Total number of taxa:
    N.taxa <- length(levels(factor(all.comparisons$taxonID[all.comparisons$comparisonID %in% c(1,2,3)])))
    N.taxa
    #Number of taxa that grew with no added glucose:
    inds <- all.comparisons$comparisonID == 1
    sum(all.comparisons$r.boot.CI.L[inds] > 0, na.rm=TRUE)  
    w1.control.growers <- factor(all.comparisons$taxonID[inds & all.comparisons$r.boot.CI.L > 0 & !is.na(all.comparisons$r.boot.CI.L)])
    #Number of taxa that grew with added glucose:
    inds <- all.comparisons$comparisonID == 2
    sum(all.comparisons$r.boot.CI.L[inds] > 0, na.rm=TRUE)  
    w1.addedC.growers <- factor(all.comparisons$taxonID[inds & all.comparisons$r.boot.CI.L > 0 & !is.na(all.comparisons$r.boot.CI.L)])
    #Number of taxa that grew on added glucose:
    inds <- all.comparisons$comparisonID == 3
    sum(all.comparisons$r.boot.CI.L[inds] > 0, na.rm=TRUE)  
    w1.glucose.growers <- factor(all.comparisons$taxonID[inds & all.comparisons$r.boot.CI.L > 0 & !is.na(all.comparisons$r.boot.CI.L)])
    #Find the members of the first group that are not in the second group:
    length(setdiff(w1.glucose.growers, w1.addedC.growers))
    setdiff(w1.glucose.growers, w1.addedC.growers)  #In week 1, there was 1 taxon that grew on glucose but did not grow in the presence of added carbon
    length(setdiff(w1.addedC.growers, w1.glucose.growers))
    setdiff(w1.addedC.growers, w1.glucose.growers)  #In week 1, there were 188 taxa that grew in the presence of added carbon but did not grow on added glucose
    length(setdiff(w1.control.growers, w1.addedC.growers))
    setdiff(w1.control.growers, w1.addedC.growers)  #In week 1, there were 8 taxa that grew in the control treatment (no added carbon) but did not grow when carbon was added
    length(setdiff(w1.addedC.growers, w1.control.growers))
    setdiff(w1.addedC.growers, w1.control.growers)  #In week 1, there were 189 taxa that grew in the presence of added carbon but did not grow in the control treatment (no added carbon)
    length(setdiff(w1.control.growers, w1.glucose.growers))
    setdiff(w1.control.growers, w1.glucose.growers)  #In week 1, there were 83 taxa that grew in the control treatment (no added carbon) but did not grow on added glucose
    length(setdiff(w1.glucose.growers, w1.control.growers))
    setdiff(w1.glucose.growers, w1.control.growers)  #In week 1, there were 77 taxa that grew on added glucose but did not grow in the control treatment (no added carbon)
  
  #Week 6: calculate the number of taxa that grew, grew with priming, grew on glucose with priming:
    #Total number of taxa:
    N.taxa <- length(levels(factor(all.comparisons$taxonID[all.comparisons$comparisonID %in% c(5,6,7)])))
    N.taxa
    #Number of taxa that grew with no added glucose:
    inds <- all.comparisons$comparisonID == 5
    sum(all.comparisons$r.boot.CI.L[inds] > 0, na.rm=TRUE)  
    w6.control.growers <- factor(all.comparisons$taxonID[inds & all.comparisons$r.boot.CI.L > 0 & !is.na(all.comparisons$r.boot.CI.L)])
    #Number of taxa that grew with added glucose:
    inds <- all.comparisons$comparisonID == 6
    sum(all.comparisons$r.boot.CI.L[inds] > 0, na.rm=TRUE)  
    w6.addedC.growers <- factor(all.comparisons$taxonID[inds & all.comparisons$r.boot.CI.L > 0 & !is.na(all.comparisons$r.boot.CI.L)])
    #Number of taxa that grew on added glucose:
    inds <- all.comparisons$comparisonID == 7
    sum(all.comparisons$r.boot.CI.L[inds] > 0, na.rm=TRUE)  
    w6.glucose.growers <- factor(all.comparisons$taxonID[inds & all.comparisons$r.boot.CI.L > 0 & !is.na(all.comparisons$r.boot.CI.L)])
    #Find the members of the first group that are not in the second group:
    length(setdiff(w6.glucose.growers, w6.addedC.growers))
    setdiff(w6.glucose.growers, w6.addedC.growers)  #In week 6, there were 11 taxa that grew on glucose but did not grow in the presence of added carbon
    length(setdiff(w6.addedC.growers, w6.glucose.growers))
    setdiff(w6.addedC.growers, w6.glucose.growers)  #In week 6, there were 96 taxa that grew in the presence of added carbon but did not grow on added glucose
    length(setdiff(w6.control.growers, w6.addedC.growers))
    setdiff(w6.control.growers, w6.addedC.growers)  #In week 6, there were 151 taxa that grew in the control treatment (no added carbon) but did not grow when carbon was added
    length(setdiff(w6.addedC.growers, w6.control.growers))
    setdiff(w6.addedC.growers, w6.control.growers)  #In week 6, there were 19 taxa that grew in the presence of added carbon but did not grow in the control treatment (no added carbon)
    length(setdiff(w6.control.growers, w6.glucose.growers))
    setdiff(w6.control.growers, w6.glucose.growers)  #In week 6, there were 229 taxa that grew in the control treatment (no added carbon) but did not grow on added glucose
    length(setdiff(w6.glucose.growers, w6.control.growers))
    setdiff(w6.glucose.growers, w6.control.growers)  #In week 6, there were 12 taxa that grew on added glucose but did not grow in the control treatment (no added carbon)
  #}


#_4a_Calculate sum of carbon flux across taxa (to see how this compares to measured rates of soil respiration):
  #{
  #(units of flux in the all.comparisons data.frame are in pg C g soil-1 day-1)
  #First, a quick and dirty way: look at sum of bootstrapped medians:
    #NOTE: flux units are converted to ngC g soil-1 day-1 below:
    Tcompare
    tapply(all.comparisons$f.boot.median, all.comparisons$comparisonID, sum, na.rm=TRUE)*(1/1000)

  #Next, construct CI's for the total C flux:
    #Read in the bootstrapped C fluxes from the all.taxa.calcs output:
    C.flux.boots <- read.table("qSIP_output/BH_bootstrapped_C_fluxes.txt", header=TRUE, sep="")
    names(C.flux.boots)
    dim(C.flux.boots)
  
    inds1 <- C.flux.boots$comparisonID == 1            #growth in control treatment (week 1)
    inds2 <- C.flux.boots$comparisonID == 2            #growth with priming (week 1)
    inds3 <- C.flux.boots$comparisonID == 3            #growth on glucose with priming (week 1)
    inds5 <- C.flux.boots$comparisonID == 5            #growth in control treatment (week 6)
    inds6 <- C.flux.boots$comparisonID == 6            #growth with priming (week 6)
    inds7 <- C.flux.boots$comparisonID == 7            #growth on glucose with priming (week 6)

    #Create a dataframe for all assemblage-wide C-fluxes (ng C g soil-1 day-1):
    C.flux.assemblage <- data.frame(comparisonID=numeric(6), percentile05=numeric(6), percentile50=numeric(6), percentile95=numeric(6))

    #Week 1, Control, 18O
    C.f.boots.1 <- as.numeric(apply(C.flux.boots[inds1,5:dim(C.flux.boots)[2]], 2, sum, na.rm=TRUE))
    C.flux.assemblage[1,2:4] <- (1/1000)*quantile(C.f.boots.1, probs=c(0.05, 0.50, 0.95), na.rm=TRUE)   #NOTE: flux units are converted to ngC g soil-1 day-1 here
    C.flux.assemblage[1,1] <- unique(C.flux.boots$comparisonID[inds1])
  
    #Week 1, added carbon, 18O
    C.f.boots.2 <- as.numeric(apply(C.flux.boots[inds2,5:dim(C.flux.boots)[2]], 2, sum, na.rm=TRUE))
    C.flux.assemblage[2,2:4] <- (1/1000)*quantile(C.f.boots.2, probs=c(0.05, 0.50, 0.95), na.rm=TRUE)   #NOTE: flux units are converted to ngC g soil-1 day-1 here
    C.flux.assemblage[2,1] <- unique(C.flux.boots$comparisonID[inds2])

    #Week 1, added carbon, 13C
    C.f.boots.3 <- as.numeric(apply(C.flux.boots[inds3,5:dim(C.flux.boots)[2]], 2, sum, na.rm=TRUE))
    C.flux.assemblage[3,2:4] <- (1/1000)*quantile(C.f.boots.3, probs=c(0.05, 0.50, 0.95), na.rm=TRUE)   #NOTE: flux units are converted to ngC g soil-1 day-1 here
    C.flux.assemblage[3,1] <- unique(C.flux.boots$comparisonID[inds3])

    #Week 6, Control, 18O
    C.f.boots.5 <- as.numeric(apply(C.flux.boots[inds5,5:dim(C.flux.boots)[2]], 2, sum, na.rm=TRUE))
    C.flux.assemblage[4,2:4] <- (1/1000)*quantile(C.f.boots.5, probs=c(0.05, 0.50, 0.95), na.rm=TRUE)   #NOTE: flux units are converted to ngC g soil-1 day-1 here
    C.flux.assemblage[4,1] <- unique(C.flux.boots$comparisonID[inds5])

    #Week 6, added carbon, 18O
    C.f.boots.6 <- as.numeric(apply(C.flux.boots[inds6,5:dim(C.flux.boots)[2]], 2, sum, na.rm=TRUE))
    C.flux.assemblage[5,2:4] <- (1/1000)*quantile(C.f.boots.6, probs=c(0.05, 0.50, 0.95), na.rm=TRUE)   #NOTE: flux units are converted to ngC g soil-1 day-1 here
    C.flux.assemblage[5,1] <- unique(C.flux.boots$comparisonID[inds6])

    #Week 6, added carbon, 13C
    C.f.boots.7 <- as.numeric(apply(C.flux.boots[inds7,5:dim(C.flux.boots)[2]], 2, sum, na.rm=TRUE))
    C.flux.assemblage[6,2:4] <- (1/1000)*quantile(C.f.boots.7, probs=c(0.05, 0.50, 0.95), na.rm=TRUE)   #NOTE: flux units are converted to ngC g soil-1 day-1 here
    C.flux.assemblage[6,1] <- unique(C.flux.boots$comparisonID[inds7])

    C.flux.assemblage <- data.frame(comparisonID=C.flux.assemblage[,1], week=c(1,1,1,6,6,6), tmt=c("control","added carbon","added carbon","control","added carbon","added carbon"), isotope=c("18O","18O","13C","18O","18O","13C"), C.flux.assemblage[,2:4])
    C.flux.assemblage

    #Write the results (C.flux.assemblage) to a text file:
    write.table(C.flux.assemblage, "qSIP_output/BH_assemblage_C_fluxes.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
  #}


#_4b_Calculate average carbon flux across all taxa (across all bootstrapped estimates):
  #{
    #(units of flux in the C.flux.boots data.frame are in pg C g soil-1 day-1)
      mean(unlist(C.flux.boots[C.flux.boots$comparisonID == 1, (5:dim(C.flux.boots)[2])]), na.rm=TRUE)
      mean(unlist(C.flux.boots[C.flux.boots$comparisonID == 2, (5:dim(C.flux.boots)[2])]), na.rm=TRUE)
      median(unlist(C.flux.boots[C.flux.boots$comparisonID == 1, (5:dim(C.flux.boots)[2])]), na.rm=TRUE)
      median(unlist(C.flux.boots[C.flux.boots$comparisonID == 2, (5:dim(C.flux.boots)[2])]), na.rm=TRUE)
  #}


#_4c_Calculate indirect effect of added carbon on C flux:
  #{
    #increase in carbon flux due to added glucose minus the 13C-based carbon flux
    #(units of flux in the C.flux.boots data.frame are in pg C g soil-1 day-1)
      #increase in carbon flux due to added glucose:
      C.flux.boots[C.flux.boots$comparisonID == 2, (5:dim(C.flux.boots)[2])] - C.flux.boots[C.flux.boots$comparisonID == 1, (5:dim(C.flux.boots)[2])]
      #13C-based carbon flux:
      C.flux.boots[C.flux.boots$comparisonID == 3, (5:dim(C.flux.boots)[2])]
      #indirect effect of added carbon on C flux:
      C.flux.indirect.taxa <- data.frame(taxonID=rep(NA, length(C.flux.boots$taxonID[C.flux.boots$comparisonID == 1])), f.mean=rep(NA, length(C.flux.boots$taxonID[C.flux.boots$comparisonID == 1])), f.median=rep(NA, length(C.flux.boots$taxonID[C.flux.boots$comparisonID == 1])), f.CI.05=rep(NA, length(C.flux.boots$taxonID[C.flux.boots$comparisonID == 1])), f.CI.95=rep(NA, length(C.flux.boots$taxonID[C.flux.boots$comparisonID == 1])))
      C.flux.indirect.taxa$taxonID <- C.flux.boots$taxonID[C.flux.boots$comparisonID == 1]
      #mean for all taxa:
      C.flux.indirect.taxa$f.mean <- as.numeric(apply((C.flux.boots[C.flux.boots$comparisonID == 2, (5:dim(C.flux.boots)[2])] - C.flux.boots[C.flux.boots$comparisonID == 1, (5:dim(C.flux.boots)[2])]) - C.flux.boots[C.flux.boots$comparisonID == 3, (5:dim(C.flux.boots)[2])], 1, mean, na.rm=TRUE))
      #median for all taxa:
      C.flux.indirect.taxa$f.median <- as.numeric(apply((C.flux.boots[C.flux.boots$comparisonID == 2, (5:dim(C.flux.boots)[2])] - C.flux.boots[C.flux.boots$comparisonID == 1, (5:dim(C.flux.boots)[2])]) - C.flux.boots[C.flux.boots$comparisonID == 3, (5:dim(C.flux.boots)[2])], 1, median, na.rm=TRUE))
      #CIs for all taxa:
      C.flux.indirect.taxa$f.CI.05 <- as.numeric(apply((C.flux.boots[C.flux.boots$comparisonID == 2, (5:dim(C.flux.boots)[2])] - C.flux.boots[C.flux.boots$comparisonID == 1, (5:dim(C.flux.boots)[2])]) - C.flux.boots[C.flux.boots$comparisonID == 3, (5:dim(C.flux.boots)[2])], 1, quantile, probs=0.05, na.rm=TRUE))
      C.flux.indirect.taxa$f.CI.95 <- as.numeric(apply((C.flux.boots[C.flux.boots$comparisonID == 2, (5:dim(C.flux.boots)[2])] - C.flux.boots[C.flux.boots$comparisonID == 1, (5:dim(C.flux.boots)[2])]) - C.flux.boots[C.flux.boots$comparisonID == 3, (5:dim(C.flux.boots)[2])], 1, quantile, probs=0.95, na.rm=TRUE))
      C.flux.indirect.taxa

      #Create phylum means & medians of taxon medians (ROUGH):
      C.flux.indirect.phyla <- data.frame(phylum=rep(NA, length(levels(data.melted$phylum))), f.mean=rep(NA, length(levels(data.melted$phylum))), f.median=rep(NA, length(levels(data.melted$phylum))))
      for (p in 1:length(levels(data.melted$phylum))){
        taxa.in.phylum <- data.melted$taxon[data.melted$phylum == levels(data.melted$phylum)[p]]
        C.flux.indirect.phyla$phylum[p] <- levels(data.melted$phylum)[p]
        C.flux.indirect.phyla$f.mean[p] <- mean(C.flux.indirect.taxa$f.median[C.flux.indirect.taxa$taxonID %in% taxa.in.phylum], na.rm=TRUE)
        C.flux.indirect.phyla$f.median[p] <- median(C.flux.indirect.taxa$f.median[C.flux.indirect.taxa$taxonID %in% taxa.in.phylum], na.rm=TRUE)
      }
      C.flux.indirect.phyla
      plot(x=c(1:length(C.flux.indirect.phyla$phylum)), y=C.flux.indirect.phyla$f.median[order(C.flux.indirect.phyla$f.median, decreasing=TRUE)])

      dev.off()
  #}


#_5a_Calculate increment in growth rate on soil organic matter over baseline growth
  #{
    #This quantity is calculated as: (18O12C vs 16O12C) - (16O13C va 16O12C) - (18O vs 16O)
      #Bootstrapped C fluxes from the all.taxa.calcs output were already read in above (in step 4a):
      #Read in the bootstrapped r's from the all.taxa.calcs output:
      r.boots <- read.table("qSIP_output/BH_bootstrapped_r.txt", header=TRUE, sep="")
      names(r.boots)
      dim(r.boots)
      Tcompare

      #Read in the bootstrapped aef's from the all.taxa.calcs output:
      ape.boots <- read.table("qSIP_output/BH_bootstrapped_ape.txt", header=TRUE, sep="")
      names(ape.boots)
      dim(ape.boots)
      Tcompare
      
      #These are the comparisonID's to do for Week 1 & Week 6:
      #  2 - 3 - 1   #Week 1
      #  6 - 7 - 5   #Week 6
      
      #First calculate the desired quantity for ape (atom excess fraction):
        #Get indices for these comparisons
        inds1 <- ape.boots$comparisonID == 1            #growth in control treatment (week 1)
        inds2 <- ape.boots$comparisonID == 2            #growth with priming (week 1)
        inds3 <- ape.boots$comparisonID == 3            #growth on glucose with priming (week 1)
        inds5 <- ape.boots$comparisonID == 5            #growth in control treatment (week 6)
        inds6 <- ape.boots$comparisonID == 6            #growth with priming (week 6)
        inds7 <- ape.boots$comparisonID == 7            #growth on glucose with priming (week 6)

        #Create data frames to store results:
        SOM.growth.incr.wk1 <- data.frame(taxonID=factor(r.boots$taxonID[inds2]), ape.boot.median=numeric(length(r.boots$taxonID[inds2])), ape.boot.CI.L=numeric(length(r.boots$taxonID[inds2])), ape.boot.CI.U=numeric(length(r.boots$taxonID[inds2])), r.boot.median=numeric(length(r.boots$taxonID[inds2])), r.boot.CI.L=numeric(length(r.boots$taxonID[inds2])), r.boot.CI.U=numeric(length(r.boots$taxonID[inds2])), f.boot.median=numeric(length(r.boots$taxonID[inds2])), f.boot.CI.L=numeric(length(r.boots$taxonID[inds2])), f.boot.CI.U=numeric(length(r.boots$taxonID[inds2])))
        SOM.growth.incr.wk6 <- data.frame(taxonID=factor(r.boots$taxonID[inds6]), ape.boot.median=numeric(length(r.boots$taxonID[inds6])), ape.boot.CI.L=numeric(length(r.boots$taxonID[inds6])), ape.boot.CI.U=numeric(length(r.boots$taxonID[inds6])), r.boot.median=numeric(length(r.boots$taxonID[inds6])), r.boot.CI.L=numeric(length(r.boots$taxonID[inds6])), r.boot.CI.U=numeric(length(r.boots$taxonID[inds6])), f.boot.median=numeric(length(r.boots$taxonID[inds6])), f.boot.CI.L=numeric(length(r.boots$taxonID[inds6])), f.boot.CI.U=numeric(length(r.boots$taxonID[inds6])))

        #Ensure that taxa are in the same order for the three comparison IDs:
        sum(ape.boots$taxonID[inds2] != ape.boots$taxonID[inds3])
        sum(ape.boots$taxonID[inds3] != ape.boots$taxonID[inds1])
        sum(ape.boots$taxonID[inds6] != ape.boots$taxonID[inds7])
        sum(ape.boots$taxonID[inds7] != ape.boots$taxonID[inds5])

        #Week 1 (increment in growth rate on soil organic matter nover baseline growth):
        SOM.growth.incr.ape.wk1.boots <- ape.boots[inds2, 5:dim(ape.boots)[2]] - ape.boots[inds3, 5:dim(ape.boots)[2]] - ape.boots[inds1, 5:dim(ape.boots)[2]]
        SOM.growth.incr.wk1$ape.boot.median <- as.numeric(apply(SOM.growth.incr.ape.wk1.boots, 1, median, na.rm=TRUE))
        SOM.growth.incr.wk1$ape.boot.CI.L <- as.numeric(apply(SOM.growth.incr.ape.wk1.boots, 1, quantile, probs=0.05, na.rm=TRUE))
        SOM.growth.incr.wk1$ape.boot.CI.U <- as.numeric(apply(SOM.growth.incr.ape.wk1.boots, 1, quantile, probs=0.95, na.rm=TRUE))

        #Week 6 (increment in growth rate on soil organic matter nover baseline growth):
        SOM.growth.incr.ape.wk6.boots <- ape.boots[inds6, 5:dim(ape.boots)[2]] - ape.boots[inds7, 5:dim(ape.boots)[2]] - ape.boots[inds5, 5:dim(ape.boots)[2]]
        SOM.growth.incr.wk6$ape.boot.median <- as.numeric(apply(SOM.growth.incr.ape.wk6.boots, 1, median, na.rm=TRUE))
        SOM.growth.incr.wk6$ape.boot.CI.L <- as.numeric(apply(SOM.growth.incr.ape.wk6.boots, 1, quantile, probs=0.05, na.rm=TRUE))
        SOM.growth.incr.wk6$ape.boot.CI.U <- as.numeric(apply(SOM.growth.incr.ape.wk6.boots, 1, quantile, probs=0.95, na.rm=TRUE))


      #Next calculate the desired quantity for r (growth rate):
        #Get indices for these comparisons
        inds1 <- r.boots$comparisonID == 1            #growth in control treatment (week 1)
        inds2 <- r.boots$comparisonID == 2            #growth with priming (week 1)
        inds3 <- r.boots$comparisonID == 3            #growth on glucose with priming (week 1)
        inds5 <- r.boots$comparisonID == 5            #growth in control treatment (week 6)
        inds6 <- r.boots$comparisonID == 6            #growth with priming (week 6)
        inds7 <- r.boots$comparisonID == 7            #growth on glucose with priming (week 6)

        #Ensure that taxa are in the same order for the three comparison IDs:
        sum(r.boots$taxonID[inds2] != r.boots$taxonID[inds3])
        sum(r.boots$taxonID[inds3] != r.boots$taxonID[inds1])
        sum(r.boots$taxonID[inds6] != r.boots$taxonID[inds7])
        sum(r.boots$taxonID[inds7] != r.boots$taxonID[inds5])

        #Week 1 (increment in growth rate on soil organic matter nover baseline growth):
        SOM.growth.incr.r.wk1.boots <- r.boots[inds2, 5:dim(r.boots)[2]] - r.boots[inds3, 5:dim(r.boots)[2]] - r.boots[inds1, 5:dim(r.boots)[2]]
        SOM.growth.incr.wk1$r.boot.median <- as.numeric(apply(SOM.growth.incr.r.wk1.boots, 1, median, na.rm=TRUE))
        SOM.growth.incr.wk1$r.boot.CI.L <- as.numeric(apply(SOM.growth.incr.r.wk1.boots, 1, quantile, probs=0.05, na.rm=TRUE))
        SOM.growth.incr.wk1$r.boot.CI.U <- as.numeric(apply(SOM.growth.incr.r.wk1.boots, 1, quantile, probs=0.95, na.rm=TRUE))

        #Week 6 (increment in growth rate on soil organic matter nover baseline growth):
        SOM.growth.incr.r.wk6.boots <- r.boots[inds6, 5:dim(r.boots)[2]] - r.boots[inds7, 5:dim(r.boots)[2]] - r.boots[inds5, 5:dim(r.boots)[2]]
        SOM.growth.incr.wk6$r.boot.median <- as.numeric(apply(SOM.growth.incr.r.wk6.boots, 1, median, na.rm=TRUE))
        SOM.growth.incr.wk6$r.boot.CI.L <- as.numeric(apply(SOM.growth.incr.r.wk6.boots, 1, quantile, probs=0.05, na.rm=TRUE))
        SOM.growth.incr.wk6$r.boot.CI.U <- as.numeric(apply(SOM.growth.incr.r.wk6.boots, 1, quantile, probs=0.95, na.rm=TRUE))


      #Next calculate the desired quantity for f (C flux):
        #Get indices for these comparisons
        inds1 <- C.flux.boots$comparisonID == 1            #growth in control treatment (week 1)
        inds2 <- C.flux.boots$comparisonID == 2            #growth with priming (week 1)
        inds3 <- C.flux.boots$comparisonID == 3            #growth on glucose with priming (week 1)
        inds5 <- C.flux.boots$comparisonID == 5            #growth in control treatment (week 6)
        inds6 <- C.flux.boots$comparisonID == 6            #growth with priming (week 6)
        inds7 <- C.flux.boots$comparisonID == 7            #growth on glucose with priming (week 6)

        #Ensure that taxa are in the same order for the three comparison IDs:
        sum(C.flux.boots$taxonID[inds2] != C.flux.boots$taxonID[inds3])
        sum(C.flux.boots$taxonID[inds3] != C.flux.boots$taxonID[inds1])
        sum(C.flux.boots$taxonID[inds6] != C.flux.boots$taxonID[inds7])
        sum(C.flux.boots$taxonID[inds7] != C.flux.boots$taxonID[inds5])

        #Week 1 (increment in growth rate on soil organic matter nover baseline growth):
        SOM.growth.incr.f.wk1.boots <- C.flux.boots[inds2, 5:dim(C.flux.boots)[2]] - C.flux.boots[inds3, 5:dim(C.flux.boots)[2]] - C.flux.boots[inds1, 5:dim(C.flux.boots)[2]]
        SOM.growth.incr.wk1$f.boot.median <- as.numeric(apply(SOM.growth.incr.f.wk1.boots, 1, median, na.rm=TRUE))
        SOM.growth.incr.wk1$f.boot.CI.L <- as.numeric(apply(SOM.growth.incr.f.wk1.boots, 1, quantile, probs=0.05, na.rm=TRUE))
        SOM.growth.incr.wk1$f.boot.CI.U <- as.numeric(apply(SOM.growth.incr.f.wk1.boots, 1, quantile, probs=0.95, na.rm=TRUE))

        #Week 6 (increment in growth rate on soil organic matter nover baseline growth):
        SOM.growth.incr.f.wk6.boots <- C.flux.boots[inds6, 5:dim(C.flux.boots)[2]] - C.flux.boots[inds7, 5:dim(C.flux.boots)[2]] - C.flux.boots[inds5, 5:dim(C.flux.boots)[2]]
        SOM.growth.incr.wk6$f.boot.median <- as.numeric(apply(SOM.growth.incr.f.wk6.boots, 1, median, na.rm=TRUE))
        SOM.growth.incr.wk6$f.boot.CI.L <- as.numeric(apply(SOM.growth.incr.f.wk6.boots, 1, quantile, probs=0.05, na.rm=TRUE))
        SOM.growth.incr.wk6$f.boot.CI.U <- as.numeric(apply(SOM.growth.incr.f.wk6.boots, 1, quantile, probs=0.95, na.rm=TRUE))

        #Write the results (SOM.growth.incr.wk1 & SOM.growth.incr.wk6) to text files:
        write.table(SOM.growth.incr.wk1, "qSIP_output/BH_SOM_growth_increment_Week1.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
        write.table(SOM.growth.incr.wk6, "qSIP_output/BH_SOM_growth_increment_Week6.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
  #}


#_5b_Graph the increment in growth rate on soil organic matter over baseline growth:
  #{
    #(NOTE THAT THERE IS ONE TAXON WHICH HAS PHYLUM == NA AND IS THUS NOT PLOTTED ON THIS FIGURE)
    #Growth rate on SOM (increment over baseline growth rate) - in units of r and C flux:
      dev.off()
      # dev.new(width=7.5, height=7.5)
      pdf(file="qSIP_output/Figures/BH_Taxa_same_order_PhylumGroups_SOM_growth_increment.pdf", width=7.5, height=7.5)
      par(mfcol=c(2,2))
      par(mai=c(0.44,0.42,0.2,0.05), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
      phyla <- levels(data.melted$phylum)
      phyla.cols <- cols.for.phyla$col
      ranks <- order(SOM.growth.incr.wk1$r.boot.median)
      tax.order <- data.frame(axis.loc=seq(1, N.taxa), ranks)
      x.min.r <- min(c(SOM.growth.incr.wk1$r.boot.CI.L, SOM.growth.incr.wk6$r.boot.CI.L), na.rm=TRUE)
      x.max.r <- max(c(SOM.growth.incr.wk1$r.boot.CI.U, SOM.growth.incr.wk6$r.boot.CI.U), na.rm=TRUE)
      x.min.f <- min(c(SOM.growth.incr.wk1$f.boot.CI.L, SOM.growth.incr.wk6$f.boot.CI.L), na.rm=TRUE)/1000
      x.max.f <- max(c(SOM.growth.incr.wk1$f.boot.CI.U, SOM.growth.incr.wk6$f.boot.CI.U), na.rm=TRUE)/1000
      
      r.add.panel <- function(DATA.WK, week){
        plot(y=1:N.taxa, x=DATA.WK$r.boot.median[tax.order$ranks], type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min.r, x.max.r), main="")
        mtext(paste("Week ", week, sep=""), side=3, line=-0.1, cex=0.75)
        counter <- 1
        for (p in 1:length(phyla)){
          curr.comp.ranked <- DATA.WK[tax.order$ranks,]
          mids <- curr.comp.ranked$r.boot.median[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
          lowers <- curr.comp.ranked$r.boot.CI.L[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
          uppers <- curr.comp.ranked$r.boot.CI.U[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
          tax.nums <- counter:(counter+length(mids)-1)
          counter <- counter+length(mids)
          arrows(x0=lowers, y0=tax.nums, x1=uppers, y1=tax.nums, length=0, angle=90, code=3, col=as.character(phyla.cols[p], lwd=0.5))
          points(x=mids, y=tax.nums, pch=21, cex=0.25, lwd=0.125, col="gray30", bg=as.character(phyla.cols[p]))
        }
        abline(v=0, col="black")
        par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
        axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
        # mtext(expression(paste("r (day"^-1, ")", sep="")), side=1, line=1.4, cex=0.75)
        mtext(expression(paste("increment over baseline in growth on SOM (day"^-1, ")", sep="")), side=1, line=1.4, cex=0.65)
        par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
        axis(side=2, at=1:N.taxa, labels=rep(NA, 103), tck=-0.015, las=1, cex.axis=0.6)
        mtext("Taxon", side=2, line=1.35, cex=0.75)
      }
      
      f.add.panel <- function(DATA.WK, week){
        plot(y=1:N.taxa, x=DATA.WK$f.boot.median[tax.order$ranks]/1000, type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min.f, x.max.f), main="")
        mtext(paste("Week ", week, sep=""), side=3, line=-0.1, cex=0.75)
        counter <- 1
        for (p in 1:length(phyla)){
          curr.comp.ranked <- DATA.WK[tax.order$ranks,]
          mids <- curr.comp.ranked$f.boot.median[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
          lowers <- curr.comp.ranked$f.boot.CI.L[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
          uppers <- curr.comp.ranked$f.boot.CI.U[as.character(curr.comp.ranked$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
          tax.nums <- counter:(counter+length(mids)-1)
          counter <- counter+length(mids)
          arrows(x0=lowers/1000, y0=tax.nums, x1=uppers/1000, y1=tax.nums, length=0, angle=90, code=3, col=as.character(phyla.cols[p]), lwd=0.5)
          points(x=mids/1000, y=tax.nums, pch=21, cex=0.25, lwd=0.125, col="gray30", bg=as.character(phyla.cols[p]))
        }
        abline(v=0, col="black")
        par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
        axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
        # mtext(expression(paste("C flux into biomass (ng C g soil"^-1, " day"^-1, ")", sep="")), side=1, line=1.4, cex=0.75)
        mtext(expression(paste("increment over baseline in growth on SOM (ng C g soil"^-1, " day"^-1, ")", sep="")), side=1, line=1.4, cex=0.65)
        par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
        axis(side=2, at=1:N.taxa, labels=rep(NA, 103), tck=-0.015, las=1, cex.axis=0.6)
        mtext("Taxon", side=2, line=1.35, cex=0.75)
      }

      r.add.panel(DATA.WK=SOM.growth.incr.wk1, week = "1")
      r.add.panel(DATA.WK=SOM.growth.incr.wk6, week = "6")
      f.add.panel(DATA.WK=SOM.growth.incr.wk1, week = "1")
      f.add.panel(DATA.WK=SOM.growth.incr.wk6, week = "6")
      par(xpd=NA)
      legend(x=7, y=0.82*N.taxa, xjust=0.5, yjust=1, legend=phyla[1:ceiling(0.5*length(phyla))], bty="n", text.col=as.character(phyla.cols)[1:ceiling(0.5*length(phyla))], cex=0.4)
      legend(x=13, y=0.82*N.taxa, xjust=0.5, yjust=1, legend=phyla[(ceiling(0.5*length(phyla))+1):length(phyla)], bty="n", text.col=as.character(phyla.cols)[(ceiling(0.5*length(phyla))+1):length(phyla)], cex=0.4)
      par(xpd=FALSE)
      rm(r.add.panel)
      rm(f.add.panel)
      par(mfrow=c(1,1))

      dev.off()
  #}


#_6a_Calculate increase in excess atom fraction 18O caused by glucose addition:
  #{
    #This quantity is calculated as: (18O12C vs 16O12C) - (18O vs 16O)
      #Bootstrapped aef's from the all.taxa.calcs output were already read in above (in step 5a):
      names(ape.boots)
      dim(ape.boots)
      Tcompare
    
      #These are the comparisonID's to do for Week 1 & Week 6:
      #  2 - 1   #Week 1
      #  6 - 5   #Week 6
    
      #First calculate the desired quantity for r (growth rate):
        #Get indices for these comparisons
        inds1 <- ape.boots$comparisonID == 1            #growth in control treatment (week 1)
        inds2 <- ape.boots$comparisonID == 2            #growth with priming (week 1)
        inds3 <- ape.boots$comparisonID == 3            #growth on glucose with priming (week 1)
        inds5 <- ape.boots$comparisonID == 5            #growth in control treatment (week 6)
        inds6 <- ape.boots$comparisonID == 6            #growth with priming (week 6)
        inds7 <- ape.boots$comparisonID == 7            #growth on glucose with priming (week 6)

        #Create data frames to store results:
        Gluc.ape.incr.wk1 <- data.frame(taxonID=factor(ape.boots$taxonID[inds2]), ape.boot.median=numeric(length(ape.boots$taxonID[inds2])), ape.boot.CI.L=numeric(length(ape.boots$taxonID[inds2])), ape.boot.CI.U=numeric(length(ape.boots$taxonID[inds2])))
        Gluc.ape.incr.wk6 <- data.frame(taxonID=factor(ape.boots$taxonID[inds6]), ape.boot.median=numeric(length(ape.boots$taxonID[inds2])), ape.boot.CI.L=numeric(length(ape.boots$taxonID[inds2])), ape.boot.CI.U=numeric(length(ape.boots$taxonID[inds2])))

        #Ensure that taxa are in the same order for the two comparison IDs:
        sum(ape.boots$taxonID[inds2] != ape.boots$taxonID[inds1])
        sum(ape.boots$taxonID[inds6] != ape.boots$taxonID[inds5])

        #Week 1 (increment in growth rate on soil organic matter nover baseline growth):
        Gluc.ape.incr.wk1.boots <- ape.boots[inds2, 5:dim(ape.boots)[2]] - ape.boots[inds1, 5:dim(ape.boots)[2]]
        Gluc.ape.incr.wk1$ape.boot.median <- as.numeric(apply(Gluc.ape.incr.wk1.boots, 1, median, na.rm=TRUE))
        Gluc.ape.incr.wk1$ape.boot.CI.L <- as.numeric(apply(Gluc.ape.incr.wk1.boots, 1, quantile, probs=0.05, na.rm=TRUE))
        Gluc.ape.incr.wk1$ape.boot.CI.U <- as.numeric(apply(Gluc.ape.incr.wk1.boots, 1, quantile, probs=0.95, na.rm=TRUE))

        #Week 6 (increment in growth rate on soil organic matter nover baseline growth):
        Gluc.ape.incr.wk6.boots <- ape.boots[inds6, 5:dim(ape.boots)[2]] - ape.boots[inds5, 5:dim(ape.boots)[2]]
        Gluc.ape.incr.wk6$ape.boot.median <- as.numeric(apply(Gluc.ape.incr.wk6.boots, 1, median, na.rm=TRUE))
        Gluc.ape.incr.wk6$ape.boot.CI.L <- as.numeric(apply(Gluc.ape.incr.wk6.boots, 1, quantile, probs=0.05, na.rm=TRUE))
        Gluc.ape.incr.wk6$ape.boot.CI.U <- as.numeric(apply(Gluc.ape.incr.wk6.boots, 1, quantile, probs=0.95, na.rm=TRUE))


        #Write the results (Gluc.ape.incr.wk1 & Gluc.ape.incr.wk6) to text files:
        write.table(Gluc.ape.incr.wk1, "qSIP_output/BH_Gluc_ape_increment_Week1.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
        write.table(Gluc.ape.incr.wk6, "qSIP_output/BH_Gluc_ape_increment_Week6.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
  #}


#_6b_Graph the increase in excess atom fraction 18O caused by glucose addition vs. excess atom fraction 13C (full axes):
  #{
    #EAF increment on glucose (18O) vs EAF 13C:
    
      #(NOTE THAT THERE IS ONE TAXON WHICH HAS PHYLUM == NA AND IS THUS NOT PLOTTED ON THIS FIGURE)
      #First get the values for excess atom fraction 13C:
      ape.13C.wk1 <- data.frame(taxonID=factor(ape.boots$taxonID[inds3]), ape.boot.median=numeric(length(ape.boots$taxonID[inds3])), ape.boot.CI.L=numeric(length(ape.boots$taxonID[inds3])), ape.boot.CI.U=numeric(length(ape.boots$taxonID[inds3])))
      ape.13C.wk1.boots <- ape.boots[inds3, 5:dim(ape.boots)[2]]
      ape.13C.wk1$ape.boot.median <- as.numeric(apply(ape.13C.wk1.boots, 1, median, na.rm=TRUE))
      ape.13C.wk1$ape.boot.CI.L <- as.numeric(apply(ape.13C.wk1.boots, 1, quantile, probs=0.05, na.rm=TRUE))
      ape.13C.wk1$ape.boot.CI.U <- as.numeric(apply(ape.13C.wk1.boots, 1, quantile, probs=0.95, na.rm=TRUE))
    
      #Make sure taxon IDs match between the 18O and 13C data:
      sum(ape.13C.wk1$taxonID == Gluc.ape.incr.wk1$taxonID)
      
      dev.off()
      dimensions <- c((8.7/2.54), (((((8.7/2.54)-0.49-0.70)+0.40+0.08)*2.54)/2.54))
      # dev.new(width=dimensions[1], height=dimensions[2])
      pdf(file="qSIP_output/Figures/BH_Gluc_ape_increment_18O_vs_ape_13C_phyla.pdf", width=dimensions[1], height=dimensions[2])
      phyla <- levels(data.melted$phylum)
      # phyla.cols <- rainbow(length(phyla))
      phyla.cols <- cols.for.phyla$col
      # x.min <- min(ape.13C.wk1$ape.boot.CI.L, na.rm=TRUE)
      # x.max <- max(ape.13C.wk1$ape.boot.CI.U, na.rm=TRUE)
      # y.min <- min(Gluc.ape.incr.wk1$ape.boot.CI.L, na.rm=TRUE)
      # y.max <- max(Gluc.ape.incr.wk1$ape.boot.CI.U, na.rm=TRUE)
      x.min <- -0.4
      x.max <- 0.62
      y.min <- -0.4
      y.max <- 0.62

      par(mai=c(0.40,0.49,0.08,0.70), cex=1, mex=0.75)
      plot(y=Gluc.ape.incr.wk1$ape.boot.CI.U, x=ape.13C.wk1$ape.boot.CI.U, type="n", bty="n", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min, x.max), ylim=c(y.min,y.max), main="")
      for (p in 1:length(phyla)){
        x.pts <- ape.13C.wk1$ape.boot.median[as.character(ape.13C.wk1$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        y.pts <- Gluc.ape.incr.wk1$ape.boot.median[as.character(Gluc.ape.incr.wk1$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        lowers.13C <- ape.13C.wk1$ape.boot.CI.L[as.character(ape.13C.wk1$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        uppers.13C <- ape.13C.wk1$ape.boot.CI.U[as.character(ape.13C.wk1$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        lowers.18O <- Gluc.ape.incr.wk1$ape.boot.CI.L[as.character(Gluc.ape.incr.wk1$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        uppers.18O <- Gluc.ape.incr.wk1$ape.boot.CI.U[as.character(Gluc.ape.incr.wk1$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        arrows(x0=lowers.13C, y0=y.pts, x1=uppers.13C, y1=y.pts, length=0, angle=90, code=3, col=as.character(phyla.cols)[p], lwd=0.85)
        arrows(x0=x.pts, y0=lowers.18O, x1=x.pts, y1=uppers.18O, length=0, angle=90, code=3, col=as.character(phyla.cols)[p], lwd=0.85)
        points(x=x.pts, y=y.pts, pch=21, cex=0.6, lwd=0.25, col="gray30", bg=as.character(phyla.cols)[p])
      }
      axis(side=1, at=c(-x.max, x.max)*1000, labels=FALSE, tck=0)
      axis(side=2, at=c(-y.max, y.max)*1000, labels=FALSE, tck=0)
      abline(a=0, b=0.33, col="black", lty=1)

      par(mgp=c(3,-0.08,0))				  #set spacing for the x-axis labels (the second value)
      axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
      par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
      axis(side=2, tck=-0.015, las=1, cex.axis=0.6)
      mtext(expression(paste("Increase in atom fraction excess "^18, "O", sep="")), side=2, line=1.85, cex=0.75)
      mtext(expression(paste("Atom fraction excess "^13, "C", sep="")), side=1, line=1.2, cex=0.75)

      par(xpd=NA)
      legend(x=x.max, y=y.max+((y.max-y.min)*0.09), xjust=0, yjust=1, legend=phyla, bty="n", text.col=as.character(phyla.cols), cex=0.4)
      par(xpd=FALSE)

      dev.off()
  #}


#_6c_Graph the increase in excess atom fraction 18O caused by glucose addition vs. excess atom fraction 13C (cropped axes):
  #{
    #EAF increment on glucose (18O) vs EAF 13C:
  
      #(NOTE THAT THERE IS ONE TAXON WHICH HAS PHYLUM == NA AND IS THUS NOT PLOTTED ON THIS FIGURE)
      #First get the values for excess atom fraction 13C:
      ape.13C.wk1 <- data.frame(taxonID=factor(ape.boots$taxonID[inds3]), ape.boot.median=numeric(length(ape.boots$taxonID[inds3])), ape.boot.CI.L=numeric(length(ape.boots$taxonID[inds3])), ape.boot.CI.U=numeric(length(ape.boots$taxonID[inds3])))
      ape.13C.wk1.boots <- ape.boots[inds3, 5:dim(ape.boots)[2]]
      ape.13C.wk1$ape.boot.median <- as.numeric(apply(ape.13C.wk1.boots, 1, median, na.rm=TRUE))
      ape.13C.wk1$ape.boot.CI.L <- as.numeric(apply(ape.13C.wk1.boots, 1, quantile, probs=0.05, na.rm=TRUE))
      ape.13C.wk1$ape.boot.CI.U <- as.numeric(apply(ape.13C.wk1.boots, 1, quantile, probs=0.95, na.rm=TRUE))
  
      #Make sure taxon IDs match between the 18O and 13C data:
      sum(ape.13C.wk1$taxonID == Gluc.ape.incr.wk1$taxonID)
    
      dev.off()
      dimensions <- c((8.7/2.54), (((((8.7/2.54)-0.49-0.70)+0.40+0.08)*2.54)/2.54))
      # dev.new(width=dimensions[1], height=dimensions[2])
      pdf(file="qSIP_output/Figures/BH_Gluc_ape_increment_18O_vs_ape_13C_phyla_cropped.pdf", width=dimensions[1], height=dimensions[2])
      phyla <- levels(data.melted$phylum)
      # phyla.cols <- rainbow(length(phyla))
      phyla.cols <- cols.for.phyla$col
      # x.min <- min(ape.13C.wk1$ape.boot.CI.L, na.rm=TRUE)
      # x.max <- max(ape.13C.wk1$ape.boot.CI.U, na.rm=TRUE)
      # y.min <- min(Gluc.ape.incr.wk1$ape.boot.CI.L, na.rm=TRUE)
      # y.max <- max(Gluc.ape.incr.wk1$ape.boot.CI.U, na.rm=TRUE)
      x.min <- 0
      x.max <- 0.62
      y.min <- 0
      y.max <- 0.62

      par(mai=c(0.40,0.49,0.08,0.70), cex=1, mex=0.75)
      plot(y=Gluc.ape.incr.wk1$ape.boot.CI.U, x=ape.13C.wk1$ape.boot.CI.U, type="n", bty="n", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min, x.max), ylim=c(y.min,y.max), main="")
      for (p in 1:length(phyla)){
        x.pts <- ape.13C.wk1$ape.boot.median[as.character(ape.13C.wk1$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        y.pts <- Gluc.ape.incr.wk1$ape.boot.median[as.character(Gluc.ape.incr.wk1$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        lowers.13C <- ape.13C.wk1$ape.boot.CI.L[as.character(ape.13C.wk1$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        uppers.13C <- ape.13C.wk1$ape.boot.CI.U[as.character(ape.13C.wk1$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        lowers.18O <- Gluc.ape.incr.wk1$ape.boot.CI.L[as.character(Gluc.ape.incr.wk1$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        uppers.18O <- Gluc.ape.incr.wk1$ape.boot.CI.U[as.character(Gluc.ape.incr.wk1$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        arrows(x0=lowers.13C, y0=y.pts, x1=uppers.13C, y1=y.pts, length=0, angle=90, code=3, col=as.character(phyla.cols)[p], lwd=0.85)
        arrows(x0=x.pts, y0=lowers.18O, x1=x.pts, y1=uppers.18O, length=0, angle=90, code=3, col=as.character(phyla.cols)[p], lwd=0.85)
        points(x=x.pts, y=y.pts, pch=21, cex=0.6, lwd=0.25, col="gray30", bg=as.character(phyla.cols)[p])
      }
      axis(side=1, at=c(-x.max, x.max)*1000, labels=FALSE, tck=0)
      axis(side=2, at=c(-y.max, y.max)*1000, labels=FALSE, tck=0)
      abline(a=0, b=0.33, col="black", lty=1)

      par(mgp=c(3,-0.08,0))				  #set spacing for the x-axis labels (the second value)
      axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
      par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
      axis(side=2, tck=-0.015, las=1, cex.axis=0.6)
      mtext(expression(paste("Increase in atom fraction excess "^18, "O", sep="")), side=2, line=1.85, cex=0.75)
      mtext(expression(paste("Atom fraction excess "^13, "C", sep="")), side=1, line=1.2, cex=0.75)

      par(xpd=NA)
      legend(x=x.max, y=y.max+((y.max-y.min)*0.09), xjust=0, yjust=1, legend=phyla, bty="n", text.col=as.character(phyla.cols), cex=0.4)
      par(xpd=FALSE)

      dev.off()
  #}


#_6d_Graph the increase in excess atom fraction 18O caused by glucose addition vs. excess atom fraction 13C (full axes) - B&W (no phyla):
  #{
    #EAF increment on glucose (18O) vs EAF 13C:
  
      #(NOTE THAT THERE IS ONE TAXON WHICH HAS PHYLUM == NA AND IS THUS NOT PLOTTED ON THIS FIGURE)
      #First get the values for excess atom fraction 13C:
      ape.13C.wk1 <- data.frame(taxonID=factor(ape.boots$taxonID[inds3]), ape.boot.median=numeric(length(ape.boots$taxonID[inds3])), ape.boot.CI.L=numeric(length(ape.boots$taxonID[inds3])), ape.boot.CI.U=numeric(length(ape.boots$taxonID[inds3])))
      ape.13C.wk1.boots <- ape.boots[inds3, 5:dim(ape.boots)[2]]
      ape.13C.wk1$ape.boot.median <- as.numeric(apply(ape.13C.wk1.boots, 1, median, na.rm=TRUE))
      ape.13C.wk1$ape.boot.CI.L <- as.numeric(apply(ape.13C.wk1.boots, 1, quantile, probs=0.05, na.rm=TRUE))
      ape.13C.wk1$ape.boot.CI.U <- as.numeric(apply(ape.13C.wk1.boots, 1, quantile, probs=0.95, na.rm=TRUE))
  
      #Make sure taxon IDs match between the 18O and 13C data:
      sum(ape.13C.wk1$taxonID == Gluc.ape.incr.wk1$taxonID)
    
      dev.off()
      dimensions <- c((8.7/2.54), (((((8.7/2.54)-0.49-0.05)+0.40+0.08)*2.54)/2.54))
      # dev.new(width=dimensions[1], height=dimensions[2])
      pdf(file="qSIP_output/Figures/BH_Gluc_ape_increment_18O_vs_ape_13C.pdf", width=dimensions[1], height=dimensions[2])
      phyla <- levels(data.melted$phylum)
      # phyla.cols <- rainbow(length(phyla))
      # phyla.cols <- cols.for.phyla$col
      # x.min <- min(ape.13C.wk1$ape.boot.CI.L, na.rm=TRUE)
      # x.max <- max(ape.13C.wk1$ape.boot.CI.U, na.rm=TRUE)
      # y.min <- min(Gluc.ape.incr.wk1$ape.boot.CI.L, na.rm=TRUE)
      # y.max <- max(Gluc.ape.incr.wk1$ape.boot.CI.U, na.rm=TRUE)
      x.min <- -0.4
      x.max <- 0.62
      y.min <- -0.4
      y.max <- 0.62

      par(mai=c(0.40,0.49,0.08,0.05), cex=1, mex=0.75)
      plot(y=Gluc.ape.incr.wk1$ape.boot.CI.U, x=ape.13C.wk1$ape.boot.CI.U, type="n", bty="n", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min, x.max), ylim=c(y.min,y.max), main="")
      for (p in 1:length(phyla)){
        x.pts <- ape.13C.wk1$ape.boot.median[as.character(ape.13C.wk1$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        y.pts <- Gluc.ape.incr.wk1$ape.boot.median[as.character(Gluc.ape.incr.wk1$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        lowers.13C <- ape.13C.wk1$ape.boot.CI.L[as.character(ape.13C.wk1$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        uppers.13C <- ape.13C.wk1$ape.boot.CI.U[as.character(ape.13C.wk1$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        lowers.18O <- Gluc.ape.incr.wk1$ape.boot.CI.L[as.character(Gluc.ape.incr.wk1$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        uppers.18O <- Gluc.ape.incr.wk1$ape.boot.CI.U[as.character(Gluc.ape.incr.wk1$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        arrows(x0=lowers.13C, y0=y.pts, x1=uppers.13C, y1=y.pts, length=0, angle=90, code=3, col="gray20", lwd=0.85)
        arrows(x0=x.pts, y0=lowers.18O, x1=x.pts, y1=uppers.18O, length=0, angle=90, code=3, col="gray20", lwd=0.85)
        points(x=x.pts, y=y.pts, pch=21, cex=0.6, lwd=0.5, col="gray20", bg="white")
      }
      axis(side=1, at=c(-x.max, x.max)*1000, labels=FALSE, tck=0)
      axis(side=2, at=c(-y.max, y.max)*1000, labels=FALSE, tck=0)
      abline(a=0, b=0.33, col="red", lty=1)

      par(mgp=c(3,-0.08,0))				  #set spacing for the x-axis labels (the second value)
      axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
      par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
      axis(side=2, tck=-0.015, las=1, cex.axis=0.6)
      mtext(expression(paste("Increase in atom fraction excess "^18, "O", sep="")), side=2, line=1.85, cex=0.75)
      mtext(expression(paste("Atom fraction excess "^13, "C", sep="")), side=1, line=1.2, cex=0.75)

      # par(xpd=NA)
      # legend(x=x.max, y=y.max+((y.max-y.min)*0.09), xjust=0, yjust=1, legend=phyla, bty="n", text.col=as.character(phyla.cols), cex=0.4)
      # par(xpd=FALSE)

      dev.off()
  #}


#_6e_Graph the increase in excess atom fraction 18O caused by glucose addition vs. excess atom fraction 13C (cropped axes) - B&W (no phyla):
  #{
    #EAF increment on glucose (18O) vs EAF 13C:

      #(NOTE THAT THERE IS ONE TAXON WHICH HAS PHYLUM == NA AND IS THUS NOT PLOTTED ON THIS FIGURE)
      #First get the values for excess atom fraction 13C:
      ape.13C.wk1 <- data.frame(taxonID=factor(ape.boots$taxonID[inds3]), ape.boot.median=numeric(length(ape.boots$taxonID[inds3])), ape.boot.CI.L=numeric(length(ape.boots$taxonID[inds3])), ape.boot.CI.U=numeric(length(ape.boots$taxonID[inds3])))
      ape.13C.wk1.boots <- ape.boots[inds3, 5:dim(ape.boots)[2]]
      ape.13C.wk1$ape.boot.median <- as.numeric(apply(ape.13C.wk1.boots, 1, median, na.rm=TRUE))
      ape.13C.wk1$ape.boot.CI.L <- as.numeric(apply(ape.13C.wk1.boots, 1, quantile, probs=0.05, na.rm=TRUE))
      ape.13C.wk1$ape.boot.CI.U <- as.numeric(apply(ape.13C.wk1.boots, 1, quantile, probs=0.95, na.rm=TRUE))

      #Make sure taxon IDs match between the 18O and 13C data:
      sum(ape.13C.wk1$taxonID == Gluc.ape.incr.wk1$taxonID)
  
      dev.off()
      dimensions <- c((8.7/2.54), (((((8.7/2.54)-0.49-0.05)+0.40+0.08)*2.54)/2.54))
      # dev.new(width=dimensions[1], height=dimensions[2])
      pdf(file="qSIP_output/Figures/BH_Gluc_ape_increment_18O_vs_ape_13C_cropped.pdf", width=dimensions[1], height=dimensions[2])
      phyla <- levels(data.melted$phylum)
      # phyla.cols <- rainbow(length(phyla))
      # phyla.cols <- cols.for.phyla$col
      # x.min <- min(ape.13C.wk1$ape.boot.CI.L, na.rm=TRUE)
      # x.max <- max(ape.13C.wk1$ape.boot.CI.U, na.rm=TRUE)
      # y.min <- min(Gluc.ape.incr.wk1$ape.boot.CI.L, na.rm=TRUE)
      # y.max <- max(Gluc.ape.incr.wk1$ape.boot.CI.U, na.rm=TRUE)
      x.min <- 0
      x.max <- 0.62
      y.min <- 0
      y.max <- 0.62

      par(mai=c(0.40,0.49,0.08,0.05), cex=1, mex=0.75)
      plot(y=Gluc.ape.incr.wk1$ape.boot.CI.U, x=ape.13C.wk1$ape.boot.CI.U, type="n", bty="n", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min, x.max), ylim=c(y.min,y.max), main="")
      for (p in 1:length(phyla)){
        x.pts <- ape.13C.wk1$ape.boot.median[as.character(ape.13C.wk1$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        y.pts <- Gluc.ape.incr.wk1$ape.boot.median[as.character(Gluc.ape.incr.wk1$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        lowers.13C <- ape.13C.wk1$ape.boot.CI.L[as.character(ape.13C.wk1$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        uppers.13C <- ape.13C.wk1$ape.boot.CI.U[as.character(ape.13C.wk1$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        lowers.18O <- Gluc.ape.incr.wk1$ape.boot.CI.L[as.character(Gluc.ape.incr.wk1$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        uppers.18O <- Gluc.ape.incr.wk1$ape.boot.CI.U[as.character(Gluc.ape.incr.wk1$taxonID) %in% as.character(taxa.id$taxon[taxa.id$phylum == phyla[p]])]
        arrows(x0=lowers.13C, y0=y.pts, x1=uppers.13C, y1=y.pts, length=0, angle=90, code=3, col="gray20", lwd=0.85)
        arrows(x0=x.pts, y0=lowers.18O, x1=x.pts, y1=uppers.18O, length=0, angle=90, code=3, col="gray20", lwd=0.85)
        points(x=x.pts, y=y.pts, pch=21, cex=0.6, lwd=0.5, col="gray20", bg="white")
      }
      axis(side=1, at=c(-x.max, x.max)*1000, labels=FALSE, tck=0)
      axis(side=2, at=c(-y.max, y.max)*1000, labels=FALSE, tck=0)
      abline(a=0, b=0.33, col="red", lty=1)

      par(mgp=c(3,-0.08,0))				  #set spacing for the x-axis labels (the second value)
      axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
      par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
      axis(side=2, tck=-0.015, las=1, cex.axis=0.6)
      mtext(expression(paste("Increase in atom fraction excess "^18, "O", sep="")), side=2, line=1.85, cex=0.75)
      mtext(expression(paste("Atom fraction excess "^13, "C", sep="")), side=1, line=1.2, cex=0.75)

      # par(xpd=NA)
      # legend(x=x.max, y=y.max+((y.max-y.min)*0.09), xjust=0, yjust=1, legend=phyla, bty="n", text.col=as.character(phyla.cols), cex=0.4)
      # par(xpd=FALSE)

      dev.off()
  #}
