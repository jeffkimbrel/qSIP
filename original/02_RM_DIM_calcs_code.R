# This code performs the basic qSIP calculations on the Dimensions week 1 qSIP data sets
# There are 20 comparisons, so we will subset the comparisons per ecosystem (5 comparisons each)
# This code only performs grassland (GL) calculations


graphics.off()	#close all graphics windows


#Set working directory:
  #ALREADY DONE FOR THIS WORKSPACE IN PREVIOUSLY RUN CODE; only reset it here if loading previously saved workspace (see below)


#Reload the saved workspace resulting from the previous script:
  setwd("/Users/bk/Research/Projects/SIP_Modeling/qSIP")
  load("qSIP_workspaces/RM_DIM_01/.RData")


#For the raw, uncorrected data: Run all wad.diff, ape, r, & flux calculations for all taxa and all comparisons:
  #NOTES: only wad and ape results are valid; r, f (C fluxes) results are not calculated here (soil data is not specified); note that incubation duration values in Tcompare are made-up to enable running the 'all.taxa.calcs' script (but the made-up durations do not affect any of the results that are returned)
  #       using set.seed ensures that subsequent function calls utilize the same random number seed, so results are replicated exactly)
  #       using system.time calculates how long it took the function to run
  #       set prop.O.from.water=0.33 according to McHugh et al. unpublished data from E. coli pure cultures incubated with 18O-H2O
  #       rename the bootstrapped output files so that they are not overwritten if 'all.taxa.calcs' is re-run
  set.seed(100)
  system.time(all.comparisons.raw <- all.taxa.calcs(X.all=data.melted, comparisons=Tcompare.GL, vars=c("taxon", "Density.g.ml", "copies.ul", "Sample", "iso.treat.eco", "DNA.ng.ul"), growth.model="exponential", prop.O.from.water=0.33, v.frac=50, copies.cell=6, pgC.cell=0.1, CI=0.90, draws=1000, tailed.test=1))
  bootstrapped.filenames <- paste("qSIP_output/", c("bootstrapped_wad1.txt", "bootstrapped_wad2.txt", "bootstrapped_wad_diff.txt", "bootstrapped_ape.txt"), sep="")
  new.bootstrapped.filenames <- paste("qSIP_output/", "RM_DIM_GL_", c("bootstrapped_wad1_raw.txt", "bootstrapped_wad2_raw.txt", "bootstrapped_wad_diff_raw.txt", "bootstrapped_ape_raw.txt"), sep="")
  file.rename(from=bootstrapped.filenames, to=new.bootstrapped.filenames)
  summary(all.comparisons.raw)
  dim(all.comparisons.raw)


#Write the results (all.comparisons) to a text file:
  dir.create(path=paste(getwd(), "/qSIP_output", sep=""), showWarnings=FALSE)
  write.table(all.comparisons.raw, "qSIP_output/RM_DIM_GL_all_comparisons_raw.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)


#Write the taxa.id dataframe to a text file:
  write.table(taxa.id, "qSIP_output/RM_DIM_GL_taxa_ID_raw.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)



#For the corrected data: Run all wad.diff, ape, r, & flux calculations for all taxa and all comparisons:
  #NOTES: only wad and ape results are valid; r, f (C fluxes) results are not calculated here (soil data is not specified); note that incubation duration values in Tcompare are made-up to enable running the 'all.taxa.calcs' script (but the made-up durations do not affect any of the results that are returned)
  #       using set.seed ensures that subsequent function calls utilize the same random number seed, so results are replicated exactly)
  #       using system.time calculates how long it took the function to run
  #       set prop.O.from.water=0.33 according to McHugh et al. unpublished data from E. coli pure cultures incubated with 18O-H2O
  #       rename the bootstrapped output files so that they are not overwritten if 'all.taxa.calcs' is re-run
  set.seed(100)
  system.time(all.comparisons <- all.taxa.calcs(X.all=data.corr.melted, comparisons=Tcompare.GL, vars=c("taxon", "Density.g.ml", "copies.ul", "Sample", "iso.treat.eco", "DNA.ng.ul"), growth.model="exponential", prop.O.from.water=0.33, v.frac=50, copies.cell=6, pgC.cell=0.1, CI=0.90, draws=1000, tailed.test=1))
  bootstrapped.filenames <- paste("qSIP_output/", c("bootstrapped_wad1.txt", "bootstrapped_wad2.txt", "bootstrapped_wad_diff.txt", "bootstrapped_ape.txt"), sep="")
  new.bootstrapped.filenames <- paste("qSIP_output/", "RM_DIM_GL_", c("bootstrapped_wad1.txt", "bootstrapped_wad2.txt", "bootstrapped_wad_diff.txt", "bootstrapped_ape.txt"), sep="")
  file.rename(from=bootstrapped.filenames, to=new.bootstrapped.filenames)
  summary(all.comparisons)
  dim(all.comparisons)


#Write the results (all.comparisons) to a text file:
  dir.create(path=paste(getwd(), "/qSIP_output", sep=""), showWarnings=FALSE)
  write.table(all.comparisons, "qSIP_output/RM_DIM_GL_all_comparisons.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)


#Write the taxa.id dataframe to a text file:
  write.table(taxa.id, "qSIP_output/RM_DIM_GL_taxa_ID.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)



#Compare EAF plots for the raw, uncorrected data and for the corrected data:
  #Raw, uncorrected data:
  #excess atom fraction (atom percent excess, expressed as a fraction rather than as a percentage):
    dev.off()
    # dev.new(width=10, height=7.5)
    pdf(file="qSIP_output/Figures/RM_DIM_GL_Taxa_resorted_NoPhyla_ape_raw.pdf", width=10, height=7.5)
    par(mfrow=c(2,3))
    par(mai=c(0.44,0.42,0.2,0.05), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
    x.min <- min(all.comparisons.raw$ape.boot.CI.L, na.rm=TRUE)
    x.max <- max(all.comparisons.raw$ape.boot.CI.U, na.rm=TRUE)
    N.taxa <- length(levels(factor(all.comparisons.raw$taxonID)))

    ape.add.panel <- function(comparisonID){
      inds <- all.comparisons.raw$comparisonID == comparisonID
      ranks <- order(all.comparisons.raw$ape.boot.median[inds])
      plot(1:N.taxa~all.comparisons.raw$ape.boot.median[inds][ranks], type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min, x.max), main="")
      mtext(paste(all.comparisons.raw$trt.code.1[inds][1], "v", all.comparisons.raw$trt.code.2[inds][1], sep=" "), side=3, line=-0.1, cex=0.75)
      points(x=all.comparisons.raw$ape.boot.median[inds][ranks], y=1:N.taxa, pch=21, cex=0.6, bg="black")
      arrows(x0=all.comparisons.raw$ape.boot.CI.L[inds][ranks], y0=1:N.taxa, x1=all.comparisons.raw$ape.boot.CI.U[inds][ranks], y1=1:N.taxa, length=0, angle=90, code=3)
      pts.with.message <- c(1:N.taxa)[all.comparisons.raw$message[inds][ranks] != "none"]
      # Mark taxa with messages in the results with red arrows (commented-out below):
      # arrows(x0=all.comparisons.raw$ape.boot.CI.U[inds][ranks][pts.with.message]+((x.max-x.min)*0.01), y0=pts.with.message, x1=all.comparisons.raw$ape.boot.CI.U[inds][ranks][pts.with.message]+((x.max-x.min)*0.1), y1=pts.with.message, length=0.06, angle=30, code=1, col="red")
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
    ape.add.panel(comparisonID=4)
    ape.add.panel(comparisonID=5)
    rm(ape.add.panel)
    par(mfrow=c(1,1))

    dev.off()


  #Corrected data:
  #excess atom fraction (atom percent excess, expressed as a fraction rather than as a percentage):
    dev.off()
    # dev.new(width=10, height=7.5)
    pdf(file="qSIP_output/Figures/RM_DIM_GL_Taxa_resorted_NoPhyla_ape_corr.pdf", width=10, height=7.5)
    par(mfrow=c(2,3))
    par(mai=c(0.44,0.42,0.2,0.05), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
    x.min <- min(all.comparisons$ape.boot.CI.L, na.rm=TRUE)
    x.max <- max(all.comparisons$ape.boot.CI.U, na.rm=TRUE)
    N.taxa <- length(levels(factor(all.comparisons$taxonID)))

    ape.add.panel <- function(comparisonID){
      inds <- all.comparisons$comparisonID == comparisonID
      ranks <- order(all.comparisons$ape.boot.median[inds])
      plot(1:N.taxa~all.comparisons$ape.boot.median[inds][ranks], type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min, x.max), main="")
      mtext(paste(all.comparisons$trt.code.1[inds][1], "v", all.comparisons$trt.code.2[inds][1], sep=" "), side=3, line=-0.1, cex=0.75)
      points(x=all.comparisons$ape.boot.median[inds][ranks], y=1:N.taxa, pch=21, cex=0.6, bg="black")
      arrows(x0=all.comparisons$ape.boot.CI.L[inds][ranks], y0=1:N.taxa, x1=all.comparisons$ape.boot.CI.U[inds][ranks], y1=1:N.taxa, length=0, angle=90, code=3)
      pts.with.message <- c(1:N.taxa)[all.comparisons$message[inds][ranks] != "none"]
      # Mark taxa with messages in the results with red arrows (commented-out below):
      # arrows(x0=all.comparisons$ape.boot.CI.U[inds][ranks][pts.with.message]+((x.max-x.min)*0.01), y0=pts.with.message, x1=all.comparisons$ape.boot.CI.U[inds][ranks][pts.with.message]+((x.max-x.min)*0.1), y1=pts.with.message, length=0.06, angle=30, code=1, col="red")
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
    ape.add.panel(comparisonID=4)
    ape.add.panel(comparisonID=5)
    rm(ape.add.panel)
    par(mfrow=c(1,1))

    dev.off()




#Save workspace:
  save(list=ls(all.names=TRUE), file="qSIP_workspaces/RM_DIM_01-02/.RData", envir=.GlobalEnv)



