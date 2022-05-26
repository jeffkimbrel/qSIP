# This code performs the basic qSIP calculations on Bri's soil type qSIP data sets
# There are 18 comparisons, and we will perform the comparisons per soil type (6 comparisons each)
# This code only performs soil type "AN" calculations


graphics.off()	#close all graphics windows


#Set working directory:
  #ALREADY DONE FOR THIS WORKSPACE IN PREVIOUSLY RUN CODE; only reset it here if loading previously saved workspace (see below)


#Reload the saved workspace resulting from the previous script:
  setwd("/Users/bk/Research/Projects/SIP_Modeling/qSIP")
  load("qSIP_workspaces/BF_01/.RData")


#For the corrected data: Run all wad.diff, ape, r, & flux calculations for all taxa and all comparisons:
  #NOTES: only wad and ape results are valid; r, f (C fluxes) results are not calculated here (soil data is not specified)
  #       using set.seed ensures that subsequent function calls utilize the same random number seed, so results are replicated exactly)
  #       using system.time calculates how long it took the function to run
  #       set prop.O.from.water=0.33 according to McHugh et al. unpublished data from E. coli pure cultures incubated with 18O-H2O
  #       rename the bootstrapped output files so that they are not overwritten if 'all.taxa.calcs' is re-run
  set.seed(100)
  system.time(all.comparisons.1 <- all.taxa.calcs(X.all=data.melted.1, comparisons=Tcompare1, vars=c("taxon", "Density_g_ml", "copies.ul", "Tube", "ComboTrt", "DNA_ng_uL"), growth.model="exponential", prop.O.from.water=0.33, v.frac=50, copies.cell=6, pgC.cell=0.1, CI=0.90, draws=1000, tailed.test=1))
  bootstrapped.filenames <- paste("qSIP_output/", c("bootstrapped_wad1.txt", "bootstrapped_wad2.txt", "bootstrapped_wad_diff.txt", "bootstrapped_ape.txt"), sep="")
  new.bootstrapped.filenames <- paste("qSIP_output/", "BF_AN_", c("bootstrapped_wad1.txt", "bootstrapped_wad2.txt", "bootstrapped_wad_diff.txt", "bootstrapped_ape.txt"), sep="")
  file.rename(from=bootstrapped.filenames, to=new.bootstrapped.filenames)
  summary(all.comparisons.1)
  dim(all.comparisons.1)


#Write the results (all.comparisons.1) to a text file:
  dir.create(path=paste(getwd(), "/qSIP_output", sep=""), showWarnings=FALSE)
  write.table(all.comparisons.1, "qSIP_output/BF_AN_all_comparisons.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)


#Write the taxa.id dataframe to a text file:
  write.table(taxa.id, "qSIP_output/BF_AN_taxa_ID.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)



#EAF plots ('AN' soil type):
  #Corrected data:
  #excess atom fraction (atom percent excess, expressed as a fraction rather than as a percentage):
    dev.off()
    # dev.new(width=10, height=7.5)
    pdf(file="qSIP_output/Figures/BF_AN_Taxa_resorted_NoPhyla_ape_corr.pdf", width=10, height=7.5)
    par(mfrow=c(2,3))
    par(mai=c(0.44,0.42,0.2,0.05), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
    x.min <- min(all.comparisons.1$ape.boot.CI.L, na.rm=TRUE)
    x.max <- max(all.comparisons.1$ape.boot.CI.U, na.rm=TRUE)
    N.taxa <- length(levels(factor(all.comparisons.1$taxonID)))

    ape.add.panel <- function(DATA, comparisonID){
      inds <- DATA$comparisonID == comparisonID
      ranks <- order(DATA$ape.boot.median[inds])
      plot(1:N.taxa~DATA$ape.boot.median[inds][ranks], type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min, x.max), main="")
      mtext(paste(DATA$trt.code.1[inds][1], "v", DATA$trt.code.2[inds][1], sep=" "), side=3, line=-0.1, cex=0.75)
      points(x=DATA$ape.boot.median[inds][ranks], y=1:N.taxa, pch=21, cex=0.6, bg="black")
      arrows(x0=DATA$ape.boot.CI.L[inds][ranks], y0=1:N.taxa, x1=DATA$ape.boot.CI.U[inds][ranks], y1=1:N.taxa, length=0, angle=90, code=3)
      pts.with.message <- c(1:N.taxa)[DATA$message[inds][ranks] != "none"]
      # Mark taxa with messages in the results with red arrows (commented-out below):
      # arrows(x0=DATA$ape.boot.CI.U[inds][ranks][pts.with.message]+((x.max-x.min)*0.01), y0=pts.with.message, x1=DATA$ape.boot.CI.U[inds][ranks][pts.with.message]+((x.max-x.min)*0.1), y1=pts.with.message, length=0.06, angle=30, code=1, col="red")
      abline(v=0, col="red")
      par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
      axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
      mtext("Excess atom fraction", side=1, line=1.4, cex=0.75)
      par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
      axis(side=2, at=1:N.taxa, labels=rep(NA, 103), tck=-0.015, las=1, cex.axis=0.6)
      mtext("Taxon", side=2, line=1.35, cex=0.75)
    }

    ape.add.panel(DATA=all.comparisons.1, comparisonID=1)
    ape.add.panel(DATA=all.comparisons.1, comparisonID=2)
    ape.add.panel(DATA=all.comparisons.1, comparisonID=3)
    ape.add.panel(DATA=all.comparisons.1, comparisonID=4)
    ape.add.panel(DATA=all.comparisons.1, comparisonID=5)
    ape.add.panel(DATA=all.comparisons.1, comparisonID=6)
    rm(ape.add.panel)
    par(mfrow=c(1,1))

    dev.off()




#Save workspace:
  save(list=ls(all.names=TRUE), file="qSIP_workspaces/BF_01-02/.RData", envir=.GlobalEnv)



