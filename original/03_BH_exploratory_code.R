# This code performs various calculations & analyses that test various parameter assumptions & explore how to improve the models:


#_Linear vs exponential growth:
  #{
  #Compare linear growth estimates with exponential growth estimates:
    set.seed(100)
    all.comparisons.lin <- all.taxa.calcs(X.all=data.melted, comparisons=Tcompare, M.soil=Sdat, vars=c("taxon", "density.g.ml", "copies", "tube", "trt.code", "g.soil"), growth.model="linear", prop.O.from.water=0.33, v.frac=50, copies.cell=6, pgC.cell=0.1, CI=0.90, draws=1000, tailed.test=1)
    bootstrapped.filenames <- paste("qSIP_output/", c("bootstrapped_wad1.txt", "bootstrapped_wad2.txt", "bootstrapped_wad_diff.txt", "bootstrapped_ape.txt", "bootstrapped_r.txt", "bootstrapped_C_fluxes.txt"), sep="")
    new.bootstrapped.filenames <- paste("qSIP_output/", "BH_linear_", c("bootstrapped_wad1.txt", "bootstrapped_wad2.txt", "bootstrapped_wad_diff.txt", "bootstrapped_ape.txt", "bootstrapped_r.txt", "bootstrapped_C_fluxes.txt"), sep="")
    file.rename(from=bootstrapped.filenames, to=new.bootstrapped.filenames)

    #Compare estimates for each treatment:
    all.comparisons.exp <- all.comparisons
    comparisons.exp.1 <- all.comparisons.exp[all.comparisons.exp$comparisonID == 1,]
    comparisons.lin.1 <- all.comparisons.lin[all.comparisons.lin$comparisonID == 1,]
    comparisons.exp.2 <- all.comparisons.exp[all.comparisons.exp$comparisonID == 2,]
    comparisons.lin.2 <- all.comparisons.lin[all.comparisons.lin$comparisonID == 2,]
    comparisons.exp.3 <- all.comparisons.exp[all.comparisons.exp$comparisonID == 3,]
    comparisons.lin.3 <- all.comparisons.lin[all.comparisons.lin$comparisonID == 3,]
    comparisons.exp.5 <- all.comparisons.exp[all.comparisons.exp$comparisonID == 5,]
    comparisons.lin.5 <- all.comparisons.lin[all.comparisons.lin$comparisonID == 5,]
    comparisons.exp.6 <- all.comparisons.exp[all.comparisons.exp$comparisonID == 6,]
    comparisons.lin.6 <- all.comparisons.lin[all.comparisons.lin$comparisonID == 6,]
    comparisons.exp.7 <- all.comparisons.exp[all.comparisons.exp$comparisonID == 7,]
    comparisons.lin.7 <- all.comparisons.lin[all.comparisons.lin$comparisonID == 7,]

    min(all.comparisons.exp$r.boot.median, na.rm=TRUE)
    max(all.comparisons.exp$r.boot.median, na.rm=TRUE)
    min(all.comparisons.lin$r.boot.median, na.rm=TRUE)
    max(all.comparisons.lin$r.boot.median, na.rm=TRUE)

    pdf(file="qSIP_output/Figures/BH_Linear_vs_ExponentialGrowth.pdf", width=10, height=7.5)
    par(mfrow=c(2,3))
    par(mai=c(1,1,0.2,0.05))        #set the margins (bottom, left, top, right) in inches
    plot(x=comparisons.lin.1$r.boot.median, y=comparisons.exp.1$r.boot.median, bty="l", type="p", pch=21, bg="black", xlim=c(-30000, 750000), ylim=c(-0.05, 0.5), xlab=expression(paste("r (copies day"^-1, "g soil"^-1, ") linear model", sep="")), ylab=expression(paste("r (day"^-1, ") exponential model", sep="")))
    mtext(paste(comparisons.exp.1$trt.code.1[1], "v", comparisons.exp.1$trt.code.2[1], sep=" "), side=3, line=-0.1, cex=0.75)
    plot(x=comparisons.lin.2$r.boot.median, y=comparisons.exp.2$r.boot.median, bty="l", type="p", pch=21, bg="black", xlim=c(-30000, 750000), ylim=c(-0.05, 0.5), xlab=expression(paste("r (copies day"^-1, "g soil"^-1, ") linear model", sep="")), ylab=expression(paste("r (day"^-1, ") exponential model", sep="")))
    mtext(paste(comparisons.exp.2$trt.code.1[1], "v", comparisons.exp.2$trt.code.2[1], sep=" "), side=3, line=-0.1, cex=0.75)
    plot(x=comparisons.lin.3$r.boot.median, y=comparisons.exp.3$r.boot.median, bty="l", type="p", pch=21, bg="black", xlim=c(-30000, 750000), ylim=c(-0.05, 0.5), xlab=expression(paste("r (copies day"^-1, "g soil"^-1, ") linear model", sep="")), ylab=expression(paste("r (day"^-1, ") exponential model", sep="")))
    mtext(paste(comparisons.exp.3$trt.code.1[1], "v", comparisons.exp.3$trt.code.2[1], sep=" "), side=3, line=-0.1, cex=0.75)
    plot(x=comparisons.lin.5$r.boot.median, y=comparisons.exp.5$r.boot.median, bty="l", type="p", pch=21, bg="black", xlim=c(-30000, 750000), ylim=c(-0.05, 0.5), xlab=expression(paste("r (copies day"^-1, "g soil"^-1, ") linear model", sep="")), ylab=expression(paste("r (day"^-1, ") exponential model", sep="")))
    mtext(paste(comparisons.exp.5$trt.code.1[1], "v", comparisons.exp.5$trt.code.2[1], sep=" "), side=3, line=-0.1, cex=0.75)
    plot(x=comparisons.lin.6$r.boot.median, y=comparisons.exp.6$r.boot.median, bty="l", type="p", pch=21, bg="black", xlim=c(-30000, 750000), ylim=c(-0.05, 0.5), xlab=expression(paste("r (copies day"^-1, "g soil"^-1, ") linear model", sep="")), ylab=expression(paste("r (day"^-1, ") exponential model", sep="")))
    mtext(paste(comparisons.exp.6$trt.code.1[1], "v", comparisons.exp.6$trt.code.2[1], sep=" "), side=3, line=-0.1, cex=0.75)
    plot(x=comparisons.lin.7$r.boot.median, y=comparisons.exp.7$r.boot.median, bty="l", type="p", pch=21, bg="black", xlim=c(-30000, 750000), ylim=c(-0.05, 0.5), xlab=expression(paste("r (copies day"^-1, "g soil"^-1, ") linear model", sep="")), ylab=expression(paste("r (day"^-1, ") exponential model", sep="")))
    mtext(paste(comparisons.exp.7$trt.code.1[1], "v", comparisons.exp.7$trt.code.2[1], sep=" "), side=3, line=-0.1, cex=0.75)
    dev.off()
  #}




#_Sensitivity to the proportion of oxygen from water parameter:
  #{
  #Compare 18O growth results to the case where prop.O.from.water=0.50:
    set.seed(100)
    all.comparisons.5 <- all.taxa.calcs(X.all=data.melted, comparisons=Tcompare, M.soil=Sdat, vars=c("taxon", "density.g.ml", "copies", "tube", "trt.code", "g.soil"), growth.model="exponential", prop.O.from.water=0.50, v.frac=50, copies.cell=6, pgC.cell=0.1, CI=0.90, draws=1000, tailed.test=1)
      #Delete the bootstrapped output files created in the call above:
      bootstrapped.filenames <- paste("qSIP_output/", c("bootstrapped_wad1.txt", "bootstrapped_wad2.txt", "bootstrapped_wad_diff.txt", "bootstrapped_ape.txt", "bootstrapped_r.txt", "bootstrapped_C_fluxes.txt"), sep="")
      unlink(bootstrapped.filenames)
  
    all.comparisons.3 <- all.comparisons
    #Get only the 18O comparisons:
    comparisons.18O.3 <- all.comparisons.3[all.comparisons.3$comparisonID %in% c(1,2,5,6),]
    comparisons.18O.5 <- all.comparisons.5[all.comparisons.5$comparisonID %in% c(1,2,5,6),]

    min(c(comparisons.18O.3$r.boot.median, comparisons.18O.5$r.boot.median), na.rm=TRUE)
    max(c(comparisons.18O.3$r.boot.median, comparisons.18O.5$r.boot.median), na.rm=TRUE)
    min(c(comparisons.18O.3$r.boot.CI.L, comparisons.18O.5$r.boot.CI.L), na.rm=TRUE)
    max(c(comparisons.18O.3$r.boot.CI.U, comparisons.18O.5$r.boot.CI.U), na.rm=TRUE)
  
    pdf(file="qSIP_output/Figures/BH_ProportionOxygenFromWater_33_vs_50.pdf", width=6, height=6)
    par(mai=c(1.02, 1.02, 0.82, 0.42))
    plot(x=comparisons.18O.5$r.boot.median, y=comparisons.18O.3$r.boot.median, bty="l", type="p", pch=21, bg="black", xlim=c(-0.08, 0.6), ylim=c(-0.08, 0.6), xlab=expression(paste("r (day"^-1, ") using 0.50", sep="")), ylab=expression(paste("r (day"^-1, ") using 0.33", sep="")), main="18O growth rates (median & CI of bootstrap estimates)")
    # arrows(x0=comparisons.18O.5$r.boot.CI.L, y0=comparisons.18O.3$r.boot.median, x1=comparisons.18O.5$r.boot.CI.U, y1=comparisons.18O.3$r.boot.median, length=0, angle=90, code=3, col="black")
    # arrows(x0=comparisons.18O.5$r.boot.median, y0=comparisons.18O.3$r.boot.CI.L, x1=comparisons.18O.5$r.boot.median, y1=comparisons.18O.3$r.boot.CI.U, length=0, angle=90, code=3, col="black")
    abline(a=0, b=1, col="red")
    dev.off()


  #Look at median and max growth results across a range of values for prop.O.from.water (for Week1 priming treatment 18O growth rate only):
    #First make sure that r.obs is reflective of r.boot.median:
    min(c(all.comparisons$r.obs, all.comparisons$r.boot.median), na.rm=TRUE)
    max(c(all.comparisons$r.obs, all.comparisons$r.boot.median), na.rm=TRUE)
    plot(x=all.comparisons$r.obs, y=all.comparisons$r.boot.median, bty="l", type="p", pch=21, col="blue", bg="blue", xlim=c(-0.06,0.60), ylim=c(-0.06,0.60))
    abline(a=0, b=1)

    POFW <- data.frame(prop.O.from.water=seq(0.1, 1, 0.1), r.median=numeric(10), r.max=numeric(10))
    for (i in 1:dim(POFW)[1]){
      set.seed(100)
      comparisons <- all.taxa.calcs(X.all=data.melted, comparisons=Tcompare[2,], M.soil=Sdat, vars=c("taxon", "density.g.ml", "copies", "tube", "trt.code", "g.soil"), growth.model="exponential", prop.O.from.water=POFW$prop.O.from.water[i], v.frac=50, copies.cell=6, pgC.cell=0.1, CI=0.90, draws=1, tailed.test=1)
        #Delete the bootstrapped output files created in the call above:
        bootstrapped.filenames <- paste("qSIP_output/", c("bootstrapped_wad1.txt", "bootstrapped_wad2.txt", "bootstrapped_wad_diff.txt", "bootstrapped_ape.txt", "bootstrapped_r.txt", "bootstrapped_C_fluxes.txt"), sep="")
        unlink(bootstrapped.filenames)
      POFW$r.median[i] <- median(comparisons$r.obs, na.rm=TRUE)
      POFW$r.max[i] <- max(comparisons$r.obs, na.rm=TRUE)
    }

    #Plot results:
    dev.off()
    min(c(POFW$r.median, POFW$r.max), na.rm=TRUE)
    max(c(POFW$r.median, POFW$r.max), na.rm=TRUE)
  
    pdf(file="qSIP_output/Figures/BH_PropOxygenFromWaterSensitivity.pdf", width=6, height=6)
    par(mai=c(1.02, 1.02, 0.82, 0.42))
    plot(x=POFW$prop.O.from.water, y=POFW$r.max, bty="l", type="n", xlim=c(0,1), ylim=c(0,1), xlab="Proportion of oxygen in DNA from water (U)", ylab=expression(paste("r (day"^-1, ")", sep="")))
    points(x=POFW$prop.O.from.water, y=POFW$r.max, pch=21, col="black", bg="green")
    points(x=POFW$prop.O.from.water, y=POFW$r.median, pch=21, col="black", bg="red")
    legend(x=1, y=1, xjust=1, yjust=1, legend=c("median observed r across all taxa", "maximum observed r across all taxa"), pch=21, col="black", pt.bg=c("red","green"), pt.cex=1)
    dev.off()
  #}  




#_18O culture incubation data and analysis:
  #{
  #Import Theresa and Ember's 18O data & look at density-excess atom fraction relationship:
    #Read in 18O data:
    data.18O <- read.table("qSIP_data/BH_18O_alltaxa.txt", header=T, sep="\t", stringsAsFactors=T)
    data.18O

    #Natural abundance 18O EAF of water:
    min(data.18O$eaf.18O.H2O)
      #same as: 
      (2005.20/(1000000+379.9+2005.20))

    #Two ways to calculate the mean proportion of oxygen in DNA from environmental water using the labeling data (E. coli only); both yield nearly identical answers:
    mean(data.18O$prop.O.from.water, na.rm=TRUE)     #average of values already computed previously using the mixing model
    eaf.18O.H2O <- (data.18O$eaf.18O.DNA - min(data.18O$eaf.18O.H2O))/(data.18O$eaf.18O.H2O - min(data.18O$eaf.18O.H2O))     #independent verification of mixing model method
    mean(eaf.18O.H2O[is.finite(eaf.18O.H2O)])
    lm(data.18O$eaf.18O.DNA~data.18O$eaf.18O.H2O)$coef[2]     #as the slope of the regression line of measured EAF 18O of DNA vs. EAF 18O of H2O in the incubation

    #Create a data frame listing the mean weighted-average densities, GC content, MW, & C content for the reference treatments for each taxon:
    wad.ref.18O <- data.frame(taxon=levels(data.18O$taxon), wad=NA, GC=NA, GC.calc=NA, MW=NA, Catoms=NA, Oatoms=NA)
    for (k in 1:dim(wad.ref.18O)[1]){
      wad.ref.18O$wad[k] <- mean(data.18O$wad.g.ml[data.18O$taxon == wad.ref.18O$taxon[k] & data.18O$eaf.18O.H2O <= min(data.18O$eaf.18O.H2O)])
      wad.ref.18O$GC[k] <- unique(data.18O$GC[data.18O$taxon == wad.ref.18O$taxon[k]])
    }
      # Calculate GC content, molecular weight, and C content:
      # wad.ref.18O$GC.calc <- (1/0.098)*(wad.ref.18O$wad-1.66)  #From Schildkraut et al. 1962
      wad.ref.18O$GC.calc <- (1/0.0835059954345993)*(wad.ref.18O$wad-1.64605745338531)  #From McHugh & Morrissey unpublished data
      wad.ref.18O$MW <- (wad.ref.18O$GC.calc*0.496) + 307.691     #Molecular weight (g/mol) varies with GC content
      wad.ref.18O$Catoms <- (-0.5*wad.ref.18O$GC.calc) + 10       #Carbon content varies with GC content
      wad.ref.18O$Oatoms <- 6                                     #Oxygen content is constant for DNA regardless of GC content

    #Calculate excess atom fraction 18O:
    ape.of.fully.labeled.18O <- 1-(2005.20/(1000000+379.9+2005.20))   #atom percent excess of a substance with 100% 18O atoms (relative to VSMOW; see IAEA 1995, Werner & Brand 2001)
    for (j in 1:dim(data.18O)[1]){
      tax.data <- wad.ref.18O[wad.ref.18O$taxon == data.18O$taxon[j],]
      obs.labeled.MW <- tax.data$MW * (data.18O$wad.g.ml[j]/tax.data$wad)
      heavy.max.MW <- tax.data$MW + (tax.data$Oatoms*2)
      data.18O$eaf.18O.DNA.calc[j] <- (obs.labeled.MW-tax.data$MW)/(heavy.max.MW-tax.data$MW) * ape.of.fully.labeled.18O
      data.18O$wad.diff[j] <- data.18O$wad.g.ml[j] - tax.data$wad
    }

    #Plot excess atom fraction 18O vs difference in WAD:
    min(c(data.18O$wad.diff, all.comparisons$wad.diff.boot.median[all.comparisons$comparisonID %in% c(1,2,5,6)]), na.rm=TRUE)
    max(c(data.18O$wad.diff, all.comparisons$wad.diff.boot.median[all.comparisons$comparisonID %in% c(1,2,5,6)]), na.rm=TRUE)
    min(c(data.18O$eaf.18O.DNA, data.18O$eaf.18O.DNA.calc), na.rm=TRUE)
    max(c(data.18O$eaf.18O.DNA, data.18O$eaf.18O.DNA.calc), na.rm=TRUE)
    par(mai=c(1.02, 1.02, 0.82, 0.42))
    plot(x=data.18O$wad.diff, y=data.18O$eaf.18O.DNA.calc, bty="l", type="p", pch=21, col="blue", bg="blue", xlim=c(-0.015,0.025), ylim=c(-0.20,0.40), xlab=expression(paste("Weighted average density difference (g ml"^-1, ")", sep="")), ylab=expression(paste("Excess atom fraction "^18, "O", sep="")))
    points(x=data.18O$wad.diff, y=data.18O$eaf.18O.DNA, pch=21, col="red", bg="red")
    points(x=all.comparisons$wad.diff.boot.median[all.comparisons$comparisonID %in% c(1,2,5,6)], y=all.comparisons$ape.boot.median[all.comparisons$comparisonID %in% c(1,2,5,6)], pch=21, col="black", bg="green", cex=0.5)
    legend(x=-0.015, y=0.40, legend=c("measured EAF (E. coli)", "calculated EAF (E. coli)", "calculated EAF (all taxa)"), pch=21, col=c("red", "blue", "black"), pt.bg=c("red", "blue", "green"), pt.cex=c(1,1,0.5))

    #Plot excess atom fraction 18O vs proportion oxygen from water
    min(all.comparisons$wad.diff.boot.median[all.comparisons$comparisonID %in% c(1,2,5,6)], na.rm=TRUE)
    max(all.comparisons$wad.diff.boot.median[all.comparisons$comparisonID %in% c(1,2,5,6)], na.rm=TRUE)
    par(mai=c(1.02, 1.02, 0.82, 0.42))
    plot(x=data.18O$prop.O.from.water, y=data.18O$eaf.18O.DNA, bty="l", type="p", pch=21, col="red", bg="red", xlim=c(-0.1,1), ylim=c(-0.10,0.30), xlab="Proportion of oxygen in DNA from water", ylab=expression(paste("Excess atom fraction "^18, "O", sep="")))
    # points(x=data.18O$prop.O.from.water, y=data.18O$eaf.18O.DNA.calc, pch=21, col="blue", bg="blue")
    summary(lm(data.18O$eaf.18O.DNA~data.18O$prop.O.from.water))
    # summary(lm(data.18O$eaf.18O.DNA.calc~data.18O$prop.O.from.water))

    dev.off()
  #}




#_Compare GC content calculated using McHugh vs Schildkraut equations
  #{
  #Calculate the GC content for all taxa (using reference treatments from Week 1):
      #(This uses the default McHugh equation for GC content)
      taxa <- levels(factor(taxa.id.present$taxon))
      N.taxa <- length(levels(factor(taxa.id.present$taxon)))
      GC.taxa <- cbind(GC=numeric(dim(taxa.id.present)[1]), taxa.id.present)
      for (n in 1:N.taxa){
        Tref <- data.melted[data.melted$taxon==taxa[n] & (data.melted$trt.code=="1_C_16O" | data.melted$trt.code=="1_P_12Cplus16O"),]
        Tref.MW.out <- MW.calc(X=Tref, vars=c("density.g.ml", "copies", "tube"))
        GC.taxa$GC[n] <- Tref.MW.out$GC
      }

    #Write the results (GC.taxa) to a text file:
      write.table(GC.taxa, "qSIP_output/BH_GC_taxa.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)


  #Create a graph to compare GC content calculated with the density-GC equations from McHugh and Schildkraut:
    #First, save the McHugh GC content data frame:
      GC.taxa.McHugh <- GC.taxa

    #Next run the 'MW.calc.Schildkraut' script (with the Schildkraut equation in effect) to calculate the GC content for all taxa (using reference treatments from Week 1):
      GC.taxa.Schildkraut <- cbind(GC=numeric(dim(taxa.id.present)[1]), taxa.id.present)
      for (n in 1:N.taxa){
        Tref <- data.melted[data.melted$taxon==taxa[n] & (data.melted$trt.code=="1_C_16O" | data.melted$trt.code=="1_P_12Cplus16O"),]
        Tref.MW.out <- MW.calc.Schildkraut(X=Tref, vars=c("density.g.ml", "copies", "tube"))
        GC.taxa.Schildkraut$GC[n] <- Tref.MW.out$GC
      }

    #Plot the comparison of GC content calculated with the density-GC equations from McHugh and Schildkraut:
      # dev.new(width=6, height=6)
      pdf(file="qSIP_output/Figures/BH_GCcontent_McHugh_vs_Schildkraut.pdf", width=6, height=6)
      par(mfrow=c(2,1))
      hist(GC.taxa.McHugh$GC, xlim=c(0,1), ylim=c(0,120), xlab="GC content", main="McHugh")
      hist(GC.taxa.Schildkraut$GC, xlim=c(0,1), ylim=c(0,120), xlab="GC content", main="Schildkraut")
      par(mfrow=c(1,1))
      dev.off()
  #}
