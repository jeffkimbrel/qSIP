### Given two treatments and the incubation period, returns an estimate (with uncertainty) of the flux of carbon into biomass (units: picograms of carbon per gram of soil per day) based on the degree of labeling by the heavy isotope
#
#     output = f.calc(X.light, X.heavy, boot.out.light, boot.out.heavy, MW.out, iso.compare, days, var="trt.code", growth.model="exponential", prop.O.from.water=0.50, copies.cell=6, pgC.cell=0.1, CI=0.90)
#
#     X.light: data frame with data for treatment 1, the unlabeled, or 'light', treatment
#     X.heavy: data frame with data for treatment 2, the labeled, or 'heavy', treatment
#     boot.out.light: list containing output from boot.TUBE.func for the 'light' treatment (corresponding to X.light)
#     boot.out.heavy: list containing output from boot.TUBE.func for the 'heavy' treatment (corresponding to X.heavy)
#     MW.out: data frame containing output from MW.calc for the unlabeled treatment(s) for a taxon
#     iso.compare: code specifying which isotopes are being compared; "13C" specifies the comparison between 13C (heavy) vs 12C (light) treatments; "18O" specifies the comparison between 18O (heavy) vs 16O (light) treatments; "15N" specifies the comparison between 15N (heavy) vs 14N (light) treatments; any other value returns an error
#     days: duration of the incubation period in days
#     var: variable name for the treatment ID columns in the data frames; default is "trt.code"
#     growth.model: model to use for calculating growth ("linear" or "exponential"); default is "exponential"
#     prop.O.from.water: proportion of oxygen atoms in DNA that come from environmental water; default is 0.50
#     copies.cell: assumed number of 16S copies per cell; default is 6
#     pgC.cell: assumed mass of carbon per cell (units: picograms); default is 0.1
#     CI: confidence interval (0-1); default is 0.90
#     -------------------------------------------------------
#     output: 
#     list of:        boot.f: vector of bootstrapped estimates of carbon flux (picograms of carbon per gram of soil per day) based on the two treatments
#                     obs.f: observed carbon flux based on the two treatments
#                     boot.f.mean: mean of the bootstrapped carbon flux values
#                     boot.f.median: median of the bootstrapped carbon flux values
#                     boot.f.CI: upper and lower confidence intervals (defined by argument 'CI') of the bootstrapped carbon flux values
#                     message: warning message listing replicates in either of the 2 treatments (does not include the reference group) that had no occurrences (i.e., all y values for a rep are 0); value for message is "none" when all replicates are present (single character value)
#
#     notes:          requires the functions 'WAD.func', 'boot.TUBE.func', 'MW.calc', & 'comparison.message'
#                     number of reps need not be equal between the two groups
#     assumptions:    assumes either a linear growth model ( Nt=N0+(rt) ) or an exponential growth model ( Nt=N0*e^(rt) )
#                     note that the units of r differ between the linear model (copies * day-1 * g soil-1) and the exponential model (1/day), but the flux calculations are adjusted according to the model specified so that the final units of carbon flux are the same (picograms of carbon per gram of soil per day)
#                     assumes a fixed GC content, molecular weight, carbon, and nitrogen content for each taxon based on the mean observed WAD
#                     oxygen content is fixed for DNA regardless of GC content
#                     the weighted average density of unlabeled DNA is a good proxy for GC content of the taxon, according to the well-known relationship between GC content and the density of DNA (Schildkraut et al. 1962)
#                       (NOTE: This is not caused by the effect of GC content on molecular weight, which is negligible. It is caused by differences in the relationship between base composition and binding with water that occurs in CsCl.)
#                     any change in weighted average density of the 16S DNA molecule is caused by incorporation of stable isotopes into the DNA
#                     the relative change in weighted average density is equal to the relative change in molecular weight
#                       (Incorporation of stable isotopes will not alter the binding properties between DNA and water in the CsCl gradient. Those properties are really important for affecting density, but they're a function of GC content, and should not change just because of heavy isotope incorporation.)
#
#     Written by Ben Koch & Natasja van Gestel


f.calc <- function(X.light, X.heavy, boot.out.light, boot.out.heavy, MW.out, iso.compare, days, var="trt.code", growth.model="exponential", prop.O.from.water=0.50, copies.cell=6, pgC.cell=0.1, CI=0.90){
  obs.labeled.MW <- MW.out$MW * (boot.out.heavy$obs.wad.mean/boot.out.light$obs.wad.mean)
  boot.labeled.MW <- MW.out$MW * (boot.out.heavy$boot.wads/boot.out.light$boot.wads)

  if (is.na(iso.compare)){                                            #if the iso.compare code is 'NA', print an error message
    stop(paste("Error: unable to calculate ape for the specified isotope comparison '", iso.compare, "'", sep=""))
  }
  else if (iso.compare == "13C"){                                     #if the treatments being compared are 13C (heavy) vs 12C (light) do the 13-carbon calculation...
    heavy.MW <- MW.out$MW + (MW.out$Catoms * (1.008665 * (1000000/(1000000+11237.2))))    #assumes unlabeled DNA already contains a minute amount of 13C (at the natural abundance level of VPDB)
  }
  else if (iso.compare == "18O"){                                     #if the treatments being compared are 18O (heavy) vs 16O (light) do the 18-oxygen calculation...
    heavy.MW <- MW.out$MW + (MW.out$Oatoms * prop.O.from.water * ((1.008665 * 2 * (1000000/(1000000+379.9+2005.20))) + (1.008665 * 1 * (379.9/(1000000+379.9+2005.20)))))    #assumes unlabeled DNA already contains minute amounts of 18O and 17O (at the natural abundance levels of those isotopes in VSMOW)
  }
  else if (iso.compare == "15N"){                                     #if the treatments being compared are 15N (heavy) vs 14N (light) do the 15-nitrogen calculation...
    heavy.MW <- MW.out$MW + (MW.out$Natoms * (1.008665 * (1000000/(1000000+(1000000/272)))))  #assumes unlabeled DNA already contains minute amounts of 15N (at the natural abundance level of AIR-N2)
  }
  else {                                                              #if the iso.compare code is anything else, print an error message
    stop(paste("Error: unable to calculate flux for the specified isotope comparison '", iso.compare, "'", sep=""))
  }
  if (growth.model=="linear"){                                         #if growth.model is "linear", calculate r using a linear growth model
    obs.copies.soil <- mean(c(boot.out.light$obs.wads$tot.copies/boot.out.light$obs.wads$g.soil, boot.out.heavy$obs.wads$tot.copies/boot.out.heavy$obs.wads$g.soil), na.rm=T)
    boot.copies.soil <- apply(cbind(boot.out.light$boot.tot.copies/boot.out.light$boot.g.soil, boot.out.heavy$boot.tot.copies/boot.out.heavy$boot.g.soil), 1, mean, na.rm=T)
    
    obs.r <- obs.copies.soil * (((obs.labeled.MW-MW.out$MW)/(heavy.MW-MW.out$MW))/days)
    boot.r <- boot.copies.soil * (((boot.labeled.MW-MW.out$MW)/(heavy.MW-MW.out$MW))/days)

    obs.f <- obs.r * (1/copies.cell) * pgC.cell
    boot.f <- boot.r * (1/copies.cell) * pgC.cell
  }
  else if (growth.model=="exponential"){                              #if growth.model is "exponential", calculate r using an exponential growth model
    obs.r <- log((heavy.MW-MW.out$MW)/(heavy.MW-obs.labeled.MW))/days
    boot.r <- log((heavy.MW-MW.out$MW)/(heavy.MW-boot.labeled.MW))/days
    
    obs.copies.soil <- mean(c(boot.out.light$obs.wads$tot.copies/boot.out.light$obs.wads$g.soil, boot.out.heavy$obs.wads$tot.copies/boot.out.heavy$obs.wads$g.soil), na.rm=T)
    boot.copies.soil <- apply(cbind(boot.out.light$boot.tot.copies/boot.out.light$boot.g.soil, boot.out.heavy$boot.tot.copies/boot.out.heavy$boot.g.soil), 1, mean, na.rm=T)

    obs.f <- obs.r * obs.copies.soil * (1/copies.cell) * pgC.cell
    boot.f <- boot.r * boot.copies.soil * (1/copies.cell) * pgC.cell
  }
  else {                                                              #if the growth.model is anything else, print an error message
    stop(paste("Error: unable to calculate r for the specified growth model '", growth.model, "'", sep=""))
  }

  boot.f.CI <- quantile(boot.f, probs=c((1-CI)/2, 1-((1-CI)/2)), na.rm=T)

  message <- comparison.message(X.light=X.light, X.heavy=X.heavy, boot.out.light=boot.out.light, boot.out.heavy=boot.out.heavy, var=var)

  # Collect all output into a list:
   return(list(boot.f=boot.f, 
               obs.f=obs.f, 
               boot.f.mean=mean(boot.f, na.rm=T), 
               boot.f.median=median(boot.f, na.rm=T), 
               boot.f.CI=boot.f.CI, 
               message=message)) 
}
