# This code performs the basic qSIP and population calculations on the TM qSIP data sets


graphics.off()	#close all graphics windows


#Set working directory:
  #ALREADY DONE FOR THIS WORKSPACE IN PREVIOUSLY RUN CODE; only reset it here if loading previously saved workspace (see below)


#Reload the saved workspace resulting from the previous script:
  setwd("/Users/bk/Research/Projects/SIP_Modeling/qSIP")
  load("qSIP_workspaces/TM_01-02/.RData")


#Load libraries & scripts:
  #(none to load here)


#Import the data.frame containing the treatment comparisons to perform:
  Tcompare <- read.table("qSIP_data/TM_TreatmentComparisons.txt", header=TRUE, sep="\t", colClasses=c("numeric","factor","character","character","character","numeric","character"))
  summary(Tcompare)
  Tcompare


#Import the data.frame containing the Time 0 treatment comparisons to perform:
  TcompareTime0 <- read.table("qSIP_data/TM_TreatmentComparisons_Time0.txt", header=TRUE, sep="\t", colClasses=c("factor","character","character","numeric","character"))
  summary(TcompareTime0)
  TcompareTime0


##########{___________Calculate wad.diff, ape, r, & flux for all taxa and all comparisons____________________________

  #NOTES: using set.seed ensures that subsequent function calls utilize the same random number seed, so results are replicated exactly)
  #       using system.time calculates how long it took the function to run
  #       use exponential growth model
  #       set prop.O.from.water = U.consensus.obs.exp according to previous analysis constraining U for an exponential growth model
  #       rename the bootstrapped output files so that they are not overwritten if 'all.taxa.calcs' is re-run
    set.seed(100)
    system.time(all.comparisons <- all.taxa.calcs(X.all=data.melted, comparisons=Tcompare, M.soil=Sdat, vars=c("taxon", "density.g.ml", "copies", "tube", "tmt", "g.soil"), growth.model="exponential", prop.O.from.water=U.consensus.obs.exp, v.frac=50, copies.cell=6, pgC.cell=0.1, CI=0.90, draws=1000, tailed.test=1))
    bootstrapped.filenames <- paste("qSIP_output/", c("bootstrapped_wad1.txt", "bootstrapped_wad2.txt", "bootstrapped_wad_diff.txt", "bootstrapped_ape.txt", "bootstrapped_r.txt", "bootstrapped_C_fluxes.txt", "bootstrapped_N1.txt", "bootstrapped_N2.txt", "bootstrapped_N.txt"), sep="")
    new.bootstrapped.filenames <- paste("qSIP_output/", "TM_", c("bootstrapped_wad1.txt", "bootstrapped_wad2.txt", "bootstrapped_wad_diff.txt", "bootstrapped_ape.txt", "bootstrapped_r.txt", "bootstrapped_C_fluxes.txt", "bootstrapped_N1.txt", "bootstrapped_N2.txt", "bootstrapped_N.txt"), sep="")
    file.rename(from=bootstrapped.filenames, to=new.bootstrapped.filenames)
    summary(all.comparisons)
    dim(all.comparisons)


  #Write the results (all.comparisons) to a text file:
    dir.create(path=paste(getwd(), "/qSIP_output", sep=""), showWarnings=FALSE)
    write.table(all.comparisons, "qSIP_output/TM_all_comparisons.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)


  #Write the taxa.id dataframe to a text file:
    write.table(taxa.id, "qSIP_output/TM_taxa_ID.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
    
##########}__________________________________________________________________________________________________________


##########{___________Calculate population-based r & flux for all taxa_______________________________________________

  #NOTES: the treatments for time t (day 10) abundances are both the 16O and 18O treatments ('Time0' is the only treatment used for time 0)
  #       using set.seed ensures that subsequent function calls utilize the same random number seed, so results are replicated exactly)
  #       using system.time calculates how long it took the function to run
  #       use exponential growth model
  #       rename the bootstrapped output files so that they are not overwritten if 'all.taxa.calcs' is re-run
    set.seed(100)
    system.time(all.comparisons.pop <- all.taxa.calcs.pop(X.all=ncopies.tube.melted, comparisons=TcompareTime0, M.soil=Sdat, vars=c("taxon", "copies", "tube", "tmt", "g.soil"), growth.model="exponential", vol=c(100, 50), copies.cell=6, pgC.cell=0.1, CI=0.90, draws=1000))
    bootstrapped.filenames <- paste("qSIP_output/", c("bootstrapped_r_pop.txt", "bootstrapped_C_fluxes_pop.txt", "bootstrapped_N_T0_pop.txt", "bootstrapped_N_Tt_pop.txt"), sep="")
    new.bootstrapped.filenames <- paste("qSIP_output/", "TM_", c("bootstrapped_r_pop.txt", "bootstrapped_C_fluxes_pop.txt", "bootstrapped_N_T0_pop.txt", "bootstrapped_N_Tt_pop.txt"), sep="")
    file.rename(from=bootstrapped.filenames, to=new.bootstrapped.filenames)
    summary(all.comparisons.pop)
    dim(all.comparisons.pop)


  #Write the results (all.comparisons.pop) to a text file:
    dir.create(path=paste(getwd(), "/qSIP_output", sep=""), showWarnings=FALSE)
    write.table(all.comparisons.pop, "qSIP_output/TM_all_comparisons_pop.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)

##########}__________________________________________________________________________________________________________


##########{___________Calculate & verify birth, death, net growth, and flux rates & ape & N__________________________

  #Read in the necessary bootstrapped growth, C flux, and abundance estimates from the all.taxa.calcs & all.taxa.calcs.pop output:
    r.gross.boots <- read.table("qSIP_output/TM_bootstrapped_r.txt", header=TRUE, sep="")
    f.gross.boots <- read.table("qSIP_output/TM_bootstrapped_C_fluxes.txt", header=TRUE, sep="")
    r.net.boots <- read.table("qSIP_output/TM_bootstrapped_r_pop.txt", header=TRUE, sep="")
    f.net.boots <- read.table("qSIP_output/TM_bootstrapped_C_fluxes_pop.txt", header=TRUE, sep="")
    N.Tt.boots <- read.table("qSIP_output/TM_bootstrapped_N_Tt_pop.txt", header=TRUE, sep="")
    N.T0.boots <- read.table("qSIP_output/TM_bootstrapped_N_T0_pop.txt", header=TRUE, sep="")

  #Verify that they are all of the same dimension:
    dim(r.gross.boots)
    dim(f.gross.boots)
    dim(r.net.boots)
    dim(f.net.boots)
    dim(N.Tt.boots)
    dim(N.T0.boots)

  #Get indices for the comparisons of interest (there is only one comparison for this dataset)
    inds1 <- r.gross.boots$comparisonID == 1            #growth with added water

  #Ensure that taxa are in the same order in the different results data.frames and bootstrapped results data.frames:
    sum(r.gross.boots$taxonID[inds1] != all.comparisons$taxonID[inds1])
    sum(r.gross.boots$taxonID[inds1] != all.comparisons.pop$taxonID[inds1])
    sum(r.gross.boots$taxonID[inds1] != f.gross.boots$taxonID[inds1])
    sum(r.gross.boots$taxonID[inds1] != r.net.boots$taxonID[inds1])
    sum(r.gross.boots$taxonID[inds1] != f.net.boots$taxonID[inds1])
    sum(r.gross.boots$taxonID[inds1] != N.Tt.boots$taxonID[inds1])
    sum(r.gross.boots$taxonID[inds1] != N.T0.boots$taxonID[inds1])

  #Set desired CI:
    CI <- 0.90
  
  #Set growth model that was used to derive these results:
    growth.model <- "exponential"

  #Establish the duration of the incubation (in days):
    days <- unique(Tcompare$days)

  #Create a data frame to store results:
    #First write relevant columns from 'all.comparisons':
      all.rates <- data.frame(  taxonID=factor(as.character(all.comparisons$taxonID[inds1])),
                                comparisonID=factor(as.character(all.comparisons$comparisonID[inds1])),
                                trt.code.1=factor(as.character(all.comparisons$trt.code.1[inds1])),
                                trt.code.2=factor(as.character(all.comparisons$trt.code.2[inds1]))    
      )
    #Next add relevant columns from 'all.comparisons.pop':
      all.rates <- data.frame(  all.rates,
                                comparisonID.pop=factor(as.character(all.comparisons.pop$comparisonID[inds1])),
                                group.code.0=factor(as.character(all.comparisons.pop$group.code.0[inds1])),
                                group.code.t=factor(as.character(all.comparisons.pop$group.code.t[inds1]))
      )

  #Calculate and write the values for all rates & abundances to the 'all.rates' data.frame:
    #{
      #Excess atom fraction (ape):
        all.rates$ape.obs <- all.comparisons$ape.obs[inds1]
        all.rates$ape.boot.median <- all.comparisons$ape.boot.median[inds1]
        all.rates$ape.boot.CI.L <- all.comparisons$ape.boot.CI.L[inds1]
        all.rates$ape.boot.CI.U <- all.comparisons$ape.boot.CI.U[inds1]
      #Net growth rate (birth + death):
        all.rates$r.net.obs <- all.comparisons.pop$net.r.obs[inds1]
        all.rates$r.net.boot.median <- all.comparisons.pop$net.r.boot.median[inds1]
        all.rates$r.net.boot.CI.L <- all.comparisons.pop$net.r.boot.CI.L[inds1]
        all.rates$r.net.boot.CI.U <- all.comparisons.pop$net.r.boot.CI.U[inds1]
      #Gross growth rate (birth rate):
        all.rates$b.obs <- all.comparisons$r.obs[inds1]
        all.rates$b.boot.median <- all.comparisons$r.boot.median[inds1]
        all.rates$b.boot.CI.L <- all.comparisons$r.boot.CI.L[inds1]
        all.rates$b.boot.CI.U <- all.comparisons$r.boot.CI.U[inds1]
      #Death rate (turnover):
        #Calculate death rate as d = r - b (convention is that death rates are negative)...observed:
          all.rates$d.obs <- all.comparisons.pop$net.r.obs[inds1] - all.comparisons$r.obs[inds1]
            #Calculate death rate as d = r - b (convention is that death rates are negative)...bootstrapped:
              r.net.boots.only <- r.net.boots[inds1, 5:dim(r.net.boots)[2]]
              b.boots.only <- r.gross.boots[inds1, 5:dim(r.gross.boots)[2]]
              d.boots.only <- r.net.boots.only - b.boots.only
          all.rates$d.boot.median <- as.numeric(apply(d.boots.only, 1, median, na.rm=TRUE))
          all.rates$d.boot.CI.L <- as.numeric(apply(d.boots.only, 1, quantile, probs=(1-CI)/2, na.rm=TRUE))
          all.rates$d.boot.CI.U <- as.numeric(apply(d.boots.only, 1, quantile, probs=1-((1-CI)/2), na.rm=TRUE))
      #Net C flux into biomass:
        all.rates$f.net.obs <- all.comparisons.pop$net.f.obs[inds1]
        all.rates$f.net.boot.median <- all.comparisons.pop$net.f.boot.median[inds1]
        all.rates$f.net.boot.CI.L <- all.comparisons.pop$net.f.boot.CI.L[inds1]
        all.rates$f.net.boot.CI.U <- all.comparisons.pop$net.f.boot.CI.U[inds1]
      #Gross C flux into biomass:
        all.rates$f.gross.obs <- all.comparisons$f.obs[inds1]
        all.rates$f.gross.boot.median <- all.comparisons$f.boot.median[inds1]
        all.rates$f.gross.boot.CI.L <- all.comparisons$f.boot.CI.L[inds1]
        all.rates$f.gross.boot.CI.U <- all.comparisons$f.boot.CI.U[inds1]
      #C flux into dead cells (turnover):
        #Calculate C flux into dead cells as f.death = f.net - f.gross (convention is that the flux into dead cells is negative)...observed:
          all.rates$f.death.obs <- all.comparisons.pop$net.f.obs[inds1] - all.comparisons$f.obs[inds1]
            #Calculate C flux into dead cells as f.death = f.net - f.gross (convention is that the flux into dead cells is negative)...bootstrapped:
              f.net.boots.only <- f.net.boots[inds1, 5:dim(f.net.boots)[2]]
              f.gross.boots.only <- f.gross.boots[inds1, 5:dim(f.gross.boots)[2]]
              f.death.boots.only <- f.net.boots.only - f.gross.boots.only
          all.rates$f.death.boot.median <- as.numeric(apply(f.death.boots.only, 1, median, na.rm=TRUE))
          all.rates$f.death.boot.CI.L <- as.numeric(apply(f.death.boots.only, 1, quantile, probs=(1-CI)/2, na.rm=TRUE))
          all.rates$f.death.boot.CI.U <- as.numeric(apply(f.death.boots.only, 1, quantile, probs=1-((1-CI)/2), na.rm=TRUE))
      #Total abundance at time 0 (total copies/g soil):
        all.rates$N.tot.T0.obs <- all.comparisons.pop$N.T0.obs.mean[inds1]
        all.rates$N.tot.T0.boot.median <- all.comparisons.pop$N.T0.boot.median[inds1]
        all.rates$N.tot.T0.boot.CI.L <- all.comparisons.pop$N.T0.boot.CI.L[inds1]
        all.rates$N.tot.T0.boot.CI.U <- all.comparisons.pop$N.T0.boot.CI.U[inds1]
      #Total abundance at time t (total copies/g soil):
        all.rates$N.tot.Tt.obs <- all.comparisons.pop$N.Tt.obs.mean[inds1]
        all.rates$N.tot.Tt.boot.median <- all.comparisons.pop$N.Tt.boot.median[inds1]
        all.rates$N.tot.Tt.boot.CI.L <- all.comparisons.pop$N.Tt.boot.CI.L[inds1]
        all.rates$N.tot.Tt.boot.CI.U <- all.comparisons.pop$N.Tt.boot.CI.U[inds1]
      #Abundance of light copies at time t (light copies/g soil):
        #Calculate abundance of light copies at time t according to the growth model used:
          N.Tt.boots.only <- N.Tt.boots[inds1, 5:dim(N.Tt.boots)[2]]
          N.T0.boots.only <- N.T0.boots[inds1, 5:dim(N.T0.boots)[2]]
          if (growth.model == "exponential"){
            all.rates$N.light.Tt.obs <- all.comparisons.pop$N.Tt.obs.mean[inds1] / (exp(all.comparisons$r.obs[inds1]*days))
            N.light.Tt.boots.only <- N.Tt.boots.only / (exp(b.boots.only*days))
          }
          if (growth.model == "linear"){
            all.rates$N.light.Tt.obs <- all.comparisons.pop$N.Tt.obs.mean[inds1] - (all.comparisons$r.obs[inds1]*days)
            N.light.Tt.boots.only <- N.Tt.boots.only - (b.boots.only*days)
          }
        all.rates$N.light.Tt.boot.median <- as.numeric(apply(N.light.Tt.boots.only, 1, median, na.rm=TRUE))
        all.rates$N.light.Tt.boot.CI.L <- as.numeric(apply(N.light.Tt.boots.only, 1, quantile, probs=(1-CI)/2, na.rm=TRUE))
        all.rates$N.light.Tt.boot.CI.U <- as.numeric(apply(N.light.Tt.boots.only, 1, quantile, probs=1-((1-CI)/2), na.rm=TRUE))
    #}

  #Double check the calculations above:
    #{
      #Create a copy of 'all.rates' to facilitate the checks:
        all.rates.check <- all.rates
      #Net r addition check:
        r.net.add.check.boots.only <- b.boots.only + d.boots.only
        all.rates.check$r.net.add.check.boot.median <- as.numeric(apply(r.net.add.check.boots.only, 1, median, na.rm=TRUE))
        all.rates.check$r.net.add.check.boot.CI.L <- as.numeric(apply(r.net.add.check.boots.only, 1, quantile, probs=(1-CI)/2, na.rm=TRUE))
        all.rates.check$r.net.add.check.boot.CI.U <- as.numeric(apply(r.net.add.check.boots.only, 1, quantile, probs=1-((1-CI)/2), na.rm=TRUE))
      #Net f addition check:
        f.net.add.check.boots.only <- f.gross.boots.only + f.death.boots.only
        all.rates.check$f.net.add.check.boot.median <- as.numeric(apply(f.net.add.check.boots.only, 1, median, na.rm=TRUE))
        all.rates.check$f.net.add.check.boot.CI.L <- as.numeric(apply(f.net.add.check.boots.only, 1, quantile, probs=(1-CI)/2, na.rm=TRUE))
        all.rates.check$f.net.add.check.boot.CI.U <- as.numeric(apply(f.net.add.check.boots.only, 1, quantile, probs=1-((1-CI)/2), na.rm=TRUE))
      #Alternative calculation of death rate:
        #Create matrices for storing bootstrapped estimates of total copies / g soil for each taxon:
        taxa.tot.copies.g.soil.boots <- Q.LIGHT.0 <- matrix(NA, nrow=length(all.rates.check$taxonID), ncol=length(5:dim(r.gross.boots)[2]))
        #Fill the matrix with the bootstrapped estimates (random generation may differ from what was done in all.taxa.calcs.pop):
        set.seed(100)
        for (i in 1:length(all.rates.check$taxonID)){
          TiT12.tube <- ncopies.tube.melted[ncopies.tube.melted$taxon==all.rates.check$taxonID[i] & (ncopies.tube.melted$tmt=="16O" | ncopies.tube.melted$tmt=="18O"),]
          TiT12.tube.copies.out <- boot.TUBE.pop(X=TiT12.tube, M.soil=Sdat, vars=c("copies", "tube", "g.soil"), vol=50, CI=0.90, draws=1000)
          TiT12.tube.copies.out.tot.copies.g.soil <- TiT12.tube.copies.out$boot.tot.copies / TiT12.tube.copies.out$boot.g.soil
          taxa.tot.copies.g.soil.boots[i,] <- apply(TiT12.tube.copies.out.tot.copies.g.soil, 1, mean, na.rm=TRUE)

          TiT0.tube <- ncopies.tube.melted[ncopies.tube.melted$taxon==all.rates.check$taxonID[i] & ncopies.tube.melted$tmt=="Time0",]
          TiT0.tube.copies.out <- boot.TUBE.pop(X=TiT0.tube, M.soil=Sdat, vars=c("copies", "tube", "g.soil"), vol=100, CI=0.90, draws=1000)
          TiT0.tube.copies.out.tot.copies.g.soil <- TiT0.tube.copies.out$boot.tot.copies / TiT0.tube.copies.out$boot.g.soil
          Q.LIGHT.0[i,] <- apply(TiT0.tube.copies.out.tot.copies.g.soil, 1, mean, na.rm=TRUE)
        }
        Q.LIGHT.10 <- (1 / exp(r.gross.boots[inds1, 5:dim(r.gross.boots)[2]] * days)) * taxa.tot.copies.g.soil.boots
        d.alt.boots.only <- log(Q.LIGHT.10/Q.LIGHT.0) * (1/days)
        all.rates.check$d.alt.boot.median <- as.numeric(apply(d.alt.boots.only, 1, median, na.rm=TRUE))
        all.rates.check$d.alt.boot.CI.L <- as.numeric(apply(d.alt.boots.only, 1, quantile, probs=(1-CI)/2, na.rm=TRUE))
        all.rates.check$d.alt.boot.CI.U <- as.numeric(apply(d.alt.boots.only, 1, quantile, probs=1-((1-CI)/2), na.rm=TRUE))
    #}

  #Graphically verify checks on the calculations:
    #Compare observed values to bootstrapped medians (good agreement indicates that the value of U used resulted in minimal problems realted to potential cases of MW.lab > MW.max in the U sensitivity analysis)
      #{
        graphics.off()
        dev.new(width=7.5, height=5.0)
        par(mfrow=c(2,3))
        obs.r.values <- c(all.rates.check$r.net.obs, all.rates.check$r.net.boot.median)
        obs.r.values.real <- obs.r.values[obs.r.values != Inf & obs.r.values != -Inf & !is.na(obs.r.values)]
        plot(all.rates.check$r.net.boot.median~all.rates.check$r.net.obs, xlim=c(min(obs.r.values.real), max(obs.r.values.real)), ylim=c(min(obs.r.values.real), max(obs.r.values.real)))
        abline(a=0, b=1)

        obs.b.values <- c(all.rates.check$b.obs, all.rates.check$b.boot.median)
        obs.b.values.real <- obs.b.values[obs.b.values != Inf & obs.b.values != -Inf & !is.na(obs.b.values)]
        plot(all.rates.check$b.boot.median~all.rates.check$b.obs, xlim=c(min(obs.b.values.real), max(obs.b.values.real)), ylim=c(min(obs.b.values.real), max(obs.b.values.real)))
        abline(a=0, b=1)
 
        obs.d.values <- c(all.rates.check$d.obs, all.rates.check$d.boot.median)
        obs.d.values.real <- obs.d.values[obs.d.values != Inf & obs.d.values != -Inf & !is.na(obs.d.values)]
        plot(all.rates.check$d.boot.median~all.rates.check$d.obs, xlim=c(min(obs.d.values.real), max(obs.d.values.real)), ylim=c(min(obs.d.values.real), max(obs.d.values.real)))
        abline(a=0, b=1)
  
        obs.f.values <- c(all.rates.check$f.net.obs, all.rates.check$f.net.boot.median)
        obs.f.values.real <- obs.f.values[obs.f.values != Inf & obs.f.values != -Inf & !is.na(obs.f.values)]
        plot(all.rates.check$f.net.boot.median~all.rates.check$f.net.obs, xlim=c(min(obs.f.values.real), max(obs.f.values.real)), ylim=c(min(obs.f.values.real), max(obs.f.values.real)))
        abline(a=0, b=1)

        obs.f.gross.values <- c(all.rates.check$f.gross.obs, all.rates.check$f.gross.boot.median)
        obs.f.gross.values.real <- obs.f.gross.values[obs.f.gross.values != Inf & obs.f.gross.values != -Inf & !is.na(obs.f.gross.values)]
        plot(all.rates.check$f.gross.boot.median~all.rates.check$f.gross.obs, xlim=c(min(obs.f.gross.values.real), max(obs.f.gross.values.real)), ylim=c(min(obs.f.gross.values.real), max(obs.f.gross.values.real)))
        abline(a=0, b=1)

        obs.f.death.values <- c(all.rates.check$f.death.obs, all.rates.check$f.death.boot.median)
        obs.f.death.values.real <- obs.f.death.values[obs.f.death.values != Inf & obs.f.death.values != -Inf & !is.na(obs.f.death.values)]
        plot(all.rates.check$f.death.boot.median~all.rates.check$f.death.obs, xlim=c(min(obs.f.death.values.real), max(obs.f.death.values.real)), ylim=c(min(obs.f.death.values.real), max(obs.f.death.values.real)))
        abline(a=0, b=1)
        par(mfrow=c(1,1))
      #}

    #Compare calculations with alternative calculations performed above:
      #{
        dev.new(width=7.5, height=7.5)
        par(mfrow=c(3,3))
        r.values <- c(all.rates.check$r.net.boot.median, all.rates.check$r.net.add.check.boot.median)
        r.values.real <- r.values[r.values != Inf & r.values != -Inf & !is.na(r.values)]
        plot(all.rates.check$r.net.add.check.boot.median~all.rates.check$r.net.boot.median, xlim=c(min(r.values.real), max(r.values.real)), ylim=c(min(r.values.real), max(r.values.real)))
        abline(a=0, b=1)

        r.CI.Ls <- c(all.rates.check$r.net.boot.CI.L, all.rates.check$r.net.add.check.boot.CI.L)
        r.CI.Ls.real <- r.CI.Ls[r.CI.Ls != Inf & r.CI.Ls != -Inf & !is.na(r.CI.Ls)]
        plot(all.rates.check$r.net.add.check.boot.CI.L~all.rates.check$r.net.boot.CI.L, xlim=c(min(r.CI.Ls.real), max(r.CI.Ls.real)), ylim=c(min(r.CI.Ls.real), max(r.CI.Ls.real)))
        abline(a=0, b=1)

        r.CI.Us <- c(all.rates.check$r.net.boot.CI.U, all.rates.check$r.net.add.check.boot.CI.U)
        r.CI.Us.real <- r.CI.Us[r.CI.Us != Inf & r.CI.Us != -Inf & !is.na(r.CI.Us)]
        plot(all.rates.check$r.net.add.check.boot.CI.U~all.rates.check$r.net.boot.CI.U, xlim=c(min(r.CI.Us.real), max(r.CI.Us.real)), ylim=c(min(r.CI.Us.real), max(r.CI.Us.real)))
        abline(a=0, b=1)

        f.values <- c(all.rates.check$f.net.boot.median, all.rates.check$f.net.add.check.boot.median)
        f.values.real <- f.values[f.values != Inf & f.values != -Inf & !is.na(f.values)]
        plot(all.rates.check$f.net.add.check.boot.median~all.rates.check$f.net.boot.median, xlim=c(min(f.values.real), max(f.values.real)), ylim=c(min(f.values.real), max(f.values.real)))
        abline(a=0, b=1)

        f.CI.Ls <- c(all.rates.check$f.net.boot.CI.L, all.rates.check$f.net.add.check.boot.CI.L)
        f.CI.Ls.real <- f.CI.Ls[f.CI.Ls != Inf & f.CI.Ls != -Inf & !is.na(f.CI.Ls)]
        plot(all.rates.check$f.net.add.check.boot.CI.L~all.rates.check$f.net.boot.CI.L, xlim=c(min(f.CI.Ls.real), max(f.CI.Ls.real)), ylim=c(min(f.CI.Ls.real), max(f.CI.Ls.real)))
        abline(a=0, b=1)

        f.CI.Us <- c(all.rates.check$f.net.boot.CI.U, all.rates.check$f.net.add.check.boot.CI.U)
        f.CI.Us.real <- f.CI.Us[f.CI.Us != Inf & f.CI.Us != -Inf & !is.na(f.CI.Us)]
        plot(all.rates.check$f.net.add.check.boot.CI.U~all.rates.check$f.net.boot.CI.U, xlim=c(min(f.CI.Us.real), max(f.CI.Us.real)), ylim=c(min(f.CI.Us.real), max(f.CI.Us.real)))
        abline(a=0, b=1)

        d.values <- c(all.rates.check$d.alt.boot.median, all.rates.check$d.boot.median)
        d.values.real <- d.values[d.values != Inf & d.values != -Inf & !is.na(d.values)]
        plot(all.rates.check$d.alt.boot.median~all.rates.check$d.boot.median, xlim=c(min(d.values.real), max(d.values.real)), ylim=c(min(d.values.real), max(d.values.real)))
        abline(a=0, b=1)

        d.CI.Ls <- c(all.rates.check$d.alt.boot.CI.L, all.rates.check$d.boot.CI.L)
        d.CI.Ls.real <- d.CI.Ls[d.CI.Ls != Inf & d.CI.Ls != -Inf & !is.na(d.CI.Ls)]
        plot(all.rates.check$d.alt.boot.CI.L~all.rates.check$d.boot.CI.L, xlim=c(min(d.CI.Ls.real), max(d.CI.Ls.real)), ylim=c(min(d.CI.Ls.real), max(d.CI.Ls.real)))
        abline(a=0, b=1)

        d.CI.Us <- c(all.rates.check$d.alt.boot.CI.U, all.rates.check$d.boot.CI.U)
        d.CI.Us.real <- d.CI.Us[d.CI.Us != Inf & d.CI.Us != -Inf & !is.na(d.CI.Us)]
        plot(all.rates.check$d.alt.boot.CI.U~all.rates.check$d.boot.CI.U, xlim=c(min(d.CI.Us.real), max(d.CI.Us.real)), ylim=c(min(d.CI.Us.real), max(d.CI.Us.real)))
        abline(a=0, b=1)
        par(mfrow=c(1,1))
      #}

  graphics.off()

##########}__________________________________________________________________________________________________________


#Write the results (all.rates) to a text file:
  write.table(all.rates, "qSIP_output/TM_all_rates.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)




#Save workspace:
  save(list=ls(all.names=TRUE), file="qSIP_workspaces/TM_01-02-03/.RData", envir=.GlobalEnv)



