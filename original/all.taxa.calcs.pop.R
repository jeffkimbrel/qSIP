
all.taxa.calcs.pop <- function(X.all, comparisons, M.soil, vars=c("taxon", "copies", "tube", "trt.code", "g.soil"), growth.model="exponential", vol=c(100, 50), copies.cell=6, pgC.cell=0.1, CI=0.90, draws=1000){
  #Create an empty data frame for the output:
  info <- data.frame(matrix(nrow=0, ncol=25))

  #Create empty data frames for the bootstrapped abundances (T0 & Tt), bootstrapped net C flux, and bootstrapped net growth rate outputs:
  N.Tt.boots <- N.T0.boots <- net.C.flux.boots <- net.r.boots <- data.frame(matrix(nrow=0, ncol=draws+4))

  #Establish the the number of comparisons and the number of taxa:
  N.comparisons <- length(levels(factor(comparisons$comparisonID)))
  N.taxa <- length(levels(X.all[,vars[1]]))

  #Establish the treatment groups and number of unique treatment groups for the T0 and Tt designations across all comparisons:
  group.t.list <- group.0.list <- list(NULL)
  group.desig <- data.frame(group.code.0=character(N.comparisons), group.code.t=character(N.comparisons), stringsAsFactors=FALSE)
  for (c in 1:N.comparisons){  #for each comparison...
    m0 <- gregexpr(pattern="([[:space:]]*(,|;)[[:space:]]*)|([[:space:]]+)", comparisons$trt.code.0[c], perl=TRUE)
    mt <- gregexpr(pattern="([[:space:]]*(,|;)[[:space:]]*)|([[:space:]]+)", comparisons$trt.code.t[c], perl=TRUE)
    group.0.list[[c]] <- sort(unlist(regmatches(comparisons$trt.code.0[c], m0, invert=TRUE)))
    group.t.list[[c]] <- sort(unlist(regmatches(comparisons$trt.code.t[c], mt, invert=TRUE)))
    group.desig$group.code.0[c] <- names(group.0.list)[[c]] <- paste(group.0.list[[c]], collapse=".")
    group.desig$group.code.t[c] <- names(group.t.list)[[c]] <- paste(group.t.list[[c]], collapse=".")
  }
  groups.list <- c(group.0.list, group.t.list)
  groups <- levels(factor(as.character(unlist(c(group.desig)))))
  N.groups <- length(groups)

  #Establish the treatments and the total number of treatments:
  treatments <- levels(factor(unlist(c(group.0.list, group.t.list))))
  N.treatments <- length(treatments)

  #Create a dataframe for building a key that relates treatment codes to the names of treatments used in the output lists:
  GROUPID <- data.frame(group.mat.name=character(N.groups), group.code=character(N.groups), trt.codes=character(N.groups), stringsAsFactors=FALSE)
  GROUPID_comparisons <- data.frame(matrix(NA, nrow=N.groups, ncol=N.comparisons))
  names(GROUPID_comparisons) <- paste("comparison", seq(1, N.comparisons, 1), sep="")
  GROUPID <- data.frame(GROUPID, GROUPID_comparisons)

  for (p in 1:N.groups){
    #Calculate the number of reps over which to resample - for each group:
    test.data <- X.all[X.all[,vars[4]] %in% groups.list[[groups[p]]],]
    N.reps <- length(levels(factor(test.data[,vars[3]])))

    #Calculate bootstrap vectors of the resampled indices by which to calculate mean abundances and mass of soil across reps:
    group.indices <- matrix(nrow=draws, ncol=N.reps)
    for (h in 1:draws){
      group.indices[h,] <- sample.vec(1:N.reps, N.reps, replace=TRUE)
    }
  
    #Name the matrix for the current treatment:
    assign(paste("GROUP_", as.character(groups[p]), sep=""), group.indices)
    #Write the current group code (and treatment codes) and the name of the group used in the output matrix to the 'key' dataframe:
    GROUPID$group.mat.name[p] <- paste("GROUP_", as.character(groups[p]), sep="")
    GROUPID$group.code[p] <- as.character(groups[p])
    GROUPID$trt.codes[p] <- paste(groups.list[[groups[p]]], collapse="; ")
    
    #Define the T0 and Tt treatment designations associated with each comparison:
    for (c in 1:N.comparisons){  #for each comparison...
      if (names(group.0.list)[c] %in% groups[p]){
        GROUPID[p, 3+c] <- "T0"
      }
      if (names(group.t.list)[c] %in% groups[p]){  
        GROUPID[p, 3+c] <- "Tt"
      }
    }
  }
  #Convert columns in the GROUPID 'key' dataframe to factor:
  GROUPID <- data.frame(unclass(GROUPID))

  for (c in 1:N.comparisons){  #for each comparison...
    #Get the duration of the incubation period for the current comparison:
    days <- comparisons$days[c]

    #Get the names of the matrices of resampled indices for the 'T0' and 'Tt' groups:
    mat.name.T0 <- as.character(GROUPID$group.mat.name[GROUPID[,3+c] == "T0" & !is.na(GROUPID[,3+c])])
    mat.name.Tt <-  as.character(GROUPID$group.mat.name[GROUPID[,3+c] == "Tt" & !is.na(GROUPID[,3+c])])

    #Reference the bootstrap samples created above to do the qPCR calculations:
    for (i in 1:N.taxa){  #for each taxon...

      #Subset the data by taxon:
      T0 <- X.all[X.all[,vars[1]]==levels(X.all[,vars[1]])[i] & X.all[,vars[4]] %in% group.0.list[[c]],]
      Tt <- X.all[X.all[,vars[1]]==levels(X.all[,vars[1]])[i] & X.all[,vars[4]] %in% group.t.list[[c]],]

      #First, for the T0 group:
        #Calculate observed abundance, total 16S copies, and mass of soil for each rep:
        obs.N <- data.frame(matrix(nrow=length(levels(factor(T0[,vars[3]]))), ncol=4))
        names(obs.N) <- c("copies", "tot.copies", "g.soil", "rep")
        for (g in 1:length(levels(factor(T0[,vars[3]])))){
          obs.N$rep[g] <- levels(factor(T0[,vars[3]]))[g]
          obs.N$copies[g] <- T0[T0[,vars[3]]==levels(factor(T0[,vars[3]]))[g], vars[2]]
          obs.N$tot.copies[g] <- T0[T0[,vars[3]]==levels(factor(T0[,vars[3]]))[g], vars[2]]*vol[1]
          obs.N$g.soil[g] <- M.soil[M.soil[,vars[3]] == levels(factor(T0[,vars[3]]))[g], vars[5]]
        }
        obs.N$rep <- factor(obs.N$rep)
  
        #Calculate a bootstrap vector of mean abundance across reps along with total 16S copies & mass of soil:
        resampled.tot.copies <- matrix(obs.N$tot.copies[eval(parse(text=mat.name.T0))], nrow=draws, ncol=dim(obs.N)[1])
        resampled.g.soil <- matrix(obs.N$g.soil[eval(parse(text=mat.name.T0))], nrow=draws, ncol=dim(obs.N)[1])
        resampled.copies <- matrix(obs.N$copies[eval(parse(text=mat.name.T0))], nrow=draws, ncol=dim(obs.N)[1])
        boot.N <- apply(resampled.copies, 1, mean, na.rm=T)
  
        boot.out.T0 <- list(
                            boot.tot.copies = resampled.tot.copies,
                            boot.g.soil = resampled.g.soil,
                            boot.N = boot.N,
                            obs.N = obs.N,
                            obs.N.mean = mean(obs.N$copies, na.rm=T), 
                            boot.N.mean=mean(boot.N, na.rm=T), 
                            boot.N.median=median(boot.N, na.rm=T), 
                            boot.N.CI=quantile(boot.N, probs=c((1-CI)/2, 1-((1-CI)/2)), na.rm=T)
        )

      #Next, for the Tt group:
        #Calculate observed abundance, total 16S copies, and mass of soil for each rep:
        obs.N <- data.frame(matrix(nrow=length(levels(factor(Tt[,vars[3]]))), ncol=4))
        names(obs.N) <- c("copies", "tot.copies", "g.soil", "rep")
        for (g in 1:length(levels(factor(Tt[,vars[3]])))){
          obs.N$rep[g] <- levels(factor(Tt[,vars[3]]))[g]
          obs.N$copies[g] <- Tt[Tt[,vars[3]]==levels(factor(Tt[,vars[3]]))[g], vars[2]]
          obs.N$tot.copies[g] <- Tt[Tt[,vars[3]]==levels(factor(Tt[,vars[3]]))[g], vars[2]]*vol[2]
          obs.N$g.soil[g] <- M.soil[M.soil[,vars[3]] == levels(factor(Tt[,vars[3]]))[g], vars[5]]
        }
        obs.N$rep <- factor(obs.N$rep)

        #Calculate a bootstrap vector of mean abundance across reps along with total 16S copies & mass of soil:
        resampled.tot.copies <- matrix(obs.N$tot.copies[eval(parse(text=mat.name.Tt))], nrow=draws, ncol=dim(obs.N)[1])
        resampled.g.soil <- matrix(obs.N$g.soil[eval(parse(text=mat.name.Tt))], nrow=draws, ncol=dim(obs.N)[1])
        resampled.copies <- matrix(obs.N$copies[eval(parse(text=mat.name.Tt))], nrow=draws, ncol=dim(obs.N)[1])
        boot.N <- apply(resampled.copies, 1, mean, na.rm=T)

        boot.out.Tt <- list(
                            boot.tot.copies = resampled.tot.copies,
                            boot.g.soil = resampled.g.soil,
                            boot.N = boot.N,
                            obs.N = obs.N,
                            obs.N.mean = mean(obs.N$copies, na.rm=T), 
                            boot.N.mean=mean(boot.N, na.rm=T), 
                            boot.N.median=median(boot.N, na.rm=T), 
                            boot.N.CI=quantile(boot.N, probs=c((1-CI)/2, 1-((1-CI)/2)), na.rm=T)
        )

      #Calculate the appropriate output values:
      r.out <- r.calc.pop(T0=T0, Tt=Tt, boot.out.T0=boot.out.T0, boot.out.Tt=boot.out.Tt, days=days, var=vars[4], growth.model=growth.model, CI=CI)
      f.out <- f.calc.pop(T0=T0, Tt=Tt, boot.out.T0=boot.out.T0, boot.out.Tt=boot.out.Tt, days=days, var=vars[4], growth.model=growth.model, copies.cell=copies.cell, pgC.cell=pgC.cell, CI=CI)

      #Write the vector of bootstrapped abundances (T0 & Tt), bootstrapped net C fluxes, and bootstrapped net growth rates to output dataframes to enable downstream calculations:
      N.Tt.to.add <- N.T0.to.add <- net.f.to.add <- net.r.to.add <- data.frame(matrix(nrow=1, ncol=draws+4))
      N.Tt.to.add[1] <- N.T0.to.add[1] <- net.f.to.add[1] <- net.r.to.add[1] <- levels(X.all[,vars[1]])[i]
      N.Tt.to.add[2] <- N.T0.to.add[2] <- net.f.to.add[2] <- net.r.to.add[2] <- as.character(comparisons$comparisonID[c])
      N.Tt.to.add[3] <- N.T0.to.add[3] <- net.f.to.add[3] <- net.r.to.add[3] <- as.character(GROUPID$group.code[GROUPID[,3+c] == "T0" & !is.na(GROUPID[,3+c])])
      N.Tt.to.add[4] <- N.T0.to.add[4] <- net.f.to.add[4] <- net.r.to.add[4] <- as.character(GROUPID$group.code[GROUPID[,3+c] == "Tt" & !is.na(GROUPID[,3+c])])
      net.r.to.add[5:dim(net.r.to.add)[2]] <- r.out$boot.r
      net.f.to.add[5:dim(net.f.to.add)[2]] <- f.out$boot.f
      N.T0.to.add[5:dim(N.T0.to.add)[2]] <- apply(boot.out.T0$boot.tot.copies/boot.out.T0$boot.g.soil, 1, mean, na.rm=TRUE)
      N.Tt.to.add[5:dim(N.Tt.to.add)[2]] <- apply(boot.out.Tt$boot.tot.copies/boot.out.Tt$boot.g.soil, 1, mean, na.rm=TRUE)
      names(net.r.to.add) <- c("taxonID", "comparisonID", "group.code.0", "group.code.t", paste("net.r", 1:draws, sep=""))
      names(net.f.to.add) <- c("taxonID", "comparisonID", "group.code.0", "group.code.t", paste("net.f", 1:draws, sep=""))
      names(N.T0.to.add) <- c("taxonID", "comparisonID", "group.code.0", "group.code.t", paste("N.T0.", 1:draws, sep=""))
      names(N.Tt.to.add) <- c("taxonID", "comparisonID", "group.code.0", "group.code.t", paste("N.Tt.", 1:draws, sep=""))
    
      net.r.boots <- rbind(net.r.boots, net.r.to.add)
      net.C.flux.boots <- rbind(net.C.flux.boots, net.f.to.add)
      N.T0.boots <- rbind(N.T0.boots, N.T0.to.add)
      N.Tt.boots <- rbind(N.Tt.boots, N.Tt.to.add)

      #Write the appropriate output values for the current taxon-comparison to a data frame:
      info.to.add <- data.frame(
                        taxonID=levels(X.all[,vars[1]])[i],
                        comparisonID=as.character(comparisons$comparisonID[c]),
                        group.code.0=as.character(GROUPID$group.code[GROUPID[,3+c] == "T0" & !is.na(GROUPID[,3+c])]),
                        group.code.t=as.character(GROUPID$group.code[GROUPID[,3+c] == "Tt" & !is.na(GROUPID[,3+c])]),

                        net.r.obs=r.out$obs.r,
                        net.r.boot.mean=r.out$boot.r.mean,
                        net.r.boot.median=r.out$boot.r.median,
                        net.r.boot.CI.L=as.numeric(r.out$boot.r.CI[1]),
                        net.r.boot.CI.U=as.numeric(r.out$boot.r.CI[2]),

                        net.f.obs=f.out$obs.f,
                        net.f.boot.mean=f.out$boot.f.mean,
                        net.f.boot.median=f.out$boot.f.median,
                        net.f.boot.CI.L=as.numeric(f.out$boot.f.CI[1]),
                        net.f.boot.CI.U=as.numeric(f.out$boot.f.CI[2]),

                        N.T0.obs.mean=mean(boot.out.T0$obs.N$tot.copies/boot.out.T0$obs.N$g.soil, na.rm=TRUE),
                        N.T0.boot.mean=mean(apply(boot.out.T0$boot.tot.copies/boot.out.T0$boot.g.soil, 1, mean, na.rm=TRUE)),
                        N.T0.boot.median=median(apply(boot.out.T0$boot.tot.copies/boot.out.T0$boot.g.soil, 1, mean, na.rm=TRUE)),
                        N.T0.boot.CI.L=quantile(apply(boot.out.T0$boot.tot.copies/boot.out.T0$boot.g.soil, 1, mean, na.rm=TRUE), probs=(1-CI)/2),
                        N.T0.boot.CI.U=quantile(apply(boot.out.T0$boot.tot.copies/boot.out.T0$boot.g.soil, 1, mean, na.rm=TRUE), probs=1-((1-CI)/2)),

                        N.Tt.obs.mean=mean(boot.out.Tt$obs.N$tot.copies/boot.out.Tt$obs.N$g.soil, na.rm=TRUE),
                        N.Tt.boot.mean=mean(apply(boot.out.Tt$boot.tot.copies/boot.out.Tt$boot.g.soil, 1, mean, na.rm=TRUE)),
                        N.Tt.boot.median=median(apply(boot.out.Tt$boot.tot.copies/boot.out.Tt$boot.g.soil, 1, mean, na.rm=TRUE)),
                        N.Tt.boot.CI.L=quantile(apply(boot.out.Tt$boot.tot.copies/boot.out.Tt$boot.g.soil, 1, mean, na.rm=TRUE), probs=(1-CI)/2),
                        N.Tt.boot.CI.U=quantile(apply(boot.out.Tt$boot.tot.copies/boot.out.Tt$boot.g.soil, 1, mean, na.rm=TRUE), probs=1-((1-CI)/2)),

                        message=r.out$message
      )
      info <- rbind(info, info.to.add)
    }
  }
  #Convert appropriate variables in 'N.Tt.boots', 'N.T0.boots', 'net.C.flux.boots', and 'net.r.boots' to factors:
  net.r.boots$taxonID <- factor(net.r.boots$taxonID)
  net.r.boots$comparisonID <- factor(net.r.boots$comparisonID)
  net.r.boots$group.code.0 <- factor(net.r.boots$group.code.0)
  net.r.boots$group.code.t <- factor(net.r.boots$group.code.t)
  net.C.flux.boots$taxonID <- factor(net.C.flux.boots$taxonID)
  net.C.flux.boots$comparisonID <- factor(net.C.flux.boots$comparisonID)
  net.C.flux.boots$group.code.0 <- factor(net.C.flux.boots$group.code.0)
  net.C.flux.boots$group.code.t <- factor(net.C.flux.boots$group.code.t)
  N.T0.boots$taxonID <- factor(N.T0.boots$taxonID)
  N.T0.boots$comparisonID <- factor(N.T0.boots$comparisonID)
  N.T0.boots$group.code.0 <- factor(N.T0.boots$group.code.0)
  N.T0.boots$group.code.t <- factor(N.T0.boots$group.code.t)
  N.Tt.boots$taxonID <- factor(N.Tt.boots$taxonID)
  N.Tt.boots$comparisonID <- factor(N.Tt.boots$comparisonID)
  N.Tt.boots$group.code.0 <- factor(N.Tt.boots$group.code.0)
  N.Tt.boots$group.code.t <- factor(N.Tt.boots$group.code.t)
  #Export the abundance (T0 & Tt) bootstrapped estimates, the net C flux bootstrapped estimates, and the net growth rate bootstrapped estimates to a text file:
  dir.create(path=paste(getwd(), "/qSIP_output", sep=""), showWarnings=FALSE)
  write.table(net.r.boots, "qSIP_output/bootstrapped_r_pop.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
  write.table(net.C.flux.boots, "qSIP_output/bootstrapped_C_fluxes_pop.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
  write.table(N.T0.boots, "qSIP_output/bootstrapped_N_T0_pop.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)
  write.table(N.Tt.boots, "qSIP_output/bootstrapped_N_Tt_pop.txt", append=F, quote=F, sep="\t", eol="\n", na="NA", dec=".", row.names=F, col.names=T)

  #Convert appropriate variables in 'info' to factors:
  info$taxonID <- factor(info$taxonID)
  info$comparisonID <- factor(info$comparisonID)
  info$group.code.0 <- factor(info$group.code.0)
  info$group.code.t <- factor(info$group.code.t)
  return(info) 
}
