#Define a function to use to filter the data so that only taxa with an appropriate level of occurrence among all fractions within replicate(s) of a treatment are retained for further analysis:

  filter.fractions.taxa <- function(DATA, trt.code=NULL, vars=c("taxon", "density.g.ml", "copies", "tube", "trt.code"), min.fracs, min.reps){

    #First, determine the specified treatment:
      trt.to.filter <- trt.code
  
    #Subset the data into only those taxon-reps with copies present in the specified treatments:
      DATA.occurrences <- DATA[!is.na(DATA[,vars[3]]) & DATA[,vars[3]] > 0 & DATA[,vars[5]] %in% trt.to.filter,]
        #Renumber row.names:
        row.names(DATA.occurrences) <- 1:dim(DATA.occurrences)[1]
        #Convert the factor columns to factor:
        DATA.occurrences <- as.data.frame(lapply(DATA.occurrences, function(x) if(is.factor(x)) factor(x) else x))
  
    #Calculate the number of unique tubes in the specified treatments with copies present for each taxon:
      tubes.per.taxon <- tapply(DATA.occurrences[,vars[4]], DATA.occurrences[,vars[1]], function(x) length(unique(x)))
      tubes.per.taxon <- sort(tubes.per.taxon)

    #Subset the data to include all rows for the specified treatment:
      DATA.treatment <- DATA[DATA[,vars[5]] %in% trt.to.filter,]
        #Renumber row.names:
        row.names(DATA.treatment) <- 1:dim(DATA.treatment)[1]
        #Convert the factor columns to factor:
        DATA.treatment <- as.data.frame(lapply(DATA.treatment, function(x) if(is.factor(x)) factor(x) else x))

    #For each replicate, calculate the number of total fractions and the number of fractions that taxon 'X' occured in:
      reps.in.trt <- unique(DATA.treatment[,vars[4]])
      fracs.in.rep <- data.frame(rep=reps.in.trt, num.fracs=numeric(length(reps.in.trt)))
      LIST <- LIST.FRAC <- TAXA.TO.KEEP <- NULL
      for (i in 1:length(reps.in.trt)){
        DATA.treatment.rep_fracs.all <- DATA.treatment[DATA.treatment[,vars[4]] == reps.in.trt[i], vars[2]]
        fracs.in.rep$num.fracs[i] <- length(unique(DATA.treatment.rep_fracs.all))
        DATA.treatment.rep_taxa <- DATA.occurrences[DATA.occurrences[,vars[4]] == reps.in.trt[i], vars[1]]
        DATA.treatment.rep_fracs <- DATA.occurrences[DATA.occurrences[,vars[4]] == reps.in.trt[i], vars[2]]
        #Calculate the number of unique fractions in the specified replicate with copies present for each taxon:
          fracs.per.taxon <- tapply(DATA.treatment.rep_fracs, DATA.treatment.rep_taxa, function(x) length(unique(x)))
          fracs.per.taxon[is.na(fracs.per.taxon)] <- 0
          LIST[[i]] <- fracs.per.taxon
          names(LIST)[i] <- paste("rep.", reps.in.trt[i], sep="")
      }

    #Calculate the fraction of total fractions that each taxon was present in ... for all replicates:
      for (i in 1:dim(fracs.in.rep)[1]){
        LIST.FRAC[[i]] <- LIST[[i]] / fracs.in.rep$num.fracs[i]
        names(LIST.FRAC)[i] <- names(LIST)[i]
      }

    #Make sure 'min.fracs' is specified correctly & use the appropriate list & threshold to screen taxa:
      if (is.factor(min.fracs)){
        stop("'min.fracs' must either be an integer specifying the minimum number of fractions that a taxon must be present within a replicate or a character string specifying the minimum percentage of total fractions (e.g., '40%')")
      }  else  if (is.numeric(min.fracs)){
        is.wholenumber <- function(x, tol=.Machine$double.eps^0.5){
          abs(x - round(x)) < tol
        }
        if (is.wholenumber(min.fracs)){
          FINAL.LIST <- LIST
          threshold <- min.fracs
        }  else  stop("'min.fracs' must either be an integer specifying the minimum number of fractions that a taxon must be present within a replicate or a character string specifying the minimum percentage of total fractions (e.g., '40%')")
      }  else  if (is.character(min.fracs)){
        if (grepl(pattern="\\%", x=min.fracs, perl=TRUE)){
          FINAL.LIST <- LIST.FRAC
          threshold <- as.numeric(gsub(pattern="^(\\d+)\\%\\.*", replacement="\\1", x=min.fracs, perl=TRUE))/100
        }  else  stop("'min.fracs' must either be an integer specifying the minimum number of fractions that a taxon must be present within a replicate or a character string specifying the minimum percentage of total fractions (e.g., '40%')")
      }
    
    #Screen taxa according to the specified threshold for the minimum number of reps:
      for (i in 1:dim(fracs.in.rep)[1]){
        TAXA.TO.KEEP[[i]] <- names(FINAL.LIST[[i]])[FINAL.LIST[[i]] >= threshold]
      }
      taxa.to.keep <- names(table(unlist(TAXA.TO.KEEP)))[table(unlist(TAXA.TO.KEEP)) >= min.reps]

    #Number of taxa filtered:
      tot.starting.taxa <- length(levels(factor(DATA[,vars[1]])))
      tot.treatment.taxa <- length(levels(factor(DATA.treatment[,vars[1]])))
      print(paste("Total number of taxa occuring in the specified data: ", tot.starting.taxa, sep=""))
      print(paste("Number of taxa occuring in the specified treatment: ", tot.treatment.taxa, sep=""))
      print(paste("Number of taxa that occurred in ≥ ", min.fracs, " fractions of ≥ ", min.reps, " replicates of the specified treatment: ", length(taxa.to.keep), sep=""))

    #Now, subset the DATA dataframe FOR REAL so that it only contains taxon-tubes for taxa occurring in at least 'min.fracs' fractions of 'min.reps' tubes across the specified treatment:
      # print(paste("Dimensions of the specified data frame (before filtering): ", paste(dim(DATA), collapse="   "), sep=""))
      # print(paste("Dimensions of the specified data frame including only those taxa that do not occur in the specified treatment: ", paste(dim(DATA[!DATA[,vars[1]] %in% DATA.treatment[,vars[1]],]), collapse="   "), sep=""))
      # print(paste("Dimensions of the specified data frame including only those taxa that do not occur in ≥ ", min.fracs, " fractions of ≥ ", min.reps, " replicates of the specified treatment: ", paste(dim(DATA[!DATA[,vars[1]] %in% taxa.to.keep,]), collapse="   "), sep=""))
      DATA <- DATA[DATA[,vars[1]] %in% taxa.to.keep,]
        #Renumber row.names:
        row.names(DATA) <- 1:dim(DATA)[1]
        #Convert the factor columns to factor:
        DATA <- as.data.frame(lapply(DATA, function(x) if(is.factor(x)) factor(x) else x))
      # print(paste("Dimensions of the specified data frame (after filtering): ", paste(dim(DATA), collapse="   "), sep=""))
      DATA
  }



