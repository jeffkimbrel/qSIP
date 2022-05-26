### Given a list of two data frames in the form of output from 'WAD.by.taxon.func', create plots of WAD metrics by taxon for each treatment and create a data.frame of these values
#
#     output = SE.WAD.by.taxon.plot(LIST, percentile=0.95)
#
#     LIST: a list of two data frames in the same form as output from 'WAD.by.taxon.func'; the first is observed WADs by replicate for all taxa; the second gives the replicates for each treatment
#     percentile: a single value specifying where to draw a line on each graph indicating the specified percentile of non-NA standard error values; default is 0.95
#     -------------------------------------------------------
#     output:  plots (one for each treatment contained within 'LIST') of the standard error (SE) of observed WAD values for each taxon
#              a data.frame of the SEs of observed WAD values for all taxa (rows) and all treatments (columns)
#
#     notes:    this function is best called in conjuntion with pdf() or similar to write graphs to a multi-page file where each page correspnds to a treatment
#               color-coding indicates the level of replication for each taxon; bars indicate taxa with 'missing' SE values because N=0 or N=1
# 
#     Written by Ben Koch


  SE.WAD.by.taxon.plot <- function(LIST, percentile=0.95){
    
    obs.wads.by.taxon <- LIST[[1]]
    reps.by.trt <- LIST[[2]]
    
    #Colors for different levels of replication:
    ink <- c("grey80", rainbow(n=12)[c(3,1,5,7,9,11)])
    ink[4] <-"#24B300FF"
    
    #Establish common xlim for all plots:
    all.se <- matrix(NA, nrow=dim(obs.wads.by.taxon)[1], ncol=dim(reps.by.trt)[1])
    for (m in 1:dim(reps.by.trt)[1]){
      curr.reps <- as.character(unlist(reps.by.trt[m,2:dim(reps.by.trt)[2]]))
      curr.reps <- curr.reps[!is.na(curr.reps)]
      curr.data <- cbind(taxon=obs.wads.by.taxon[,1], obs.wads.by.taxon[,names(obs.wads.by.taxon) %in% curr.reps])
      if (length(2:dim(curr.data)[2]) > 1){
        all.se[,m] <- apply(curr.data[,2:dim(curr.data)[2]], 1, function(x) sd(x, na.rm=TRUE)/sqrt(sum(!is.na(x))))
      }  else  all.se[,m] <- rep(NA, dim(curr.data)[1])
    }
    x.min <- 0
    x.max <- max(all.se, na.rm=TRUE)
    all.se <- data.frame(obs.wads.by.taxon[,1], all.se)
    names(all.se) <- c("taxon", as.character(reps.by.trt[,1]))
    all.se <- all.se[order(as.numeric(as.character(all.se[,1]))),]
    row.names(all.se) <- 1:dim(all.se)[1]
    
    #Establish y-axis tick-mark locations:
    AT.y <- seq(0, max(as.numeric(as.character(curr.data$taxon))), 50)
    AT.y.minor <- seq(0, max(as.numeric(as.character(curr.data$taxon))), 10)
    AT.y.minor2 <- seq(0, max(as.numeric(as.character(curr.data$taxon))), 5)
  
    #Plot of SE of observed WADS for all taxa by treatment:
    for (m in 1:dim(reps.by.trt)[1]){
      curr.reps <- as.character(unlist(reps.by.trt[m,2:dim(reps.by.trt)[2]]))
      curr.reps <- curr.reps[!is.na(curr.reps)]
      curr.data <- data.frame(taxon=obs.wads.by.taxon[,1], obs.wads.by.taxon[,names(obs.wads.by.taxon) %in% curr.reps])
      names(curr.data)[2:dim(curr.data)[2]] <- names(obs.wads.by.taxon)[names(obs.wads.by.taxon) %in% curr.reps]
      if (length(2:dim(curr.data)[2]) > 1){
        curr.data$se <- apply(curr.data[,2:dim(curr.data)[2]], 1, function(x) sd(x, na.rm=TRUE)/sqrt(sum(!is.na(x))))
        curr.data$n <- apply(curr.data[,2:(dim(curr.data)[2]-1)], 1, function(x) sum(!is.na(x)))
      }  else  if (length(2:dim(curr.data)[2]) <= 1){
        curr.data$se <- rep(NA, dim(curr.data)[1])
        curr.data$n <- as.numeric(!is.na(curr.data[,2:(dim(curr.data)[2]-1)]))
      }
      curr.data <- curr.data[order(as.numeric(as.character(curr.data[,1]))),]
      row.names(curr.data) <- 1:dim(curr.data)[1]
      par(mai=c(0.44,0.57,0.22,1.60), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
      plot(x=curr.data$se, y=as.numeric(as.character(curr.data$taxon)), type="n", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(x.min, x.max))
        if (sum(curr.data$n == 0) > 0){
          rect(xleft=x.max*-100, ybottom=as.numeric(as.character(curr.data$taxon[curr.data$n == 0]))-0.5, xright=x.max*100, ytop=as.numeric(as.character(curr.data$taxon[curr.data$n == 0]))+0.5, col=ink[1], border="grey60", lwd=1)
        }
        if (sum(curr.data$n == 1) > 0){
          rect(xleft=x.max*-100, ybottom=as.numeric(as.character(curr.data$taxon[curr.data$n == 1]))-0.5, xright=x.max*100, ytop=as.numeric(as.character(curr.data$taxon[curr.data$n == 1]))+0.5, col=ink[2], border="grey60", lwd=1)
        }
        abline(v=quantile(curr.data$se, probs=percentile, na.rm=TRUE), lwd=2.5, lty=2, col="grey70")
        text(x=quantile(curr.data$se, probs=percentile, na.rm=TRUE), y=min(as.numeric(as.character(curr.data$taxon)))-(0.04*(max(as.numeric(as.character(curr.data$taxon)))-1)), labels=paste(eval(percentile*100), "%", sep=""), adj=c(1.1,-0.25), col="grey70", cex=0.75)
        text(x=quantile(curr.data$se, probs=percentile, na.rm=TRUE), y=max(as.numeric(as.character(curr.data$taxon)))+(0.04*(max(as.numeric(as.character(curr.data$taxon)))-1)), labels=paste(eval(percentile*100), "%", sep=""), adj=c(1.1,1.25), col="grey70", cex=0.75)
        points(x=curr.data$se[curr.data$n == 2], y=as.numeric(as.character(curr.data$taxon[curr.data$n == 2])), pch=21, cex=0.7, col="black", bg=ink[3])
        points(x=curr.data$se[curr.data$n == 3], y=as.numeric(as.character(curr.data$taxon[curr.data$n == 3])), pch=21, cex=0.7, col="black", bg=ink[4])
        points(x=curr.data$se[curr.data$n == 4], y=as.numeric(as.character(curr.data$taxon[curr.data$n == 4])), pch=21, cex=0.7, col="black", bg=ink[5])
        points(x=curr.data$se[curr.data$n == 5], y=as.numeric(as.character(curr.data$taxon[curr.data$n == 5])), pch=21, cex=0.7, col="black", bg=ink[6])
        points(x=curr.data$se[curr.data$n >= 6], y=as.numeric(as.character(curr.data$taxon[curr.data$n >= 6])), pch=21, cex=0.7, col="black", bg=ink[7])
        par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
        axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
        mtext("SE of observed WADs", side=1, line=1.35, cex=0.75)
        par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
        axis(side=2, at=AT.y, labels=AT.y, tck=-0.015, las=1, cex.axis=0.6)
        axis(side=2, at=AT.y.minor, labels=FALSE, tck=-0.0075)
        axis(side=2, at=AT.y.minor2, labels=FALSE, tck=-0.0075/2)
        axis(side=4, at=AT.y, labels=AT.y, tck=-0.015, las=1, cex.axis=0.6)
        axis(side=4, at=AT.y.minor, labels=FALSE, tck=-0.0075)
        axis(side=4, at=AT.y.minor2, labels=FALSE, tck=-0.0075/2)
        mtext("Taxon", side=2, line=2.4, cex=0.75)
        mtext(reps.by.trt[m,1], side=3, line=0, cex=1)
        par(xpd=NA)
        legend(x=x.max*1.20, y=max(as.numeric(as.character(curr.data$taxon)))+(0.04*(max(as.numeric(as.character(curr.data$taxon)))-1)), legend=c("N = 0", "N = 1", "N = 2", "N = 3", "N = 4", "N = 5", expression("N">=6)), col=c("grey60","grey60","white","white","white","white","white"), pt.bg=ink, pch=c(NA,NA,21,21,21,21,21), pt.cex=0, lwd=c(4.5,4.5,0,0,0,0,0), bty="o", cex=0.6, xjust=0, yjust=1)
        legend(x=x.max*1.20, y=max(as.numeric(as.character(curr.data$taxon)))+(0.04*(max(as.numeric(as.character(curr.data$taxon)))-1)), legend=c("N = 0", "N = 1", "N = 2", "N = 3", "N = 4", "N = 5", expression("N">=6)), col=c(ink[1],ink[2],"white","white","white","white","white"), pt.bg=ink, pch=c(NA,NA,21,21,21,21,21), pt.cex=0, lwd=c(2.0,2.0,0,0,0,0,0), bty="o", cex=0.6, xjust=0, yjust=1)
        legend(x=x.max*1.20, y=max(as.numeric(as.character(curr.data$taxon)))+(0.04*(max(as.numeric(as.character(curr.data$taxon)))-1)), legend=c("N = 0", "N = 1", "N = 2", "N = 3", "N = 4", "N = 5", expression("N">=6)), col="black", pt.bg=ink, pch=c(NA,NA,21,21,21,21,21), pt.cex=0.7, lty=0, lwd=c(0,0,1,1,1,1,1), bty="n", cex=0.6, xjust=0, yjust=1)
        par(xpd=FALSE)
    }
    return(all.se)

  }
