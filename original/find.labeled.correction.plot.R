#Function to create plots for assessing the fit of each iterated subset of taxa for identifying non-growers:

#NOTE: this plotting code assumes the units of density are g/mL
#NOTE: these plots can only display up to 120 iterations
find.labeled.correction.plot <- function(find.labeled.correction.list, filename="", path="qSIP_output/Figures/", highlight=NULL, highlight.col=NULL){
  putative.nongrower.metrics <- find.labeled.correction.list$putative.nongrower.metrics
  curr.data.lowSD <- find.labeled.correction.list$data.low.SD
  lab.replicate <- find.labeled.correction.list$lab.replicate
  unlab.corr.names <- find.labeled.correction.list$unlab.corr.names
  lab.names <- find.labeled.correction.list$lab.names
  method <- find.labeled.correction.list$method
  unlab.SD.percentile <- find.labeled.correction.list$unlab.SD.percentile
  lab.SD.percentile <- find.labeled.correction.list$lab.SD.percentile
  min.num.nongrowers <- find.labeled.correction.list$min.num.nongrowers
  delta.WAD <- curr.data.lowSD[, names(curr.data.lowSD) == lab.replicate] - curr.data.lowSD$unlab.mean
  #Set the proper number of colors to match the highlighted points:
    if(length(highlight) > length(highlight.col)){
      highlight.col <- rep(highlight.col, length.out=length(highlight))
    } else if(length(highlight) < length(highlight.col)){
      highlight.col <- highlight.col[1:length(highlight)]
    }
  ###PLOT OF WADS FROM LABELED TUBE VS. UNLABELED MEAN WAD; FROM L TO R AND TOP TO BOTTOM:
  #Black line = regression line with a fixed slope of 1 and a fitted intercept term
  #Blue line = 1:1 line (the expected relationship after correcting for the offset - i.e., the intercept term in the fixed slope=1 regression above - assuming all growers have been removed from the dataset)
  #Red dotted line = regression line (with fitted intercept and slope terms) of the offset-corrected data
    # dev.new(width=24, height=18.5)
    pdf(file=paste(path, filename, ".pdf", sep=""), width=24, height=18.5)
    par(mfrow=c(10, 12))
    par(mai=c(0.44,0.57,0.22,0.11), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
    for (j in 1:dim(putative.nongrower.metrics)[1]){
      if (method == "bu.abs.resid"){
        iteration.j.points <- c(as.character(curr.data.lowSD$taxon[delta.WAD %in% sort(delta.WAD)[1:min.num.nongrowers]]), putative.nongrower.metrics$next.taxon[1:j])
        X.plot <- curr.data.lowSD$unlab.mean[curr.data.lowSD$taxon %in% iteration.j.points]
        Y.plot <- curr.data.lowSD[curr.data.lowSD$taxon %in% iteration.j.points, names(curr.data.lowSD) == lab.replicate]
        XLIM <- YLIM <- c(min(curr.data.lowSD$unlab.mean, curr.data.lowSD[, names(curr.data.lowSD) == lab.replicate]), max(curr.data.lowSD$unlab.mean, curr.data.lowSD[, names(curr.data.lowSD) == lab.replicate]))
      }
      else{   #for method="td.abs.resid" & "method=td.pos.resid"
        X.plot <- curr.data.lowSD$unlab.mean[!(curr.data.lowSD$taxon %in% putative.nongrower.metrics$next.taxon[1:j])]
        Y.plot <- curr.data.lowSD[!(curr.data.lowSD$taxon %in% putative.nongrower.metrics$next.taxon[1:j]), names(curr.data.lowSD) == lab.replicate]
        if (j == 1){
          XLIM <- YLIM <- c(min(c(X.plot, Y.plot)), max(c(X.plot, Y.plot)))
        }
      }
      plot(x=X.plot, y=Y.plot, xlim=XLIM, ylim=YLIM, cex=0.5, pch=21, col="black", bg="black", xlab="", ylab="", xaxt="n", yaxt="n")
      par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
      axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
      mtext("Unlabeled mean WAD", side=1, line=1.35, cex=0.75)
      par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
      axis(side=2, tck=-0.015, las=1, cex.axis=0.6)
      mtext("Labeled WAD", side=2, line=2.4, cex=0.75)
      mtext(paste("iteration ", j, "; taxon ", putative.nongrower.metrics$next.taxon[j], sep=""), side=3, line=0.1, cex=0.75)
      if (j %in% highlight){
        mtext(paste("iteration ", j, "; taxon ", putative.nongrower.metrics$next.taxon[j], sep=""), side=3, line=0.1, cex=0.75, col=highlight.col[j == highlight])
      }
      abline(a=0, b=1, col="blue", lwd=2)
      abline(a=putative.nongrower.metrics$fixed.intercept[j], b=putative.nongrower.metrics$fixed.slope[j], col="black", lwd=2)
      abline(a=putative.nongrower.metrics$corr.free.intercept[j], b=putative.nongrower.metrics$corr.free.slope[j], col="red", lwd=2, lty=3)
      if (method == "bu.abs.resid" & j == 1){
        legend(x=XLIM[1]-0.03*(XLIM[2]-XLIM[1]), y=YLIM[2]+0.03*(YLIM[2]-YLIM[1]), legend=c(expression(paste("corr, ", italic("a"), "&", italic("b"), "=free", sep="")), expression(paste(italic("b"), "=1, ", italic("a"), "=free", sep="")), "1:1"), bty="n", lty=c(3,1,1), lwd=c(2,2,2), col=c("red", "black", "blue"), cex=0.6)
      }
      if (method != "bu.abs.resid" & j == dim(putative.nongrower.metrics)[1]){
        legend(x=XLIM[1]-0.03*(XLIM[2]-XLIM[1]), y=YLIM[2]+0.03*(YLIM[2]-YLIM[1]), legend=c(expression(paste("corr, ", italic("a"), "&", italic("b"), "=free", sep="")), expression(paste(italic("b"), "=1, ", italic("a"), "=free", sep="")), "1:1"), bty="n", lty=c(3,1,1), lwd=c(2,2,2), col=c("red", "black", "blue"), cex=0.6)
      }
    }
    par(mfrow=c(1, 1))

  #NEXT PAGE:
  ###PLOT OF RESIDUALS OF WADS FROM LABELED TUBE VS. UNLABELED MEAN WAD; FROM L TO R AND TOP TO BOTTOM:
  #Black line = predicted mean line (residual=0)
    # dev.new(width=24, height=18.5)
    par(mfrow=c(10, 12))
    par(mai=c(0.44,0.57,0.22,0.11), cex=1, mex=0.75)				     #set the margins (bottom, left, top, right) in inches
    for (j in 1:dim(putative.nongrower.metrics)[1]){
      if (method == "bu.abs.resid"){
        X.all <- curr.data.lowSD$unlab.mean
        Y.all <- curr.data.lowSD[, names(curr.data.lowSD) == lab.replicate]
        iteration.j.points <- c(as.character(curr.data.lowSD$taxon[delta.WAD %in% sort(delta.WAD)[1:min.num.nongrowers]]), putative.nongrower.metrics$next.taxon[1:j])
        X.plot <- curr.data.lowSD$unlab.mean[curr.data.lowSD$taxon %in% iteration.j.points]
        Y.plot <- curr.data.lowSD[curr.data.lowSD$taxon %in% iteration.j.points, names(curr.data.lowSD) == lab.replicate]
        mod.all.temp <- lm(I(Y.all-1*X.all)~1)                       #fixed slope of 1; note that the function "I()" is used to inhibit the interpretation of operators such as "+", "-", "*" and "^" as formula operators, so they are instead used as arithmetical operators.
      }
      else{   #for method="td.abs.resid" & "method=td.pos.resid"
        X.plot <- curr.data.lowSD$unlab.mean[!(curr.data.lowSD$taxon %in% putative.nongrower.metrics$next.taxon[1:j])]
        Y.plot <- curr.data.lowSD[!(curr.data.lowSD$taxon %in% putative.nongrower.metrics$next.taxon[1:j]), names(curr.data.lowSD) == lab.replicate]
      }
      mod.temp <- lm(I(Y.plot-1*X.plot)~1)                       #fixed slope of 1; note that the function "I()" is used to inhibit the interpretation of operators such as "+", "-", "*" and "^" as formula operators, so they are instead used as arithmetical operators.
      resid.mod.temp <- lm(as.numeric(resid(mod.temp))~X.plot)   #regression of residuals on x-values (unabeled mean WAD)
      if (method == "bu.abs.resid"){
        XLIM <- c(min(X.all), max(X.all))    
        YLIM <- c(min(as.numeric(resid(mod.all.temp))), max(as.numeric(resid(mod.all.temp))))
      }
      else{   #for method="td.abs.resid" & "method=td.pos.resid"
        if (j == 1){
          XLIM <- c(min(X.plot), max(X.plot))    
          YLIM <- c(min(as.numeric(resid(mod.temp))), max(as.numeric(resid(mod.temp))))
        }
      }
      plot(x=X.plot, y=as.numeric(resid(mod.temp)), xlim=XLIM, ylim=YLIM, cex=0.5, pch=21, col="black", bg="black", xlab="", ylab="", xaxt="n", yaxt="n")
      par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
      axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
      mtext("Unlabeled mean WAD", side=1, line=1.35, cex=0.75)
      par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
      axis(side=2, tck=-0.015, las=1, cex.axis=0.6)
      mtext("Residual of labeled WAD", side=2, line=2.4, cex=0.75)
      mtext(paste("iteration ", j, "; taxon ", putative.nongrower.metrics$next.taxon[j], sep=""), side=3, line=0.1, cex=0.75)
      if (j %in% highlight){
        mtext(paste("iteration ", j, "; taxon ", putative.nongrower.metrics$next.taxon[j], sep=""), side=3, line=0.1, cex=0.75, col=highlight.col[j == highlight])
      }
      abline(h=0, col="black", lwd=2)
      abline(a=putative.nongrower.metrics$fixed.resid.intercept[j], b=putative.nongrower.metrics$fixed.resid.slope[j], col="red", lwd=2, lty=3)
      if (method == "bu.abs.resid" & j == 1){
        legend(x=XLIM[1]-0.03*(XLIM[2]-XLIM[1]), y=YLIM[2]+0.03*(YLIM[2]-YLIM[1]), legend=c("residuals vs. x", expression(paste("pred mean with ", italic("b"), "=1", sep=""))), bty="n", lty=c(3,1), lwd=c(2,2), col=c("red", "black"), cex=0.6)
      }
      if (method != "bu.abs.resid" & j == dim(putative.nongrower.metrics)[1]){
        legend(x=XLIM[1]-0.03*(XLIM[2]-XLIM[1]), y=YLIM[2]+0.03*(YLIM[2]-YLIM[1]), legend=c("residuals vs. x", expression(paste("pred mean with ", italic("b"), "=1", sep=""))), bty="n", lty=c(3,1), lwd=c(2,2), col=c("red", "black"), cex=0.6)
      }
    }
    par(mfrow=c(1, 1))

  #NEXT PAGE:
  par(mfrow=c(3,4))
  par(mai=c(0.75,0.80,0.40,0.15), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches

  #PLOT OF N VS ROW (ITERATION):
  plot(x=1:dim(putative.nongrower.metrics)[1], y=putative.nongrower.metrics$N, pch=21, col="black", bg="skyblue", xlab="", ylab="", xaxt="n", yaxt="n")
  points(x=as.numeric(row.names(putative.nongrower.metrics))[highlight], y=putative.nongrower.metrics$N[highlight], pch=21, col="black", bg=highlight.col)
  par(mgp=c(3,0.60,0))			#set spacing for the x-axis labels (the second value)
  axis(side=1, tck=-0.010, las=1, cex.axis=1.25)
  mtext("Iteration", side=1, line=3.0, cex=1.6)
  par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
  axis(side=2, tck=-0.010, cex.axis=1.25)
  mtext("N", side=2, line=2.55, cex=1.6)

  #NORMAL-FITTED TUBE-LEVEL CORRECTION FACTOR VS ROW (ITERATION):
  intercept.correction <- putative.nongrower.metrics$fixed.intercept
  norm.correction <- (putative.nongrower.metrics$fixed.norm.mean.Y - putative.nongrower.metrics$norm.mean.X)
  YLIM <- c(min(c(intercept.correction, norm.correction)), max(c(intercept.correction, norm.correction)))
  plot(x=1:dim(putative.nongrower.metrics)[1], y=norm.correction, pch=21, col="black", bg="skyblue", ylim=YLIM, xlab="", ylab="", xaxt="n", yaxt="n")  
  points(x=as.numeric(row.names(putative.nongrower.metrics))[highlight], y=norm.correction[highlight], pch=21, col="black", bg=highlight.col)
  par(mgp=c(3,0.60,0))			#set spacing for the x-axis labels (the second value)
  axis(side=1, tck=-0.010, las=1, cex.axis=1.25)
  mtext("Iteration", side=1, line=3.0, cex=1.6)
  par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
  axis(side=2, tck=-0.010, cex.axis=1.25)
  mtext(expression(paste("Normal-fitted WAD correction (g ml"^-1, ")", sep="")), side=2, line=2.55, cex=1.6)

  #COMPARING THE TWO WAYS TO CORRECT WADS IN LABELED TUBES:
  #Note that (putative.nongrower.metrics$fixed.mean.Y - putative.nongrower.metrics$mean.X) ...is equal to... putative.nongrower.metrics$fixed.intercept
  intercept.correction <- putative.nongrower.metrics$fixed.intercept
  norm.correction <- (putative.nongrower.metrics$fixed.norm.mean.Y - putative.nongrower.metrics$norm.mean.X)
  XLIM <- YLIM <- c(min(c(intercept.correction, norm.correction)), max(c(intercept.correction, norm.correction)))
  plot(x=intercept.correction, y=norm.correction, pch=21, col="black", bg="skyblue", xlim=XLIM, ylim=YLIM, xlab="", ylab="", xaxt="n", yaxt="n")
  points(x=intercept.correction[highlight], y=norm.correction[highlight], pch=21, col="black", bg=highlight.col)
  arrows(x0=intercept.correction[1:(length(intercept.correction)-1)], y0=norm.correction[1:(length(norm.correction)-1)], x1=intercept.correction[2:length(intercept.correction)], y1=norm.correction[2:length(norm.correction)], lwd=1, length=0.05, angle=30, code=2, col="blue")
  abline(a=0, b=1)
  par(mgp=c(3,0.60,0))			#set spacing for the x-axis labels (the second value)
  axis(side=1, tck=-0.010, las=1, cex.axis=1.25)
  mtext(expression(paste("Intercept-calculated WAD correction (g ml"^-1, ")", sep="")), side=1, line=3.4, padj=0, cex=1.6)  #padj needed to adjust vertical alignment when string contains superscript characters
  par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
  axis(side=2, tck=-0.010, cex.axis=1.25)
  mtext(expression(paste("Normal-fitted WAD correction (g ml"^-1, ")", sep="")), side=2, line=2.55, cex=1.6)

  #NORMAL-FITTED AICC (ATTEMPTS TO CORRECT FOR SAMPLE SIZE) VS ROW (ITERATION):
  k <- 2
  AICc.X <- (2*k) + (2*putative.nongrower.metrics$norm.NLL.X) + ((2*(k+1)*(k+2))/(putative.nongrower.metrics$N-k-2))
  AICc.Y <- (2*k) + (2*putative.nongrower.metrics$corr.free.norm.NLL.Y) + ((2*(k+1)*(k+2))/(putative.nongrower.metrics$N-k-2))
  YLIM <- c(min(c(AICc.X, AICc.Y)), max(c(AICc.X, AICc.Y)))
  plot(x=0, y=1, type="n", xlim=c(1, dim(putative.nongrower.metrics)[1]), ylim=YLIM, xlab="", ylab="", xaxt="n", yaxt="n")
  points(x=1:dim(putative.nongrower.metrics)[1], y=AICc.X, pch=21, col="black", bg="white")
  points(x=1:dim(putative.nongrower.metrics)[1], y=AICc.Y, pch=21, col="black", bg="green")
  points(x=as.numeric(row.names(putative.nongrower.metrics))[highlight], y=AICc.Y[highlight], pch=21, col="black", bg=highlight.col)
  par(mgp=c(3,0.60,0))			#set spacing for the x-axis labels (the second value)
  axis(side=1, tck=-0.010, las=1, cex.axis=1.25)
  mtext("Iteration", side=1, line=3.0, cex=1.6)
  par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
  axis(side=2, tck=-0.010, cex.axis=1.25)
  mtext("Normal-fitted AICc", side=2, line=2.55, cex=1.6)

  #CORRELATION VS ROW (ITERATION):
  YLIM.corr <- c(min(putative.nongrower.metrics$fixed.correl, putative.nongrower.metrics$corr.free.correl), max(putative.nongrower.metrics$fixed.correl, putative.nongrower.metrics$corr.free.correl))
  YLIM.r2 <- c(min(putative.nongrower.metrics$fixed.r2.adj, putative.nongrower.metrics$corr.free.r2.adj), max(putative.nongrower.metrics$fixed.r2.adj, putative.nongrower.metrics$corr.free.r2.adj))
  YLIM <- c(min(0, YLIM.corr, YLIM.r2), max(1, YLIM.corr, YLIM.r2))
  plot(x=1:dim(putative.nongrower.metrics)[1], y=putative.nongrower.metrics$fixed.correl, pch=21, col="black", bg="grey80", ylim=YLIM, xlab="", ylab="", xaxt="n", yaxt="n")
  points(x=1:dim(putative.nongrower.metrics)[1], y=putative.nongrower.metrics$corr.free.correl, pch=21, col="black", bg="black")
  points(x=as.numeric(row.names(putative.nongrower.metrics))[highlight], y=putative.nongrower.metrics$corr.free.correl[highlight], pch=21, col="black", bg=highlight.col)
  par(mgp=c(3,0.60,0))			#set spacing for the x-axis labels (the second value)
  axis(side=1, tck=-0.010, las=1, cex.axis=1.25)
  mtext("Iteration", side=1, line=3.0, cex=1.6)
  par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
  axis(side=2, tck=-0.010, cex.axis=1.25)
  mtext("Correlation coefficient", side=2, line=2.55, cex=1.6)
  # fixed.correl is the same as corr.free.correl

  #ADJUSTED R2 VS ROW (ITERATION):
  YLIM.corr <- c(min(putative.nongrower.metrics$fixed.correl, putative.nongrower.metrics$corr.free.correl), max(putative.nongrower.metrics$fixed.correl, putative.nongrower.metrics$corr.free.correl))
  YLIM.r2 <- c(min(putative.nongrower.metrics$fixed.r2.adj, putative.nongrower.metrics$corr.free.r2.adj), max(putative.nongrower.metrics$fixed.r2.adj, putative.nongrower.metrics$corr.free.r2.adj))
  YLIM <- c(min(0, YLIM.corr, YLIM.r2), max(1, YLIM.corr, YLIM.r2))
  plot(x=1:dim(putative.nongrower.metrics)[1], y=putative.nongrower.metrics$fixed.r2.adj, pch=21, col="black", bg="grey80", ylim=YLIM, xlab="", ylab="", xaxt="n", yaxt="n")
  points(x=1:dim(putative.nongrower.metrics)[1], y=putative.nongrower.metrics$corr.free.r2.adj, pch=21, col="black", bg="black")
  points(x=as.numeric(row.names(putative.nongrower.metrics))[highlight], y=putative.nongrower.metrics$corr.free.r2.adj[highlight], pch=21, col="black", bg=highlight.col)
  par(mgp=c(3,0.60,0))			#set spacing for the x-axis labels (the second value)
  axis(side=1, tck=-0.010, las=1, cex.axis=1.25)
  mtext("Iteration", side=1, line=3.0, cex=1.6)
  par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
  axis(side=2, tck=-0.010, cex.axis=1.25)
  mtext(expression(paste("Adjusted ", italic("r")^2, sep="")), side=2, line=2.55, cex=1.6)

  #SLOPE VS ROW (ITERATION):
  YLIM <- c(min(c(1, putative.nongrower.metrics$fixed.slope, putative.nongrower.metrics$corr.free.slope)), max(c(1, putative.nongrower.metrics$fixed.slope, putative.nongrower.metrics$corr.free.slope)))
  plot(x=0, y=0, type="n", xlim=c(1, dim(putative.nongrower.metrics)[1]), ylim=YLIM, xlab="", ylab="", xaxt="n", yaxt="n")
  abline(h=1)
  points(x=1:dim(putative.nongrower.metrics)[1], y=putative.nongrower.metrics$fixed.slope, pch=21, col="black", bg="grey80")
  points(x=1:dim(putative.nongrower.metrics)[1], y=putative.nongrower.metrics$corr.free.slope, pch=21, col="black", bg="black")
  points(x=as.numeric(row.names(putative.nongrower.metrics))[highlight], y=putative.nongrower.metrics$corr.free.slope[highlight], pch=21, col="black", bg=highlight.col)
  par(mgp=c(3,0.60,0))			#set spacing for the x-axis labels (the second value)
  axis(side=1, tck=-0.010, las=1, cex.axis=1.25)
  mtext("Iteration", side=1, line=3.0, cex=1.6)
  par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
  axis(side=2, tck=-0.010, cex.axis=1.25)
  mtext(expression(paste("Slope (g ml"^-1, "/ g ml"^-1, ")", sep="")), side=2, line=2.55, cex=1.6)

  #INTERCEPT VS ROW (ITERATION):
  YLIM <- c(min(c(0, putative.nongrower.metrics$fixed.intercept, putative.nongrower.metrics$corr.free.intercept)), max(c(0, putative.nongrower.metrics$fixed.intercept, putative.nongrower.metrics$corr.free.intercept)))
  plot(x=0, y=0, type="n", xlim=c(1, dim(putative.nongrower.metrics)[1]), ylim=YLIM, xlab="", ylab="", xaxt="n", yaxt="n")
  abline(h=0)
  points(x=1:dim(putative.nongrower.metrics)[1], y=putative.nongrower.metrics$fixed.intercept, pch=21, col="black", bg="grey80")
  points(x=1:dim(putative.nongrower.metrics)[1], y=putative.nongrower.metrics$corr.free.intercept, pch=21, col="black", bg="black")
  points(x=as.numeric(row.names(putative.nongrower.metrics))[highlight], y=putative.nongrower.metrics$corr.free.intercept[highlight], pch=21, col="black", bg=highlight.col)
  par(mgp=c(3,0.60,0))			#set spacing for the x-axis labels (the second value)
  axis(side=1, tck=-0.010, las=1, cex.axis=1.25)
  mtext("Iteration", side=1, line=3.0, cex=1.6)
  par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
  axis(side=2, tck=-0.010, cex.axis=1.25)
  mtext(expression(paste("Intercept (g ml"^-1, ")", sep="")), side=2, line=2.55, cex=1.6)

  #(TWO EMPTY PLOTS):
  XLIM <- YLIM <- c(0,1)
  plot(x=-100, y=-100, xlim=XLIM, ylim=YLIM, xlab="", ylab="", bty="n", xaxt="n", yaxt="n")
  par(xpd=NA)
  legend(x=XLIM[1]-0.03*(XLIM[2]-XLIM[1]), y=YLIM[2]+0.03*(YLIM[2]-YLIM[1]), legend=c("from fixed slope (b=1, a=free) regression of labeled WADs vs. unlabeled mean WADs", "from free parameter (a&b=free) regression of labeled WADs vs. unlabeled mean WADs", "from unlabeled mean WADs", "from normal-fit-corrected labeled WADs", "target slope & intercept"), bty="o", pch=c(21,21,21,21,3), col=c("black", "black", "black", "black", "black"), pt.bg=c("grey80", "black", "white", "green", NULL), pt.cex=c(1,1,1,1,1.5), cex=1, xjust=0, yjust=1)
  text(x=XLIM[1]-0.03*(XLIM[2]-XLIM[1]), y=0.76, labels=paste("lab.replicate: ", eval(lab.replicate), sep=""), adj=c(0,1))
  text(x=XLIM[1]-0.03*(XLIM[2]-XLIM[1]), y=0.72, labels=paste("method: ", eval(method), sep=""), adj=c(0,1))
  text(x=XLIM[1]-0.03*(XLIM[2]-XLIM[1]), y=0.68, labels=paste("unlab.SD.percentile: ", eval(unlab.SD.percentile), sep=""), adj=c(0,1))
  text(x=XLIM[1]-0.03*(XLIM[2]-XLIM[1]), y=0.64, labels=paste("lab.SD.percentile: ", eval(lab.SD.percentile), sep=""), adj=c(0,1))
  text(x=XLIM[1]-0.03*(XLIM[2]-XLIM[1]), y=0.60, labels=paste("min.num.nongrowers: ", eval(min.num.nongrowers), sep=""), adj=c(0,1))
  text(x=XLIM[1]-0.03*(XLIM[2]-XLIM[1]), y=0.56, labels=paste("lab.names: ", paste(eval(lab.names), collapse=", "), sep=""), adj=c(0,1))
  text(x=XLIM[1]-0.03*(XLIM[2]-XLIM[1]), y=0.52, labels=paste("unlab.corr.names: ", paste(eval(unlab.corr.names), collapse="\n                             "), sep=""), adj=c(0,1))   #note: spaces in the 'collapse' argument here ensure proper spacing of text on the graph
  par(xpd=FALSE)
  plot(x=-100, y=-100, xlim=XLIM, ylim=YLIM, xlab="", ylab="", bty="n", xaxt="n", yaxt="n")

  #DIFFERENCE IN SLOPE (CORRECTED-FREE - FIXED) VS ROW (ITERATION):
  YLIM <- c(min(0, putative.nongrower.metrics$corr.free.slope - putative.nongrower.metrics$fixed.slope), max(0, putative.nongrower.metrics$corr.free.slope - putative.nongrower.metrics$fixed.slope))
  plot(x=0, y=0, type="n", xlim=c(1, dim(putative.nongrower.metrics)[1]), ylim=YLIM, xlab="", ylab="", xaxt="n", yaxt="n")
  abline(h=0)
  points(x=1:dim(putative.nongrower.metrics)[1], y=putative.nongrower.metrics$corr.free.slope - putative.nongrower.metrics$fixed.slope, pch=21, col="black", bg="skyblue")
  points(x=as.numeric(row.names(putative.nongrower.metrics))[highlight], y=putative.nongrower.metrics$corr.free.slope[highlight] - putative.nongrower.metrics$fixed.slope[highlight], pch=21, col="black", bg=highlight.col)
  par(mgp=c(3,0.60,0))			#set spacing for the x-axis labels (the second value)
  axis(side=1, tck=-0.010, las=1, cex.axis=1.25)
  mtext("Iteration", side=1, line=3.0, cex=1.6)
  par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
  axis(side=2, tck=-0.010, cex.axis=1.25)
  mtext(expression(paste("Difference in slope (g ml"^-1, "/ g ml"^-1, ")", sep="")), side=2, line=2.55, cex=1.6)

  #SLOPE VS INTERCEPT:
  XLIM <- c(min(c(0, putative.nongrower.metrics$fixed.intercept, putative.nongrower.metrics$corr.free.intercept)), max(c(0, putative.nongrower.metrics$fixed.intercept, putative.nongrower.metrics$corr.free.intercept)))
  YLIM <- c(min(c(1, putative.nongrower.metrics$fixed.slope, putative.nongrower.metrics$corr.free.slope)), max(c(1, putative.nongrower.metrics$fixed.slope, putative.nongrower.metrics$corr.free.slope)))
  plot(x=0, y=1, pch=3, col="black", xlim=XLIM, ylim=YLIM, xlab="", ylab="", xaxt="n", yaxt="n", cex=1.5)
  points(x=putative.nongrower.metrics$fixed.intercept, y=putative.nongrower.metrics$fixed.slope, pch=21, col="black", bg="grey80")
  points(x=putative.nongrower.metrics$corr.free.intercept, y=putative.nongrower.metrics$corr.free.slope, pch=21, col="black", bg="black")
  points(x=putative.nongrower.metrics$corr.free.intercept[highlight], y=putative.nongrower.metrics$corr.free.slope[highlight], pch=21, col="black", bg=highlight.col)
  par(mgp=c(3,0.60,0))			#set spacing for the x-axis labels (the second value)
  axis(side=1, tck=-0.010, las=1, cex.axis=1.25)
  mtext(expression(paste("Intercept (g ml"^-1, ")", sep="")), side=1, line=3.4, padj=0, cex=1.6)  #padj needed to adjust vertical alignment when string contains superscript characters
  par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
  axis(side=2, tck=-0.010, cex.axis=1.25)
  mtext(expression(paste("Slope (g ml"^-1, "/ g ml"^-1, ")", sep="")), side=2, line=2.55, cex=1.6)

  #NEXT PAGE:
  layout(mat=matrix(c(0,1,1,0,0,6,6,0,2,2,3,3,7,7,8,8,4,4,5,5,9,9,10,10), 3, 8, byrow=TRUE), widths=rep(0.5,8), heights=rep(1,3))
  par(mai=c(0.75,0.80,0.40,0.15), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches

  #EMPIRICAL SD VS MEAN (UNCORRECTED X = UNLABELED; CORRECTED FREE Y = LABELED):
  XLIM.emp.mean <- c(min(c(putative.nongrower.metrics$mean.X, putative.nongrower.metrics$corr.free.mean.Y)), max(c(putative.nongrower.metrics$mean.X, putative.nongrower.metrics$corr.free.mean.Y)))
  XLIM.norm.mean <- c(min(c(putative.nongrower.metrics$norm.mean.X, putative.nongrower.metrics$corr.free.norm.mean.Y)), max(c(putative.nongrower.metrics$norm.mean.X, putative.nongrower.metrics$corr.free.norm.mean.Y)))
  XLIM <- c(min(XLIM.emp.mean, XLIM.norm.mean), max(XLIM.emp.mean, XLIM.norm.mean))
  YLIM.emp.sd <- c(min(c(putative.nongrower.metrics$sd.X, putative.nongrower.metrics$corr.free.sd.Y)), max(c(putative.nongrower.metrics$sd.X, putative.nongrower.metrics$corr.free.sd.Y)))
  YLIM.norm.sd <- c(min(c(putative.nongrower.metrics$norm.sd.X, putative.nongrower.metrics$corr.free.norm.sd.Y)), max(c(putative.nongrower.metrics$norm.sd.X, putative.nongrower.metrics$corr.free.norm.sd.Y)))
  YLIM <- c(min(YLIM.emp.sd, YLIM.norm.sd), max(YLIM.emp.sd, YLIM.norm.sd))
  plot(x=0, y=0, type="n", xlim=XLIM, ylim=YLIM, xlab="", ylab="", xaxt="n", yaxt="n")
  points(x=putative.nongrower.metrics$mean.X, y=putative.nongrower.metrics$sd.X, pch=21, col="black", bg="white")
  arrows(x0=putative.nongrower.metrics$mean.X[1:(dim(putative.nongrower.metrics)[1]-1)], y0=putative.nongrower.metrics$sd.X[1:(dim(putative.nongrower.metrics)[1]-1)], x1=putative.nongrower.metrics$mean.X[2:dim(putative.nongrower.metrics)[1]], y1=putative.nongrower.metrics$sd.X[2:dim(putative.nongrower.metrics)[1]], lwd=1, length=0.05, angle=30, code=2, col="black")
  points(x=putative.nongrower.metrics$corr.free.mean.Y, y=putative.nongrower.metrics$corr.free.sd.Y, pch=21, col="black", bg="green")
  points(x=putative.nongrower.metrics$corr.free.mean.Y[highlight], y=putative.nongrower.metrics$corr.free.sd.Y[highlight], pch=21, col="black", bg=highlight.col)
  arrows(x0=putative.nongrower.metrics$corr.free.mean.Y[1:(dim(putative.nongrower.metrics)[1]-1)], y0=putative.nongrower.metrics$corr.free.sd.Y[1:(dim(putative.nongrower.metrics)[1]-1)], x1=putative.nongrower.metrics$corr.free.mean.Y[2:dim(putative.nongrower.metrics)[1]], y1=putative.nongrower.metrics$corr.free.sd.Y[2:dim(putative.nongrower.metrics)[1]], lwd=1, length=0.05, angle=30, code=2, col="blue")
  par(mgp=c(3,0.60,0))			#set spacing for the x-axis labels (the second value)
  axis(side=1, tck=-0.010, las=1, cex.axis=1.25)
  mtext(expression(paste("Empirical mean WAD (g ml"^-1, ")", sep="")), side=1, line=3.4, padj=0, cex=1.6)  #padj needed to adjust vertical alignment when string contains superscript characters
  par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
  axis(side=2, tck=-0.010, cex.axis=1.25)
  mtext(expression(paste("Empirical SD of WAD (g ml"^-1, ")", sep="")), side=2, line=2.55, cex=1.6)
  mtext("Empirical parameters", side=3, line=0.4, cex=2.0)

  #EMPIRICAL MEANS VS ROW (ITERATION):
  YLIM.emp.mean <- c(min(c(putative.nongrower.metrics$mean.X, putative.nongrower.metrics$corr.free.mean.Y)), max(c(putative.nongrower.metrics$mean.X, putative.nongrower.metrics$corr.free.mean.Y)))
  YLIM.norm.mean <- c(min(c(putative.nongrower.metrics$norm.mean.X, putative.nongrower.metrics$corr.free.norm.mean.Y)), max(c(putative.nongrower.metrics$norm.mean.X, putative.nongrower.metrics$corr.free.norm.mean.Y)))
  YLIM <- c(min(YLIM.emp.mean, YLIM.norm.mean), max(YLIM.emp.mean, YLIM.norm.mean))
  plot(x=0, y=1, type="n", xlim=c(1, dim(putative.nongrower.metrics)[1]), ylim=YLIM, xlab="", ylab="", xaxt="n", yaxt="n")
  points(x=1:dim(putative.nongrower.metrics)[1], y=putative.nongrower.metrics$mean.X, pch=21, col="black", bg="white")
  points(x=1:dim(putative.nongrower.metrics)[1], y=putative.nongrower.metrics$corr.free.mean.Y, pch=21, col="black", bg="green")
  par(mgp=c(3,0.60,0))			#set spacing for the x-axis labels (the second value)
  axis(side=1, tck=-0.010, las=1, cex.axis=1.25)
  mtext("Iteration", side=1, line=3.0, cex=1.6)
  par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
  axis(side=2, tck=-0.010, cex.axis=1.25)
  mtext(expression(paste("Empirical mean WAD (g ml"^-1, ")", sep="")), side=2, line=2.55, cex=1.6)
  # mean.X is the same as corr.free.mean.Y if intercept correction is used

  #EMPIRICAL SD VS ROW (ITERATION):
  YLIM.emp.sd <- c(min(c(putative.nongrower.metrics$sd.X, putative.nongrower.metrics$corr.free.sd.Y)), max(c(putative.nongrower.metrics$sd.X, putative.nongrower.metrics$corr.free.sd.Y)))
  YLIM.norm.sd <- c(min(c(putative.nongrower.metrics$norm.sd.X, putative.nongrower.metrics$corr.free.norm.sd.Y)), max(c(putative.nongrower.metrics$norm.sd.X, putative.nongrower.metrics$corr.free.norm.sd.Y)))
  YLIM <- c(min(YLIM.emp.sd, YLIM.norm.sd), max(YLIM.emp.sd, YLIM.norm.sd))
  plot(x=0, y=1, type="n", xlim=c(1, dim(putative.nongrower.metrics)[1]), ylim=YLIM, xlab="", ylab="", xaxt="n", yaxt="n")
  points(x=1:dim(putative.nongrower.metrics)[1], y=putative.nongrower.metrics$sd.X, pch=21, col="black", bg="white")
  points(x=1:dim(putative.nongrower.metrics)[1], y=putative.nongrower.metrics$corr.free.sd.Y, pch=21, col="black", bg="green")
  par(mgp=c(3,0.60,0))			#set spacing for the x-axis labels (the second value)
  axis(side=1, tck=-0.010, las=1, cex.axis=1.25)
  mtext("Iteration", side=1, line=3.0, cex=1.6)
  par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
  axis(side=2, tck=-0.010, cex.axis=1.25)
  mtext(expression(paste("Empirical SD of WAD (g ml"^-1, ")", sep="")), side=2, line=2.55, cex=1.6)

  #DIFF IN EMPIRICAL MEAN (LABELED - UNLABELED) VS ROW (ITERATION):
  YLIM.diff.emp.mean <- c(min(putative.nongrower.metrics$corr.free.mean.Y - putative.nongrower.metrics$mean.X), max(putative.nongrower.metrics$corr.free.mean.Y - putative.nongrower.metrics$mean.X))
  YLIM.diff.norm.mean <- c(min(putative.nongrower.metrics$corr.free.norm.mean.Y - putative.nongrower.metrics$norm.mean.X), max(putative.nongrower.metrics$corr.free.norm.mean.Y - putative.nongrower.metrics$norm.mean.X))
  YLIM <- c(min(0, YLIM.diff.emp.mean, YLIM.diff.norm.mean), max(0, YLIM.diff.emp.mean, YLIM.diff.norm.mean))
  plot(x=0, y=1, type="n", xlim=c(1, dim(putative.nongrower.metrics)[1]), ylim=YLIM, xlab="", ylab="", xaxt="n", yaxt="n")
  points(x=1:dim(putative.nongrower.metrics)[1], y=putative.nongrower.metrics$corr.free.mean.Y-putative.nongrower.metrics$mean.X, pch=21, col="black", bg="skyblue")
  points(x=as.numeric(row.names(putative.nongrower.metrics))[highlight], y=putative.nongrower.metrics$corr.free.mean.Y[highlight] - putative.nongrower.metrics$mean.X[highlight], pch=21, col="black", bg=highlight.col)
  abline(h=0)
  par(mgp=c(3,0.60,0))			#set spacing for the x-axis labels (the second value)
  axis(side=1, tck=-0.010, las=1, cex.axis=1.25)
  mtext("Iteration", side=1, line=3.0, cex=1.6)
  par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
  axis(side=2, tck=-0.010, cex.axis=1.25)
  mtext(expression(paste("Difference in empirical mean WAD (g ml"^-1, ")", sep="")), side=2, line=2.55, cex=1.6)
  #difference in empirical means are essentially all zero if intercept correction is used

  #DIFF IN EMPIRICAL SD (LABELED - UNLABELED) VS ROW (ITERATION):
  YLIM.diff.emp.sd <- c(min(putative.nongrower.metrics$corr.free.sd.Y - putative.nongrower.metrics$sd.X), max(putative.nongrower.metrics$corr.free.sd.Y - putative.nongrower.metrics$sd.X))
  YLIM.diff.norm.sd <- c(min(putative.nongrower.metrics$corr.free.norm.sd.Y - putative.nongrower.metrics$norm.sd.X), max(putative.nongrower.metrics$corr.free.norm.sd.Y - putative.nongrower.metrics$norm.sd.X))
  YLIM <- c(min(0, YLIM.diff.emp.sd, YLIM.diff.norm.sd), max(0, YLIM.diff.emp.sd, YLIM.diff.norm.sd))
  plot(x=0, y=1, type="n", xlim=c(1, dim(putative.nongrower.metrics)[1]), ylim=YLIM, xlab="", ylab="", xaxt="n", yaxt="n")
  points(x=1:dim(putative.nongrower.metrics)[1], y=putative.nongrower.metrics$corr.free.sd.Y-putative.nongrower.metrics$sd.X, pch=21, col="black", bg="skyblue")
  points(x=as.numeric(row.names(putative.nongrower.metrics))[highlight], y=putative.nongrower.metrics$corr.free.sd.Y[highlight] - putative.nongrower.metrics$sd.X[highlight], pch=21, col="black", bg=highlight.col)
  abline(h=0)
  par(mgp=c(3,0.60,0))			#set spacing for the x-axis labels (the second value)
  axis(side=1, tck=-0.010, las=1, cex.axis=1.25)
  mtext("Iteration", side=1, line=3.0, cex=1.6)
  par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
  axis(side=2, tck=-0.010, cex.axis=1.25)
  mtext(expression(paste("Difference in empirical SD of WAD (g ml"^-1, ")", sep="")), side=2, line=2.55, cex=1.6)

  #NORMAL-FITTED SD VS MEAN (UNCORRECTED X = UNLABELED; CORRECTED FREE Y = LABELED):
  XLIM.emp.mean <- c(min(c(putative.nongrower.metrics$mean.X, putative.nongrower.metrics$corr.free.mean.Y)), max(c(putative.nongrower.metrics$mean.X, putative.nongrower.metrics$corr.free.mean.Y)))
  XLIM.norm.mean <- c(min(c(putative.nongrower.metrics$norm.mean.X, putative.nongrower.metrics$corr.free.norm.mean.Y)), max(c(putative.nongrower.metrics$norm.mean.X, putative.nongrower.metrics$corr.free.norm.mean.Y)))
  XLIM <- c(min(XLIM.emp.mean, XLIM.norm.mean), max(XLIM.emp.mean, XLIM.norm.mean))
  YLIM.emp.sd <- c(min(c(putative.nongrower.metrics$sd.X, putative.nongrower.metrics$corr.free.sd.Y)), max(c(putative.nongrower.metrics$sd.X, putative.nongrower.metrics$corr.free.sd.Y)))
  YLIM.norm.sd <- c(min(c(putative.nongrower.metrics$norm.sd.X, putative.nongrower.metrics$corr.free.norm.sd.Y)), max(c(putative.nongrower.metrics$norm.sd.X, putative.nongrower.metrics$corr.free.norm.sd.Y)))
  YLIM <- c(min(YLIM.emp.sd, YLIM.norm.sd), max(YLIM.emp.sd, YLIM.norm.sd))
  plot(x=0, y=0, type="n", xlim=XLIM, ylim=YLIM, xlab="", ylab="", xaxt="n", yaxt="n")
  points(x=putative.nongrower.metrics$norm.mean.X, y=putative.nongrower.metrics$norm.sd.X, pch=21, col="black", bg="white")
  arrows(x0=putative.nongrower.metrics$norm.mean.X[1:(dim(putative.nongrower.metrics)[1]-1)], y0=putative.nongrower.metrics$norm.sd.X[1:(dim(putative.nongrower.metrics)[1]-1)], x1=putative.nongrower.metrics$norm.mean.X[2:dim(putative.nongrower.metrics)[1]], y1=putative.nongrower.metrics$norm.sd.X[2:dim(putative.nongrower.metrics)[1]], lwd=1, length=0.05, angle=30, code=2, col="black")
  points(x=putative.nongrower.metrics$corr.free.norm.mean.Y, y=putative.nongrower.metrics$corr.free.norm.sd.Y, pch=21, col="black", bg="green")
  points(x=putative.nongrower.metrics$corr.free.norm.mean.Y[highlight], y=putative.nongrower.metrics$corr.free.norm.sd.Y[highlight], pch=21, col="black", bg=highlight.col)
  arrows(x0=putative.nongrower.metrics$corr.free.norm.mean.Y[1:(dim(putative.nongrower.metrics)[1]-1)], y0=putative.nongrower.metrics$corr.free.norm.sd.Y[1:(dim(putative.nongrower.metrics)[1]-1)], x1=putative.nongrower.metrics$corr.free.norm.mean.Y[2:dim(putative.nongrower.metrics)[1]], y1=putative.nongrower.metrics$corr.free.norm.sd.Y[2:dim(putative.nongrower.metrics)[1]], lwd=1, length=0.05, angle=30, code=2, col="blue")
  par(mgp=c(3,0.60,0))			#set spacing for the x-axis labels (the second value)
  axis(side=1, tck=-0.010, las=1, cex.axis=1.25)
  mtext(expression(paste("Normal-fitted mean WAD (g ml"^-1, ")", sep="")), side=1, line=3.4, padj=0, cex=1.6)  #padj needed to adjust vertical alignment when string contains superscript characters
  par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
  axis(side=2, tck=-0.010, cex.axis=1.25)
  mtext(expression(paste("Normal-fitted SD of WAD (g ml"^-1, ")", sep="")), side=2, line=2.55, cex=1.6)
  mtext("Normal-fitted parameters", side=3, line=0.4, cex=2.0)

  #NORMAL-FITTED MEANS VS ROW (ITERATION):
  YLIM.emp.mean <- c(min(c(putative.nongrower.metrics$mean.X, putative.nongrower.metrics$corr.free.mean.Y)), max(c(putative.nongrower.metrics$mean.X, putative.nongrower.metrics$corr.free.mean.Y)))
  YLIM.norm.mean <- c(min(c(putative.nongrower.metrics$norm.mean.X, putative.nongrower.metrics$corr.free.norm.mean.Y)), max(c(putative.nongrower.metrics$norm.mean.X, putative.nongrower.metrics$corr.free.norm.mean.Y)))
  YLIM <- c(min(YLIM.emp.mean, YLIM.norm.mean), max(YLIM.emp.mean, YLIM.norm.mean))
  plot(x=0, y=1, type="n", xlim=c(1, dim(putative.nongrower.metrics)[1]), ylim=YLIM, xlab="", ylab="", xaxt="n", yaxt="n")
  points(x=1:dim(putative.nongrower.metrics)[1], y=putative.nongrower.metrics$norm.mean.X, pch=21, col="black", bg="white")
  points(x=1:dim(putative.nongrower.metrics)[1], y=putative.nongrower.metrics$corr.free.norm.mean.Y, pch=21, col="black", bg="green")
  par(mgp=c(3,0.60,0))			#set spacing for the x-axis labels (the second value)
  axis(side=1, tck=-0.010, las=1, cex.axis=1.25)
  mtext("Iteration", side=1, line=3.0, cex=1.6)
  par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
  axis(side=2, tck=-0.010, cex.axis=1.25)
  mtext(expression(paste("Normal-fitted mean WAD (g ml"^-1, ")", sep="")), side=2, line=2.55, cex=1.6)
  # norm.mean.X is the same as corr.free.norm.mean.Y if norm-fit correction is used

  #NORMAL-FITTED SD VS ROW (ITERATION):
  YLIM.emp.sd <- c(min(c(putative.nongrower.metrics$sd.X, putative.nongrower.metrics$corr.free.sd.Y)), max(c(putative.nongrower.metrics$sd.X, putative.nongrower.metrics$corr.free.sd.Y)))
  YLIM.norm.sd <- c(min(c(putative.nongrower.metrics$norm.sd.X, putative.nongrower.metrics$corr.free.norm.sd.Y)), max(c(putative.nongrower.metrics$norm.sd.X, putative.nongrower.metrics$corr.free.norm.sd.Y)))
  YLIM <- c(min(YLIM.emp.sd, YLIM.norm.sd), max(YLIM.emp.sd, YLIM.norm.sd))
  plot(x=0, y=1, type="n", xlim=c(1, dim(putative.nongrower.metrics)[1]), ylim=YLIM, xlab="", ylab="", xaxt="n", yaxt="n")
  points(x=1:dim(putative.nongrower.metrics)[1], y=putative.nongrower.metrics$norm.sd.X, pch=21, col="black", bg="white")
  points(x=1:dim(putative.nongrower.metrics)[1], y=putative.nongrower.metrics$corr.free.norm.sd.Y, pch=21, col="black", bg="green")
  par(mgp=c(3,0.60,0))			#set spacing for the x-axis labels (the second value)
  axis(side=1, tck=-0.010, las=1, cex.axis=1.25)
  mtext("Iteration", side=1, line=3.0, cex=1.6)
  par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
  axis(side=2, tck=-0.010, cex.axis=1.25)
  mtext(expression(paste("Normal-fitted SD of WAD (g ml"^-1, ")", sep="")), side=2, line=2.55, cex=1.6)

  #DIFF IN NORMAL-FITTED MEAN (LABELED - UNLABELED) VS ROW (ITERATION):
  YLIM.diff.emp.mean <- c(min(putative.nongrower.metrics$corr.free.mean.Y - putative.nongrower.metrics$mean.X), max(putative.nongrower.metrics$corr.free.mean.Y - putative.nongrower.metrics$mean.X))
  YLIM.diff.norm.mean <- c(min(putative.nongrower.metrics$corr.free.norm.mean.Y - putative.nongrower.metrics$norm.mean.X), max(putative.nongrower.metrics$corr.free.norm.mean.Y - putative.nongrower.metrics$norm.mean.X))
  YLIM <- c(min(0, YLIM.diff.emp.mean, YLIM.diff.norm.mean), max(0, YLIM.diff.emp.mean, YLIM.diff.norm.mean))
  plot(x=0, y=1, type="n", xlim=c(1, dim(putative.nongrower.metrics)[1]), ylim=YLIM, xlab="", ylab="", xaxt="n", yaxt="n")
  points(x=1:dim(putative.nongrower.metrics)[1], y=putative.nongrower.metrics$corr.free.norm.mean.Y-putative.nongrower.metrics$norm.mean.X, pch=21, col="black", bg="skyblue")
  points(x=as.numeric(row.names(putative.nongrower.metrics))[highlight], y=putative.nongrower.metrics$corr.free.norm.mean.Y[highlight] - putative.nongrower.metrics$norm.mean.X[highlight], pch=21, col="black", bg=highlight.col)
  abline(h=0)
  par(mgp=c(3,0.60,0))			#set spacing for the x-axis labels (the second value)
  axis(side=1, tck=-0.010, las=1, cex.axis=1.25)
  mtext("Iteration", side=1, line=3.0, cex=1.6)
  par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
  axis(side=2, tck=-0.010, cex.axis=1.25)
  mtext(expression(paste("Difference in normal-fitted mean WAD (g ml"^-1, ")", sep="")), side=2, line=2.55, cex=1.6)
  #difference in normal-fitted means are essentially all zero if normal-fit correction is used

  #DIFF IN NORMAL-FITTED SD (LABELED - UNLABELED) VS ROW (ITERATION):
  YLIM.diff.emp.sd <- c(min(putative.nongrower.metrics$corr.free.sd.Y - putative.nongrower.metrics$sd.X), max(putative.nongrower.metrics$corr.free.sd.Y - putative.nongrower.metrics$sd.X))
  YLIM.diff.norm.sd <- c(min(putative.nongrower.metrics$corr.free.norm.sd.Y - putative.nongrower.metrics$norm.sd.X), max(putative.nongrower.metrics$corr.free.norm.sd.Y - putative.nongrower.metrics$norm.sd.X))
  YLIM <- c(min(0, YLIM.diff.emp.sd, YLIM.diff.norm.sd), max(0, YLIM.diff.emp.sd, YLIM.diff.norm.sd))
  plot(x=0, y=1, type="n", xlim=c(1, dim(putative.nongrower.metrics)[1]), ylim=YLIM, xlab="", ylab="", xaxt="n", yaxt="n")
  points(x=1:dim(putative.nongrower.metrics)[1], y=putative.nongrower.metrics$corr.free.norm.sd.Y-putative.nongrower.metrics$norm.sd.X, pch=21, col="black", bg="skyblue")
  points(x=as.numeric(row.names(putative.nongrower.metrics))[highlight], y=putative.nongrower.metrics$corr.free.norm.sd.Y[highlight] - putative.nongrower.metrics$norm.sd.X[highlight], pch=21, col="black", bg=highlight.col)
  abline(h=0)
  par(mgp=c(3,0.60,0))			#set spacing for the x-axis labels (the second value)
  axis(side=1, tck=-0.010, las=1, cex.axis=1.25)
  mtext("Iteration", side=1, line=3.0, cex=1.6)
  par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
  axis(side=2, tck=-0.010, cex.axis=1.25)
  mtext(expression(paste("Difference in normal-fitted SD of WAD (g ml"^-1, ")", sep="")), side=2, line=2.55, cex=1.6)

  dev.off()   #close the pdf device
}



