#Define a function to plot graphs of relevant variables across a range of values for U, given the appropriate data for a soecific taxon:

U.sens.plot.func <- function(LIST, taxonID, U, growth.model="exponential", selected.Us=seq(0.1,1,length.out=10)){
  U.sens <- LIST[[1]]
  U.sens.MW <- LIST[[2]]
  layout( rbind(c(1,1,1,1,1,2,2,2,2,2),
                c(3,3,4,4,5,5,6,6,7,7),
                c(8,8,9,9,10,10,11,11,12,12),
                c(13,13,13,13,13,14,14,14,14,14)),
          widths=lcm(2.54*rep(0.85, 10)),
          heights=lcm(2.54*c(4,1.5,1.5,4)),
          respect=T)
  # layout.show(n=14)
  inds <- round(U.sens.MW$U, 5) %in% round(selected.Us, 5)
  #Plot 1: growth rate sensitivity to U
  par(mai=c(0.44,0.63,0.2,0.05), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
  y.values <- c(U.sens$r.net.obs, U.sens$r.gross.obs, U.sens$d.obs, U.sens$r.net.boot.median, U.sens$r.net.boot.CI.L, U.sens$r.net.boot.CI.U, U.sens$r.gross.boot.median, U.sens$r.gross.boot.CI.L, U.sens$r.gross.boot.CI.U, U.sens$d.boot.median, U.sens$d.boot.CI.L, U.sens$d.boot.CI.U)
  y.values <- y.values[y.values != -Inf & y.values != Inf & !is.na(y.values)]
  y.min <- min(y.values)
  y.max <- max(y.values)
  plot(y=U.sens$r.net.obs, x=U.sens$U, type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(0, 1), ylim=c(y.min, y.max))
    abline(v=U, lty=2, col="black")
    abline(h=0, col="black")
    arrows(x0=U.sens$U, y0=U.sens$d.boot.CI.L, x1=U.sens$U, y1=U.sens$d.boot.CI.U, length=0, angle=90, code=3, col="dodgerblue", lwd=1.5)
    # arrows(x0=U.sens$U[inds], y0=U.sens$d.boot.CI.L[inds], x1=U.sens$U[inds], y1=U.sens$d.boot.CI.U[inds], length=0, angle=90, code=3, col="purple", lwd=1.5)
    arrows(x0=U.sens$U, y0=U.sens$r.net.boot.CI.L, x1=U.sens$U, y1=U.sens$r.net.boot.CI.U, length=0, angle=90, code=3, col="tomato", lwd=0.5)
    arrows(x0=U.sens$U, y0=U.sens$r.gross.boot.CI.L, x1=U.sens$U, y1=U.sens$r.gross.boot.CI.U, length=0, angle=90, code=3, col="limegreen", lwd=1)
    points(y=U.sens$r.net.boot.median, x=U.sens$U, pch=21, cex=0.6, col="black", bg="tomato")
    lines(y=U.sens$r.net.obs, x=U.sens$U, lwd=2.5, col="black")
    lines(y=U.sens$r.net.obs, x=U.sens$U, lwd=2, col="yellow")
    points(y=U.sens$r.gross.boot.median, x=U.sens$U, pch=21, cex=0.6, col="black", bg="limegreen")
    lines(y=U.sens$r.gross.obs, x=U.sens$U, lwd=2.5, col="black")
    lines(y=U.sens$r.gross.obs, x=U.sens$U, lwd=2, col="limegreen")
    points(y=U.sens$d.boot.median, x=U.sens$U, pch=21, cex=0.6, col="black", bg="dodgerblue")
    # points(y=U.sens$d.boot.median[inds], x=U.sens$U[inds], pch=21, cex=0.6, col="black", bg="purple")
    lines(y=U.sens$d.obs, x=U.sens$U, lwd=2.5, col="black")
    lines(y=U.sens$d.obs, x=U.sens$U, lwd=2, col="dodgerblue")
    par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
    axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
    mtext("U", side=1, line=1.35, cex=0.75)
    par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
    axis(side=2, tck=-0.015, las=1, cex.axis=0.6)
    if (growth.model == "exponential"){
      mtext(expression(paste("day"^-1, sep="")), side=2, line=2.4, cex=0.75)
    }
    if (growth.model == "linear"){
      mtext(expression(paste("copies g soil"^-1, " day"^-1, sep="")), side=2, line=2.4, cex=0.75)
    }
    par(xpd=NA)
    legend(x=1, y=y.max+(0.09*(y.max-y.min)), legend=c("net r (obs)", "net r (boot)", "gross r (obs)", "gross r (boot)", "death (obs)", "death (boot)"), col=c("black","white","black","white","black","white"), pt.bg=c("yellow","white","limegreen","white","dodgerblue","white"), pch=c(NA,21,NA,21,NA,21), pt.cex=0, lwd=c(2.5,0,2.5,0,2.5,0), bty="o", cex=0.6, xjust=1, yjust=1)
    legend(x=1, y=y.max+(0.09*(y.max-y.min)), legend=c("net r (obs)", "net r (boot)", "gross r (obs)", "gross r (boot)", "death (obs)", "death (boot)"), col=c("yellow","white","limegreen","white","dodgerblue","white"), pt.bg=c("yellow","white","limegreen","white","dodgerblue","white"), pch=c(NA,21,NA,21,NA,21), pt.cex=0, lwd=c(2,0,2,0,2,0), bty="o", cex=0.6, xjust=1, yjust=1)
    legend(x=1, y=y.max+(0.09*(y.max-y.min)), legend=c("net r (obs)", "net r (boot)", "gross r (obs)", "gross r (boot)", "death (obs)", "death (boot)"), col="black", pt.bg=c("yellow","tomato","limegreen","limegreen","dodgerblue","dodgerblue"), pch=c(NA,21,NA,21,NA,21), pt.cex=0.6,  bty="n", cex=0.6, xjust=1, yjust=1)
    par(xpd=FALSE)
    mtext(paste("Taxon ", taxonID, sep=""), side=3, line=0, cex=1)

  #Plot 2: abundance sensitivity to U
  par(mai=c(0.44,0.63,0.2,0.05), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
  y.values <- c(U.sens$tot.copies.g.soil.t0.obs, U.sens$tot.copies.g.soil.tT.obs, U.sens$light.copies.g.soil.tT.obs, U.sens$tot.copies.g.soil.t0.boot.median, U.sens$tot.copies.g.soil.t0.boot.CI.L, U.sens$tot.copies.g.soil.t0.boot.CI.U, U.sens$tot.copies.g.soil.tT.boot.median, U.sens$tot.copies.g.soil.tT.boot.CI.L, U.sens$tot.copies.g.soil.tT.boot.CI.U, U.sens$light.copies.g.soil.tT.boot.median, U.sens$light.copies.g.soil.tT.boot.CI.L, U.sens$light.copies.g.soil.tT.boot.CI.U)
  y.values <- y.values[y.values != -Inf & y.values != Inf & !is.na(y.values)]
  y.min <- min(y.values)
  # y.max <- max(y.values)
  y.max <- max(2*unique(U.sens$tot.copies.g.soil.tT.obs), U.sens$tot.copies.g.soil.t0.boot.CI.U)
  plot(y=U.sens$tot.copies.g.soil.tT.obs, x=U.sens$U, type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(0, 1), ylim=c(y.min, y.max))
    abline(v=U, lty=2, col="black")
    abline(h=0, col="black")
    arrows(x0=U.sens$U, y0=U.sens$light.copies.g.soil.tT.boot.CI.L, x1=U.sens$U, y1=U.sens$light.copies.g.soil.tT.boot.CI.U, length=0, angle=90, code=3, col="dodgerblue", lwd=1.5)
    # arrows(x0=U.sens$U[inds], y0=U.sens$light.copies.g.soil.tT.boot.CI.L[inds], x1=U.sens$U[inds], y1=U.sens$light.copies.g.soil.tT.boot.CI.U[inds], length=0, angle=90, code=3, col="purple", lwd=1.5)
    arrows(x0=U.sens$U, y0=U.sens$tot.copies.g.soil.tT.boot.CI.L, x1=U.sens$U, y1=U.sens$tot.copies.g.soil.tT.boot.CI.U, length=0, angle=90, code=3, col="tomato", lwd=0.5)
    arrows(x0=U.sens$U, y0=U.sens$tot.copies.g.soil.t0.boot.CI.L, x1=U.sens$U, y1=U.sens$tot.copies.g.soil.t0.boot.CI.U, length=0, angle=90, code=3, col="limegreen", lwd=1)
    points(y=U.sens$tot.copies.g.soil.tT.boot.median, x=U.sens$U, pch=21, cex=0.6, col="black", bg="tomato")
    lines(y=U.sens$tot.copies.g.soil.tT.obs, x=U.sens$U, lwd=2.5, col="black")
    lines(y=U.sens$tot.copies.g.soil.tT.obs, x=U.sens$U, lwd=2, col="yellow")
    points(y=U.sens$tot.copies.g.soil.t0.boot.median, x=U.sens$U, pch=21, cex=0.6, col="black", bg="limegreen")
    lines(y=U.sens$tot.copies.g.soil.t0.obs, x=U.sens$U, lwd=2.5, col="black")
    lines(y=U.sens$tot.copies.g.soil.t0.obs, x=U.sens$U, lwd=2, col="limegreen")
    points(y=U.sens$light.copies.g.soil.tT.boot.median, x=U.sens$U, pch=21, cex=0.6, col="black", bg="dodgerblue")
    # points(y=U.sens$light.copies.g.soil.tT.boot.median[inds], x=U.sens$U[inds], pch=21, cex=0.6, col="black", bg="purple")
    lines(y=U.sens$light.copies.g.soil.tT.obs, x=U.sens$U, lwd=2.5, col="black")
    lines(y=U.sens$light.copies.g.soil.tT.obs, x=U.sens$U, lwd=2, col="dodgerblue")
    par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
    axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
    mtext("U", side=1, line=1.35, cex=0.75)
    par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
    axis(side=2, tck=-0.015, las=1, cex.axis=0.6)
    par(xpd=NA)
    mtext(expression(paste("Abundance (16S copies g soil"^-1, ")", sep="")), side=2, line=2.9, cex=0.75)
    # legend(x=1, y=y.max+(0.09*(y.max-y.min)), legend=c("total, time 10 (obs)", "total, time 10 (boot)", "total (unlabeled), time 0 (obs)", "total (unlabeled), time 0 (boot)", "unlabeled, time 10 (obs)", "unlabeled, time 10 (boot)"), col="black", pt.bg=c("yellow","tomato","limegreen","limegreen","dodgerblue","dodgerblue"), pch=c(NA,21,NA,21,NA,21), pt.cex=0, lwd=c(2.5,0,2.5,0,2.5,0), bty="o", cex=0.6, xjust=1, yjust=1)
    # legend(x=1, y=y.max+(0.09*(y.max-y.min)), legend=c("total, time 10 (obs)", "total, time 10 (boot)", "total (unlabeled), time 0 (obs)", "total (unlabeled), time 0 (boot)", "unlabeled, time 10 (obs)", "unlabeled, time 10 (boot)"), col=c("yellow","tomato","limegreen","limegreen","dodgerblue","dodgerblue"), pt.bg=c("yellow","tomato","limegreen","limegreen","dodgerblue","dodgerblue"), pch=c(NA,21,NA,21,NA,21), pt.cex=0, lwd=c(2,0,2,0,2,0), bty="o", cex=0.6, xjust=1, yjust=1)
    # legend(x=1, y=y.max+(0.09*(y.max-y.min)), legend=c("total, time 10 (obs)", "total, time 10 (boot)", "total (unlabeled), time 0 (obs)", "total (unlabeled), time 0 (boot)", "unlabeled, time 10 (obs)", "unlabeled, time 10 (boot)"), col="black", pt.bg=c("yellow","tomato","limegreen","limegreen","dodgerblue","dodgerblue"), pch=c(NA,21,NA,21,NA,21), pt.cex=0.6,  bty="n", cex=0.6, xjust=1, yjust=1)
    legend(x=1, y=y.max+(0.09*(y.max-y.min)), legend=c("total, time 10 (obs)", "total, time 10 (boot)", "total (unlabeled), time 0 (obs)", "total (unlabeled), time 0 (boot)", "unlabeled, time 10 (obs)", "unlabeled, time 10 (boot)"), col=c("black","white","black","white","black","white"), pt.bg=c("yellow","white","limegreen","white","dodgerblue","white"), pch=c(NA,21,NA,21,NA,21), pt.cex=0, lwd=c(2.5,0,2.5,0,2.5,0), bty="o", cex=0.6, xjust=1, yjust=1)
    legend(x=1, y=y.max+(0.09*(y.max-y.min)), legend=c("total, time 10 (obs)", "total, time 10 (boot)", "total (unlabeled), time 0 (obs)", "total (unlabeled), time 0 (boot)", "unlabeled, time 10 (obs)", "unlabeled, time 10 (boot)"), col=c("yellow","white","limegreen","white","dodgerblue","white"), pt.bg=c("yellow","white","limegreen","white","dodgerblue","white"), pch=c(NA,21,NA,21,NA,21), pt.cex=0, lwd=c(2,0,2,0,2,0), bty="o", cex=0.6, xjust=1, yjust=1)
    legend(x=1, y=y.max+(0.09*(y.max-y.min)), legend=c("total, time 10 (obs)", "total, time 10 (boot)", "total (unlabeled), time 0 (obs)", "total (unlabeled), time 0 (boot)", "unlabeled, time 10 (obs)", "unlabeled, time 10 (boot)"), col="black", pt.bg=c("yellow","tomato","limegreen","limegreen","dodgerblue","dodgerblue"), pch=c(NA,21,NA,21,NA,21), pt.cex=0.6,  bty="n", cex=0.6, xjust=1, yjust=1)
    par(xpd=FALSE)

  #Plots 3-12: histograms of observed molecular weights in 18O treatment at different values of U
  par(mai=c(0.44,0.63,0.2,0.05), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
  for (i in 1:length(selected.Us)){
    MW.boots <- as.numeric(U.sens.MW[round(U.sens.MW$U, 5) == round(selected.Us[i], 5), 5:1004])
    hist(MW.boots, xlab="", ylab="", xaxt="n", yaxt="n", main="")
    abline(v=U.sens.MW$MW.light[round(U.sens.MW$U, 5) == round(selected.Us[i], 5)], lwd=1, col="blue")
    abline(v=U.sens.MW$MW.heavy[round(U.sens.MW$U, 5) == round(selected.Us[i], 5)], lwd=1, col="red")
    rect(xleft=U.sens.MW$MW.heavy[round(U.sens.MW$U, 5) == round(selected.Us[i], 5)], ybottom=-5000, xright=5000, ytop=5000, border=NA, col=rgb(255, 0, 0, alpha=0.4*255, maxColorValue=255))
    # abline(v=U.sens.MW$MW.obs.labeled[round(U.sens.MW$U, 5) == round(selected.Us[i], 5)], lwd=1, lty=2, col="red")
    arrows(x0=quantile(MW.boots, probs=0.05, na.rm=TRUE), y0=0.125*max(hist(MW.boots, plot=FALSE)$counts), x1=quantile(MW.boots, probs=0.95, na.rm=TRUE), y1=0.125*max(hist(MW.boots, plot=FALSE)$counts), length=0, angle=90, code=3, col="limegreen", lwd=2.5)
    points(x=median(MW.boots, na.rm=TRUE), y=0.125*max(hist(MW.boots, plot=FALSE)$counts), pch=21, cex=0.7, col="black", bg="limegreen")
    MW.boots.remaining <- MW.boots[MW.boots < U.sens.MW$MW.heavy[round(U.sens.MW$U, 5) == round(selected.Us[i], 5)]]
    arrows(x0=quantile(MW.boots.remaining, probs=0.05, na.rm=TRUE), y0=0.25*max(hist(MW.boots, plot=FALSE)$counts), x1=quantile(MW.boots.remaining, probs=0.95, na.rm=TRUE), y1=0.25*max(hist(MW.boots, plot=FALSE)$counts), length=0, angle=90, code=3, col="purple", lwd=2.5)
    points(x=median(MW.boots.remaining, na.rm=TRUE), y=0.25*max(hist(MW.boots, plot=FALSE)$counts), pch=21, cex=0.7, col="black", bg="purple")
    par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
    axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
    mtext(expression(paste("MW (g mol"^-1, ")", sep="")), side=1, line=1.35, cex=0.75)
    par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
    axis(side=2, tck=-0.015, las=1, cex.axis=0.6)
    mtext("Frequency", side=2, line=2.4, cex=0.75)
    mtext(paste("U = ", format(selected.Us, digits=2)[i], sep=""), side=3, line=0, cex=0.75)
  }

  #Plot 13: bias in MW versus U
  par(mai=c(0.44,0.63,0.2,0.05), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
  y.min <- min(as.numeric(U.sens.MW[1,5:1004]), na.rm=TRUE) - U.sens.MW$MW.obs.labeled[1]
  y.max <- max(as.numeric(U.sens.MW[(dim(U.sens.MW)[1]),5:1004]), na.rm=TRUE) - U.sens.MW$MW.obs.labeled[(dim(U.sens.MW)[1])]
  plot(y=0, x=U.sens.MW$U[1], type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(0, 1), ylim=c(y.min, y.max))
  abline(v=U, lty=2, col="black")
  abline(h=0, col="black")
  for (i in 1:length(U.sens.MW$U)){
    MW.boots <- as.numeric(U.sens.MW[U.sens.MW$U == U.sens.MW$U[i], 5:1004])
    MW.boots.remaining <- MW.boots[MW.boots < U.sens.MW$MW.heavy[U.sens.MW$U == U.sens.MW$U[i]]]
    MW.bias <- MW.boots.remaining - U.sens.MW$MW.obs.labeled[U.sens.MW$U == U.sens.MW$U[i]]
    arrows(x0=U.sens.MW$U[i], y0=quantile(MW.bias, probs=0.05, na.rm=TRUE), x1=U.sens.MW$U[i], y1=quantile(MW.bias, probs=0.95, na.rm=TRUE), length=0, angle=90, code=3, col="purple", lwd=1.5)
    points(y=median(MW.bias,na.rm=TRUE), x=U.sens.MW$U[i], pch=21, cex=0.6, col="black", bg="purple")
  }
    par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
    axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
    mtext("U", side=1, line=1.35, cex=0.75)
    par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
    axis(side=2, tck=-0.015, las=1, cex.axis=0.6)
    mtext(expression(paste("Bias in molecular weight (g mol"^-1, ")", sep="")), side=2, line=2.4, cex=0.75)

  #Plot 14: bias in abundance of unlabeled copies at time 10 versus U
  par(mai=c(0.44,0.63,0.2,0.05), cex=1, mex=0.75)				#set the margins (bottom, left, top, right) in inches
  Q.L.10.bias <- U.sens$light.copies.g.soil.tT.boot.median - U.sens$light.copies.g.soil.tT.obs
  if(length(Q.L.10.bias[!is.na(Q.L.10.bias)]) > 0){
    y.min <- min(Q.L.10.bias, na.rm=TRUE)
    y.max <- max(Q.L.10.bias, na.rm=TRUE)
  }
  else {
    y.min <- -5000000
    y.max <- 5000000   
  }
  plot(y=Q.L.10.bias, x=U.sens$U, type="n", bty="l", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(0, 1), ylim=c(y.min, y.max))
    abline(v=U, lty=2, col="black")
    abline(h=0, col="black")
    points(y=Q.L.10.bias, x=U.sens$U, pch=21, cex=0.6, col="black", bg="purple")
    par(mgp=c(3,0,0))				  #set spacing for the x-axis labels (the second value)
    axis(side=1, tck=-0.015, las=1, cex.axis=0.6)
    mtext("U", side=1, line=1.35, cex=0.75)
    par(mgp=c(3,0.45,0))			#set spacing for the y-axis labels (the second value)
    axis(side=2, tck=-0.015, las=1, cex.axis=0.6)
    mtext(expression(paste("Bias in abundance (16S copies g soil"^-1, ")", sep="")), side=2, line=2.4, cex=0.75)
}
