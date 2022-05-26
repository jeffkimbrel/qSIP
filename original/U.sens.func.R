#Define a function to create two dataframes of relevant variables across a range of values for U:

U.sens.func <- function(DATA, DATA.POP, taxonID, time0, tmt1, tmt2, M.soil, iso.compare, days, vars=c("taxon", "density", "copies", "tube", "trt.code", "g.soil"), growth.model="exponential", vol=c(100, 50), CI=0.90, draws=1000, sens.seq=seq(0.01, 1, 0.01)){
  TAXONT1 <- DATA[DATA[,vars[1]]==taxonID & DATA[,vars[5]]==tmt1,]
  TAXONT2 <- DATA[DATA[,vars[1]]==taxonID & DATA[,vars[5]]==tmt2,]
  TAXONref <- DATA[DATA[,vars[1]]==taxonID & DATA[,vars[5]]==tmt1,]
  TAXONT0.tube <- DATA.POP[DATA.POP[,vars[1]]==taxonID & DATA.POP[,vars[5]]==time0,]
  TAXONT1.tube <- DATA.POP[DATA.POP[,vars[1]]==taxonID & DATA.POP[,vars[5]]==tmt1,]
  TAXONT2.tube <- DATA.POP[DATA.POP[,vars[1]]==taxonID & DATA.POP[,vars[5]]==tmt2,]
  TAXONT1.boot.out <- boot.WAD.func(X=TAXONT1, vars=vars[2:4], CI=CI, draws=draws)
  TAXONT2.boot.out <- boot.WAD.func(X=TAXONT2, vars=vars[2:4], CI=CI, draws=draws)
  TAXONref.MW.out <- MW.calc(X=TAXONref, vars=vars[2:4])
  TAXONT0v12.r.out <- boot.r.pop(T0=TAXONT0.tube, Tt=rbind(TAXONT1.tube, TAXONT2.tube), M.soil=M.soil, days=days, vars=vars[3:6], growth.model=growth.model, vol=vol, CI=CI, draws=draws)
  TAXONT0.tube.copies.out <- boot.TUBE.pop(X=TAXONT0.tube, M.soil=M.soil, vars=vars[c(3,4,6)], vol=vol[1], CI=CI, draws=draws)
  TAXONT12.tube.copies.out <- boot.TUBE.pop(X=rbind(TAXONT1.tube, TAXONT2.tube), M.soil=M.soil, vars=vars[c(3,4,6)], vol=vol[2], CI=CI, draws=draws)
  U.sens.MW.TAXON <- data.frame(U=sens.seq)
  U.sens.MW.TAXON$MW.heavy <- rep(NA, length(U.sens.MW.TAXON$U))
  U.sens.MW.TAXON$MW.light <- TAXONref.MW.out$MW
  U.sens.MW.TAXON$MW.obs.labeled <- TAXONref.MW.out$MW * (TAXONT2.boot.out$obs.wad.mean/TAXONT1.boot.out$obs.wad.mean)
  U.sens.MW.TAXON <- data.frame(U.sens.MW.TAXON, matrix(NA, nrow=length(U.sens.MW.TAXON$U), ncol=draws))
  names(U.sens.MW.TAXON)[5:(draws+4)] <- paste("MW.boot.labeled.", 1:draws, sep="")
  U.sens.TAXON <- data.frame(U=U.sens.MW.TAXON$U)
  U.sens.TAXON$tot.copies.g.soil.t0.obs <- mean(TAXONT0.tube.copies.out$obs.N$tot.copies / TAXONT0.tube.copies.out$obs.N$g.soil)
  U.sens.TAXON$tot.copies.g.soil.tT.obs <- mean(TAXONT12.tube.copies.out$obs.N$tot.copies / TAXONT12.tube.copies.out$obs.N$g.soil)
  U.sens.TAXON$light.copies.g.soil.tT.obs <- rep(NA, length(U.sens.TAXON$U))
  U.sens.TAXON$r.net.obs <- TAXONT0v12.r.out$obs.r
  U.sens.TAXON$r.gross.obs <- rep(NA, length(U.sens.TAXON$U))
  U.sens.TAXON$d.obs <- rep(NA, length(U.sens.TAXON$U))
  U.sens.TAXON$tot.copies.g.soil.t0.boot.median <- median(apply(TAXONT0.tube.copies.out$boot.tot.copies / TAXONT0.tube.copies.out$boot.g.soil, 1, mean, na.rm=TRUE), na.rm=TRUE)
  U.sens.TAXON$tot.copies.g.soil.t0.boot.CI.L <- quantile(apply(TAXONT0.tube.copies.out$boot.tot.copies / TAXONT0.tube.copies.out$boot.g.soil, 1, mean, na.rm=TRUE), probs=(1-CI)/2, na.rm=TRUE)
  U.sens.TAXON$tot.copies.g.soil.t0.boot.CI.U <- quantile(apply(TAXONT0.tube.copies.out$boot.tot.copies / TAXONT0.tube.copies.out$boot.g.soil, 1, mean, na.rm=TRUE), probs=1-((1-CI)/2), na.rm=TRUE)
  U.sens.TAXON$tot.copies.g.soil.tT.boot.median <- median(apply(TAXONT12.tube.copies.out$boot.tot.copies / TAXONT12.tube.copies.out$boot.g.soil, 1, mean, na.rm=TRUE), na.rm=TRUE)
  U.sens.TAXON$tot.copies.g.soil.tT.boot.CI.L <- quantile(apply(TAXONT12.tube.copies.out$boot.tot.copies / TAXONT12.tube.copies.out$boot.g.soil, 1, mean, na.rm=TRUE), probs=(1-CI)/2, na.rm=TRUE)
  U.sens.TAXON$tot.copies.g.soil.tT.boot.CI.U <- quantile(apply(TAXONT12.tube.copies.out$boot.tot.copies / TAXONT12.tube.copies.out$boot.g.soil, 1, mean, na.rm=TRUE), probs=1-((1-CI)/2), na.rm=TRUE)
  U.sens.TAXON$light.copies.g.soil.tT.boot.median <- rep(NA, length(U.sens.TAXON$U))
  U.sens.TAXON$light.copies.g.soil.tT.boot.CI.L <- rep(NA, length(U.sens.TAXON$U))
  U.sens.TAXON$light.copies.g.soil.tT.boot.CI.U <- rep(NA, length(U.sens.TAXON$U))
  U.sens.TAXON$r.net.boot.median <- TAXONT0v12.r.out$boot.r.median
  U.sens.TAXON$r.net.boot.CI.L <- TAXONT0v12.r.out$boot.r.CI[1]
  U.sens.TAXON$r.net.boot.CI.U <- TAXONT0v12.r.out$boot.r.CI[2]
  for (i in 1:length(U.sens.TAXON$U)){
    U.sens.MW.TAXON$MW.heavy[i] <- TAXONref.MW.out$MW + (TAXONref.MW.out$Oatoms*2*U.sens.TAXON$U[i])
    U.sens.MW.TAXON[i,5:(draws+4)] <- TAXONref.MW.out$MW * (TAXONT2.boot.out$boot.wads/TAXONT1.boot.out$boot.wads)
    TAXONT1v2.r.out <- boot.diff.r(X.light=TAXONT1, X.heavy=TAXONT2, X.reference=TAXONref, M.soil=M.soil, iso.compare=iso.compare, days=days, vars=vars[2:6], growth.model=growth.model, prop.O.from.water=U.sens.TAXON$U[i], v.frac=vol[2], CI=CI, draws=draws)
    U.sens.TAXON$r.gross.obs[i] <- TAXONT1v2.r.out$obs.r
    U.sens.TAXON$d.obs[i] <- TAXONT0v12.r.out$obs.r - TAXONT1v2.r.out$obs.r
    if (growth.model == "exponential"){
      U.sens.TAXON$light.copies.g.soil.tT.obs[i] <- U.sens.TAXON$tot.copies.g.soil.tT.obs[i] * (1 / exp(TAXONT1v2.r.out$obs.r * days))
      U.sens.TAXON$light.copies.g.soil.tT.boot.median[i] <- median(U.sens.TAXON$tot.copies.g.soil.tT.obs[i] * (1 / exp(TAXONT1v2.r.out$boot.r * days)), na.rm=TRUE)
      U.sens.TAXON$light.copies.g.soil.tT.boot.CI.L[i] <- quantile(U.sens.TAXON$tot.copies.g.soil.tT.obs[i] * (1 / exp(TAXONT1v2.r.out$boot.r * days)), probs=(1-CI)/2, na.rm=TRUE)
      U.sens.TAXON$light.copies.g.soil.tT.boot.CI.U[i] <- quantile(U.sens.TAXON$tot.copies.g.soil.tT.obs[i] * (1 / exp(TAXONT1v2.r.out$boot.r * days)), probs=1-((1-CI)/2), na.rm=TRUE)
    }
    if (growth.model == "linear"){
      U.sens.TAXON$light.copies.g.soil.tT.obs[i] <- U.sens.TAXON$tot.copies.g.soil.tT.obs[i] - (TAXONT1v2.r.out$obs.r * days)
      U.sens.TAXON$light.copies.g.soil.tT.boot.median[i] <- median(U.sens.TAXON$tot.copies.g.soil.tT.obs[i] - (TAXONT1v2.r.out$boot.r * days), na.rm=TRUE)
      U.sens.TAXON$light.copies.g.soil.tT.boot.CI.L[i] <- quantile(U.sens.TAXON$tot.copies.g.soil.tT.obs[i] - (TAXONT1v2.r.out$boot.r * days), probs=(1-CI)/2, na.rm=TRUE)
      U.sens.TAXON$light.copies.g.soil.tT.boot.CI.U[i] <- quantile(U.sens.TAXON$tot.copies.g.soil.tT.obs[i] - (TAXONT1v2.r.out$boot.r * days), probs=1-((1-CI)/2), na.rm=TRUE)
    }
    U.sens.TAXON$r.gross.boot.median[i] <- TAXONT1v2.r.out$boot.r.median
    U.sens.TAXON$r.gross.boot.CI.L[i] <- TAXONT1v2.r.out$boot.r.CI[1]
    U.sens.TAXON$r.gross.boot.CI.U[i] <- TAXONT1v2.r.out$boot.r.CI[2]
    U.sens.TAXON$d.boot.median[i] <- median(TAXONT0v12.r.out$boot.r - TAXONT1v2.r.out$boot.r, na.rm=TRUE)
    U.sens.TAXON$d.boot.CI.L[i] <- quantile(TAXONT0v12.r.out$boot.r - TAXONT1v2.r.out$boot.r, probs=(1-CI)/2, na.rm=TRUE)
    U.sens.TAXON$d.boot.CI.U[i] <- quantile(TAXONT0v12.r.out$boot.r - TAXONT1v2.r.out$boot.r, probs=1-((1-CI)/2), na.rm=TRUE)
  }
  return(list(U.sens=U.sens.TAXON, U.sens.MW=U.sens.MW.TAXON))
}
