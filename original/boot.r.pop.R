
boot.r.pop <- function(T0, Tt, M.soil, days, vars=c("copies", "tube", "trt.code", "g.soil"), growth.model="exponential", vol=c(100, 50), CI=0.90, draws=1000){
  group1 <- boot.TUBE.pop(X=T0, M.soil=M.soil, vars=vars[c(1:2,4)], vol=vol[1], CI=CI, draws=draws)
  group2 <- boot.TUBE.pop(X=Tt, M.soil=M.soil, vars=vars[c(1:2,4)], vol=vol[2], CI=CI, draws=draws)

  r.out <- r.calc.pop(T0=T0, Tt=Tt, boot.out.T0=group1, boot.out.Tt=group2, days=days, var=vars[3], growth.model=growth.model, CI=CI)
  r.out
}
