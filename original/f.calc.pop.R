#     boot.out.T0: list containing output from boot.TUBE.pop for the 'time 0' treatment (corresponding to T0)
#     boot.out.Tt: list containing output from boot.TUBE.pop for the 'time t' treatment (corresponding to Tt)
#     (this function requires output from boot.TUBE.pop; output from boot.pop won't work)


f.calc.pop <- function(T0, Tt, boot.out.T0, boot.out.Tt, days, var="trt.code", growth.model="exponential", copies.cell=6, pgC.cell=0.1, CI=0.90){

  if (growth.model=="linear"){                                        #if growth.model is "linear", calculate r using a linear growth model
    # Growth rate in units of total copies per g soil per day:
    obs.r <- (mean(boot.out.Tt$obs.N$tot.copies/boot.out.Tt$obs.N$g.soil, na.rm=TRUE) - mean(boot.out.T0$obs.N$tot.copies/boot.out.T0$obs.N$g.soil, na.rm=TRUE))/days
    boot.r <- (apply(boot.out.Tt$boot.tot.copies/boot.out.Tt$boot.g.soil, 1, mean, na.rm=TRUE) - apply(boot.out.T0$boot.tot.copies/boot.out.T0$boot.g.soil, 1, mean, na.rm=TRUE))/days

    # Flux of carbon in pg C per g soil per day:
    obs.f <- obs.r * (1/copies.cell) * pgC.cell
    boot.f <- boot.r * (1/copies.cell) * pgC.cell
  }
  else if (growth.model=="exponential"){                              #if growth.model is "exponential", calculate r using an exponential growth model
    # Growth rate in units of per day:
    obs.r <- log(mean(boot.out.Tt$obs.N$tot.copies/boot.out.Tt$obs.N$g.soil, na.rm=TRUE) / mean(boot.out.T0$obs.N$tot.copies/boot.out.T0$obs.N$g.soil, na.rm=TRUE))/days
    boot.r <- log(apply(boot.out.Tt$boot.tot.copies/boot.out.Tt$boot.g.soil, 1, mean, na.rm=TRUE) / apply(boot.out.T0$boot.tot.copies/boot.out.T0$boot.g.soil, 1, mean, na.rm=TRUE))/10

    obs.copies.soil <- mean(c(boot.out.Tt$obs.N$tot.copies/boot.out.Tt$obs.N$g.soil, boot.out.T0$obs.N$tot.copies/boot.out.T0$obs.N$g.soil), na.rm=T)
    boot.copies.soil <- apply(cbind(boot.out.Tt$boot.tot.copies/boot.out.Tt$boot.g.soil, boot.out.T0$boot.tot.copies/boot.out.T0$boot.g.soil), 1, mean, na.rm=T)

    # Flux of carbon in pg C per g soil per day:
    obs.f <- obs.r * obs.copies.soil * (1/copies.cell) * pgC.cell
    boot.f <- boot.r * boot.copies.soil * (1/copies.cell) * pgC.cell
  }
  else {                                                              #if the growth.model is anything else, print an error message
    stop(paste("Error: unable to calculate r for the specified growth model '", growth.model, "'", sep=""))
  }

  boot.f.CI <- quantile(boot.f, probs=c((1-CI)/2, 1-((1-CI)/2)), na.rm=T)

  message <- comparison.message.pop(T0=T0, Tt=Tt, boot.out.T0=boot.out.T0, boot.out.Tt=boot.out.Tt, var=var)

  # Collect all output into a list:
   return(list(boot.f=boot.f, 
               obs.f=obs.f, 
               boot.f.mean=mean(boot.f, na.rm=T), 
               boot.f.median=median(boot.f, na.rm=T), 
               boot.f.CI=boot.f.CI, 
               message=message)) 
}
