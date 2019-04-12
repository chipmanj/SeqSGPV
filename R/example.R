devtools::install_github("chipmanj/sgpvAM")
library(sgpvAM)

# Simulate AM trial
# -- n set to 2000 to motivate location shift example
#    Will be updated so needn't be so large.
am <-  sgpvAM(nreps = 100,
              maxAlertSteps = 100, lookSteps = 1, waitWidths = seq(0.15, 0.6, by = 0.05),
              dataGeneration = rnorm,   dataGenArgs = list(n=2000),
              effectGeneration = 0,
              modelFit = lmCI,
              pointNull = 0, deltaL2 = -0.5, deltaL1=-0.2, deltaG1=0.2, deltaG2=0.5,
              monitoringIntervalLevel=0.05)


# Names of results
names(am)
names(am[["mcmcEndOfStudy"]])
names(am[["mcmcEndOfStudy"]][["width_0.35"]])


# Example of end of study averages for particular wait time
head(am[["mcmcEndOfStudy"]][["width_0.35"]][["mcmcEndOfStudyAve"]])


# Example of ecdf of sample sizes
plot(am[["mcmcEndOfStudy"]][["width_0.35"]][["mcmcEndOfStudyEcdfSize"]][["alertK_50"]],
     xlab="sample size", main="ECDF of sample size with treatment effect of 0",las=1)



# Example of location shift to get power curve
amShifted <- locationShift(am,shiftedThetas = seq(-0.5,0.5,by=0.2))

names(amShifted)

plotPower <- function(am, amShifted, waitWidth, alertK){
  plot(x=0,y=0,xlim=c(-0.5,0.5),ylim=c(0,1),type="n",las=1,xlab="",ylab="",
       main=paste0("P( reject PN )\n wait time ", waitWidth, " and requiring ", alertK, " affirmation steps"))

  o     <- am[["mcmcEndOfStudy"]][[paste0("width_",waitWidth)]][["mcmcEndOfStudyAve"]]
  rejPN <- o[o[,"alertK"]==alertK,"rejPN"]
  x     <- o[o[,"alertK"]==alertK,"theta"]
  points(x=x,y=rejPN)

  for(i in names(amShifted)){
    o     <- amShifted[[i]][["mcmcEndOfStudy"]][[paste0("width_",waitWidth)]][["mcmcEndOfStudyAve"]]
    rejPN <- o[o[,"alertK"]==alertK,"rejPN"]
    x     <- as.numeric(strsplit(i,split="_")[[1]][2])
    points(x=x,y=rejPN)
  }

}

plotPower(am,amShifted, waitWidth = 0.35, alertK = 50)




