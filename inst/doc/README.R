## ---- echo=FALSE---------------------------------------------------------
knitr::opts_chunk$set(fig.width=7, fig.height=5) 

## ------------------------------------------------------------------------
# devtools::install_github("chipmanj/sgpvAM")
library(sgpvAM)

system.time(am <-  sgpvAM(nreps            = 1000,
                          maxAlertSteps    = 100,       lookSteps = 10,  kSteps = 10,
                          waitWidths       = seq(0.2, 0.35, length.out = 5),
                          dataGeneration   = rnorm,   dataGenArgs = list(n=200),
                          effectGeneration = 0,
                          modelFit         = lmCI,
                          pointNull = 0, deltaL2 = NA, deltaL1=NA, deltaG1=0.2, deltaG2=0.5,
                          monitoringIntervalLevel = 0.05,
                          printProgress = TRUE,
                          maxN = 200, lagOutcomeN = 50, 
                          cores = detectCores()))

system.time(amShifted <- locationShift(am, shiftedThetas = seq(-0.5, 1, by = 0.025)))

# Unrestricted Sample Size
# Explore wait width (first with alert k specified at 50)
plot(amShifted, "n",                 alertK = 20,      xlim=c(-0.5, 1),   ylim=c(0,500))
plot(amShifted, "rejPN",             alertK = 20,      xlim=c(-0.5, 1),   ylim=c(0,1))
plot(amShifted, "bias",              alertK = 20,      xlim=c(-0.5, 1),   ylim=c(-0.05,0.05))
plot(amShifted, "mse",               alertK = 20,      xlim=c(-0.5, 1),   ylim=c(0,0.04))
plot(amShifted, "cover",             alertK = 20,      xlim=c(-0.5, 1),   ylim=c(0.75,1))
plot(amShifted, "stopNotImpactful",  alertK = 20,      xlim=c(-0.5, 1),   ylim=c(0,1))
plot(amShifted, "stopNotTrivial",    alertK = 20,      xlim=c(-0.5, 1),   ylim=c(0,1))

# Explore number of steps in affirmation step
plot(amShifted, "n",                 waitWidth = 0.275, alertK = c(0,20,50,100), xlim=c(-1, 1), ylim=c(0,500))
plot(amShifted, "bias",              waitWidth = 0.275, alertK = c(0,20,50,100), xlim=c(-1, 1), ylim=c(-0.05,0.05))
plot(amShifted, "mse",               waitWidth = 0.275, alertK = c(0,20,50,100), xlim=c(-1, 1), ylim=c(0,0.04))
plot(amShifted, "cover",             waitWidth = 0.275, alertK = c(0,20,50,100), xlim=c(-1, 1), ylim=c(0.90,1))
plot(amShifted, "stopNotImpactful",  waitWidth = 0.275, alertK = c(0,20,50,100), xlim=c(-1, 1), ylim=c(0,1))
plot(amShifted, "stopNotTrivial",    waitWidth = 0.275, alertK = c(0,20,50,100), xlim=c(-1, 1), ylim=c(0,1))


# Unrestricted Sample Size with lag of remaining outcomes
plot(amShifted, "n",                 alertK = 20,     xlim=c(-0.5,1),   ylim=c(0,500), sizeRestrictions = "lag")
plot(amShifted, "stopInconclusive",  alertK = 20,     xlim=c(-0.5,1),   ylim=c(0,1),   sizeRestrictions = "lag")
plot(amShifted, "stopInconclusive",  waitWidth = 0.275, alertK = c(0, 20, 50, 100), xlim=c(-0.5,1), ylim=c(0,1),   sizeRestrictions = "lag")

# Max N with immediate outcomes
plot(amShifted, "n",                 alertK = 20,     xlim=c(-0.5,1),   ylim=c(0,250), sizeRestrictions = "maxN")
plot(amShifted, "stopInconclusive",  alertK = 20,     xlim=c(-0.5,1),   ylim=c(0,1),   sizeRestrictions = "maxN")


summary(amShifted, alertK = 10, waitTime = 0.275, treatEffect = 0)
# summary(amShifted)
# summary(am)


## ------------------------------------------------------------------------
system.time(am2 <-  sgpvAM(nreps           = 1000,
                           maxAlertSteps   = 100,       
                           lookSteps = 1, kSteps=10, waitWidths = seq(0.4,0.7,length.out=5),
                           dataGeneration   = rnorm,   dataGenArgs = list(n=200, sd=2),
                           effectGeneration = 0,
                           modelFit         = lmCI,
                           pointNull = 0, deltaL2 = -1, deltaL1=-0.2, deltaG1=0.2, deltaG2=1,
                           monitoringIntervalLevel = 0.05,
                           lagOutcomeN = 100,
                           cores=detectCores()))

system.time(amShifted2 <- locationShift(am2, shiftedThetas = seq(-2,2,by=0.4)))



# Unrestricted Sample Size
# Explore wait width (first with alert k specified at 50)
plot(amShifted2, "n",                 alertK = 50,      xlim=c(-2,2),   ylim=c(0,500))
plot(amShifted2, "rejPN",             alertK = 50,      xlim=c(-2,2),   ylim=c(0,1))
plot(amShifted2, "bias",              alertK = 50,      xlim=c(-2,2),   ylim=c(-0.1,0.1))
plot(amShifted2, "mse",               alertK = 50,      xlim=c(-2,2),   ylim=c(0,0.15))
plot(amShifted2, "cover",             alertK = 50,      xlim=c(-2,2),   ylim=c(0.85,1))
plot(amShifted2, "stopNotImpactful",  alertK = 50,      xlim=c(-2,2),   ylim=c(0,1))
plot(amShifted2, "stopNotTrivial",    alertK = 50,      xlim=c(-2,2),   ylim=c(0,1))

# Explore number of steps in affirmation s
plot(amShifted2, "n",                waitWidth = 0.7,  alertK = c(0, 20, 50, 100), xlim=c(-2,2),   ylim=c(0,500))
plot(amShifted2, "bias",             waitWidth = 0.7,  alertK = c(0, 20, 50, 100), xlim=c(-2,2),   ylim=c(-0.1,0.1))
plot(amShifted2, "mse",              waitWidth = 0.7,  alertK = c(0, 20, 50, 100), xlim=c(-2,2),   ylim=c(0,0.20))
plot(amShifted2, "cover",            waitWidth = 0.7,  alertK = c(0, 20, 50, 100), xlim=c(-2,2),   ylim=c(0.85,1))
plot(amShifted2, "stopNotImpactful", waitWidth = 0.7,  alertK = c(0, 20, 50, 100), xlim=c(-2,2),   ylim=c(0,1))  
plot(amShifted2, "stopNotTrivial",   waitWidth = 0.7,  alertK = c(0, 20, 50, 100), xlim=c(-2,2),   ylim=c(0,1))


# Unrestricted Sample Size with lag of remaining outcomes
plot(amShifted2, "stopInconclusive",  alertK = 50,     xlim=c(-2,1),   ylim=c(0,0.4),  sizeRestrictions = "lag")
plot(amShifted2, "stopInconclusive",  waitWidth = 0.7, alertK = c(0, 20, 50, 100), xlim=c(-2,1),   ylim=c(0,.4),  sizeRestrictions = "lag")


summary(amShifted2, alertK = 50, waitTime = 0.7, treatEffect = 0)
# summary(amShifted2)
# summary(am2)


# Explore ECDF of Bias and Sample Size
ooo <- ecdf.sgpvAM(am = amShifted2, stat = "Bias",alertK = 50,treatEffect = 0,xlim = c(-1,2))
ooo <- ecdf.sgpvAM(am = amShifted2, stat = "Bias",waitWidth = 0.55,alertK = 50,treatEffect = 0,xlim = c(-1,2))
ooo <- ecdf.sgpvAM(am = amShifted2, stat = "Size",alertK = 50,treatEffect = 0,xlim = c(0,400))

# Explore quantiles of sample size for fully specified study design
# See probability of stopping by certain sample size under fully specified design
ooo <- ecdf.sgpvAM(am = amShifted2, stat = "Size",waitWidth = 0.7,alertK = 50,treatEffect = 0,xlim = c(0,400))
quantile(ooo, probs=c(0,0.05,0.1,0.2,0.5,0.8,0.9,0.95,1))
ooo(250)




