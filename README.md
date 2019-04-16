# sgpvAM
Adaptive Monitoring Design Features using the second generation p-value


# install.packages("devtools")
devtools::install_github("chipmanj/sgpvAM")


# Examples

Beware that running 1000 reps still takes a good amount of time.  Start small with 100.

```r
devtools::install_github("chipmanj/sgpvAM")
library(sgpvAM)
?sgpvAM::sgpvAM


# Simulate AM trial
# Two-sided deltas
# defaults to sd = 1
am <-  sgpvAM(nreps            = 10,
              maxAlertSteps    = 100,       lookSteps = 1, waitWidths = seq(0.15, 0.6, by = 0.05),
              dataGeneration   = rnorm,   dataGenArgs = list(n=800),
              effectGeneration = 0,
              modelFit         = sgpvAM::lmCI,
              pointNull = 0, deltaL2 = -0.5, deltaL1=-0.2, deltaG1=0.2, deltaG2=0.5,
              monitoringIntervalLevel = 0.05,
              maxN = 200, lagOutcomeN = 50)
amShifted <- locationShift(am, shiftedThetas = seq(-1,0,by=0.1))

# Unrestricted Sample Size
# Explore wait width (first with alert k specified at 50)
plot(amShifted, "n",                 alertK = 50,      xlim=c(-1,0),   ylim=c(0,801),    maxVary = 10)
plot(amShifted, "rejPN",             alertK = 50,      xlim=c(-1,0),   ylim=c(0,1),      maxVary = 10)
plot(amShifted, "rejPN",             alertK = 50,      xlim=c(-1,0),   ylim=c(0,1),      maxVary = 5)
plot(amShifted, "bias",              alertK = 50,      xlim=c(-1,0),   ylim=c(-0.1,0.1), maxVary = 5)
plot(amShifted, "mse",               alertK = 50,      xlim=c(-1,0),   ylim=c(0,0.1),    maxVary = 5)
plot(amShifted, "cover",             alertK = 50,      xlim=c(-1,0),   ylim=c(0.75,1),   maxVary = 5)
plot(amShifted, "stopInconclusive",  alertK = 50,      xlim=c(-1,0),   ylim=c(0,1),      maxVary = 5)
plot(amShifted, "stopNotImpactful",  alertK = 50,      xlim=c(-1,0),   ylim=c(0,1),      maxVary = 5)
plot(amShifted, "stopNotTrivial",    alertK = 50,      xlim=c(-1,0),   ylim=c(0,1),      maxVary = 5)

# Explore number of steps in affirmation step
plot(amShifted, "n",                 waitWidth = 0.3,  xlim=c(-1,0),   ylim=c(0,801),    maxVary = 5)
plot(amShifted, "bias",              waitWidth = 0.3,  xlim=c(-1,0),   ylim=c(-0.1,0.1), maxVary = 5)
plot(amShifted, "mse",               waitWidth = 0.3,  xlim=c(-1,0),   ylim=c(0,0.1),    maxVary = 5)
plot(amShifted, "cover",             waitWidth = 0.3,  xlim=c(-1,0),   ylim=c(0.75,1),   maxVary = 5)
plot(amShifted, "stopInconclusive",  waitWidth = 0.3,  xlim=c(-1,0),   ylim=c(0,1),      maxVary = 5)
plot(amShifted, "stopNotImpactful",  waitWidth = 0.3,  xlim=c(-1,0),   ylim=c(0,1),      maxVary = 5)
plot(amShifted, "stopNotTrivial",    waitWidth = 0.3,  xlim=c(-1,0),   ylim=c(0,1),      maxVary = 5)

# Unrestricted Sample Size with lag of remaining outcomes
plot(amShifted, "n",                 alertK = 50,     xlim=c(-1,0),   ylim=c(0,801), maxVary = 10, sizeRestrictions = "lag")
plot(amShifted, "stopInconclusive",  alertK = 50,     xlim=c(-1,0),   ylim=c(0,1),   maxVary = 10, sizeRestrictions = "lag")
plot(amShifted, "stopInconclusive",  waitWidth = 0.3, xlim=c(-1,0),   ylim=c(0,1),   maxVary = 10, sizeRestrictions = "lag")

# Max N with immediate outcomes
plot(amShifted, "n",     alertK = 50,     xlim=c(-1,0),   ylim=c(0,801), maxVary = 10, sizeRestrictions = "maxN")
plot(amShifted, "stopInconclusive",     alertK = 50,     xlim=c(-1,0),   ylim=c(0,1), maxVary = 10, sizeRestrictions = "lag")
summary(amShifted, alertK = 50, waitTime = 0.3, treatEffect = -0.5)
summary(amShifted)
summary(am)



# Simulate AM trial
# One-sided deltas
# Outcome has sd = 2
am2 <-  sgpvAM(nreps            = 100,
               maxAlertSteps    = 100,       lookSteps = 1, waitWidths = sort(c(seq(0.2, 1.2, by = 0.2),c(0.5,0.7))),
               dataGeneration   = rnorm,   dataGenArgs = list(n=800, sd=2),
               effectGeneration = 0,
               modelFit         = lmCI,
               pointNull = 0, deltaL2 = -1, deltaL1=-0.2, deltaG1=NA, deltaG2=NA,
               monitoringIntervalLevel = 0.05,
               maxN = 200, lagOutcomeN = 50)
amShifted2 <- locationShift(am2, shiftedThetas = seq(-2,0.8,by=0.4))

# Unrestricted Sample Size
# Explore wait width (first with alert k specified at 50)
plot(amShifted2, "n",                 alertK = 0,       xlim=c(-2,1),    ylim=c(0,500))
plot(amShifted2, "n",                 alertK = 50,      xlim=c(-2,1),    ylim=c(0,500))
plot(amShifted2, "rejPN",             alertK = 50,      xlim=c(-2,1),   ylim=c(0,1))
plot(amShifted2, "bias",              alertK = 50,      xlim=c(-2,1),   ylim=c(-0.2,0.2))
plot(amShifted2, "mse",               alertK = 50,      xlim=c(-2,1),   ylim=c(0,0.3))
plot(amShifted2, "cover",             alertK = 50,      xlim=c(-2,1),   ylim=c(0.75,1))
plot(amShifted2, "stopInconclusive",  alertK = 50,      xlim=c(-2,1),   ylim=c(0,1))
plot(amShifted2, "stopNotImpactful",  alertK = 50,      xlim=c(-2,1),   ylim=c(0,1))
plot(amShifted2, "stopNotTrivial",    alertK = 50,      xlim=c(-2,1),   ylim=c(0,1))

# Explore number of steps in affirmation step
plot(amShifted2, "n",                waitWidth = 0.4,  xlim=c(-2,1),   ylim=c(0,500),    maxVary = 5)
plot(amShifted2, "bias",             waitWidth = 0.4,  xlim=c(-2,1),   ylim=c(-0.1,0.1), maxVary = 5)
plot(amShifted2, "mse",              waitWidth = 0.4,  xlim=c(-2,1),   ylim=c(0,0.3),    maxVary = 5)
plot(amShifted2, "cover",            waitWidth = 0.4,  xlim=c(-2,1),   ylim=c(0.75,1),   maxVary = 5)
plot(amShifted2, "stopInconclusive", waitWidth = 0.4,  xlim=c(-2,1),   ylim=c(0,1),      maxVary = 5)
plot(amShifted2, "stopNotImpactful", waitWidth = 0.4,  xlim=c(-2,1),   ylim=c(0,1),      maxVary = 5)
plot(amShifted2, "stopNotTrivial",   waitWidth = 0.4,  xlim=c(-2,1),   ylim=c(0,1),      maxVary = 5)

# Unrestricted Sample Size with lag of remaining outcomes
plot(amShifted2, "stopInconclusive",  alertK = 50,     xlim=c(-2,1),   ylim=c(0,1),   maxVary = 10, sizeRestrictions = "lag")
plot(amShifted2, "stopInconclusive",  waitWidth = 0.4, xlim=c(-2,1),   ylim=c(0,1),   maxVary = 10, sizeRestrictions = "lag")
plot(amShifted2, "stopInconclusive",  waitWidth = 0.6, xlim=c(-2,1),   ylim=c(0,1),   maxVary = 10, sizeRestrictions = "lag")
# Max N with immediate outcomes
plot(amShifted2, "n",                    alertK = 50, xlim=c(-2,1),   ylim=c(0,500), maxVary = 10, sizeRestrictions = "maxN")
plot(amShifted2, "stopInconclusive",     alertK = 50, xlim=c(-2,1),   ylim=c(0,1),   maxVary = 10, sizeRestrictions = "lag")
summary(amShifted2, alertK = 50, waitTime = 0.6, treatEffect = -0.4)
summary(amShifted2)
summary(am2)

# Explore ECDF of Bias and Sample Size
ooo <- ecdf.sgpvAM(am = amShifted2, stat = "Bias",alertK = 50,treatEffect = 0,xlim = c(-1,2))
ooo <- ecdf.sgpvAM(am = amShifted2, stat = "Bias",waitWidth = 0.4,alertK = 50,treatEffect = 0,xlim = c(-1,2))
ooo <- ecdf.sgpvAM(am = amShifted2, stat = "Size",alertK = 50,treatEffect = 0,xlim = c(0,400))

# Explore quantiles of sample size for fully specified study design
# See probability of stopping by certain sample size under fully specified design
ooo <- ecdf.sgpvAM(am = amShifted2, stat = "Size",waitWidth = 0.4,alertK = 50,treatEffect = 0,xlim = c(0,400))
quantile(ooo, probs=c(0,0.05,0.1,0.2,0.5,0.8,0.9,0.95,1))
ooo(200)

```
