# sgpvAM package

The sgpvAM package allows the user to obtain study design operating characteristics under a variety of settings for adaptive monitoring using the second generation p-value.  

MCMC Replicates: The user may use the sgpvAM function to generate mcmc replicates of outcomes and intervention assignments along with an estimate of the effect and a lower- and upper- interval bound; replicates are generated using parallel computing.  Alternatively, the user may provide their own generated data together with an estimated effect and interval bounds.  When using the sgpvAM function, the user specifies the data generation function (any of the r[dist] such as rnorm) along with arguments to the function.  Similarly, the user specifies effect generation.  At this point, only fixed effects have been thoroughly tested.  However, by specified a distribution for the effects, the user may explore False Discovery Probabilities and other operating characterstics dependent on distributional assumptions of the effect.

One- vs Two-Sided Hypotheses: Clinical Guideposts defining regions of Trivial and Meaningful Effects must be provided though may be one- or two-sided.  The point null must be within the Trivial Region and cannot be a boundary of the region.  For general nomenclature, inputs to define the regions are: deltaL2 (the Clinically Meaningful Boundary less than the point null), deltaL1 (the Trivial Region Boundary less than the point null), deltaG1 (the Trivial Region Boundary greater than the point null), and deltaG2 (the Clinically Meaningful Boundary greater than the point null).

Tuning study parameters: To maximize performance of operating characteristics under a given sample size, the sgpvAM function allows the user to specifiy multiple wait time settings, frequency of looks, and number of steps before affirming a stopping rule.  (The wait time is the time until the expected Confidence Interval Width achieves a certain length or less).

Operating characteristics under normal outcomes: After generating the operating characterstics under a fixed normal outcome, the user may use the locationShift function to obtain operating characterstics under a range of fixed treatment effects.  The function uses the saved mcmc replicates and adds to them if needed for additional monitoring.

ECDF of sample size and bias: Once a study design has been selected based on average performance (sample size, bias, and error probabilities), the user may use the ecdf.sgpv function to see the empirical cumulative distribution for sample size and bias under a specific design.  The user may see the estimated probability of the sample size exceeding a certain maximum sample size.

General suggestions: Computations may be time consuming.  It is recommended to start with 1000 replicates to get a general sense of average sample size and error probabilities under a variety of investigated wait times and affirmation steps.  Investigating many wait times increases the computational burden.  It is also recommended to generate mcmc replicates in the (or one of the) mid point(s) between the Clinically Trivial and Meaningful Regions.  This is the region with greatest expected sample size and reduces the burden of the locationShift function to generate more data.


# Required install of sgpv package

The sgpvAM package requires the sgpv package, which is currently available on github.  If this has not yet been installed, call:  

devtools::install_github("weltybiostat/sgpv")

# Github install of sgpvAM

devtools::install_github("chipmanj/sgpvAM")


# Example 1: A one-sided hypothesis with standard normal data

Suppose a study has a standard normal continuous outcome with a one-sided hypothesis where Trivial Effects are defined as (-∞, 0.2] and Meaningful Effects are defined as [0.5, ∞).  Practically, we are limited to collecting at most 250 observations, yet outcomes are not observed immediately.  There is a lag time such that at the end of enrollment the study will have observed 200 observations (with 50 yet to be observed).  The study has potential to receive additional funding, so we are also interested to see the operating characteristics when unrestricted in sample size.  We are still determining the best wait time and number of affirmation steps, so we will investigate wait times between 0.2 and 0.35, including 0.275.  Since patients are enrolled in batches of ten patients at a time and outcomes are also observed roughly every ten patients, we will monitor the study at every tenth patient.  

We’ll use the sgpvAM package, making use of parallel clustering (which on this machine is 8 cores) and start with only 1000 replicates to explore operating characteristics for different Confidence Interval Width wait times and different number of steps before affirming a rule to stop.  Using the sgpvAM the following call is made.


Suggest to start small with 1000 reps to explore study designs.

```r
devtools::install_github("chipmanj/sgpvAM")
library(sgpvAM)
#
# Simulate AM trial
# Two-sided deltas
# defaults to sd = 1
system.time(am <-  sgpvAM(nreps            = 100,
                          maxAlertSteps    = 100,       lookSteps = 10,  kSteps = 10,
                          waitWidths       = seq(0.2, 0.35, by = 0.025),
                          dataGeneration   = rnorm,   dataGenArgs = list(n=200),
                          effectGeneration = 0,
                          modelFit         = lmCI,
                          pointNull = 0, deltaL2 = NA, deltaL1=NA, deltaG1=0.2, deltaG2=0.5,
                          monitoringIntervalLevel = 0.05,
                          printProgress = TRUE,
                          maxN = 200, lagOutcomeN = 50, 
                          cores = detectCores()))

system.time(amShifted <- locationShift(am, shiftedThetas = seq(-1, 1, by = 0.1)))

# Unrestricted Sample Size
# Explore wait width (first with alert k specified at 50)
plot(amShifted, "n",                 alertK = 20,      xlim=c(-1, 1),   ylim=c(0,500),    maxVary = 10)
plot(amShifted, "n",                 alertK = 20,      xlim=c(-1, 1),   ylim=c(0,500),    maxVary = 5)
plot(amShifted, "rejPN",             alertK = 20,      xlim=c(-1, 1),   ylim=c(0,1),      maxVary = 10)
plot(amShifted, "bias",              alertK = 20,      xlim=c(-1, 1),   ylim=c(-0.1,0.1), maxVary = 5)
plot(amShifted, "mse",               alertK = 20,      xlim=c(-1, 1),   ylim=c(0,0.1),    maxVary = 5)
plot(amShifted, "cover",             alertK = 20,      xlim=c(-1, 1),   ylim=c(0.75,1),   maxVary = 5)
plot(amShifted, "stopInconclusive",  alertK = 20,      xlim=c(-1, 1),   ylim=c(0,1),      maxVary = 5)
plot(amShifted, "stopNotImpactful",  alertK = 20,      xlim=c(-1, 1),   ylim=c(0,1),      maxVary = 5)
plot(amShifted, "stopNotTrivial",    alertK = 20,      xlim=c(-1, 1),   ylim=c(0,1),      maxVary = 5)

# Explore number of steps in affirmation step
plot(amShifted, "n",                 waitWidth = 0.275,  xlim=c(-1, 1),   ylim=c(0,500),    maxVary = 6)
plot(amShifted, "bias",              waitWidth = 0.275,  xlim=c(-1, 1),   ylim=c(-0.1,0.1), maxVary = 5)
plot(amShifted, "mse",               waitWidth = 0.275,  xlim=c(-1, 1),   ylim=c(0,0.1),    maxVary = 5)
plot(amShifted, "cover",             waitWidth = 0.275,  xlim=c(-1, 1),   ylim=c(0.75,1),   maxVary = 5)
plot(amShifted, "stopInconclusive",  waitWidth = 0.275,  xlim=c(-1, 1),   ylim=c(0,1),      maxVary = 5)
plot(amShifted, "stopNotImpactful",  waitWidth = 0.275,  xlim=c(-1, 1),   ylim=c(0,1),      maxVary = 5)
plot(amShifted, "stopNotTrivial",    waitWidth = 0.275,  xlim=c(-1, 1),   ylim=c(0,1),      maxVary = 5)


# Unrestricted Sample Size with lag of remaining outcomes
plot(amShifted, "n",                 alertK = 10,     xlim=c(-1,1),   ylim=c(0,500), maxVary = 10, sizeRestrictions = "lag")
plot(amShifted, "stopInconclusive",  alertK = 10,     xlim=c(-1,1),   ylim=c(0,1),   maxVary = 10, sizeRestrictions = "lag")
plot(amShifted, "stopInconclusive",  waitWidth = 0.275, xlim=c(-1,1),   ylim=c(0,1),   maxVary = 10, sizeRestrictions = "lag")

# Max N with immediate outcomes
plot(amShifted, "n",     alertK = 10,     xlim=c(-1,1),   ylim=c(0,500), maxVary = 10, sizeRestrictions = "maxN")
plot(amShifted, "stopInconclusive",     alertK = 10,     xlim=c(-1,1),   ylim=c(0,1), maxVary = 10, sizeRestrictions = "lag")


summary(amShifted, alertK = 10, waitTime = 0.275, treatEffect = 0.2)
summary(amShifted)
summary(am)

```


# Example 2: A two-sided hypothesis with normal data

Suppose a study has a normal continuous outcome with standard deviation equal to two.  This study finds it beneficial to study both hypotheses of benefit and harm.  Clinically Trivial Effects are defined as [-0.4, 0.4] and clinically meaningful effects are defined as (-∞, -1] and [1, ∞).  Enrolling and observing outcomes is inexpensive and immediate, to the point that the sample size may be unrestricted and the study may be monitored fully sequentially.  We are still determining the best wait time and number of affirmation steps, so we will investigate wait times of Confidence Interval Widths between 0.4 and 0.7, including 0.55.

We’ll use the sgpvAM package, making use of parallel clustering (which on this machine is 8 cores) and start with only 1000 replicates to explore operating characteristics for different Confidence Interval Width wait times and different number of steps before affirming a rule to stop.  Using the sgpvAM the following call is made.


Suggest to start small with 1000 reps to explore study designs.

```{r}
# Simulate AM trial
# One-sided deltas
# Outcome has sd = 2
system.time(am2 <-  sgpvAM(nreps           = 100,
                          maxAlertSteps    = 100,       lookSteps = 1, kSteps=10, waitWidths = seq(0.4,0.7,length.out=5),
                          dataGeneration   = rnorm,   dataGenArgs = list(n=200, sd=2),
                          effectGeneration = 0,
                          modelFit         = lmCI,
                          pointNull = 0, deltaL2 = -1, deltaL1=-0.2, deltaG1=0.2, deltaG2=1,
                          monitoringIntervalLevel = 0.05,
                          cores=detectCores()))

system.time(amShifted2 <- locationShift(am2, shiftedThetas = seq(-2,2,by=0.4)))



# Unrestricted Sample Size
# Explore wait width (first with alert k specified at 50)
plot(amShifted2, "n",                 alertK = 0,       xlim=c(-2,2),   ylim=c(0,500))
plot(amShifted2, "rejPN",             alertK = 50,      xlim=c(-2,2),   ylim=c(0,1))
plot(amShifted2, "bias",              alertK = 50,      xlim=c(-2,2),   ylim=c(-0.2,0.2))
plot(amShifted2, "mse",               alertK = 50,      xlim=c(-2,2),   ylim=c(0,0.3))
plot(amShifted2, "cover",             alertK = 50,      xlim=c(-2,2),   ylim=c(0.75,1))
plot(amShifted2, "stopInconclusive",  alertK = 50,      xlim=c(-2,2),   ylim=c(0,1))
plot(amShifted2, "stopNotImpactful",  alertK = 50,      xlim=c(-2,2),   ylim=c(0,1))
plot(amShifted2, "stopNotTrivial",    alertK = 50,      xlim=c(-2,2),   ylim=c(0,1))

# Explore number of steps in affirmation step
plot(amShifted2, "n",                waitWidth = 0.55,  xlim=c(-2,2),   ylim=c(0,500),    maxVary = 5)
plot(amShifted2, "bias",             waitWidth = 0.55,  xlim=c(-2,2),   ylim=c(-0.1,0.1), maxVary = 5)
plot(amShifted2, "mse",              waitWidth = 0.55,  xlim=c(-2,2),   ylim=c(0,0.3),    maxVary = 5)
plot(amShifted2, "cover",            waitWidth = 0.55,  xlim=c(-2,2),   ylim=c(0.75,1),   maxVary = 5)
plot(amShifted2, "stopInconclusive", waitWidth = 0.55,  xlim=c(-2,2),   ylim=c(0,1),      maxVary = 5)
plot(amShifted2, "stopNotImpactful", waitWidth = 0.55,  xlim=c(-2,2),   ylim=c(0,1),      maxVary = 5)
plot(amShifted2, "stopNotTrivial",   waitWidth = 0.55,  xlim=c(-2,2),   ylim=c(0,1),      maxVary = 5)


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
