one arm, bernoulli outcomes
================

``` r
Sys.info()[c("sysname","release","version","machine")]
```

                                                                                                     sysname 
                                                                                                    "Darwin" 
                                                                                                     release 
                                                                                                    "22.6.0" 
                                                                                                     version 
    "Darwin Kernel Version 22.6.0: Mon Feb 19 19:48:53 PST 2024; root:xnu-8796.141.3.704.6~1/RELEASE_X86_64" 
                                                                                                     machine 
                                                                                                    "x86_64" 

In this example, 2000 replicates are used to get an initial sense of
operating characteristics. SeqSGPV is more time intensive for bernoulli
outcomes with multiple combinations of wait, step, N, and affirm. It is
recommended to start with a small number of replicates, even under 1000,
to get an initial sense of operating characteristics and then increase
for more precision in estimated operating characteristics.

``` r
library(SeqSGPV)
nreps <- 200
```

## Context

A prostate cancer trial is designed to assess an immediate outcome
following a novel surgery technique. Let $\pi$ be the probability of
success. The standard of care has a success rate of $\pi=0.20$. A
minimally scientifically meaningful effect is $\pi = 0.40$.

Without incorporating scientific relevance, a traditional hypothesis
could be:

H0: $\pi$ $\le$ 0.20  
H1: $\pi$ \> 0.40

For this early phase study, the clinician deems success probabilities
less than 0.225 as essentially equivalent or worse than the standard of
care (i.e. ROWPE). The PRISM is defined by ROE$`_{(0.225, 0.40)}`$.

The investigators says the study can afford up to 40 participants and
wants to know the design-based average sample size, Type I error, and
power across a range of treatment effects. The investigator wishes to
use repeated 95% credible intervals with a Beta(0.005, 0.005) prior on
the success probability.

As a secondary outcome, the investigator is interested in 3m quality of
life. Before measuring a single patient’s outcome, an additional 5-10
participants may be placed on the trial. The impact of delayed outcomes
will be assessed in [Section: Delayed Outcomes](#Delayed-outcomes) yet
simulated in the initial call to SeqSGPV assuming immediate outcomes.

The investigator wants a Type I error $\le$ 0.05 and is willing to
monitor the trial at every outcome.

## Simon’s two-stage design

For context, the operating characteristics of Simon’s two-stage design
are:

``` r
clinfun::ph2simon(pu = 0.2, pa = 0.4, ep1 = 0.05, ep2 = 0.2,nmax = 40)
```


     Simon 2-stage Phase II design 

    Unacceptable response rate:  0.2 
    Desirable response rate:  0.4 
    Error rates: alpha =  0.05 ; beta =  0.2 

            r1 n1  r  n EN(p0) PET(p0)   qLo   qHi
    Minimax  4 18 10 33  22.25  0.7164 0.168 1.000
    Optimal  3 14 11 38  21.24  0.6982 0.000 0.168

## SeqSGPV

``` r
# PRISM: deltaG1 = 0.225, deltaG2 = 0.4
# Wait time until first evaluation = 5 - 10
# Monitoring every 1 up to 3 outcomes
# maximum sample size = 35 - 40
# possible number of lag/delayed outcomes: 0, 5, 10
system.time(PRISM <-  SeqSGPV(nreps            = nreps,
                              dataGeneration   = rbinom, dataGenArgs = list(n=40, size=1, prob = .2),
                              effectGeneration = 0, effectGenArgs=NULL,  effectScale  = "identity",
                              allocation       = 1,
                              effectPN         = 0.2,
                              null             = "less",
                              PRISM            = list(deltaL2 = NA, deltaL1 = NA, 
                                                      deltaG1 = .225, deltaG2 = .4),
                              modelFit         = binomCI,
                              modelFitArgs     = list(conf.level=.95, 
                                                      prior.shape1=0.005, prior.shape2=0.005,
                                                      methods="bayes", type="central"),
                              wait             = c(15,20,25),
                              steps            = c(3,5),
                              affirm           = c(0,3,5),
                              lag              = c(0,5,10),
                              N                = c(38,40,45),
                              printProgress    = FALSE))
```

       user  system elapsed 
     35.168   3.413   7.436 

``` r
# Note: This step is typically done after evaluating operating characteristics
# under the point null. It will be shown again later.
# This step is done here for the sake of saving an Rmd cache with
# minimal retained data (after removing the simulated date).

# Obtain design under range of effects
se <- seq(-0.05, 0.3, by = 0.125)
system.time(PRISMse <- fixedDesignEffects(PRISM, shift = se))
```

    [1] "effect: -0.05"
    [1] "effect: 0.075"
    [1] "effect: 0.2"

       user  system elapsed 
    113.450  10.643  22.610 

``` r
# This next step is not required but is done for reducing the size of the Rmd cache.
PRISM$mcmcMonitoring <- NULL
```

Type I error under different monitoring frequencies. Increasing the
number of observations between assessments (steps from 1 to 5) and
requiring a stopping rule to be affirmed decreases Type I error.

``` r
par(mfrow=c(2,2))
plot(PRISM,stat = "rejH0", affirm=0, steps=3,lag=0,ylim=c(0.02, 0.10))
abline(h=.05)
plot(PRISM,stat = "rejH0", affirm=0, steps=3,lag=0,ylim=c(0.02, 0.10))
abline(h=.05)
plot(PRISM,stat = "rejH0", affirm=5, steps=5,lag=0,ylim=c(0.02, 0.10))
abline(h=.05)
plot(PRISM,stat = "rejH0", affirm=5, steps=5,lag=0,ylim=c(0.02, 0.10))
abline(h=.05)
```

<img src="README_files/figure-gfm/unnamed-chunk-5-1.png" style="display: block; margin: auto;" />

``` r
plot(PRISM,stat = "cover", affirm=5, steps=5,lag=0)
abline(h=.05)
```

<img src="README_files/figure-gfm/unnamed-chunk-5-2.png" style="display: block; margin: auto;" />

We may be interested in the impact upon average sample size for these
same changes to monitoring frequency.

``` r
par(mfrow=c(2,2))
plot(PRISM,stat = "n", affirm=0, steps=3,lag=0,ylim=c(15, 30))
abline(h=.05)
plot(PRISM,stat = "n", affirm=0, steps=5,lag=0,ylim=c(15, 30))
abline(h=.05)
plot(PRISM,stat = "n", affirm=5, steps=3,lag=0,ylim=c(15, 30))
abline(h=.05)
plot(PRISM,stat = "n", affirm=5, steps=5,lag=0,ylim=c(15, 30))
abline(h=.05)
```

<img src="README_files/figure-gfm/unnamed-chunk-6-1.png" style="display: block; margin: auto;" />

With a minimal increase in average sample size, $S=1, A=1$ controls Type
I error for all $W$ considered.

For comparison, consider SGPV interval monitoring the ROE with a
point-null boundary (ROE-PN-BOUNDED) – determining the trial successful
when there is evidence for an effect \> 0.2 and futile when there is
evidence the effect is \< 0.4.

``` r
# Change to monitoring null-bound ROE
inputs <- PRISM$inputs
inputs$PRISM$deltaG1 <- 0.20
inputs$PRISM$outData <- FALSE
system.time(ROE_PN_BOUNDED <-  do.call(SeqSGPV, inputs))
```

       user  system elapsed 
    398.824   9.506  72.357 

``` r
par(mfrow=c(2,2))
# Compare different designs
plot(PRISM,        stat = "rejH0", affirm=0, steps=3,lag=0,ylim=c(0.02, 0.10))
title(sub="ROE-PRISM", adj=0)
abline(h=.05)
plot(PRISM,        stat = "rejH0", affirm=0, steps=5,lag=0,ylim=c(0.02, 0.10))
title(sub="ROE-PRISM", adj=0)
abline(h=.05)
plot(ROE_PN_BOUNDED,stat = "rejH0", affirm=0, steps=3,lag=0,ylim=c(0.02, 0.10))
title(sub="ROE-PN-BOUNDED", adj=0)
abline(h=.05)
plot(ROE_PN_BOUNDED,stat = "rejH0", affirm=0, steps=5,lag=0,ylim=c(0.02, 0.10))
title(sub="ROE-PN-BOUNDED", adj=0)
abline(h=.05)
```

<img src="README_files/figure-gfm/unnamed-chunk-8-1.png" style="display: block; margin: auto;" />

Back to PRISM SeqSGPV, we would want to know the operating
characteristics of PRISM design across a range of effects.

``` r
# Obtain design under range of effects
se <- seq(-0.05, 0.3, by = 0.025)
system.time(PRISMse <- fixedDesignEffects(PRISM, shift = se))
```

``` r
par(mfrow=c(3,2))
plot(PRISMse, stat = "rejH0",            steps = 3, wait=25, affirm = 0, lag = 0)
plot(PRISMse, stat = "stopNotROPE",      steps = 3, wait=25, affirm = 0, lag = 0)
plot(PRISMse, stat = "stopNotROME",      steps = 3, wait=25, affirm = 0, lag = 0)
plot(PRISMse, stat = "stopInconclusive", steps = 3, wait=25, affirm = 0, lag = 0)
plot(PRISMse, stat = "n",                steps = 3, wait=25, affirm = 0, lag = 0)
plot(PRISMse, stat = "cover",            steps = 3, wait=25, affirm = 0, lag = 0)
```

<img src="README_files/figure-gfm/unnamed-chunk-10-1.png" style="display: block; margin: auto;" />

The interval coverage suffers

## Example interpretations following SeqSGPV monitoring of PRISM:

1.  The estimated success probability was 0.54 (95% credible interval:
    0.26, 0.81) which is evidence that the treatment effect is at least
    trivially better than the null hypothesis (p$`_{ROWPE}`$ = 0) and
    the evidence for being scientifically meaningful (p$`_{ROME}`$ =
    0.74).

2.  The estimated success probability was 0.11 (95% credible interval:
    ~0.00, 0.36) which is evidence that the treatment effect is not
    scientifically meaningful (p$`_{ROME}`$ = 0) and the evidence for
    being practically equivalent or worse than the point null is
    p$`_{ROWPE}`$=0.56.

3.  The estimated success probability was 0.35 (95% credible interval:
    0.21, 0.50), which is suggestive though inconclusive evidence to
    rule out at essentially null effects (p$`_{ROWPE}`$ = 0.04) and
    insufficient to rule out scientifically meaningful effects
    (p$`_{ROME}`$=0.35).

## Delayed outcomes

The same investigator in wants to study short term (3m) quality of life
as an endpoint. The below operating characteristics reflect a lag (or
delay) of five outcomes.

``` r
par(mfrow=c(2,2))
plot(PRISM,        stat = "lag.rejH0",         affirm=3, steps=3,lag=5,N=38,ylim=c(0.02, 0.10))
abline(h=0.05)
plot(PRISM,        stat = "lag.n",             affirm=3, steps=3,lag=5,N=38,ylim=c(10, 30))
plot(PRISMse,      stat = "lag.rejH0",         affirm=3, steps=3,lag=5,N=38,ylim=c(0, 1))
abline(h=c(0.05, 0.8, 0.9))
plot(PRISMse,      stat = "lag.n",             affirm=3, steps=3,lag=5,N=38,ylim=c(10, 30))
```
