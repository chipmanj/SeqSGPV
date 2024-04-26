one arm, bernoulli outcomes, prior effect generation
================

This example uses a small number of replicates to demonstrate how to
obtain end of study replicates where the effect generation follows a
random distribution. For each replicate, a success probability is drawn
uniformly between 0.1 and 0.5.

``` r
system.time(PRISM1c <-  SeqSGPV(nreps            = nreps,
                               dataGeneration   = rbinom, dataGenArgs = list(n=40, size=1, prob = .2),
                               effectGeneration = runif, effectGenArgs=list(n=1, min=-.1, max = .3),  effectScale  = "identity",
                               allocation       = 1,
                               effectPN         = 0.2,
                               null             = "less",
                               PRISM            = list(deltaL2 = NA, deltaL1 = NA, 
                                                       deltaG1 = .225, deltaG2 = .4),
                               modelFit         = binomCI,
                               modelFitArgs     = list(conf.level=.95, 
                                                       prior.shape1=0.005, prior.shape2=0.005,
                                                       methods="bayes", type="central"),
                               wait             = 5:10,
                               steps            = 1:3,
                               affirm           = 0:1,
                               lag              = c(0,5,10),
                               N                = 35:40,
                               printProgress    = FALSE))
```

       user  system elapsed 
      1.184   0.600   0.525 

``` r
 head(PRISM1c$mcmcEOS$W5_S1_A0_L0_N35)
```

         theta0     effect0  n       est         bias rejH0 cover stopNotROPE
    [1,]    0.2 -0.03855292 13 0.1541122 -0.007334855     0     1           0
    [2,]    0.2  0.19400437 35 0.3429020 -0.051102342     0     1           0
    [3,]    0.2  0.06706704 35 0.2572122 -0.009854818     0     1           0
    [4,]    0.2  0.06114171 24 0.2084548 -0.052686898     0     1           0
    [5,]    0.2 -0.03432873  9 0.1115427 -0.054128543     0     1           0
    [6,]    0.2  0.23191886 24 0.4167014 -0.015217486     1     1           1
         stopNotROME stopInconclusive lag.n   lag.est     lag.bias lag.rejH0
    [1,]           1                0    13 0.1541122 -0.007334855         0
    [2,]           0                1    35 0.3429020 -0.051102342         0
    [3,]           0                1    35 0.2572122 -0.009854818         0
    [4,]           1                0    24 0.2084548 -0.052686898         0
    [5,]           1                0     9 0.1115427 -0.054128543         0
    [6,]           0                0    24 0.4167014 -0.015217486         1
         lag.cover lag.stopNotROPE lag.stopNotROME lag.stopInconclusive
    [1,]         1               0               1                    0
    [2,]         1               0               0                    1
    [3,]         1               0               0                    1
    [4,]         1               0               1                    0
    [5,]         1               0               1                    0
    [6,]         1               1               0                    0
         lag.stopInconsistent lag.stopRejH0_YN lag.stopRejH0_NY wait steps affirm
    [1,]                    0                0                0    5     1      0
    [2,]                    0                0                0    5     1      0
    [3,]                    0                0                0    5     1      0
    [4,]                    0                0                0    5     1      0
    [5,]                    0                0                0    5     1      0
    [6,]                    0                0                0    5     1      0
         lag  N
    [1,]   0 35
    [2,]   0 35
    [3,]   0 35
    [4,]   0 35
    [5,]   0 35
    [6,]   0 35
