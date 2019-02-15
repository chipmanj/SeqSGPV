# sgpvAM
Adaptive Monitoring Design Features using the second generation p-value


# install.packages("devtools")
devtools::install_github("chipmanj/sgpvAM")


# Example call
am <- sgpvAM(nreps = 10, maxAlertSteps = 100, lookSteps = 5, waitWidths = seq(0.15, 0.6, by = 0.05),
             dataGeneration = rnorm,   dataGenArgs = list(n=800),
             effectGeneration = 0.3,
             deltaL2 = -0.4, deltaL1=-0.3, deltaG1=0.3, deltaG2=0.4,
             monitoringIntervalLevel=0.05)
