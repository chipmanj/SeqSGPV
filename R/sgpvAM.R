# sgpvAM.R
# J Chipman


sgpvAM <- function(mcmcData=NULL, nreps, maxAlertSteps=100, lookSteps=1,
                   waitWidths = c(0.15, 0.20, 0.30, 0.35, 0.40, 0.45, 0.50, 0.60),
                   monitoringIntervalLevel = 0.05, returnSimData = FALSE, cores, ...){


  # 1 collect list of simulated data
  if(is.null(mcmcData)){
         mcmcMonitoring <- sgpvAMdata(nreps = nreps,
                                      monitoringIntervalLevel = monitoringIntervalLevel,
                                      cores = cores, ... )
         # mcmcMonitoring <- sgpvAMdata(nreps = nreps, monitoringIntervalLevel = monitoringIntervalLevel,
         #                              maxAlertSteps = 100, lookSteps = 1, waitWidths = seq(0.15 ,0.6, by = 0.05),
         #                              dataGenArgs = list(n=800), dataGeneration = rnorm,
         #                              effectGeneration = 0,
         #                              deltaL2 = -0.4, deltaL1=-0.3, deltaG1=0.3, deltaG2=0.4)

  } else mcmcMonitoring <- mcmcData



  # 2 Make sure all simulations will continue until completion
  #   Look for stability of sgpv for last set of maxAlert patients
  getMore      <- unlist(lapply(mcmcMonitoring, mcmcMonitoringEnoughCheck,
                                maxAlertSteps = maxAlertSteps,
                                minWW         = min(waitWidths)))
  getMoreWhich <- which(getMore > 0)
  getMoreWhich

  while(sum(getMore) > 0){
    print(paste("Adding another", max(getMore), "observations to ensure simulations go to completion."))

    for (i in getMoreWhich){

      mcmcMonitoring[[i]] <- sgpvAMdataSingle(monitoringIntervalLevel = 0.05,
                                              existingData = mcmcMonitoring[[i]],
                                              ... )

    }

    getMore      <- unlist(lapply(mcmcMonitoring, mcmcMonitoringEnoughCheck,
                                  maxAlertSteps = maxAlertSteps,
                                  minWW         = min(waitWidths)))
    getMoreWhich <- which(getMore > 0)

  }


  # 3 adaptively monitor simulated data across multiple burn ins
  # 4 aggregate simulated data
  #   average performance and mse
  #   ecdf of n and bias

  mcmcEndOfStudyAve <- list()
  mcmcEndOfStudyVar <- list()
  mcmcEndOfStudy    <- list()
  for (w in 1:length(waitWidths)){

    ww <- waitWidths[w]

    # 3 adaptively monitor simulated data across multiple burn ins
    mcmcEOS <- simplify2array(lapply(mcmcMonitoring, sgpvAMrules, waitWidth = ww,
                                     lookSteps = lookSteps,  maxAlertSteps = maxAlertSteps,
                                     monitoringIntervalLevel = monitoringIntervalLevel))


    # 4 aggregate simulated data
    #   average performance and mse
    #   ecdf of n and bias
    ooAve <- plyr::aaply(mcmcEOS, .margins = c(1,2), .fun = mean)
    ooVar <- plyr::aaply(mcmcEOS, .margins = c(1,2), .fun = var )

    mcmcEndOfStudyAve                <- cbind(ooAve, mse = ooVar[,"bias"] + ooAve[,"bias"]^2)
    mcmcEndOfStudyEcdfSize           <- apply(mcmcEOS[,"n",], 1, ecdf)
    mcmcEndOfStudyEcdfBias           <- apply(mcmcEOS[,"bias",], 1, ecdf)
    names(mcmcEndOfStudyEcdfSize)    <- paste0("alertK_",mcmcEOS[,"alertK",1])
    names(mcmcEndOfStudyEcdfBias)    <- paste0("alertK_",mcmcEOS[,"alertK",1])


    mcmcEndOfStudy[[paste0("width_",ww)]] <-
      list(mcmcEndOfStudyAve      = mcmcEndOfStudyAve,
           mcmcEndOfStudyEcdfSize = mcmcEndOfStudyEcdfSize,
           mcmcEndOfStudyEcdfBias = mcmcEndOfStudyEcdfBias)

  }


  if(returnSimData){
    out <- list(mcmcMonitoring = mcmcMonitoring,
                mcmcEndOfStudy = mcmcEndOfStudy)
  } else {
    out <- mcmcEndOfStudy
  }

  return(out)

}


# Examples
# Generate Data
# am <- sgpvAMdata(nreps = 20, maxAlertSteps = 100, lookSteps = 1, waitWidths = seq(0.25, 0.6, by = 0.05),
#                  dataGeneration = rnorm,   dataGenArgs = list(n=800),
#                  modelFit = lmCI,
#                  effectGeneration = 0,
#                  deltaL2 = -0.4, deltaL1=-0.3, deltaG1=0.3, deltaG2=0.4,
#                  monitoringIntervalLevel=0.05)


# No previously generated data
# am1 <- sgpvAM(nreps = 20, maxAlertSteps = 100, lookSteps = 1, waitWidths = seq(0.25, 0.6, by = 0.05),
#              dataGeneration = rnorm,   dataGenArgs = list(n=800),
#              modelFit = lmCI,
#              effectGeneration = 0,
#              deltaL2 = -0.4, deltaL1=-0.3, deltaG1=0.3, deltaG2=0.4,
#              monitoringIntervalLevel=0.05)
# save(am1,file="~/Dropbox/test/am1.RData")

# Testing and development
# mcmcData = NULL
# nreps = 20; maxAlertSteps = 100; lookSteps = 1; waitWidths = seq(0.15 ,0.6, by = 0.05);
# dataGenArgs = list(n=800);
# effectGeneration = 0;
# deltaL2 = -0.4; deltaL1=-0.3; deltaG1=0.3; deltaG2=0.4;
# monitoringIntervalLevel=0.05




# Previously generated data
# amData <- sgpvAMdata(dataGeneration = rnorm,   dataGenArgs = list(n=800),
#                      effectGeneration = 0.3,
#                      deltaL2 = -0.4, deltaL1=-0.3, deltaG1=0.3, deltaG2=0.4,
#                      monitoringIntervalLevel=0.05)
#
#
# am <- sgpvAM(mcmcData = amData, maxAlertSteps = 100, lookSteps = 1, waitWidths = seq(0.15, 0.6, by = 0.05))
# save(am,file="~/Dropbox/test/am1.RData")
