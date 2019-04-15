# sgpvAM.R
# J Chipman
#
# Simulate operating characteristics of adaptive monitoring with SGPV


sgpvAM <- function(mcmcData=NULL, nreps, maxAlertSteps=100, lookSteps=1,
                   waitWidths = c(0.15, 0.20, 0.30, 0.35, 0.40, 0.45, 0.50, 0.60),
                   dataGeneration=NULL,   dataGenArgs,
                   effectGeneration=NULL, effectGenArgs,
                   modelFit,
                   pointNull, deltaL2, deltaL1, deltaG1, deltaG2,
                   maxN=NULL, lagOutcomeN=NULL,
                   monitoringIntervalLevel = 0.05, outData = TRUE){


  # 1 collect list of simulated data
  if(is.null(mcmcData)){
         mcmcMonitoring <- amData(nreps = nreps,
                                  monitoringIntervalLevel = monitoringIntervalLevel,
                                  dataGeneration   = dataGeneration,   dataGenArgs   = dataGenArgs,
                                  effectGeneration = effectGeneration, effectGenArgs = effectGenArgs,
                                  modelFit         = modelFit)

  } else mcmcMonitoring <- mcmcData



  # 2 Add stats (bias, rejPN, cover, sgpvNonTrival, sgpvFutility)
  mcmcMonitoring <- lapply(mcmcMonitoring,
                           addStats,
                           pointNull = pointNull,
                           deltaL2 = deltaL2, deltaL1=deltaL1, deltaG1=deltaG1, deltaG2=deltaG2)



  # 3 Make sure all generated simulations will continue until completion
  #   - Look for stability of sgpv for last set of maxAlert patients
  getMore      <- unlist(lapply(mcmcMonitoring, sgpvAM::mcmcMonitoringEnoughCheck,
                                maxAlertSteps = maxAlertSteps,
                                minWW         = min(waitWidths)))
  getMoreWhich <- which(getMore > 0)
  getMoreWhich

  if( !is.null(mcmcData) & sum(getMore) > 0 & is.null(dataGeneration) & is.null(effectGeneration)){

    stop("Provided mcmcData needs more observations to ensure study completes with unrestricted n\n
         Generate on own or provide data generation inputs.")

  } else {

    while(sum(getMore) > 0){

      print(paste("Adding another", max(getMore), "observations to ensure simulations go to completion."))

      for (i in getMoreWhich){

        # Add more data then statistics to generated dataset with insufficient n
        mcmcMonitoring[[i]] <- amDataSingle(monitoringIntervalLevel = monitoringIntervalLevel,
                                            existingData     = mcmcMonitoring[[i]],
                                            dataGeneration   = dataGeneration,   dataGenArgs   = dataGenArgs,
                                            effectGeneration = effectGeneration, effectGenArgs = effectGenArgs,
                                            modelFit         = modelFit)

        mcmcMonitoring[[i]] <- addStats(o = mcmcMonitoring[[i]],
                                        pointNull = pointNull,
                                        deltaL2   = deltaL2, deltaL1 = deltaL1,
                                        deltaG1   = deltaG1, deltaG2 = deltaG2)

      }

      # Continue checking until all datasets have sufficient n
      getMore      <- unlist(lapply(mcmcMonitoring, mcmcMonitoringEnoughCheck,
                                    maxAlertSteps = maxAlertSteps,
                                    minWW         = min(waitWidths)))
      getMoreWhich <- which(getMore > 0)

    }
  }


  # 4 adaptively monitor simulated data across multiple burn ins
  # 5 aggregate simulated data
  #   average performance and mse
  #   ecdf of n and bias

  mcmcEndOfStudyAve <- list()
  mcmcEndOfStudyVar <- list()
  mcmcEndOfStudy    <- list()
  for (w in 1:length(waitWidths)){

    ww <- waitWidths[w]

    # 3 adaptively monitor simulated data across multiple burn ins
    mcmcEOS <- simplify2array(lapply(mcmcMonitoring, sgpvAMrules,
                                     waitWidth               = ww,
                                     lookSteps               = lookSteps,
                                     maxAlertSteps           = maxAlertSteps,
                                     monitoringIntervalLevel = monitoringIntervalLevel,
                                     maxN = maxN, lagOutcomeN = lagOutcomeN))


    # 6 aggregate simulated data
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


  # Indicate whether to keep generated data
  if(outData==FALSE) mcmcMonitoring=NULL

  out <- list(mcmcMonitoring = mcmcMonitoring,
              mcmcEndOfStudy = mcmcEndOfStudy,
              inputs = list(mcmcData         = mcmcData,           nreps = nreps,
                            maxAlertSteps    = maxAlertSteps,  lookSteps = lookSteps,
                            waitWidths       = waitWidths,
                            dataGeneration   = dataGeneration,     dataGenArgs = dataGenArgs,
                            effectGeneration = effectGeneration,
                            modelFit         = modelFit,
                            pointNull        = pointNull,
                            deltaL2          = deltaL2, deltaL1 = deltaL1,
                            deltaG1          = deltaG1, deltaG2 = deltaG2,
                            maxN             = maxN,
                            lagOutcomeN      = lagOutcomeN,
                            monitoringIntervalLevel = monitoringIntervalLevel,
                            outData          = outData))

  return(out)

}


# Examples

# No previously generated data
# am1 <- sgpvAM(nreps = 100,
#               maxAlertSteps = 100, lookSteps = 1, waitWidths = seq(0.15, 0.6, by = 0.05),
#               dataGeneration = rnorm,   dataGenArgs = list(n=2000),
#               effectGeneration = 0,
#               modelFit = lmCI,
#               pointNull = 0, deltaL2 = -0.5, deltaL1=-0.2, deltaG1=0.2, deltaG2=0.5,
#               monitoringIntervalLevel=0.05)



# save(am1,file="~/Dropbox/test/am1.RData")

# Testing and development
# mcmcData = NULL
# nreps = 20; maxAlertSteps = 100; lookSteps = 1; waitWidths = seq(0.15 ,0.6, by = 0.05);
# dataGenArgs = list(n=800);
# effectGeneration = 0;
# deltaL2 = -0.4; deltaL1=-0.3; deltaG1=0.3; deltaG2=0.4;
# monitoringIntervalLevel=0.05
# modelFit = lmCI




# Previously generated data
# amData <- sgpvAMdata(dataGeneration = rnorm,   dataGenArgs = list(n=800),
#                      effectGeneration = 0.3,
#                      deltaL2 = -0.4, deltaL1=-0.3, deltaG1=0.3, deltaG2=0.4,
#                      monitoringIntervalLevel=0.05)
#
#
# am <- sgpvAM(mcmcData = amData, maxAlertSteps = 100, lookSteps = 1, waitWidths = seq(0.15, 0.6, by = 0.05))
# save(am,file="~/Dropbox/test/am1.RData")
