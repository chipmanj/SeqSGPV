# sgpvAM.R
# J Chipman
#
# Adaptive monitoring design
#
#
# library(plyr)


sgpvAM <- function(mcmcData=NULL, nreps, maxAlertSteps=100, lookSteps=1,
                   waitWidths = c(0.15, 0.20, 0.30, 0.35, 0.40, 0.45, 0.50, 0.60),
                   monitoringIntervalLevel = 0.05, ...){


  # 1 collect array of simulated data
  if(is.null(mcmcData)){
    mcmcMonitoring <- plyr::rlply(.n = nreps, .expr = {
      # sgpvAMdata(rnorm, dataGenArgs = list(n=800), waitWidths = seq(0.15, 0.6, by = 0.05),
      #            lookSteps = 5, effectGeneration = 0.3,
      #            deltaL2 = -0.4, deltaL1 = -0.3, deltaG1 = 0.3, deltaG2 = 0.4,
      #            monitoringIntervalLevel = 0.05)})
      sgpvAMdata( monitoringIntervalLevel = monitoringIntervalLevel, ... )
    })
  }



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

      mcmcMonitoring[[i]] <- sgpvAMdata(rnorm, dataGenArgs = list(n=800), waitWidths = seq(0.15, 0.6, by = 0.05),
                                              lookSteps = 5, effectGeneration = 0.3,
                                              deltaL2 = -0.4, deltaL1 = -0.3, deltaG1 = 0.3, deltaG2 = 0.4,
                                              monitoringIntervalLevel = 0.05)



      # sgpvAMdata(dataGeneration = rnorm, dataGenArgs = list(n=getMore[i]),
      #                                    effectGeneration = 0.5,
      #                                    deltaL2 = -1, deltaL1 = -0.15, deltaG1 = 0.15, deltaG2 = 1,
      #                                    monitoringIntervalLevel = 0.05, existingData = mcmcMonitoring[[i]])

      # mcmcMonitoring[[i]] <- sgpvAMdata(monitoringIntervalLevel = monitoringIntervalLevel, ...)

    }

    getMore      <- unlist(lapply(mcmcMonitoring, mcmcMonitoringEnoughCheck,
                                  maxAlertSteps = maxAlertSteps,
                                  minWW         = min(waitWidths)))
    getMoreWhich <- which(getMore > 0)

  }




  # waitWidths        <- c(0.15, 0.20, 0.30, 0.35, 0.40, 0.45, 0.50, 0.60)
  mcmcEndOfStudyAve <- list()
  mcmcEndOfStudyVar <- list()
  mcmcEndOfStudy    <- list()
  for (w in 1:length(waitWidths)){

    ww <- waitWidths[w]

    # 3 adaptively monitor simulated data across multiple burn ins
    mcmcEOS <- simplify2array(lapply(mcmcMonitoring, sgpvAMrules, waitWidth = ww,
                                     lookSteps = lookSteps,  maxAlertSteps = maxAlertSteps,
                                     monitoringIntervalLevel = monitoringIntervalLevel))


    # i = 1
    # i = i+1
    # sgpvAMrules(mcmcMonitoring[[i]], waitWidth = ww, lookSteps = lookSteps, maxAlertSteps = maxAlertSteps,
    #             monitoringIntervalLevel = monitoringIntervalLevel)



    # 4 aggregate simulated data
    #   average performance and mse
    #   ecdf of n and bias
    ooAve <- plyr::aaply(mcmcEOS, .margins = c(1,2), .fun = mean)
    ooVar <- plyr::aaply(mcmcEOS, .margins = c(1,2), .fun = var )

    mcmcEndOfStudyAve[[paste0("width_",ww)]] <- cbind(ooAve, mse = ooVar[,"bias"] + ooAve[,"bias"]^2)

    mcmcEndOfStudyEcdfSize        <- apply(mcmcEOS[,"n",], 1, ecdf)
    names(mcmcEndOfStudyEcdfSize) <- paste0("alertK_",mcmcEOS[,"alertK",1])

    mcmcEndOfStudyEcdfBias        <- apply(mcmcEOS[,"bias",], 1, ecdf)
    names(mcmcEndOfStudyEcdfBias) <- paste0("alertK_",mcmcEOS[,"alertK",1])


    mcmcEndOfStudy[[paste0("width_",ww)]] <-
      list(mcmcEndOfStudyAve      = mcmcEndOfStudyAve,
           mcmcEndOfStudyEcdfSize = mcmcEndOfStudyEcdfSize,
           mcmcEndOfStudyEcdfBias = mcmcEndOfStudyEcdfBias)

  }



  return(mcmcEndOfStudy)

}


# # With no previously generated data
# a <- sgpvAM(nreps = 10, maxAlertSteps = 100, lookSteps = 5, waitWidths = seq(0.15, 0.6, by = 0.05),
#             dataGeneration = rnorm,   dataGenArgs = list(n=800),
#             effectGeneration = 0.3,
#             deltaL2 = -0.4, deltaL1=-0.3, deltaG1=0.3, deltaG2=0.4,
#             monitoringIntervalLevel=0.05)






# 2 adaptively monitor simulated data across multiple burn ins
# a <- aaply(mcmcMonitoringFixed, .margins = 3, .fun = function(x) {
#   t(sgpvAMrules(data = x, waitWidth = 0.3,
#                 lookSteps = 1, maxAlertSteps = 100,
#                 monitoringIntervalLevel = 0.05))
# })
#
#
# # 3 aggregate simulated data
# aa <- t(aaply(a, .margins = c(2,3), .fun = mean))
#
