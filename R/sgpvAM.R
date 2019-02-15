# sgpvAM.R
# J Chipman

library(plyr)


sgpvAM <- function(mcmcMonitoring, nreps, maxAlertSteps=100, lookSteps=1,
                   waitWidths = c(0.15, 0.20, 0.30, 0.35, 0.40, 0.45, 0.50, 0.60),
                   monitoringIntervalLevel = 0.05, ...){


  # 1 collect array of simulated data
  if(!missing(mcmcMonitoring)){
    mcmcMonitoring <- rlply(.n = nreps, .expr = {
      # sgpvAMdata(rnorm, dataGenArgs = list(n=900), effectGeneration = 0.5,
      #            deltaL2 = -0.5, deltaL1 = -0.15, deltaG1 = 0.15, deltaG2 = 0.5,
      #            monitoringIntervalLevel = 0.05)})
      sgpvAMdata( ... )
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

      # mcmcMonitoring[[i]] <- sgpvAMdata(dataGeneration = rnorm, dataGenArgs = list(n=getMore[i]),
      #                                    effectGeneration = 0.5,
      #                                    deltaL2 = -1, deltaL1 = -0.15, deltaG1 = 0.15, deltaG2 = 1,
      #                                    monitoringIntervalLevel = 0.05, existingData = mcmcMonitoring[[i]])

      mcmcMonitoring[[i]] <- sgpvAMdata(...)

    }

    getMore      <- unlist(lapply(mcmcMonitoring, mcmcMonitoringEnoughCheck, maxAlertSteps = 50))
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



    # 4 aggregate simulated data
    #   average performance and mse
    #   ecdf of n and bias
    ooAve <- aaply(mcmcEOS, .margins = c(1,2), .fun = mean)
    ooVar <- aaply(mcmcEOS, .margins = c(1,2), .fun = var )

    mcmcEndOfStudyAve[[paste0("width_",ww)]] <- cbind(ooAve, mse = ooVar[,"bias"] + ooAve[,"bias"]^2)

    mcmcEndOfStudyEcdfSize        <- apply(mcmcEOS[,"n",], 1, ecdf)
    names(mcmcEndOfStudyEcdfSize) <- paste0("alertK_",mcmcEOS[,"alertK",1])

    mcmcEndOfStudyEcdfBias        <- apply(mcmcEOS[,"bias",], 1, ecdf)
    names(mcmcEndOfEcdfBias)      <- paste0("alertK_",mcmcEOS[,"alertK",1])


    mcmcEndOfStudy[[paste0("width_",ww)]] <-
      list(mcmcEndOfStudyAve      = mcmcEndOfStudyAve,
           mcmcEndOfStudyEcdfSize = mcmcEndOfStudyEcdfSize,
           mcmcEndOfStudyEcdfBias = mcmcEndOfStudyEcdfBias)

  }



  return(mcmcEndOfStudy)

}







# 2 adaptively monitor simulated data across multiple burn ins
a <- aaply(mcmcMonitoringFixed, .margins = 3, .fun = function(x) {
  t(sgpvAMrules(data = x, waitWidth = 0.3,
                lookSteps = 1, maxAlertSteps = 100,
                monitoringIntervalLevel = 0.05))
})


# 3 aggregate simulated data
aa <- t(aaply(a, .margins = c(2,3), .fun = mean))

