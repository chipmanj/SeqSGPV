#' amDataSingleGetMore
#'
#' Generate additional observations for single mcmc simulation
#'
#' @export
amDataSingleGetMore <- function( existingDataList, iInsufficient, getMore,
                                 monitoringIntervalLevel,
                                 dataGeneration,   dataGenArgs,
                                 effectGeneration, effectGenArgs,
                                 pointNull, deltaL2, deltaL1, deltaG1, deltaG2,
                                 modelFit ){

  dataGenArgs$n     <- getMore[iInsufficient]
  effectGeneration  <- existingDataList[[iInsufficient]][1,"theta"]

  amDataSingleMore  <- amDataSingle( monitoringIntervalLevel = monitoringIntervalLevel,
                                     existingData     = existingDataList[[iInsufficient]],
                                     dataGeneration   = dataGeneration,   dataGenArgs   = dataGenArgs,
                                     effectGeneration = effectGeneration, effectGenArgs = effectGenArgs,
                                     modelFit         = modelFit)

  amDataSingleMore  <- addStats(o = amDataSingleMore,
                                pointNull = pointNull,
                                deltaL2   = deltaL2, deltaL1 = deltaL1,
                                deltaG1   = deltaG1, deltaG2 = deltaG2)

  return(amDataSingleMore)
}


# test <- amDataSingleGetMore(mcmcMonitoring, getMoreWhich[1], getMore)
# test <- lapply(getMoreWhich, amDataSingleGetMore, existingDataList = mcmcMonitoring, getMore = getMore,
#                monitoringIntervalLevel = monitoringIntervalLevel,
#                dataGeneration   = dataGeneration,   dataGenArgs   = dataGenArgs,
#                effectGeneration = effectGeneration, effectGenArgs = effectGenArgs,
#                modelFit         = modelFit)
#
# names(mcmcMonitoring) <- 1:nreps
# names(test) <- getMoreWhich
#
# test2 <- modifyList(mcmcMonitoring,test,keep.null = TRUE)
