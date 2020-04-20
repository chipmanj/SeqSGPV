#' amDataSingleGetMore
#'
#' Generate additional observations for single mcmc simulation
#'
#' @export
amDataSingleGetMore <- function( existingDataList, iInsufficient, getMore,
                                 monitoringIntervalLevel,
                                 dataGeneration,   dataGenArgs,
                                 effectGeneration, effectGenArgs,
                                 effectScale,
                                 randomize,
                                 pointNull, deltaL2, deltaL1, deltaG1, deltaG2,
                                 modelFit ){

  dataGenArgs$n     <- getMore[iInsufficient]
  effectGeneration  <- existingDataList[[iInsufficient]][1,"theta"]

  amDataSingleMore  <- amDataSingle( monitoringIntervalLevel = monitoringIntervalLevel,
                                     existingData     = existingDataList[[iInsufficient]],
                                     dataGeneration   = dataGeneration,   dataGenArgs   = dataGenArgs,
                                     effectGeneration = effectGeneration, effectGenArgs = effectGenArgs,
                                     effectScale      = effectScale,
                                     randomize        = randomize,
                                     modelFit         = modelFit)

  amDataSingleMore  <- addStats(o = amDataSingleMore,
                                pointNull = pointNull,
                                deltaL2   = deltaL2, deltaL1 = deltaL1,
                                deltaG1   = deltaG1, deltaG2 = deltaG2)

  return(amDataSingleMore)
}
