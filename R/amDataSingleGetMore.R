#' amDataSingleGetMore
#'
#' Generate additional observations for single mcmc simulation
#'
#' @export
amDataSingleGetMore <- function( existingDataList, iInsufficient, getMore,
                                 miLevel,
                                 dataGeneration,   dataGenArgs,
                                 effectGeneration, effectGenArgs,
                                 effectScale,
                                 randomize,
                                 effectPN, deltaL2, deltaL1, deltaG1, deltaG2,
                                 modelFit ){

  dataGenArgs$n     <- getMore[iInsufficient]
  effectGeneration  <- existingDataList[[iInsufficient]][1,"theta"]

  amDataSingleMore  <- amDataSingle( miLevel          = miLevel,
                                     existingData     = existingDataList[[iInsufficient]],
                                     dataGeneration   = dataGeneration,   dataGenArgs   = dataGenArgs,
                                     effectGeneration = effectGeneration, effectGenArgs = effectGenArgs,
                                     effectScale      = effectScale,
                                     randomize        = randomize,
                                     modelFit         = modelFit)

  amDataSingleMore  <- addStats(o = amDataSingleMore,
                                effectPN = effectPN,
                                deltaL2   = deltaL2, deltaL1 = deltaL1,
                                deltaG1   = deltaG1, deltaG2 = deltaG2)

  return(amDataSingleMore)
}
