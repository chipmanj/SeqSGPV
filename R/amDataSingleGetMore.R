#' @title amDataSingleGetMore
#'
#' @description Generate additional observations for single mcmc simulation. amDataSingleGetMore is called within SeqSGPV.
#'
#' @param existingDataList Previously provided data
#' @param iInsufficient Single index of list element with insufficient data for SeqSGPV to be conclusion
#' @param getMore List of additional data to generate
#' @param dataGeneration See SeqSGPV
#' @param dataGenArgs See SeqSGPV
#' @param effectGeneration See SeqSGPV
#' @param effectGenArgs See SeqSGPV
#' @param effectScale See SeqSGPV
#' @param allocation See SeqSGPV
#' @param randomize TRUE if length(allocation) > 1
#' @param effectPN See SeqSGPV
#' @param null See SeqSGPV
#' @param deltaL2 See SeqSGPV
#' @param deltaL1 See SeqSGPV
#' @param deltaG1 See SeqSGPV
#' @param deltaG2 See SeqSGPV
#' @param modelFit See SeqSGPV
#' @param modelFitArgs See SeqSGPV
#'
#' @export
amDataSingleGetMore <- function( existingDataList, iInsufficient, getMore,
                                 dataGeneration,   dataGenArgs,
                                 effectGeneration, effectGenArgs,
                                 effectScale,
                                 allocation,
                                 randomize,
                                 effectPN, null, deltaL2, deltaL1, deltaG1, deltaG2,
                                 modelFit, modelFitArgs ){

  if(randomize==FALSE){
    effectX <- "effect0"
  } else {
    effectX <- "effect1"
  }

  dataGenArgs$n     <- getMore[iInsufficient]
  effectGeneration  <- existingDataList[[iInsufficient]][1,effectX]

  amDataSingleMore  <- amDataSingle( existingData     = existingDataList[[iInsufficient]],
                                     dataGeneration   = dataGeneration,   dataGenArgs   = dataGenArgs,
                                     effectGeneration = effectGeneration, effectGenArgs = effectGenArgs,
                                     effectScale      = effectScale,
                                     allocation       = allocation,
                                     randomize        = randomize,
                                     modelFit         = modelFit,
                                     modelFitArgs     = modelFitArgs)

  amDataSingleMore  <- addStats(o         = amDataSingleMore,
                                randomize = randomize,
                                effectPN  = effectPN,
                                null      = null,
                                deltaL2   = deltaL2, deltaL1 = deltaL1,
                                deltaG1   = deltaG1, deltaG2 = deltaG2)

  return(amDataSingleMore)
}
