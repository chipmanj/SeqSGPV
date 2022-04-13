#' @title amDataGetMore
#'
#' @description Generate additional data if previously generated data were insufficient to run to conclusion with unrestricted n. amDataGetMore is called within SeqSGPV.
#'
#' @param insufficients Indexes of list with insufficient data for drawing SeqSGPV conclusions
#' @param existingDataList List of mcmc replicates
#' @param os See SeqSGPV
#' @param fork See SeqSGPV
#' @param socket See SeqSGPV
#' @param cores See SeqSGPV
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
amDataGetMore <- function(insufficients, existingDataList, os, fork=TRUE, socket = TRUE, cores = detectCores(),
                          getMore, dataGeneration,   dataGenArgs, effectGeneration, effectGenArgs,
                          effectScale, allocation, randomize, effectPN, null, deltaL2, deltaL1, deltaG1, deltaG2,
                          modelFit, modelFitArgs){

  if(fork==TRUE & os!="Windows"){
    # Only works on POSIX systems (Mac, Linux, Unix, BSD) and not Windows.
    mcmcMonitoringGetMore <- parallel::mclapply(insufficients,     amDataSingleGetMore,
                                                existingDataList = existingDataList, getMore       = getMore,
                                                dataGeneration   = dataGeneration,   dataGenArgs   = dataGenArgs,
                                                effectGeneration = effectGeneration, effectGenArgs = effectGenArgs,
                                                effectScale      = effectScale,
                                                allocation       = allocation,
                                                randomize        = randomize,
                                                effectPN         = effectPN,
                                                null             = null,
                                                deltaL2 = deltaL2, deltaL1 = deltaL1, deltaG1 = deltaG1, deltaG2 = deltaG2,
                                                modelFit         = modelFit,
                                                modelFitArgs     = modelFitArgs,
                                                mc.cores         = cores)

  } else if(socket==TRUE){
    # Works on Mac and Windows; only slightly slower
    cl             <- parallel::makeCluster(cores)
    clusterCall(cl, function() library(sgpvAM))
    on.exit(parallel::stopCluster(cl), add = TRUE)
    mcmcMonitoringGetMore <- parallel::parLapply(cl, insufficients, amDataSingleGetMore,
                                                 existingDataList = existingDataList, getMore       = getMore,
                                                 dataGeneration   = dataGeneration,   dataGenArgs   = dataGenArgs,
                                                 effectGeneration = effectGeneration, effectGenArgs = effectGenArgs,
                                                 effectScale      = effectScale,
                                                 allocation       = allocation,
                                                 randomize        = randomize,
                                                 effectPN         = effectPN,
                                                 null             = null,
                                                 deltaL2 = deltaL2, deltaL1 = deltaL1, deltaG1 = deltaG1, deltaG2 = deltaG2,
                                                 modelFit         = modelFit,
                                                 modelFitArgs     = modelFitArgs)
  } else {

    mcmcMonitoringGetMore <- lapply(insufficients, amDataSingleGetMore,
                                    existingDataList = existingDataList, getMore       = getMore,
                                    dataGeneration   = dataGeneration,   dataGenArgs   = dataGenArgs,
                                    effectGeneration = effectGeneration, effectGenArgs = effectGenArgs,
                                    effectScale      = effectScale,
                                    allocation       = allocation,
                                    randomize        = randomize,
                                    effectPN         = effectPN,
                                    null             = null,
                                    deltaL2 = deltaL2, deltaL1 = deltaL1, deltaG1 = deltaG1, deltaG2 = deltaG2,
                                    modelFit         = modelFit,
                                    modelFitArgs     = modelFitArgs)

  }


  # update existing data
  names(existingDataList)      <- 1:length(existingDataList)
  names(mcmcMonitoringGetMore) <- insufficients
  mcmcMonitoring               <- modifyList(existingDataList, mcmcMonitoringGetMore, keep.null = TRUE)


  # Return updated
  return(mcmcMonitoring)
}
