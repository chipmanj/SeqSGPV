#' @title amData
#'
#' @description Generate mcmc simulations of adaptive monitoring with parallel computing. amData is called within SeqSGPV.
#'
#' @param nreps See SeqSGPV
#' @param os See SeqSGPV
#' @param fork See SeqSGPV
#' @param socket See SeqSGPV
#' @param cores See SeqSGPV
#' @param dataGeneration See SeqSGPV
#' @param dataGenArgs See SeqSGPV
#' @param effectGeneration See SeqSGPV
#' @param effectGenArgs See SeqSGPV
#' @param effectScale See SeqSGPV
#' @param allocation See SeqSGPV
#' @param randomize TRUE if length(allocation) > 1
#' @param modelFit See SeqSGPV
#' @param modelFitArgs See SeqSGPV
#' @param existingData Previously provided data
#'
#' @export
amData <- function(nreps, os, fork=TRUE, socket = TRUE, cores = detectCores(),
                   dataGeneration,   dataGenArgs, effectGeneration, effectGenArgs,
                   effectScale, allocation, randomize, modelFit, modelFitArgs,existingData){

  if(fork==TRUE & os!="Windows"){
    # Only works on POSIX systems (Mac, Linux, Unix, BSD) and not Windows.
    mcmcMonitoring <- parallel::mclapply(1:nreps, FUN = function(x) {
      amDataSingle(dataGeneration   = dataGeneration,   dataGenArgs   = dataGenArgs,
                   effectGeneration = effectGeneration, effectGenArgs = effectGenArgs,
                   effectScale      = effectScale,
                   allocation       = allocation,
                   randomize        = randomize,
                   modelFit         = modelFit,
                   modelFitArgs     = modelFitArgs,
                   existingData     = existingData) })
  } else if(socket==TRUE){
    # Works on Mac and Windows; only slightly slower
    cl             <- parallel::makeCluster(cores)
    clusterCall(cl, function() library(SeqSGPVAM))
    on.exit(parallel::stopCluster(cl), add = TRUE)
    mcmcMonitoring <- parallel::parLapply(cl, 1:nreps, amDataSingle,
                                          dataGeneration   = dataGeneration,   dataGenArgs   = dataGenArgs,
                                          effectGeneration = effectGeneration, effectGenArgs = effectGenArgs,
                                          effectScale      = effectScale,
                                          allocation       = allocation,
                                          randomize        = randomize,
                                          modelFit         = modelFit,
                                          modelFitArgs     = modelFitArgs,
                                          existingData     = existingData)
  } else {
     mcmcMonitoring <- plyr::rlply(.n = nreps, .expr = { amDataSingle(
      dataGeneration   = dataGeneration,   dataGenArgs   = dataGenArgs,
      effectGeneration = effectGeneration, effectGenArgs = effectGenArgs,
      effectScale      = effectScale,
      allocation       = allocation,
      randomize        = randomize,
      modelFit         = modelFit,
      modelFitArgs     = modelFitArgs,
      existingData     = existingData ) })


  }

  return(mcmcMonitoring)
}
