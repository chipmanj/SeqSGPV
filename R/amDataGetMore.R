#' amDataGetMore
#'
#' Generate additional data if previously generated data were insufficient to run to conclusion with unrestricted n.
#'
#' @param insufficients Column vector of indices corresponding to existingDataList that needs additional data
#' @param exisitingDataList List of generated mcmc data in each element.
#' @param fork Fork clustering, works on POSIX systems (Mac, Linux, Unix, BSD) and not Windows.  Defaults to TRUE.
#' @param socket Socket clustering.  Defaults to TRUE yet only applies if FORK = FALSE.
#' @param ... Inputs to amDataSingleGetMore
#'
#' @export
amDataGetMore <- function(insufficients, existingDataList, os, fork=TRUE, socket = TRUE, cores = detectCores(), ...){

  if(fork==TRUE & os!="Windows"){
    # Only works on POSIX systems (Mac, Linux, Unix, BSD) and not Windows.
    mcmcMonitoringGetMore <- parallel::mclapply(insufficients,     amDataSingleGetMore,
                                                existingDataList = existingDataList,
                                                ... ,   mc.cores = cores)
  } else if(socket==TRUE){
    # Works on Mac and Windows; only slightly slower
    cl             <- parallel::makeCluster(cores)
    clusterCall(cl, function() library(sgpvAM))
    on.exit(parallel::stopCluster(cl), add = TRUE)
    mcmcMonitoringGetMore <- parallel::parLapply(cl, insufficients, amDataSingleGetMore,
                                                 existingDataList = existingDataList, ...)
  } else {

    mcmcMonitoringGetMore <- lapply(insufficients, amDataSingleGetMore,
                                    existingDataList = mcmcMonitoring, getMore = getMore,
                                    monitoringIntervalLevel = monitoringIntervalLevel,
                                    dataGeneration   = dataGeneration,   dataGenArgs   = dataGenArgs,
                                    effectGeneration = effectGeneration, effectGenArgs = effectGenArgs,
                                    randomize        = randomize,
                                    pointNull        = pointNull,
                                    deltaL2 = deltaL2, deltaL1 = deltaL1, deltaG1 = deltaG1, deltaG2 = deltaG2,
                                    modelFit         = modelFit)
  }


  # update existing data
  names(existingDataList)      <- 1:length(existingDataList)
  names(mcmcMonitoringGetMore) <- insufficients
  mcmcMonitoring               <- modifyList(existingDataList, mcmcMonitoringGetMore, keep.null = TRUE)


  # Return updated
  return(mcmcMonitoring)
}
