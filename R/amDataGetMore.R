#' amDataGetMore
#'
#' Generate mcmc simulations of adaptive monitoring with parallel computing
#'
#' @export
amDataGetMore <- function(insufficients, existingDataList, fork=TRUE, socket = TRUE, cores = detectCores(), ...){

  if(fork==TRUE){
    # Only works on POSIX systems (Mac, Linux, Unix, BSD) and not Windows.
    mcmcMonitoringGetMore <- parallel::mclapply(insufficients,     amDataSingleGetMore,
                                                existingDataList = existingDataList,
                                                ... ,   mc.cores = cores)
  } else if(socket==TRUE){
    # Works on Mac and Windows; only slightly slower
    cl             <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    mcmcMonitoringGetMore <- parallel::parLapply(cl, insufficients, amDataSingleGetMore,
                                                 existingDataList = existingDataList, ...)
  } else {

    mcmcMonitoringGetMore <- lapply(getMoreWhich, amDataSingleGetMore,
                                    existingDataList = mcmcMonitoring, getMore = getMore,
                                    monitoringIntervalLevel = monitoringIntervalLevel,
                                    dataGeneration   = dataGeneration,   dataGenArgs   = dataGenArgs,
                                    effectGeneration = effectGeneration, effectGenArgs = effectGenArgs,
                                    pointNull = pointNull,
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

# test3 <- amDataGetMore(insufficients = getMoreWhich, existingDataList = mcmcMonitoring, getMore = getMore,
#               monitoringIntervalLevel = monitoringIntervalLevel,
#               dataGeneration   = dataGeneration,   dataGenArgs   = dataGenArgs,
#               effectGeneration = effectGeneration, effectGenArgs = effectGenArgs,
#               modelFit         = modelFit, fork=FALSE,socket=FALSE)
#
# test <- lapply(getMoreWhich, amDataGetMore, existingDataList = mcmcMonitoring, getMore = getMore,
#                monitoringIntervalLevel = monitoringIntervalLevel,
#                dataGeneration   = dataGeneration,   dataGenArgs   = dataGenArgs,
#                effectGeneration = effectGeneration, effectGenArgs = effectGenArgs,
#                modelFit         = modelFit)
