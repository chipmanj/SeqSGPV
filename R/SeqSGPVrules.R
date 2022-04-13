#' @title SeqSGPVrules
#'
#' @description Applies monitoring frequencies to PRISM monitoring and returns the final observed observation. SeqSGPVrules is called within SeqSGPV.
#'
#' @param mcmcMonitoring List of mcmc replicates
#' @param os See SeqSGPV
#' @param fork See SeqSGPV
#' @param socket See SeqSGPV
#' @param cores See SeqSGPV
#' @param monitoringFrequency Matrix of monitoring frequency combinations
#' @param randomize TRUE if length(allocations) > 1
#'
#' @export
SeqSGPVrules <- function(mcmcMonitoring, os, fork=TRUE, socket = TRUE, cores = parallel::detectCores(), monitoringFrequency, randomize){


  if(fork==TRUE & os!="Windows"){
    # Only works on POSIX systems (Mac, Linux, Unix, BSD) and not Windows.
    mcmcEOS <- parallel::mclapply(mcmcMonitoring,     SeqSGPVrulesSingle,
                                  monitoringFrequency = monitoringFrequency,
                                  randomize           = randomize,
                                  mc.cores            = cores)

  } else if(socket==TRUE){
    # Works on Mac and Windows; only slightly slower
    cl             <- parallel::makeCluster(cores)
    clusterCall(cl, function() library(sgpvAM))
    on.exit(parallel::stopCluster(cl), add = TRUE)
    mcmcEOS <- parallel::parLapply(cl, mcmcMonitoring,     SeqSGPVrulesSingle,
                                       monitoringFrequency = monitoringFrequency,
                                       randomize           = randomize)
  } else {

    mcmcEOS <- lapply(mcmcMonitoring,       SeqSGPVrulesSingle,
                      monitoringFrequency = monitoringFrequency,
                      randomize           = randomize)
  }

  # Return End of Study
  return(simplify2array(mcmcEOS))
}
