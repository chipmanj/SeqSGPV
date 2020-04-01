#' amSgpvRules
#'
#' Generate mcmc simulations of adaptive monitoring with parallel computing
#'
#' @export
#'
#'
sgpvAMrules <- function(mcmcMonitoring, os, fork=TRUE, socket = TRUE, cores = detectCores(),
                        waitTime,
                        lookSteps,
                        kSteps,
                        maxAlertSteps,
                        monitoringIntervalLevel,
                        getUnrestricted, maxN, lagOutcomeN){


  if(fork==TRUE & os!="Windows"){
    # Only works on POSIX systems (Mac, Linux, Unix, BSD) and not Windows.
    mcmcEOS <- parallel::mclapply(mcmcMonitoring, sgpvAMrulesSingle,
                                  waitTime                 = waitTime,
                                  lookSteps                = lookSteps,
                                  kSteps                   = kSteps,
                                  maxAlertSteps            = maxAlertSteps,
                                  monitoringIntervalLevel  = monitoringIntervalLevel,
                                  getUnrestricted = getUnrestricted, maxN = maxN, lagOutcomeN = lagOutcomeN,
                                  mc.cores = cores)

  } else if(socket==TRUE){
    # Works on Mac and Windows; only slightly slower
    cl             <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    mcmcEOS <- parallel::parLapply(cl, mcmcMonitoring, sgpvAMrulesSingle,
                                   waitTime                 = waitTime,
                                   lookSteps                = lookSteps,
                                   kSteps                   = kSteps,
                                   maxAlertSteps            = maxAlertSteps,
                                   monitoringIntervalLevel  = monitoringIntervalLevel,
                                   getUnrestricted = getUnrestricted, maxN = maxN, lagOutcomeN = lagOutcomeN)
  } else {

    mcmcEOS <- lapply(mcmcMonitoring, sgpvAMrulesSingle,
                      waitTime                 = waitTime,
                      lookSteps                = lookSteps,
                      kSteps                   = kSteps,
                      maxAlertSteps            = maxAlertSteps,
                      monitoringIntervalLevel  = monitoringIntervalLevel,
                      getUnrestricted = getUnrestricted, maxN = maxN, lagOutcomeN = lagOutcomeN)
  }

  # Return End of Study
  return(simplify2array(mcmcEOS))
}
