#' amSgpvRules
#'
#' Generate mcmc simulations of adaptive monitoring with parallel computing
#'
#' @export
#'
#'
sgpvAMrules <- function(mcmcMonitoring, fork=TRUE, socket = TRUE, cores = detectCores(),
                        waitWidth,
                        sdY,
                        waitEmpirical,
                        minWaitN,
                        periphLo,
                        periphUp,
                        lookSteps,
                        kSteps,
                        maxAlertSteps,
                        monitoringIntervalLevel,
                        maxN, lagOutcomeN){

  if(fork==TRUE){
    # Only works on POSIX systems (Mac, Linux, Unix, BSD) and not Windows.
    mcmcEOS <- parallel::mclapply(mcmcMonitoring, sgpvAMrulesSingle,
                                  waitWidth                = waitWidth,
                                  sdY                      = sdY,
                                  waitEmpirical            = waitEmpirical,
                                  minWaitN                 = minWaitN,
                                  periphLo                 = periphLo,
                                  periphUp                 = periphUp,
                                  lookSteps                = lookSteps,
                                  kSteps                   = kSteps,
                                  maxAlertSteps            = maxAlertSteps,
                                  monitoringIntervalLevel  = monitoringIntervalLevel,
                                  maxN = maxN, lagOutcomeN = lagOutcomeN , mc.cores = cores)

  } else if(socket==TRUE){
    # Works on Mac and Windows; only slightly slower
    cl             <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    mcmcEOS <- parallel::parLapply(cl, mcmcMonitoring, sgpvAMrulesSingle,
                                   waitWidth                = waitWidth,
                                   sdY                      = sdY,
                                   waitEmpirical            = waitEmpirical,
                                   minWaitN                 = minWaitN,
                                   periphLo                 = periphLo,
                                   periphUp                 = periphUp,
                                   lookSteps                = lookSteps,
                                   kSteps                   = kSteps,
                                   maxAlertSteps            = maxAlertSteps,
                                   monitoringIntervalLevel  = monitoringIntervalLevel,
                                   maxN = maxN, lagOutcomeN = lagOutcomeN)
  } else {

    mcmcEOS <- lapply(mcmcMonitoring, sgpvAMrulesSingle,
                      waitWidth                = waitWidth,
                      sdY                      = sdY,
                      waitEmpirical            = waitEmpirical,
                      minWaitN                 = minWaitN,
                      periphLo                 = periphLo,
                      periphUp                 = periphUp,
                      lookSteps                = lookSteps,
                      kSteps                   = kSteps,
                      maxAlertSteps            = maxAlertSteps,
                      monitoringIntervalLevel  = monitoringIntervalLevel,
                      maxN = maxN, lagOutcomeN = lagOutcomeN)
  }

  # Return End of Study
  return(simplify2array(mcmcEOS))
}
