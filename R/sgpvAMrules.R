#' amSgpvRules
#'
#' Generate mcmc simulations of adaptive monitoring with parallel computing
#'
#' @export
#'
#'
sgpvAMrules <- function(mcmcMonitoring, fork=TRUE, socket = TRUE, cores = detectCores(),
                        waitWidth,
                        sd,
                        waitEmpirical,
                        minWaitN,
                        lookSteps,
                        kSteps,
                        maxAlertSteps,
                        monitoringIntervalLevel,
                        maxN, lagOutcomeN){

  if(fork==TRUE){
    # Only works on POSIX systems (Mac, Linux, Unix, BSD) and not Windows.
    mcmcEOS <- parallel::mclapply(mcmcMonitoring, sgpvAMrulesSingle,
                                  waitWidth                = waitWidth,
                                  sd                       = sd,
                                  waitEmpirical            = waitEmpirical,
                                  minWaitN                 = minWaitN,
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
                                   sd                       = sd,
                                   waitEmpirical            = waitEmpirical,
                                   minWaitN                 = minWaitN,
                                   lookSteps                = lookSteps,
                                   kSteps                   = kSteps,
                                   maxAlertSteps            = maxAlertSteps,
                                   monitoringIntervalLevel  = monitoringIntervalLevel,
                                   maxN = maxN, lagOutcomeN = lagOutcomeN)
  } else {

    mcmcEOS <- lapply(mcmcMonitoring, sgpvAMrulesSingle,
                      waitWidth                = waitWidth,
                      sd                       = sd,
                      waitEmpirical            = waitEmpirical,
                      minWaitN                 = minWaitN,
                      lookSteps                = lookSteps,
                      kSteps                   = kSteps,
                      maxAlertSteps            = maxAlertSteps,
                      monitoringIntervalLevel  = monitoringIntervalLevel,
                      maxN = maxN, lagOutcomeN = lagOutcomeN)
  }

  # Return End of Study
  return(simplify2array(mcmcEOS))
}

# system.time(test3 <- sgpvAMrules(mcmcMonitoring=mcmcMonitoring,
#                      waitWidth                = ww,
#                      lookSteps                = lookSteps,
#                      kSteps                   = kSteps,
#                      maxAlertSteps            = maxAlertSteps,
#                      monitoringIntervalLevel  = monitoringIntervalLevel,
#                      maxN = maxN, lagOutcomeN = lagOutcomeN))
#
# system.time(mcmcEOS <- simplify2array(lapply(mcmcMonitoring, sgpvAMrulesSingle,
#                                  waitWidth               = ww,
#                                  lookSteps               = lookSteps,
#                                  kSteps                  = kSteps,
#                                  maxAlertSteps           = maxAlertSteps,
#                                  monitoringIntervalLevel = monitoringIntervalLevel,
#                                  maxN = maxN, lagOutcomeN = lagOutcomeN)))
