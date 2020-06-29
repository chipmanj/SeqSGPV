#' amSgpvRules
#'
#' Generate mcmc simulations of adaptive monitoring with parallel computing
#'
#' @export
#'
#'
sgpvAMrules <- function(mcmcMonitoring, os, fork=TRUE, socket = TRUE, cores = detectCores(),
                        waitJ,
                        lookS,
                        affirmK,
                        miLevel,
                        getUnrestricted, maxN, lagN){


  if(fork==TRUE & os!="Windows"){
    # Only works on POSIX systems (Mac, Linux, Unix, BSD) and not Windows.
    mcmcEOS <- parallel::mclapply(mcmcMonitoring, sgpvAMrulesSingle,
                                  waitJ           = waitJ,
                                  lookS           = lookS,
                                  affirmK         = affirmK,
                                  getUnrestricted = getUnrestricted, maxN = maxN, lagN = lagN,
                                  mc.cores        = cores)

  } else if(socket==TRUE){
    # Works on Mac and Windows; only slightly slower
    cl             <- parallel::makeCluster(cores)
    clusterCall(cl, function() library(sgpvAM))
    on.exit(parallel::stopCluster(cl), add = TRUE)
    mcmcEOS <- parallel::parLapply(cl, mcmcMonitoring, sgpvAMrulesSingle,
                                   waitJ           = waitJ,
                                   lookS           = lookS,
                                   affirmK         = affirmK,
                                   getUnrestricted = getUnrestricted, maxN = maxN, lagN = lagN)
  } else {

    mcmcEOS <- lapply(mcmcMonitoring, sgpvAMrulesSingle,
                      waitJ    = waitJ,
                      lookS    = lookS,
                      affirmK  = affirmK,
                      getUnrestricted = getUnrestricted, maxN = maxN, lagN = lagN)
  }

  # Return End of Study
  return(simplify2array(mcmcEOS))
}
