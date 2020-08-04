#' amSgpvRules
#'
#' Generate mcmc simulations of adaptive monitoring with parallel computing
#'
#' @export
#'
#'
sgpvAMrules <- function(mcmcMonitoring, os, fork=TRUE, socket = TRUE, cores = detectCores(),
                        designLooks,
                        miLevel,
                        getUnrestricted){


  if(fork==TRUE & os!="Windows"){
    # Only works on POSIX systems (Mac, Linux, Unix, BSD) and not Windows.
    mcmcEOS <- parallel::mclapply(mcmcMonitoring, sgpvAMrulesSingle,
                                  designLooks     = designLooks,
                                  getUnrestricted = getUnrestricted,
                                  mc.cores        = cores)

  } else if(socket==TRUE){
    # Works on Mac and Windows; only slightly slower
    cl             <- parallel::makeCluster(cores)
    clusterCall(cl, function() library(sgpvAM))
    on.exit(parallel::stopCluster(cl), add = TRUE)
    mcmcEOS <- parallel::parLapply(cl, mcmcMonitoring, sgpvAMrulesSingle,
                                       designLooks     = designLooks,
                                       getUnrestricted = getUnrestricted)
  } else {

    mcmcEOS <- lapply(mcmcMonitoring, sgpvAMrulesSingle,
                      designLooks     = designLooks,
                      getUnrestricted = getUnrestricted)
  }

  # Return End of Study
  return(simplify2array(mcmcEOS))
}
