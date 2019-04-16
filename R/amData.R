#' amData
#'
#' Generate mcmc simulations of adaptive monitoring with parallel computing
#'
#' @export
amData <- function(nreps, fork=TRUE, socket = TRUE, cores = detectCores(), ...){

  if(fork==TRUE){
    # Only works on POSIX systems (Mac, Linux, Unix, BSD) and not Windows.
    mcmcMonitoring <- parallel::mclapply(1:nreps, amDataSingle, ... , mc.cores = cores)
  } else if(socket==TRUE){
    # Works on Mac and Windows; only slightly slower
    cl             <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    mcmcMonitoring <- parallel::parLapply(cl, 1:nreps, amDataSingle, ...)
  } else {
    mcmcMonitoring <- plyr::rlply(.n = nreps, .expr = { amDataSingle( ... ) })
  }

  return(mcmcMonitoring)
}
