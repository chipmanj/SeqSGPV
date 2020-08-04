#' Adaptive Monitoring with Second Generation p-Value
#'
#' The Second Generation p-Value (SGPV) incorporates pre-specified clinical context
#' to monitor studies until either ruling out trivial or impactful effects.  Operating
#' characteristics may be determined via.  This function allows the user to use either
#' their own generated data or to generate their own data from any of R's random data
#' generators.
#'
#' Data are either provided from one's own simulations or via monte-carlo simulations of data, using \code{dataGeneration},
#' \code{dataGenArgs}, \code{effectGeneration}, and \code{effectGenArgs}.  When one's own data are provided, it must be in matrix
#' format with at least the columns: 'theta', 'est', 'lo' (lower bound), and 'up' (upper bound).  Otherwise, data are generated
#' until the final \code{maxAlertSteps} observations affirm a stopping rule under unrestricted sample size.  Without \code{dataGeneration},
#' \code{dataGenArgs}, \code{effectGeneration}, and \code{effectGenArgs} no further data may be generated.  If pre-generated data
#' are insufficient for unrestricted monitoring, the user is notified to generate more data.
#'
#' Under working knowledge of the data / clinical context, indicate the point null and clinical guideposts (2 guideposts for one-sided
#' investigations and 4 guideposts for two-sided investigations).  For a one-sided investigation, leave both of the lower (or upper)
#' deltas (ie \code{deltaL2} \code{deltaL1} or \code{deltaG1} \code{deltaG2}) set to NA.
#'
#' To reduce false discoveries, the study waits until the interval of interest (currently coded as the Confidence Interval)
#' achieves an expected maximum width.  The expected width assumes the standard deviation provided in \code{dataGenArgs}.
#' Assuming larger variability than the truth yields a longer wait time and more conservative error.
#'
#' Acruing data are monitored by intervals as specified with \code{modelFit} which may be user defined.
#'
#' Monitoring rules: An alert is raised when the effect, theta, is ruled to be not trivial or not meaningful (using sgpv).
#' To reduce bias, increase coverage, and improve operating characteristics in general, the study may require the same alert
#' be raised k observations later.  To investigate the ideal choice of k, user may specify \code{lookSteps} up until
#' \code{maxAlertSteps}.
#'
#' While data are simulated to allow unrestricted sample size, practically the study may feasibly allow up to \code{maxN} units
#' to accrue.  Further, outcomes are simulated to occur simultaneously.  However, there may be a lag time in Outcomes.
#' \code{lagN} indicates the number of units that are yet being followed.  A study may stop with unrestricted sample size
#' yet may still have a lag on observing additionally accrued units.
#'
#' Note: User may choose to save the generated or pre-specified data.  This can be substantive in size.  However, when
#' generating data from a location-shift family, it may be used to fully explore a range of shifted intervention effects (theta).
#'
#'
#'
#' @param mcmcData Previously generated data.  Default (NULL) uses mcmc generation inputs to generate new or additional data.
#' @param nreps Number of mcmc replicates to generate
#' @param waitWidths Wait time, in terms of Confidence Interval width, before monitoring data.  Assuming the data generation standard deviation, waits until appropriate n to achieve confidence interval width.
#' @param dataGeneration Function (such as rnorm) to generate outcomes.
#' @param dataGenArgs Arguments for dataGeneration function.  This includes, in the least, 'n' observations to generate.  If 'n' is insufficient for unrestricted adaptive monitoring, addtional data will be generated.
#' @param effectGeneration Function (such as rnorm) or fixed value to generate intervention effect (theta).
#' @param effectGenArgs Arguments for effectGeneration function (if any)
#' @param modelFit A existing or user-defined function to obtain intervals.  Two existing functions are provided: 1) lmCI which obtains a confidence interval from a linear model and has class 'normal' indicating normal data and 2) lrCI obtains Wald Confidence Interval from logistic regression model and has class 'binomial' indicating binomial data.
#' @param pointNull Point null.
#' @param deltaL2 Clinical guidepost less than and furthest from point null.
#' @param deltaL1 Clinical guidepost less than and closest to point null.
#' @param deltaG1 Clinical guidepost greater than and closest to point null.
#' @param deltaG2 Clinical guidepost greater than and furthest from point null.
#' @param lookSteps The frequency data are observed (defaults to 1 -- fully sequential).
#' @param kSteps Affirmation steps to consider range from 0 to maxAlertSteps by kSteps.
#' @param maxAlertSteps Maximum number of steps before affirming an alert.
#' @param maxN Total enrolled patients equals maxN observed patients plus lagN.
#' @param lagN Total enrolled patients equals MaxN observed patients plus lagN.  lagN are number of observations enrolled but awaiting to observe outcome.
#' @param monitoringIntervalLevel Traditional (1-alpha) used in monitoring intervals.
#' @param printProgress Prints when adding more data for mcmc replicates to have sufficient observations to monitor until a conclusion.  Defaults to TRUE.
#' @param outData Returns the mcmc generated data.  This can result in an out object with large memory.  Yet, with location shift data, can be re-used to obtain operating characteristics of shifted effects.
#' @param getECDF Returns the ECDF of sample size and bias for each wait width and number of steps before affirming end of study.
#' @param cores Number of cores used in parallel computing.  The default (NULL) does not run on parallel cores.
#' @param fork Fork clustering, works on POSIX systems (Mac, Linux, Unix, BSD) and not Windows.  Defaults to TRUE.
#' @param socket Socket clustering.  Defaults to TRUE yet only applies if FORK = FALSE.
#'
#' @examples
#' #' devtools::install_github("chipmanj/sgpvAM")
#' library(sgpvAM)
#'
#' @export

sgpvAM <- function(mcmcData = NULL,
                   nreps,
                   dataGeneration   = rnorm,   dataGenArgs = list(n=200),
                   effectGeneration = 0, effectGenArgs=NULL,  effectScale  = "identity",
                   randomize        = FALSE,
                   effectPN         = 0,
                   PRISM,
                   modelFit,
                   miLevel          = 0.05,
                   designLooks      = NULL,
                   wait             = 4,
                   steps            = 1,
                   affirm           = 0,
                   lag              = 0,
                   maxN             = NA,
                   getUnrestricted  = TRUE,
                   printProgress    = TRUE,
                   outData          = TRUE,
                   getECDF          = TRUE,
                   cores            = NULL,
                   fork             = TRUE,
                   socket           = TRUE){



  # 0.1 PRISM parameters
  if(missing(PRISM)){
    stop("Must specify PRISM")
  } else {
    deltaL2 <- PRISM[["deltaL2"]]
    deltaL1 <- PRISM[["deltaL1"]]
    deltaG1 <- PRISM[["deltaG1"]]
    deltaG2 <- PRISM[["deltaG2"]]
  }

  # 0.2 Design parameters
  if(is.null(designLooks)){
    designLooks        <- expand.grid(wait,steps,affirm,lag,maxN)
    names(designLooks) <- c("W","S","A","L","N")
  }
  designLooksLabels    <- paste0("W",designLooks[,"W"],
                                 "_S",designLooks[,"S"],
                                 "_A",designLooks[,"A"],
                                 "_L",designLooks[,"L"],
                                 "_N",designLooks[,"N"])


  # 0.3 set cores for parallel computing
  if(is.null(cores)) cores <- parallel::detectCores()


  # 0.4 if POSIX systems (Mac, Linux, Unix, BSD) use mcapply.  For windows use parLapply
  os <- Sys.info()["sysname"]



  # 1 collect list of simulated data
  if(is.null(mcmcData)){
    if(printProgress) cat("\rGenerating simulated data and sequential estimates")
         mcmcMonitoring <- amData(nreps            = nreps,
                                  miLevel          = miLevel,
                                  dataGeneration   = dataGeneration,   dataGenArgs   = dataGenArgs,
                                  effectGeneration = effectGeneration, effectGenArgs = effectGenArgs,
                                  effectScale      = effectScale,
                                  randomize        = randomize,
                                  modelFit         = modelFit,
                                  cores            = cores,
                                  os               = os,
                                  fork             = fork,
                                  socket           = socket)

  } else mcmcMonitoring <- mcmcData

  # 2 Add stats (bias, rejPN, cover, sgpvNotROPE, sgpvNotROME)
  if(printProgress) cat("\rAdding monitoring statistics: bias, rejPN, cover, sgpvROPE, sgpvROME")
  mcmcMonitoring <- lapply(mcmcMonitoring,
                           addStats,
                           effectPN = effectPN,
                           deltaL2 = deltaL2, deltaL1=deltaL1, deltaG1=deltaG1, deltaG2=deltaG2)


  # 3 Make sure all generated simulations will continue until completion
  #   - Look for stability of sgpv for last set of maxAlert patients
  if( any(designLooks[,"N"]==Inf) | dataGenArgs$n < max(designLooks[designLooks!=Inf,"N"])){

    if(printProgress) cat("\rEnsuring simulations with unrestricted sample size each continue to completion")

    getMore      <- unlist(lapply(mcmcMonitoring, mcmcMonitoringEnoughCheck, designLooks=designLooks))
    getMoreWhich <- which(getMore > 0)
    getMoreWhich


    if( !is.null(mcmcData)      & sum(getMore) > 0          &
        is.null(dataGeneration) & is.null(effectGeneration) ){

      stop("Provided mcmcData needs more observations to ensure study completes with unrestricted n\n
           Generate on own or provide data generation inputs.")

    } else {

      iter <- 1
      while(sum(getMore) > 0){

        if(printProgress==TRUE) cat(paste("\r",iter,". Adding up to", max(getMore), "observations for unrestricted sample size monitoring.       "))

        mcmcMonitoring <- amDataGetMore(insufficients    = getMoreWhich,
                                        existingDataList = mcmcMonitoring,   getMore       = getMore,
                                        miLevel          = miLevel,
                                        dataGeneration   = dataGeneration,   dataGenArgs   = dataGenArgs,
                                        effectGeneration = effectGeneration, effectGenArgs = effectGenArgs,
                                        effectScale      = effectScale,
                                        randomize        = randomize,
                                        effectPN         = effectPN,
                                        deltaL2 = deltaL2, deltaL1 = deltaL1, deltaG1 = deltaG1, deltaG2 = deltaG2,
                                        modelFit         = modelFit,
                                        cores            = cores,
                                        os               = os,
                                        fork             = fork,
                                        socket           = socket)

        # Continue checking until all datasets have sufficient n
        getMore      <- unlist(lapply(mcmcMonitoring, mcmcMonitoringEnoughCheck, designLooks=designLooks))
        getMoreWhich <- which(getMore > 0)

        iter <- iter + 1
      }
    }
  }


  # 4 adaptively monitor simulated data across designLooks parameters
  mcmcEndOfStudyAve <- list()
  mcmcEndOfStudyVar <- list()
  mcmcEndOfStudy    <- list()
  mcmcECDFs         <- list()

  if(printProgress) cat("\rAdaptively monitoring for each wait time                                                   ")

  mcmcEOS <- sgpvAMrules(mcmcMonitoring  = mcmcMonitoring,
                        designLooks     = designLooks,
                        os              = os,
                        fork            = fork,
                        socket          = socket,
                        getUnrestricted = getUnrestricted)


  # 5 aggregate simulated data
  #   average performance and mse
  #   ecdf of n and bias
  ooAve <- plyr::aaply(mcmcEOS, .margins = c(1,2), .fun = mean)
  ooVar <- plyr::aaply(mcmcEOS, .margins = c(1,2), .fun = var )

  mcmcEndOfStudyAve   <- ooAve

  if(getECDF==TRUE){

    # Without lag
    mcmcECDFs$mcmcEndOfStudyEcdfSize           <- apply(mcmcEOS[,"n",],    1, ecdfDataReduction)
    mcmcECDFs$mcmcEndOfStudyEcdfBias           <- apply(mcmcEOS[,"bias",], 1, ecdfDataReduction)
    names(mcmcECDFs$mcmcEndOfStudyEcdfSize)    <- designLooksLabels
    names(mcmcECDFs$mcmcEndOfStudyEcdfBias)    <- designLooksLabels

    # With lag
    mcmcECDFs$mcmcEndOfStudyEcdfSizeLag        <- apply(mcmcEOS[,"lag.n",],    1, ecdfDataReduction)
    mcmcECDFs$mcmcEndOfStudyEcdfBiasLag        <- apply(mcmcEOS[,"lag.bias",], 1, ecdfDataReduction)
    names(mcmcECDFs$mcmcEndOfStudyEcdfSizeLag) <- designLooksLabels
    names(mcmcECDFs$mcmcEndOfStudyEcdfBiasLag) <- designLooksLabels

  }




  # 6 output object

  # Indicate whether to keep generated data
  if(outData==FALSE) {
    mcmcMonitoring <- NULL
    mcmcData       <- NULL
  }

  out <- list(mcmcMonitoring = mcmcMonitoring,
              mcmcOC         = mcmcEndOfStudyAve,
              mcmcECDFs      = mcmcECDFs,
              inputs         = lapply(match.call.defaults()[-1], eval))


  # Set class
  class(out) <- append("sgpvAM",class(out))

  # Clear print progress
  if(printProgress) {
    cat("\r                                                                                    ")
    cat("\n")
  }

  return(out)

}
