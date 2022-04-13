#' @title Sequential Monitoring of Second Generation P-Value (SeqSGPV)
#'
#' @description The Second Generation p-Value (SGPV) incorporates pre-specified clinical context
#' to monitor studies until either ruling out trivial or impactful effects.  Operating
#' characteristics may be determined via simulation.  This function allows the user to use either
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
#' To reduce false discoveries, the study waits until the interval of interest  achieves an expected maximum width.  The expected width
#' assumes the standard deviation provided in \code{dataGenArgs}. Assuming larger variability than the truth yields a longer wait time and
#' more conservative error.
#'
#' Accruing data are monitored by intervals as specified with \code{modelFit} and \code{modelFitArgs} which may be user defined.
#'
#' Monitoring rules: An alert is raised when the effect, theta, is ruled to be not trivial or not meaningful (using sgpv).
#' To reduce bias, increase coverage, and improve operating characteristics in general, the study may require the same alert
#' be raised A observations later (or prior).  To investigate the ideal choice of A, user may specify \code{lookSteps} up until
#' \code{maxAlertSteps}.
#'
#'
#' Note: User may choose to save the generated or pre-specified data.  This can be substantive in size.  However, when
#' generating data from a location-shift family, it may be used to fully explore a range of shifted intervention effects (theta).
#'
#'
#'
#' @param mcmcData Previously generated data.  Default (NULL) uses mcmc generation inputs to generate new or additional data.
#' @param nreps Number of mcmc replicates to generate
#' @param dataGeneration Function (such as rnorm) to generate outcomes.
#' @param dataGenArgs Arguments for dataGeneration function.  This includes, in the least, 'n' observations to generate.  If 'n' is insufficient for unrestricted adaptive monitoring, additional data will be generated.
#' @param effectGeneration Function (such as rnorm) or fixed value to generate intervention effect (theta).
#' @param effectGenArgs Arguments for effectGeneration function (if any)
#' @param effectScale Required transformation effect to generate data under dataGeneration. Can be any of 'identity', 'log', 'oddsratio', or 'or'.
#' @param effectPN Point null of null hypothesis.  For 1-sided hypotheses, it is the boundary of the null of greatest interest (ex: H0: mu <= 0, effectPN = 0).
#' @param null Any of "two.sided", "greater", "less" to specify the type of null hypothesis.
#' @param allocation A vector of the allocation ratio.  For single arm trials, allocation is 1.  For 2:3 allocation, set to c(2,3).
#' @param PRISM A list with elements deltaL2, deltaL1, deltaG1, deltaG2.  The lower, or upper, deltas may be set to NA to assess one- or two-sided hypotheses
#' @param modelFit Function to obtain interval from data. Available functions for 1 and 2 armed trials include 'lmCI' and 'lrCI' which are respectively confidence intervals of a linear model and of a logistic regression (CI of odds ratio).  Also, binomCI, uses binom::binom.confit for credible and confidence intervals of 1-armed trials with bernoulli outcomes.
#' @param modelFitArgs List of inputs for modelFit function
#' @param wait Vector of possible wait times (W) before monitoring.
#' @param steps Vector of the number of observations (S) in between monitoring assessments.
#' @param affirm Vector of the number of observations required for affirming a stopping rule (A) in between raising an alert of stopping.
#' @param lag Vector of the number of delayed outcomes between enrolling a single observation and observing its outcome.
#' @param N Vector of maximum sample size. Can be set as Inf (indefinite).
#' @param printProgress Prints when adding more data for mcmc replicates to have sufficient observations to monitor until a conclusion.  Defaults to TRUE.
#' @param outData Returns the mcmc generated data.  This can result in an out object with large memory.  Yet, with location shift data, can be re-used to obtain operating characteristics of shifted effects.
#' @param getECDF Returns the ECDF of sample size and bias for each wait width and number of steps before affirming end of study.
#' @param cores Number of cores used in parallel computing.  The default (NULL) does not run on parallel cores.
#' @param fork Fork clustering, works on POSIX systems (Mac, Linux, Unix, BSD) and not Windows.  Defaults to TRUE.
#' @param socket Socket clustering.  Defaults to TRUE yet only applies if FORK = FALSE.
#'
#' @return If effectGeneration is single value, returns list of class SeqSGPV which includes elements: mcmcMonitoring replicates (all generated data), mcmcOC (average operating characterstics for each combination of monitoring frequencies), mcmcECDFs (ecdf functions of sample size and bias for each combination of monitoring frequencies), and inputs (inputs into SeqSGPV).
#' @return If effectGeneration is a function, returns a list of class SeqSGPVre which includes elements: mcmcEOS (the final observation of each replicate for each combination of monitoring frequencies) and inputs (inputs into SeqSGPV).
#'
#' @export
SeqSGPV <- function(mcmcData = NULL,
                    nreps,
                    dataGeneration   = rnorm,   dataGenArgs = list(n=200),
                    effectGeneration = 0, effectGenArgs=NULL,  effectScale  = "identity",
                    effectPN         = 0,
                    allocation,
                    null             = "two.sided",
                    PRISM,
                    modelFit,
                    modelFitArgs     = NULL,
                    wait             = 4,
                    steps            = 1,
                    affirm           = 0,
                    lag              = 0,
                    N                = NA,
                    printProgress    = TRUE,
                    outData          = TRUE,
                    getECDF          = TRUE,
                    cores            = NULL,
                    fork             = TRUE,
                    socket           = TRUE){


  # 0.0 Checks
  if(missing(null)){
    stop("set 'null' to one of: greater, less, or two.sided")
  } else if(null == "less" & any(is.na(PRISM[["deltaG1"]]), is.na(PRISM[["deltaG2"]]))){
    stop("if 'null == greater' set deltaL2 = NA and deltaL1 = NA")
  } else if(null == "greater" & any(is.na(PRISM[["deltaL1"]]), is.na(PRISM[["deltaL2"]]))){
    stop("if 'null == less' set deltaG1 = NA and deltaG2 = NA")
  }


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
  monitoringFrequency        <- expand.grid(wait,steps,affirm,lag,N)
  names(monitoringFrequency) <- c("W","S","A","L","N")
  monitoringFrequencyLabels  <- paste0("W", monitoringFrequency[,"W"],
                                       "_S",monitoringFrequency[,"S"],
                                       "_A",monitoringFrequency[,"A"],
                                       "_L",monitoringFrequency[,"L"],
                                       "_N",monitoringFrequency[,"N"])


  # 0.3 set cores for parallel computing
  if(is.null(cores)) cores <- parallel::detectCores()


  # 0.4 if POSIX systems (Mac, Linux, Unix, BSD) use mcapply.  For windows use parLapply
  os <- Sys.info()["sysname"]


  # 0.5 Denote randomization if multiple arms
  if(missing(allocation)){
    allocation <- 1
    cat("\nSingle arm trial is assumed. Change allocation input for randomized trial.\n")
  }

  if(length(allocation)>1) randomize <- TRUE else randomize <- FALSE


  # 0.6 Make sure A is coded as looking back retrospectively
  if(any(affirm<0)){
    stop("Set A as non-negative integers.  Affirmations look back A (non-negative) observations to affirm SGPV evidence.")
  }


  # 0.7 Neither effect generation nor effect pn can be 0 if effect scale is odds ratio
  if(toupper(effectScale) %in% c("OR","ODDSRATIO") & !is.function(effectGeneration)){
    if((effectGeneration == 0 | effectPN == 0)){
      stop("effectScale is set to 'OR' or 'ODDSRATIO' and either effectGeneration or effectPN = 0.")
    }
  }


  # 1 collect list of simulated data
  if(is.null(mcmcData)){
    if(printProgress) cat("\rGenerating simulated data and sequential estimates")
         mcmcMonitoring <- amData(nreps            = nreps,
                                  dataGeneration   = dataGeneration,   dataGenArgs   = dataGenArgs,
                                  effectGeneration = effectGeneration, effectGenArgs = effectGenArgs,
                                  effectScale      = effectScale,
                                  allocation       = allocation,
                                  randomize        = randomize,
                                  modelFit         = modelFit,
                                  modelFitArgs     = modelFitArgs,
                                  existingData     = NULL,
                                  cores            = cores,
                                  os               = os,
                                  fork             = fork,
                                  socket           = socket)


  } else mcmcMonitoring <- mcmcData

  # 2 Add stats (bias, rejPN, cover, sgpvNotROPE, sgpvNotROME)
  if(printProgress) cat("\rAdding monitoring statistics: bias, rejH0, cover, sgpvROPE, sgpvROME")
  mcmcMonitoring <- lapply(mcmcMonitoring,
                           addStats,
                           randomize = randomize,
                           effectPN  = effectPN,
                           null      = null,
                           deltaL2   = deltaL2, deltaL1=deltaL1, deltaG1=deltaG1, deltaG2=deltaG2)


  # 3 Make sure all generated simulations will continue until completion
  #   - Look for stability of sgpv for last set of maxAlert patients
  if( dataGenArgs$n < max(monitoringFrequency[,"N"]) ){

    if(printProgress) cat("\rEnsuring simulations with unrestricted sample size each continue to completion")

    getMore      <- unlist(lapply(mcmcMonitoring, mcmcMonitoringEnoughCheck, monitoringFrequency=monitoringFrequency))
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
                                        dataGeneration   = dataGeneration,   dataGenArgs   = dataGenArgs,
                                        effectGeneration = effectGeneration, effectGenArgs = effectGenArgs,
                                        effectScale      = effectScale,
                                        allocation       = allocation,
                                        randomize        = randomize,
                                        effectPN         = effectPN,
                                        null             = null,
                                        deltaL2 = deltaL2, deltaL1 = deltaL1, deltaG1 = deltaG1, deltaG2 = deltaG2,
                                        modelFit         = modelFit,
                                        modelFitArgs     = modelFitArgs,
                                        cores            = cores,
                                        os               = os,
                                        fork             = fork,
                                        socket           = socket)

        # Continue checking until all datasets have sufficient n
        getMore      <- unlist(lapply(mcmcMonitoring, mcmcMonitoringEnoughCheck, monitoringFrequency=monitoringFrequency))
        getMoreWhich <- which(getMore > 0)

        iter <- iter + 1
      }
    }
  }


  # 4 adaptively monitor simulated data across monitoringFrequency parameters
  mcmcEndOfStudyAve <- list()
  mcmcEndOfStudyVar <- list()
  mcmcEndOfStudy    <- list()
  mcmcECDFs         <- list()

  if(printProgress) cat("\rAdaptively monitoring for each wait time                                                   ")

  mcmcEOS <- SeqSGPVrules(mcmcMonitoring      = mcmcMonitoring,
                          monitoringFrequency = monitoringFrequency,
                          randomize           = randomize,
                          os                  = os,
                          fork                = fork,
                          socket              = socket)


  # 5 summarize data:
  # If conditioning on single effect
  #   - average performance
  #   - ecdf of n and bias
  # If effect generated randomly
  #   - obtain debiasing function (if requested)
  #   - return end of study observation (not averaged)
  if(!is.function(effectGeneration)){

    ooAve <- plyr::aaply(mcmcEOS, .margins = c(1,2), .fun = mean)
    ooVar <- plyr::aaply(mcmcEOS, .margins = c(1,2), .fun = var )


    if(length(monitoringFrequencyLabels)==1){

      ooAve <- t(as.matrix(ooAve, nrow=1))
      ooVar <- t(as.matrix(ooVar, nrow=1))

      mcmcEndOfStudyAve   <- ooAve

      if(getECDF==TRUE){

        mcmcECDFs$mcmcEndOfStudyEcdfN[[monitoringFrequencyLabels]]       <- ecdfDataReduction(mcmcEOS[,"n",])
        mcmcECDFs$mcmcEndOfStudyEcdfBias[[monitoringFrequencyLabels]]    <- ecdfDataReduction(mcmcEOS[,"bias",])

        mcmcECDFs$mcmcEndOfStudyEcdfNLag[[monitoringFrequencyLabels]]    <- ecdfDataReduction(mcmcEOS[,"lag.n",])
        mcmcECDFs$mcmcEndOfStudyEcdfBiasLag[[monitoringFrequencyLabels]] <- ecdfDataReduction(mcmcEOS[,"lag.bias",])

      }


    } else {

      mcmcEndOfStudyAve   <- ooAve

      if(getECDF==TRUE){

        mcmcECDFs$mcmcEndOfStudyEcdfN              <- apply(mcmcEOS[,"n",],    1, ecdfDataReduction)
        mcmcECDFs$mcmcEndOfStudyEcdfBias           <- apply(mcmcEOS[,"bias",], 1, ecdfDataReduction)
        names(mcmcECDFs$mcmcEndOfStudyEcdfN)       <- monitoringFrequencyLabels
        names(mcmcECDFs$mcmcEndOfStudyEcdfBias)    <- monitoringFrequencyLabels

        mcmcECDFs$mcmcEndOfStudyEcdfNLag           <- apply(mcmcEOS[,"lag.n",],    1, ecdfDataReduction)
        mcmcECDFs$mcmcEndOfStudyEcdfBiasLag        <- apply(mcmcEOS[,"lag.bias",], 1, ecdfDataReduction)
        names(mcmcECDFs$mcmcEndOfStudyEcdfNLag)    <- monitoringFrequencyLabels
        names(mcmcECDFs$mcmcEndOfStudyEcdfBiasLag) <- monitoringFrequencyLabels

      }

    }


  } else {

    # Save end of study simulations to list by monitoring frequency settings
    mcmcEOS        <- plyr::aaply(mcmcEOS, .margins = 1, .fun = function(x){ list(t(x)) } )
    names(mcmcEOS) <- monitoringFrequencyLabels

  }



  # 6 output object

  # Indicate whether to keep generated data
  if(outData==FALSE) {
    mcmcMonitoring <- NULL
    mcmcData       <- NULL
  }



  if(!is.function(effectGeneration)){
    out <- list(mcmcMonitoring = mcmcMonitoring,
                mcmcOC         = mcmcEndOfStudyAve,
                mcmcECDFs      = mcmcECDFs,
                inputs         = lapply(match.call.defaults()[-1], eval))

    # Set class denoting fixed effects
    class(out) <- append("SeqSGPV",class(out))
  } else {
    out <- list(mcmcEOS        = mcmcEOS,
                inputs         = lapply(match.call.defaults()[-1], eval))

    # Set class denoting random effects
    class(out) <- append("SeqSGPVre",class(out))
  }


  # Clear print progress
  if(printProgress) {
    cat("\r                                                                                    ")
    cat("\n")
  }

  return(out)

}
