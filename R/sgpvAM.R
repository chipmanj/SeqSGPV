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
#' \code{lagOutcomeN} indicates the number of units that are yet being followed.  A study may stop with unrestricted sample size
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
#' @param maxN Total enrolled patients equals maxN observed patients plus lagOutcomeN.
#' @param lagOutcomeN Total enrolled patients equals MaxN observed patients plus lagOutcomeN.  lagOutcomeN are number of observations enrolled but awaiting to observe outcome.
#' @param monitoringIntervalLevel Traditional (1-alpha) used in monitoring intervals.
#' @param printProgress Prints when adding more data for mcmc replicates to have sufficient observations to monitor until a conclusion.  Defaults to TRUE.
#' @param outData Returns the mcmc generated data.  This can result in an out object with large memory.  Yet, with location shift data, can be re-used to obtain operating characteristics of shifted effects.
#' @param getECDF Returns the ECDF of sample size and bias for each wait width and number of steps before affirming end of study.
#' @param cores Number of cores used in parallel computing.  The default (NULL) does not run on parallel cores.
#' @param fork Fork clustering, works on POSIX systems (Mac, Linux, Unix, BSD) and not Windows.  Defaults to TRUE.
#' @param socket Socket clustering.  Defaults to TRUE yet only applies if FORK = FALSE.
#'
#' @return \enumerate{
#'             \item [mcmcMonitoring] The monitored generated / pre-generated datasets
#'             \item }
#' mcmcMonitoring the monitored generated / pre-generated datasets, (2) mcmcEndOfStudy: For each waitWidth and alertK, the average sample size, bias, mse, and errors plus the ecdf of sample size and bias, and (3) the function inputs
#' @examples

sgpvAM <- function(mcmcData=NULL, nreps,
                   waitWidths,
                   dataGeneration=NULL,   dataGenArgs,
                   effectGeneration=NULL, effectGenArgs,
                   modelFit,
                   pointNull, deltaL2, deltaL1, deltaG1, deltaG2,
                   lookSteps=1,
                   kSteps=10,
                   maxAlertSteps=100,
                   maxN=NULL, lagOutcomeN=0,
                   monitoringIntervalLevel = 0.05, printProgress=TRUE, outData = TRUE,
                   getECDF=TRUE,
                   cores            = NULL,
                   fork             = TRUE,
                   socket           = TRUE){


  # 0 Checks
  # If normal data and no standard deviation provided, set sd to 1
  if(class(modelFit)=="normal"){
    if(is.null(dataGenArgs$sd)) dataGenArgs$sd <- 1
  }

  if(is.null(cores)) cores <- detectCores()

  # 1 collect list of simulated data
  if(is.null(mcmcData)){
         mcmcMonitoring <- amData(nreps = nreps,
                                  monitoringIntervalLevel = monitoringIntervalLevel,
                                  dataGeneration   = dataGeneration,   dataGenArgs   = dataGenArgs,
                                  effectGeneration = effectGeneration, effectGenArgs = effectGenArgs,
                                  modelFit         = modelFit,
                                  cores            = cores,
                                  fork             = fork,
                                  socket           = socket)

  } else mcmcMonitoring <- mcmcData

  # 2 Add stats (bias, rejPN, cover, sgpvNonTrival, sgpvFutility)
  mcmcMonitoring <- lapply(mcmcMonitoring,
                           addStats,
                           pointNull = pointNull,
                           deltaL2 = deltaL2, deltaL1=deltaL1, deltaG1=deltaG1, deltaG2=deltaG2)


  # 3 Make sure all generated simulations will continue until completion
  #   - Look for stability of sgpv for last set of maxAlert patients
  getMore      <- unlist(lapply(mcmcMonitoring, sgpvAM::mcmcMonitoringEnoughCheck,
                                maxAlertSteps = maxAlertSteps,
                                lagOutcomeN   = lagOutcomeN,
                                minWW         = min(waitWidths),
                                sd            = dataGenArgs$sd))
  getMoreWhich <- which(getMore > 0)
  getMoreWhich

  if( !is.null(mcmcData) & sum(getMore) > 0 & is.null(dataGeneration) & is.null(effectGeneration)){

    stop("Provided mcmcData needs more observations to ensure study completes with unrestricted n\n
         Generate on own or provide data generation inputs.")

  } else {

    while(sum(getMore) > 0){

      if(printProgress==TRUE) print(paste("Adding up to", max(getMore), "observations for unrestricted sample size monitoring."))

      mcmcMonitoring <- amDataGetMore(insufficients    = getMoreWhich,
                                      existingDataList = mcmcMonitoring, getMore = getMore,
                                      monitoringIntervalLevel = monitoringIntervalLevel,
                                      dataGeneration   = dataGeneration,   dataGenArgs   = dataGenArgs,
                                      effectGeneration = effectGeneration, effectGenArgs = effectGenArgs,
                                      pointNull = pointNull,
                                      deltaL2 = deltaL2, deltaL1 = deltaL1, deltaG1 = deltaG1, deltaG2 = deltaG2,
                                      modelFit         = modelFit,
                                      cores            = cores,
                                      fork             = fork,
                                      socket           = socket)

      # Continue checking until all datasets have sufficient n
      getMore      <- unlist(lapply(mcmcMonitoring, mcmcMonitoringEnoughCheck,
                                    maxAlertSteps = maxAlertSteps,
                                    lagOutcomeN   = lagOutcomeN,
                                    minWW         = min(waitWidths),
                                    sd            = dataGenArgs$sd))
      getMoreWhich <- which(getMore > 0)

    }
  }


  # 4 adaptively monitor simulated data across multiple burn ins
  # 5 aggregate simulated data
  #   average performance and mse
  #   ecdf of n and bias
  mcmcEndOfStudyAve <- list()
  mcmcEndOfStudyVar <- list()
  mcmcEndOfStudy    <- list()
  for (w in 1:length(waitWidths)){

    ww <- waitWidths[w]

    # 4 adaptively monitor simulated data across multiple burn ins
    mcmcEOS <- sgpvAMrules(mcmcMonitoring=mcmcMonitoring,
                           waitWidth                = ww,
                           sd                       = dataGenArgs$sd,
                           lookSteps                = lookSteps,
                           kSteps                   = kSteps,
                           maxAlertSteps            = maxAlertSteps,
                           monitoringIntervalLevel  = monitoringIntervalLevel,
                           maxN = maxN, lagOutcomeN = lagOutcomeN)


    # 5 aggregate simulated data
    #   average performance and mse
    #   ecdf of n and bias
    ooAve <- plyr::aaply(mcmcEOS, .margins = c(1,2), .fun = mean)
    ooVar <- plyr::aaply(mcmcEOS, .margins = c(1,2), .fun = var )

    mcmcEndOfStudyAve                <- cbind(ooAve, mse = ooVar[,"bias"] + ooAve[,"bias"]^2)
    if(getECDF==TRUE){
      mcmcEndOfStudyEcdfSize           <- apply(mcmcEOS[,"n",], 1, ecdf)
      mcmcEndOfStudyEcdfBias           <- apply(mcmcEOS[,"bias",], 1, ecdf)
      names(mcmcEndOfStudyEcdfSize)    <- paste0("alertK_",mcmcEOS[,"alertK",1])
      names(mcmcEndOfStudyEcdfBias)    <- paste0("alertK_",mcmcEOS[,"alertK",1])
    } else {
      mcmcEndOfStudyEcdfSize <- NULL
      mcmcEndOfStudyEcdfBias <- NULL
    }

    mcmcEndOfStudy[[paste0("width_",ww)]] <-
      list(mcmcEndOfStudyAve      = mcmcEndOfStudyAve,
           mcmcEndOfStudyEcdfSize = mcmcEndOfStudyEcdfSize,
           mcmcEndOfStudyEcdfBias = mcmcEndOfStudyEcdfBias)

  }

  # Indicate whether to keep generated data
  if(outData==FALSE) {
    mcmcMonitoring <- NULL
    mcmcData       <- NULL
  }

  out <- list(mcmcMonitoring = mcmcMonitoring,
              mcmcEndOfStudy = mcmcEndOfStudy,
              inputs = list(mcmcData         = mcmcData,           nreps = nreps,
                            maxAlertSteps    = maxAlertSteps,  lookSteps = lookSteps, kSteps=kSteps,
                            waitWidths       = waitWidths,
                            dataGeneration   = dataGeneration,     dataGenArgs = dataGenArgs,
                            effectGeneration = effectGeneration,
                            modelFit         = modelFit,
                            pointNull        = pointNull,
                            deltaL2          = deltaL2, deltaL1 = deltaL1,
                            deltaG1          = deltaG1, deltaG2 = deltaG2,
                            maxN             = maxN,
                            lagOutcomeN      = lagOutcomeN,
                            monitoringIntervalLevel = monitoringIntervalLevel,
                            outData          = outData,
                            getECDF          = getECDF,
                            cores            = cores,
                            fork             = fork,
                            socket           = socket))


  # Set class
  class(out) <- append(class(out),"sgpvAM")

  return(out)

}
