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
#' @param mcmcData Previously generated data
#' @param nreps Number of mcmc replicates to generate
#' @param waitWidths Width of Confidence Intervals, matching the standard deviation scale provided in dataGeneration
#' @param dataGeneration Function (such as rnorm) to generate outcomes
#' @param dataGenArgs Arguments for dataGeneration function.  This includes, in the least, 'n' observations to generate.
#' @param effectGeneration Function (such as rnorm) or fixed value to generate intervention effect (theta)
#' @param effectGenArgs Arguments for effectGeneration function (if any)
#' @param modelFit Function written to obtain intervals.  Two are provided: 1) lmCI which obtains a confidence interval from a linear model and has class 'normal' indicating normal data and 2) lrCI obtains Wald Confidence Interval from logistic regression model and has class 'binomial' indicating binomial data.
#' @param pointNull Point null
#' @param deltaL2 Clinical guidepost less than and furthest from point null
#' @param deltaL1 Clinical guidepost less than and closest to point null
#' @param deltaG1 Clinical guidepost greater than and closest to point null
#' @param deltaG2 Clinical guidepost greater than and furthest from point null
#' @param maxAlertSteps Maximum number of steps before affirming an alert
#' @param lookSteps The frequency data are observed (defaults to fully sequential)
#' @param maxN Maximum number of observations
#' @param lagOutcomeN Number of units collected while awaiting an outcome
#' @param monitoringIntervalLevel Traditional (1-alpha) used in monitoring intervals
#' @param outData T/F for returning the mcmcMonitoringData at end of call
#'
#' @return (1) mcmcMonitoring: the monitored generated / pre-generated datasets, (2) mcmcEndOfStudy: For each waitWidth and alertK, the average sample size, bias, mse, and errors plus the ecdf of sample size and bias, and (3) the function inputs
#' @examples
#' devtools::install_github("chipmanj/sgpvAM")
#' library(sgpvAM)

#' # Simulate AM trial
#' # Two-sided deltas
#' # defaults to sd = 1
#' am <-  sgpvAM(nreps            = 100,
#'               maxAlertSteps    = 100,       lookSteps = 1, waitWidths = seq(0.15, 0.6, by = 0.05),
#'               dataGeneration   = rnorm,   dataGenArgs = list(n=800),
#'               effectGeneration = 0,
#'               modelFit         = lmCI,
#'               pointNull = 0, deltaL2 = -0.5, deltaL1=-0.2, deltaG1=0.2, deltaG2=0.5,
#'               monitoringIntervalLevel = 0.05,
#'               maxN = 200, lagOutcomeN = 50)

#' amShifted <- locationShift(am, shiftedThetas = seq(-1,0,by=0.1))


#' # Unrestricted Sample Size
#' # Explore wait width (first with alert k specified at 50)
#' plot(amShifted, "n",                 alertK = 50,      xlim=c(-1,0),   ylim=c(0,801),    maxVary = 10)
#' plot(amShifted, "rejPN",             alertK = 50,      xlim=c(-1,0),   ylim=c(0,1),      maxVary = 10)
#' plot(amShifted, "rejPN",             alertK = 50,      xlim=c(-1,0),   ylim=c(0,1),      maxVary = 5)
#' plot(amShifted, "bias",              alertK = 50,      xlim=c(-1,0),   ylim=c(-0.1,0.1), maxVary = 5)
#' plot(amShifted, "mse",               alertK = 50,      xlim=c(-1,0),   ylim=c(0,0.1),    maxVary = 5)
#' plot(amShifted, "cover",             alertK = 50,      xlim=c(-1,0),   ylim=c(0.75,1),   maxVary = 5)
#' plot(amShifted, "stopInconclusive",  alertK = 50,      xlim=c(-1,0),   ylim=c(0,1),      maxVary = 5)
#' plot(amShifted, "stopNotImpactful",  alertK = 50,      xlim=c(-1,0),   ylim=c(0,1),      maxVary = 5)
#' plot(amShifted, "stopNotTrivial",    alertK = 50,      xlim=c(-1,0),   ylim=c(0,1),      maxVary = 5)

#' # Explore number of steps in affirmation step
#' plot(amShifted, "n",                 waitWidth = 0.3,  xlim=c(-1,0),   ylim=c(0,801),    maxVary = 5)
#' plot(amShifted, "bias",              waitWidth = 0.3,  xlim=c(-1,0),   ylim=c(-0.1,0.1), maxVary = 5)
#' plot(amShifted, "mse",               waitWidth = 0.3,  xlim=c(-1,0),   ylim=c(0,0.1),    maxVary = 5)
#' plot(amShifted, "cover",             waitWidth = 0.3,  xlim=c(-1,0),   ylim=c(0.75,1),   maxVary = 5)
#' plot(amShifted, "stopInconclusive",  waitWidth = 0.3,  xlim=c(-1,0),   ylim=c(0,1),      maxVary = 5)
#' plot(amShifted, "stopNotImpactful",  waitWidth = 0.3,  xlim=c(-1,0),   ylim=c(0,1),      maxVary = 5)
#' plot(amShifted, "stopNotTrivial",    waitWidth = 0.3,  xlim=c(-1,0),   ylim=c(0,1),      maxVary = 5)


#' # Unrestricted Sample Size with lag of remaining outcomes
#' plot(amShifted, "n",                 alertK = 50,     xlim=c(-1,0),   ylim=c(0,801), maxVary = 10, sizeRestrictions = "lag")
#' plot(amShifted, "stopInconclusive",  alertK = 50,     xlim=c(-1,0),   ylim=c(0,1),   maxVary = 10, sizeRestrictions = "lag")
#' plot(amShifted, "stopInconclusive",  waitWidth = 0.3, xlim=c(-1,0),   ylim=c(0,1),   maxVary = 10, sizeRestrictions = "lag")

#' # Max N with immediate outcomes
#' plot(amShifted, "n",     alertK = 50,     xlim=c(-1,0),   ylim=c(0,801), maxVary = 10, sizeRestrictions = "maxN")
#' plot(amShifted, "stopInconclusive",     alertK = 50,     xlim=c(-1,0),   ylim=c(0,1), maxVary = 10, sizeRestrictions = "lag")


#' summary(amShifted, alertK = 50, waitTime = 0.3, treatEffect = -0.5)
#' summary(amShifted)
#' summary(am)





#' # Simulate AM trial
#' # One-sided deltas
#' # Outcome has sd = 2
#' am2 <-  sgpvAM(nreps           = 100,
#'                maxAlertSteps    = 100,       lookSteps = 1, waitWidths = sort(c(seq(0.2, 1.2, by = 0.2),c(0.5,0.7))),
#'                dataGeneration   = rnorm,   dataGenArgs = list(n=800, sd=2),
#'                effectGeneration = 0,
#'                modelFit         = lmCI,
#'                pointNull = 0, deltaL2 = -1, deltaL1=-0.2, deltaG1=NA, deltaG2=NA,
#'                monitoringIntervalLevel = 0.05,
#'                maxN = 200, lagOutcomeN = 50)

#' amShifted2 <- locationShift(am2, shiftedThetas = seq(-2,0.8,by=0.4))



#' # Unrestricted Sample Size
#' # Explore wait width (first with alert k specified at 50)
#' plot(amShifted2, "n",                 alertK = 0,       xlim=c(-2,1),    ylim=c(0,500))
#' plot(amShifted2, "n",                 alertK = 50,      xlim=c(-2,1),    ylim=c(0,500))
#' plot(amShifted2, "rejPN",             alertK = 50,      xlim=c(-2,1),   ylim=c(0,1))
#' plot(amShifted2, "bias",              alertK = 50,      xlim=c(-2,1),   ylim=c(-0.2,0.2))
#' plot(amShifted2, "mse",               alertK = 50,      xlim=c(-2,1),   ylim=c(0,0.3))
#' plot(amShifted2, "cover",             alertK = 50,      xlim=c(-2,1),   ylim=c(0.75,1))
#' plot(amShifted2, "stopInconclusive",  alertK = 50,      xlim=c(-2,1),   ylim=c(0,1))
#' plot(amShifted2, "stopNotImpactful",  alertK = 50,      xlim=c(-2,1),   ylim=c(0,1))
#' plot(amShifted2, "stopNotTrivial",    alertK = 50,      xlim=c(-2,1),   ylim=c(0,1))

#' # Explore number of steps in affirmation step
#' plot(amShifted2, "n",                waitWidth = 0.4,  xlim=c(-2,1),   ylim=c(0,500),    maxVary = 5)
#' plot(amShifted2, "bias",             waitWidth = 0.4,  xlim=c(-2,1),   ylim=c(-0.1,0.1), maxVary = 5)
#' plot(amShifted2, "mse",              waitWidth = 0.4,  xlim=c(-2,1),   ylim=c(0,0.3),    maxVary = 5)
#' plot(amShifted2, "cover",            waitWidth = 0.4,  xlim=c(-2,1),   ylim=c(0.75,1),   maxVary = 5)
#' plot(amShifted2, "stopInconclusive", waitWidth = 0.4,  xlim=c(-2,1),   ylim=c(0,1),      maxVary = 5)
#' plot(amShifted2, "stopNotImpactful", waitWidth = 0.4,  xlim=c(-2,1),   ylim=c(0,1),      maxVary = 5)
#' plot(amShifted2, "stopNotTrivial",   waitWidth = 0.4,  xlim=c(-2,1),   ylim=c(0,1),      maxVary = 5)


#' # Unrestricted Sample Size with lag of remaining outcomes
#' plot(amShifted2, "stopInconclusive",  alertK = 50,     xlim=c(-2,1),   ylim=c(0,1),   maxVary = 10, sizeRestrictions = "lag")
#' plot(amShifted2, "stopInconclusive",  waitWidth = 0.4, xlim=c(-2,1),   ylim=c(0,1),   maxVary = 10, sizeRestrictions = "lag")
#' plot(amShifted2, "stopInconclusive",  waitWidth = 0.6, xlim=c(-2,1),   ylim=c(0,1),   maxVary = 10, sizeRestrictions = "lag")

#' # Max N with immediate outcomes
#' plot(amShifted2, "n",                    alertK = 50, xlim=c(-2,1),   ylim=c(0,500), maxVary = 10, sizeRestrictions = "maxN")
#' plot(amShifted2, "stopInconclusive",     alertK = 50, xlim=c(-2,1),   ylim=c(0,1),   maxVary = 10, sizeRestrictions = "lag")


#' summary(amShifted2, alertK = 50, waitTime = 0.6, treatEffect = -0.4)
#' summary(amShifted2)
#' summary(am2)


#' # Explore ECDF of Bias and Sample Size
#' ooo <- ecdf.sgpvAM(am = amShifted2, stat = "Bias",alertK = 50,treatEffect = 0,xlim = c(-1,2))
#' ooo <- ecdf.sgpvAM(am = amShifted2, stat = "Bias",waitWidth = 0.4,alertK = 50,treatEffect = 0,xlim = c(-1,2))
#' ooo <- ecdf.sgpvAM(am = amShifted2, stat = "Size",alertK = 50,treatEffect = 0,xlim = c(0,400))

#' # Explore quantiles of sample size for fully specified study design
#' # See probability of stopping by certain sample size under fully specified design
#' ooo <- ecdf.sgpvAM(am = amShifted2, stat = "Size",waitWidth = 0.4,alertK = 50,treatEffect = 0,xlim = c(0,400))
#' quantile(ooo, probs=c(0,0.05,0.1,0.2,0.5,0.8,0.9,0.95,1))
#' ooo(200)
#'
#' @export
sgpvAM <- function(mcmcData=NULL, nreps,
                   waitWidths,
                   dataGeneration=NULL,   dataGenArgs,
                   effectGeneration=NULL, effectGenArgs,
                   modelFit,
                   pointNull, deltaL2, deltaL1, deltaG1, deltaG2,
                   maxAlertSteps=100, lookSteps=1,
                   maxN=NULL, lagOutcomeN=NULL,
                   monitoringIntervalLevel = 0.05, outData = TRUE){


  # 0 Checks
  # If normal data and no standard deviation provided, set sd to 1
  if(class(modelFit)=="normal"){
    if(is.null(dataGenArgs$sd)) dataGenArgs$sd <- 1
  }


  # 1 collect list of simulated data
  if(is.null(mcmcData)){
         mcmcMonitoring <- amData(nreps = nreps,
                                  monitoringIntervalLevel = monitoringIntervalLevel,
                                  dataGeneration   = dataGeneration,   dataGenArgs   = dataGenArgs,
                                  effectGeneration = effectGeneration, effectGenArgs = effectGenArgs,
                                  modelFit         = modelFit)

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
                                minWW         = min(waitWidths),
                                sd            = dataGenArgs$sd))
  getMoreWhich <- which(getMore > 0)
  getMoreWhich

  if( !is.null(mcmcData) & sum(getMore) > 0 & is.null(dataGeneration) & is.null(effectGeneration)){

    stop("Provided mcmcData needs more observations to ensure study completes with unrestricted n\n
         Generate on own or provide data generation inputs.")

  } else {

    while(sum(getMore) > 0){

      print(paste("Adding another", max(getMore), "observations for unrestricted sample size monitoring."))

      for (i in getMoreWhich){

        # Add more data then statistics to generated dataset with insufficient n
        mcmcMonitoring[[i]] <- amDataSingle(monitoringIntervalLevel = monitoringIntervalLevel,
                                            existingData     = mcmcMonitoring[[i]],
                                            dataGeneration   = dataGeneration,   dataGenArgs   = dataGenArgs,
                                            effectGeneration = effectGeneration, effectGenArgs = effectGenArgs,
                                            modelFit         = modelFit)

        mcmcMonitoring[[i]] <- addStats(o = mcmcMonitoring[[i]],
                                        pointNull = pointNull,
                                        deltaL2   = deltaL2, deltaL1 = deltaL1,
                                        deltaG1   = deltaG1, deltaG2 = deltaG2)

      }

      # Continue checking until all datasets have sufficient n
      getMore      <- unlist(lapply(mcmcMonitoring, mcmcMonitoringEnoughCheck,
                                    maxAlertSteps = maxAlertSteps,
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

    # 3 adaptively monitor simulated data across multiple burn ins
    mcmcEOS <- simplify2array(lapply(mcmcMonitoring, sgpvAMrules,
                                     waitWidth               = ww,
                                     lookSteps               = lookSteps,
                                     maxAlertSteps           = maxAlertSteps,
                                     monitoringIntervalLevel = monitoringIntervalLevel,
                                     maxN = maxN, lagOutcomeN = lagOutcomeN))


    # 6 aggregate simulated data
    #   average performance and mse
    #   ecdf of n and bias
    ooAve <- plyr::aaply(mcmcEOS, .margins = c(1,2), .fun = mean)
    ooVar <- plyr::aaply(mcmcEOS, .margins = c(1,2), .fun = var )

    mcmcEndOfStudyAve                <- cbind(ooAve, mse = ooVar[,"bias"] + ooAve[,"bias"]^2)
    mcmcEndOfStudyEcdfSize           <- apply(mcmcEOS[,"n",], 1, ecdf)
    mcmcEndOfStudyEcdfBias           <- apply(mcmcEOS[,"bias",], 1, ecdf)
    names(mcmcEndOfStudyEcdfSize)    <- paste0("alertK_",mcmcEOS[,"alertK",1])
    names(mcmcEndOfStudyEcdfBias)    <- paste0("alertK_",mcmcEOS[,"alertK",1])


    mcmcEndOfStudy[[paste0("width_",ww)]] <-
      list(mcmcEndOfStudyAve      = mcmcEndOfStudyAve,
           mcmcEndOfStudyEcdfSize = mcmcEndOfStudyEcdfSize,
           mcmcEndOfStudyEcdfBias = mcmcEndOfStudyEcdfBias)

  }


  # Indicate whether to keep generated data
  if(outData==FALSE) mcmcMonitoring=NULL

  out <- list(mcmcMonitoring = mcmcMonitoring,
              mcmcEndOfStudy = mcmcEndOfStudy,
              inputs = list(mcmcData         = mcmcData,           nreps = nreps,
                            maxAlertSteps    = maxAlertSteps,  lookSteps = lookSteps,
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
                            outData          = outData))


  # Set class
  class(out) <- append(class(out),"sgpvAM")

  return(out)

}
