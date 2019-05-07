#' @export
sgpvAMrulesSingle <- function(data,
                              waitWidth,
                              sdY,
                              waitEmpirical, minWaitN, midPointLess, midPointGreater,
                              monitoringIntervalLevel,
                              lookSteps, kSteps, maxAlertSteps=100,
                              maxN,      lagOutcomeN=0){


  # If not using assumed sdY for determining wait Time
  if(waitEmpirical==TRUE){

    # Standard deviation on beta coefficient for trt effect = | CI Width | * sqrt ( n ) / ( 2 * 1.96 )
    # sdBeta <- abs(data[,"lo"] - data[,"up"]) * sqrt(1:nrow(data)) / (2 * qnorm(1 - monitoringIntervalLevel / 2))


    # Relax wait time (ME width) if estimated effect outside peripheral of clinical region boundaries
    # addWidth          <- rep(0,nrow(data))
    # relaxLo           <- which(data[,"est"]      < periphLo)
    # relaxUp           <- which(data[,"est"]      > periphUp)
    # addWidth[relaxLo] <- abs(data[relaxLo,"est"] - periphLo)
    # addWidth[relaxUp] <-     data[relaxUp,"est"] - periphUp

    # Add to ME distance from estimate to closest midPoint
    # Divide by 2 because margin of error is doubled for CI width
    addWidth <- pmin( abs( data[,"est"] - midPointLess ),
                      abs( data[,"est"] - midPointGreater ),
                     na.rm = TRUE) / 6

    meWaitWidth <- waitWidth + addWidth


    # Start monitoring once Margin of Error is less than wait width
    # Require minimum wait time
    # Require start affirmation of [0] steps
    possibleStarts <- (1:nrow(data))[abs(data[,"lo"] - data[,"up"]) / 2 < meWaitWidth]
    possibleStarts <- possibleStarts[possibleStarts > minWaitN & !is.na(possibleStarts)]
    waitTime <- min(possibleStarts[possibleStarts %in% c(possibleStarts + 30)])

  } else {

    # Relax wait time (ME Width) if estimated effect outside peripheral of clinical region boundaries
    theta       <- data[1,"theta"]
    addWidth    <- 0 +
                   as.numeric(theta < periphLo) * abs(theta - periphLo) +
                   as.numeric(theta > periphUp) *    (theta - periphUp)

    meWaitWidth <- waitWidth + addWidth

    # Note: Coefficient ME = 2 * 1.96 * S_Y * ( 1 / S_trt ) * ( 1 / sqrt(n) )
    # With 1:1 allocation, S_trt = 0.5
    waitTime    <- ceiling( ( 2 * qnorm(1 - monitoringIntervalLevel / 2) * sdY / meWaitWidth)^2 )

}


  # 1 Establish observations that surpass wait time and occur at lookSteps
  looks     <- waitTime + (0:nrow(data)) * lookSteps
  monitor   <- data[,"n"] %in% looks



  # 2 Raise Alert
  alertNotTrivial   <- which(data[,"sgpvTrivial"]   == 0 & monitor==TRUE)
  alertNotImpactful <- which(data[,"sgpvImpactful"] == 0 & monitor==TRUE)

  alertNotTrivialAny   <- length(alertNotTrivial)   > 0
  alertNotImpactfulAny <- length(alertNotImpactful) > 0



  # 3 Affirm alert and report end of study (eos) operating characteristics (oc)
  #   - To be applied across various k (affirmation steps)
  affirmEndOfStudy <- function(alertK){


    # Stop times for being Not Trivial and being Not Impactful
    stopNotTrivial   <- alertNotTrivial[  alertNotTrivial   %in% (alertNotTrivial   + alertK)]
    stopNotImpactful <- alertNotImpactful[alertNotImpactful %in% (alertNotImpactful + alertK)]

    # Unconstrained N stop
    stop                     <- min(stopNotTrivial, stopNotImpactful, na.rm = TRUE)
    eos                      <- data[data[,"n"]==stop,]
    eos["stopNotTrivial"]    <- as.numeric(eos["sgpvTrivial"]==0)
    eos["stopNotImpactful"]  <- as.numeric(eos["sgpvImpactful"]==0)
    eos["stopInconclusive"]  <- as.numeric(eos["stopNotTrivial"]==0 & eos["stopNotImpactful"]==0)

    # Unconstrained with lag on outcome
    eosLag                     <- data[data[,"n"]==stop + lagOutcomeN,]
    eosLag["stopNotTrivial"]   <- as.numeric(eosLag["sgpvTrivial"]==0)
    eosLag["stopNotImpactful"] <- as.numeric(eosLag["sgpvImpactful"]==0)
    eosLag["stopInconclusive"] <- as.numeric(eosLag["stopNotTrivial"]==0 &
                                             eosLag["stopNotImpactful"]==0)

    # MaxN stop
    if(!is.null(maxN)){
      eosMaxN                     <- data[data[,"n"]==min(stop,maxN),]
      eosMaxN["stopNotTrivial"]   <- as.numeric(eosMaxN["sgpvTrivial"]==0)
      eosMaxN["stopNotImpactful"] <- as.numeric(eosMaxN["sgpvImpactful"]==0)
      eosMaxN["stopInconclusive"] <- as.numeric(eosMaxN["stopNotTrivial"]==0 &
                                                eosMaxN["stopNotImpactful"]==0)

      eosLagMaxN                     <- data[data[,"n"]==min(stop,maxN)+lagOutcomeN,]
      eosLagMaxN["stopNotTrivial"]   <- as.numeric(eosLagMaxN["sgpvTrivial"]==0)
      eosLagMaxN["stopNotImpactful"] <- as.numeric(eosLagMaxN["sgpvImpactful"]==0)
      eosLagMaxN["stopInconclusive"] <- as.numeric(eosLagMaxN["stopNotTrivial"]==0 &
                                                   eosLagMaxN["stopNotImpactful"]==0)

    } else eosMaxN <- eosLagMaxN <- NULL


    # Keep theta and operating characteristics (drop y, trt, and sgpvs)
    keepStats <- c("n","est","bias","rejPN","cover", "stopNotTrivial","stopNotImpactful","stopInconclusive")


    if(!is.na(stop)){
      oc <- c(eos[c("theta", keepStats)],
              lag     = eosLag[keepStats],
              maxN    = eosMaxN[keepStats],
              lagMaxN = eosLagMaxN[keepStats])

      return(c(oc,alertK=alertK))
    } else {
      return(NULL)
    }
  }

  k           <- seq(0,maxAlertSteps,by=kSteps)
  affirmedEnd <- t(sapply(k,affirmEndOfStudy))


  return(affirmedEnd)

}
