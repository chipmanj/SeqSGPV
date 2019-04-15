# sgpvAMrules.R
# J Chipman
#
# From fully sequential monitoring data, applies a burn in and monitors
#  every desired set of steps thereafter.  Raises an alert if sgpv for
#  non-trivial effects or futility equals one.  Returns the affirmation
#  end of study results for affirming an alert starting from 0 to
#  maxAlertSteps by how frequently data are monitored.
#
# Data: matrix from sgpvAMdata.  Or a similar matrix with columns:
#       theta, n, y, lo, hi, sgpvTrivial, sgpvImpactful,
# waitWidth:               The width of the confidence interval under an assumed variance
# monitoringIntervalLevel: The traditional alpha
# lookSteps:               How frequently will data be monitored after burn in
# maxAlertSteps:           The maximum number of looks before affirming an alert.  Default is 100.



sgpvAMrules <- function(data,      waitWidth,
                        monitoringIntervalLevel,
                        lookSteps, maxAlertSteps=100,
                        maxN,      lagOutcomeN){

  # 1 Establish observations that surpass wait time and occur at lookSteps
  waitTime  <- ceiling((2 * qnorm(1 - monitoringIntervalLevel / 2) * 1 / waitWidth)^2)
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


    if(alertNotTrivialAny) {
      sNT    <- alertNotTrivial[alertNotTrivial %in% (alertNotTrivial + alertK)]
      if(length(sNT) > 0) {
             stopNotTrivial <- min(sNT)
      } else stopNotTrivial <- NA
    } else   stopNotTrivial <- NA

    if(alertNotImpactfulAny) {
      sNI     <- alertNotImpactful[alertNotImpactful %in% (alertNotImpactful + alertK)]
      if(length(sNI) > 0) {
             stopNotImpactful <- min(sNI)
      } else stopNotImpactful <- NA
    } else   stopNotImpactful <- NA



    # Unconstrained N stop
    stop                     <- min(stopNotTrivial, stopNotImpactful, na.rm = TRUE)
    eos                      <- data[data[,"n"]==stop,]
    eos["stopNotTrivial"]    <- as.numeric(eos["sgpvTrivial"]==0)
    eos["stopNotImpactful"]  <- as.numeric(eos["sgpvImpactful"]==0)
    eos["stopInconclusive"]  <- as.numeric(eos["stopNotTrivial"]==0 & eos["stopNotImpactful"]==0)

    # Unconstrained with lag on outcome
    if(!is.null(lagOutcomeN)){
      eosLag                     <- data[data[,"n"]==stop + lagOutcomeN,]
      eosLag["stopNotTrivial"]   <- as.numeric(eosLag["sgpvTrivial"]==0)
      eosLag["stopNotImpactful"] <- as.numeric(eosLag["sgpvImpactful"]==0)
      eosLag["stopInconclusive"] <- as.numeric(eosLag["stopNotTrivial"]==0 &
                                               eosLag["stopNotImpactful"]==0)
    } else eosLag <- NULL

    # MaxN stop
    if(!is.null(maxN)){
      eosMaxN                     <- data[data[,"n"]==min(stop,maxN),]
      eosMaxN["stopNotTrivial"]   <- as.numeric(eosMaxN["sgpvTrivial"]==0)
      eosMaxN["stopNotImpactful"] <- as.numeric(eosMaxN["sgpvImpactful"]==0)
      eosMaxN["stopInconclusive"] <- as.numeric(eosMaxN["stopNotTrivial"]==0 &
                                                eosMaxN["stopNotImpactful"]==0)
    } else eosMaxN <- NULL

    # lag on outcome, having enrolled MaxN stop
    if(!is.null(lagOutcomeN) & !is.null(maxN)){
      eosLagMaxN                     <- data[data[,"n"]==min(stop,maxN)+lagOutcomeN,]
      eosLagMaxN["stopNotTrivial"]   <- as.numeric(eosLagMaxN["sgpvTrivial"]==0)
      eosLagMaxN["stopNotImpactful"] <- as.numeric(eosLagMaxN["sgpvImpactful"]==0)
      eosLagMaxN["stopInconclusive"] <- as.numeric(eosLagMaxN["stopNotTrivial"]==0 &
                                                   eosLagMaxN["stopNotImpactful"]==0)
    } else eosLagMaxN <- NULL

    # Keep theta and operating characteristics (drop y, trt, and sgpvs)
    keepStats <- c("n","est","bias","rejPN","cover", "stopNotTrivial","stopNotImpactful","stopInconclusive")


    if(!is.na(stop)){
      oc <- c(eos[c("theta", keepStats)],
              lag     = eosLag[keepStats],
              maxN    = eosMaxN[keepStats],
              lagMaxN = eosLagMaxN[keepStats])
      # oc <- c(eos[c("theta", keepStats)],
      #         lag     = eosLag[keepStats],
      #         maxN    = eosMaxN[keepStats])

      return(c(oc,alertK=alertK))
    } else {
      return(NULL)
    }
  }

  k           <- seq(0,maxAlertSteps,by=lookSteps)
  affirmedEnd <- t(sapply(k,affirmEndOfStudy))


  return(affirmedEnd)

}
