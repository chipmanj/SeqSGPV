# sgpvAMrules.R
# J Chipman
#
# From fully sequential monitoring data, applies a burn in and monitors
#  every desired set of steps thereafter.  Raises an alert if sgpv for
#  non-trivial effects or futility equals one.  Returns the affirmation
#  end of study results for affirming an alert starting from 0 to
#  maxAlertSteps by how frequently data are monitored.
#
# Data:                    matrix from sgpvAMdata.  Or a similar matrix with columns:
#                          n, y, lo, hi, z
# waitWidth:               The width of the confidence interval under an assumed variance
# monitoringIntervalLevel: The traditional alpha
# lookSteps:               How frequently will data be monitored after burn in
# maxAlertSteps:           The maximum number of looks before affirming an alert.  Default is 100.



sgpvAMrules <- function(data, waitWidth, monitoringIntervalLevel, lookSteps, maxAlertSteps=100){




  # 1 Establish observations that surpass wait time and occur at lookSteps
  waitTime  <- ceiling((2 * qnorm(1 - monitoringIntervalLevel / 2) * 1 / waitWidth)^2)
  looks     <- waitTime + (0:nrow(data)) * lookSteps
  monitor   <- data[,"n"] %in% looks



  # 2 Raise Alert
  alertNonTrivial <- which(data[,"sgpvNonTrivial"] == 1 & monitor==TRUE)
  alertFutility   <- which(data[,"sgpvFutility"]   == 1 & monitor==TRUE)

  alertNonTrivialAny <- length(alertNonTrivial)>0
  alertFutilityAny   <- length(alertFutility)  >0



  # 3 Affirm alert and report end of study (eos)
  affirmEndOfStudy <- function(alertK){

    if(alertNonTrivialAny) {
      sNT    <- alertNonTrivial[alertNonTrivial %in% (alertNonTrivial + alertK)]
      if(length(sNT) > 0) {
             stopNonTrivial <- min(sNT)
      } else stopNonTrivial <- NA
    } else   stopNonTrivial <- NA

    if(alertFutilityAny) {
      sF     <- alertFutility[alertFutility %in% (alertFutility + alertK)]
      if(length(sF) > 0) {
             stopFutility <- min(sF)
      } else stopFutility <- NA
    } else   stopFutility <- NA

    stop   <- min(stopNonTrivial,stopFutility,na.rm = TRUE)
    eos    <- data[data[,"n"]==stop,]
    if(!is.na(stop)){
      return(c(eos,alertK=alertK))
    } else {
      return(NULL)
    }
  }

  k           <- seq(0,maxAlertSteps,by=lookSteps)
  affirmedEnd <- t(sapply(k,affirmEndOfStudy))


  return(affirmedEnd)

}
