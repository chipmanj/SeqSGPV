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
# deltaL1:                 The delta less than and closest to the point null.
# deltaL2:                 The delta less than and furthest from the point null.
# deltaG1:                 The delta greater than and closest to the point null.
# deltaG2:                 The delta greater than and furthest from the point null.
# waitWidth:               The width of the confidence interval under an assumed variance
# monitoringIntervalLevel: The traditional alpha
# lookSteps:               How frequently will data be monitored after burn in
# maxAlertSteps:           The maximum number of looks before affirming an alert.  Default is 100.



sgpvAMrules <- function(data,
                        deltaL2, deltaL1, deltaG1, deltaG2,
                        waitWidth, monitoringIntervalLevel, lookSteps, maxAlertSteps=100){




  # 1 Obtain sgpv
  if(!anyNA(c(deltaL2, deltaL1, deltaG1, deltaG2))){

    # Two sided
    sgpvNonTrivial <- 1 - sgpv::sgpvalue(est.lo = data[,"lo"], est.hi = data[,"up"], null.lo = deltaL1, null.hi = deltaG1)$p.delta
    sgpvFutility   <-     sgpv::sgpvalue(est.lo = data[,"lo"], est.hi = data[,"up"], null.lo = deltaL2, null.hi = deltaG2)$p.delta
  } else if(!anyNA(c(deltaL2, deltaL1))){

    # One sided: efficacy when less than null
    sgpvNonTrivial <- 1 - sgpv::sgpvalue(est.lo = data[,"lo"], est.hi = data[,"up"], null.lo =  -10^10, null.hi = deltaL1)$p.delta
    sgpvFutility   <-     sgpv::sgpvalue(est.lo = data[,"lo"], est.hi = data[,"up"], null.lo = deltaL2, null.hi = 10^10)$p.delta
  } else if(!anyNA(c(deltaG1, deltaG2))){

    # One sided: efficacy when less than null
    sgpvNonTrivial <- 1 - sgpv::sgpvalue(est.lo = data[,"lo"], est.hi = data[,"up"], null.lo = deltaG1, null.hi = 10^10)$p.delta
    sgpvFutility   <-     sgpv::sgpvalue(est.lo = data[,"lo"], est.hi = data[,"up"], null.lo =  -10^10, null.hi = deltaG2)$p.delta
  } else{

    stop("A one sided study requires both deltas to be strictly greater or lower than point null")
  }

  data <- cbind(data,sgpvNonTrivial,sgpvFutility)



  # 2 Establish observations that surpass wait time and occur at lookSteps
  waitTime  <- ceiling((2 * qnorm(1 - monitoringIntervalLevel / 2) * 1 / waitWidth)^2)
  looks     <- waitTime + (0:nrow(data)) * lookSteps
  monitor   <- data[,"n"] %in% looks



  # 3 Raise Alert
  alertNonTrivial <- which(sgpvNonTrivial == 1 & monitor==TRUE)
  alertFutility   <- which(sgpvFutility   == 1 & monitor==TRUE)

  alertNonTrivialAny <- length(alertNonTrivial)>0
  alertFutilityAny   <- length(alertFutility)  >0



  # 4 Affirm alert and report end of study (eos)
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
