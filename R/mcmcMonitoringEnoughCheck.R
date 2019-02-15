# mcmcMonitoringEnoughCheck.R
# J Chipman
#
# Makes sure enough observations collected to ensure end of study
# when applying monitoring rules



mcmcMonitoringEnoughCheck <- function(o, maxAlertSteps, minWW){

  obs <- nrow(o)

  # Check that sufficient n for burn in CI width
  minN <- ceiling((2 * 1.96 / minWW)^2)
  if(obs - minN <= 0) {
         waitMoreN <- minN - obs
  } else waitMoreN <- NA

  # Check for stability of sgpv
  stabilityNonTrivial <- sum(o[(obs-maxAlertSteps + 1):obs,"sgpvNonTrivial"]==1)
  stabilityFutility   <- sum(o[(obs-maxAlertSteps + 1):obs,"sgpvFutility"]==1)

  # Return T/F regarding stability
  # Get a sense of how many additional observations needed and multiple by arbitraty factor of 4
  getMore0 <- min(maxAlertSteps - stabilityNonTrivial, maxAlertSteps - stabilityFutility) * 5
  getMore  <- min(getMore0, waitMoreN, na.rm = TRUE)

  getMore
}

