# mcmcMonitoringEnoughCheck.R
# J Chipman
#
# Makes sure enough observations collected to ensure end of study
# when applying monitoring rules



mcmcMonitoringEnoughCheck <- function(o, maxAlertSteps){

  # Check for stability of sgpv
  stabilityNonTrivial <- sum(o[(nrow(o)-maxAlertSteps + 1):nrow(o),"sgpvNonTrivial"]==1)
  stabilityFutility   <- sum(o[(nrow(o)-maxAlertSteps + 1):nrow(o),"sgpvFutility"]==1)

  # Return T/F regarding stability
  # Get a sense of how many additional observations needed and multiple by arbitraty factor of 4
  getMore <- min(maxAlertSteps - stabilityNonTrivial, maxAlertSteps - stabilityFutility) * 4

  getMore
}

