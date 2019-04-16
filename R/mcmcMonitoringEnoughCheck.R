#' @export

# mcmcMonitoringEnoughCheck.R
# J Chipman
#
# Makes sure enough observations collected to ensure end of study
# when applying monitoring rules



mcmcMonitoringEnoughCheck <- function(o, maxAlertSteps, minWW, sd){

  obs <- nrow(o)

  # Check that sufficient n for burn in CI width
  minN <- ceiling((2 * 1.96 * sd / minWW)^2)
  if(obs - minN <= 0) {
         waitMoreN <- minN - obs + maxAlertSteps
  } else waitMoreN <- 0

  # Check for stability of sgpv
  #  - Require at least as many observations as maxAlertSteps
  if(obs > maxAlertSteps){
    stabilityTrivial   <- sum(o[(obs-maxAlertSteps + 1):obs,"sgpvTrivial"]==0)
    stabilityImpactful <- sum(o[(obs-maxAlertSteps + 1):obs,"sgpvImpactful"]==0)

    # Get a sense of how many additional observations needed and multiple by arbitraty factor of 4
    getMore0 <- min(maxAlertSteps - stabilityTrivial, maxAlertSteps - stabilityImpactful) * 5
  } else {
    getMore0 <- maxAlertSteps * 5
  }


  # Use burn in time if greater than the additional number of obs needed for stability
  # Else use n required to achieve greater stability
  if(waitMoreN > getMore0){
    getMore <- waitMoreN
  } else {
    getMore <- getMore0
  }

  getMore
}
