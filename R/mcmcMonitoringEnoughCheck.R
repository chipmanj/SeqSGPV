#' @export

# mcmcMonitoringEnoughCheck.R
# J Chipman
#
# Makes sure enough observations collected to ensure end of study
# when applying monitoring rules



mcmcMonitoringEnoughCheck <- function(o, maxAlertSteps, lagOutcomeN, minWW, sd){

  obs <- nrow(o)

  # Check that sufficient observed n for burn in CI width + the lag time on outcome for unobserved outcomes
  waitN <- ceiling((2 * 1.96 * sd / minWW)^2)
  minN  <- waitN + maxAlertSteps + lagOutcomeN
  if(obs - minN <= 0) {
         waitMoreN <- minN - obs + maxAlertSteps
  } else waitMoreN <- 0

  # Check for stability of sgpv
  #  - Require at least as many observations as maxAlertSteps
  stabilityN <- maxAlertSteps + lagOutcomeN
  if(obs > stabilityN ){
    stabilityTrivial   <- sum(o[(obs-stabilityN + 1):obs,"sgpvTrivial"]==0)
    stabilityImpactful <- sum(o[(obs-stabilityN + 1):obs,"sgpvImpactful"]==0)

    # Get a sense of how many additional observations needed and multiple by arbitraty factor of 4
    getMore0 <- min(stabilityN - stabilityTrivial, stabilityN - stabilityImpactful) * 4

  } else {
    getMore0 <- stabilityN * 4
  }


  # Use burn in time if greater than the additional number of obs needed for stability
  # Else use n required to achieve greater stability
  if(waitMoreN > getMore0){
    getMore <- waitMoreN
  } else {
    getMore <- getMore0
  }

  # If getting more, get at least 50 more observations (arbitrary)
  if(getMore > 0 & getMore < 50) getMore <- max(getMore, 50)

  getMore
}
