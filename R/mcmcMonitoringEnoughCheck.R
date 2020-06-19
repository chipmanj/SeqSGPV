#' @export

# mcmcMonitoringEnoughCheck.R
# J Chipman
#
# Makes sure enough observations collected to ensure end of study
# when applying monitoring rules



mcmcMonitoringEnoughCheck <- function(o, waitN, lookSteps, maxAlertSteps, lagOutcomeN){

  obs <- nrow(o)

  # Check that sufficient observed n
  minN  <- waitN + maxAlertSteps + lagOutcomeN
  if(obs - minN <= 0) {
         waitMoreN <- minN - obs + maxAlertSteps
  } else waitMoreN <- 0

  # Check for stability of sgpv
  #  - Additional observations needed to have SGPVs through the end of stopping rule + lag time
  addedStabilityN <- maxAlertSteps + lagOutcomeN + lookSteps
  if(obs > addedStabilityN ){

    # How many of the last set of observations would indicate to stop, +1 for indexing
    stabilityROPE <- sum(o[(obs-addedStabilityN+1):obs,"sgpvROPE"]==0)
    stabilityROME <- sum(o[(obs-addedStabilityN+1):obs,"sgpvROME"]==0)

    # Get a sense of how many additional observations needed and multiple by arbitraty factor of 4
    getMore0 <- min(addedStabilityN - stabilityROPE, addedStabilityN - stabilityROME) * 4

  } else {
    getMore0 <- addedStabilityN * 4
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
