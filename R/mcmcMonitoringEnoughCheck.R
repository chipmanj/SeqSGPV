#' @export

# mcmcMonitoringEnoughCheck.R
# J Chipman
#
# Makes sure enough observations collected to ensure end of study
# when applying monitoring rules



mcmcMonitoringEnoughCheck <- function(o, waitN, maxAlertSteps, lagOutcomeN){

  obs <- nrow(o)

  # Check that sufficient observed n
  minN  <- waitN + maxAlertSteps + lagOutcomeN
  if(obs - minN <= 0) {
         waitMoreN <- minN - obs + maxAlertSteps
  } else waitMoreN <- 0

  # Check for stability of sgpv
  #  - Additional observations needed to be consistent with stopping rule
  #  - +1 for element wise indexing: an observation with 0 alerts will stop at observing one stopping rule
  addedStabilityN <- maxAlertSteps + lagOutcomeN + 1
  if(obs > addedStabilityN ){

    # How many of the last set of observations would indicate to stop
    stabilityROPE <- sum(o[(obs-addedStabilityN):obs,"sgpvROPE"]==0)
    stabilityROME <- sum(o[(obs-addedStabilityN):obs,"sgpvROME"]==0)

    # Get a sense of how many additional observations needed and multiple by arbitraty factor of 4
    getMore0 <- min(addedStabilityN - stabilityROPE, 1 + addedStabilityN - stabilityROME) * 4

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
