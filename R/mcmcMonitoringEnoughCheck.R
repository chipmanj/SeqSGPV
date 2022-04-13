#' @title mcmcMontiroingEnoughCheck
#'
#' @description Makes sure enough observations collected to ensure end of study when applying monitoring rules. mcmcMontiroingEnoughCheck is called within SeqSPGV.
#'
#' @param o List of mcmc replicates
#' @param monitoringFrequency Matrix of monitoring frequency combinations
#'
#' @export
mcmcMonitoringEnoughCheck <- function(o, monitoringFrequency){

  obs  <- nrow(o)
  maxW <- max(monitoringFrequency[,"W"])
  maxS <- max(monitoringFrequency[,"S"])
  maxA <- max(monitoringFrequency[,"A"])
  maxL <- max(monitoringFrequency[,"L"])

  # Check that sufficient observed n
  minN  <- maxW + maxA + maxL
  if(obs - minN <= 0) {
         waitMoreN <- minN - obs + maxA
  } else waitMoreN <- 0

  # Check for stability of sgpv
  #  - Additional observations needed to have SGPVs through the end of stopping rule + lag time
  addedStabilityN <- maxA + maxS
  if( (obs - maxL) > addedStabilityN ){

    # How many of the last set of observations would indicate to stop, +1 for indexing
    stabilityROPE <- sum(o[(obs - maxL - addedStabilityN+1):(obs - maxL),"sgpvROPE"]==0)
    stabilityROME <- sum(o[(obs - maxL - addedStabilityN+1):(obs - maxL),"sgpvROME"]==0)

    # Get a sense of how many additional observations needed and multiple by arbitrary factor of 4
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
