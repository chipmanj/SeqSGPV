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


pooledVariance <- function(i){
  dd  <- d[1:i,]
  y0  <- dd[dd[,"trt"]==0,"y"]
  y1  <- dd[dd[,"trt"]==1,"y"]
  vp <- ((length(y0) - 1) * var(y0) + (length(y1) - 1) * var(y1)) / (i - 2)
  vp
}


plot(x=4:nrow(d),y=sapply(4:nrow(d),pooledVariance),type="l",xlim=c(0,200))
abline(v=15,col="red")
abline(v=30,col="red")

var(sapply(4:14,pooledVariance))
var(sapply(15:24,pooledVariance))
var(sapply(25:34,pooledVariance))
