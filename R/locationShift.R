# locationShift.R
# J Chipman
#
# For continuous data, location shift to estimate power function


locationShift <- function(o, shiftedThetas){

  # Keep only columns:
  #  theta, n, y, trt, est, lo, up
  keepColumns <- function(simData){
    simData <- simData[,c("theta", "n", "y", "trt", "est", "lo", "up")]
    simData
  }
  mcmcMonitoring <- lapply(X=o[["mcmcMonitoring"]], keepColumns)

  # Function to location shift treatment effect
  shiftTheta <- function(simData, shift){

    simData[simData[,"trt"]==1,"y"] <- simData[simData[,"trt"]==1,"y"] + shift
    simData[,c("est","lo","up")]    <- simData[,c("est","lo","up")] + shift
    simData

  }


  # Store results for each location shifted data
  mcmcEndOfStudyShifted <- list()

  for(shift in shiftedThetas){

    mcmcMonitoringShifted <- lapply(X = mcmcMonitoring, shiftTheta, shift = shift)

    mcmcEndOfStudyShifted[[paste0("theta_",shift)]] <-
      sgpvAM(mcmcData  = mcmcMonitoringShifted,
             maxAlertSteps = o$maxAlertSteps, lookSteps=o$lookSteps,
             waitWidths = o$waitWidths,
             pointNull  = o$pointNull,
             deltaL2    = o$deltaL2, deltaL1 = o$deltaL1,
             deltaG1    = o$deltaG1, deltaG2 = o$deltaG2,
             monitoringIntervalLevel = o$monitoringIntervalLevel,
             outData=FALSE)

  }

  return(mcmcEndOfStudyShifted)

}

