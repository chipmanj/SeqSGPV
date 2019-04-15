# locationShift.R
# J Chipman
#
# For continuous data, location shift to estimate power function


locationShift <- function(o, shiftedThetas, printProgress=TRUE){

  # Only perform given data generated under fixed effect and of location-shift family
  if(  is.function(o$inputs$effectGeneration) |
      !class(o$inputs$modelFit) %in% c("normal","uniform","t")){
    stop("Location shift only performed for data generated under fixed treatment effect and location-shift family.")
  }


  # Keep only columns:
  #  theta, n, y, trt, est, lo, up
  keepColumns <- function(simData){
    simData <- simData[,c("theta", "n", "y", "trt", "est", "lo", "up")]
    simData
  }
  mcmcMonitoring <- lapply(X=o[["mcmcMonitoring"]], keepColumns)

  # Function to location shift treatment effect
  shiftTheta <- function(simData, shift){

    simData[simData[,"trt"]==1,"y"]      <- simData[simData[,"trt"]==1,"y"]      + shift
    simData[,c("theta","est","lo","up")] <- simData[,c("theta","est","lo","up")] + shift
    simData

  }


  # Store results for each location shifted data
  mcmcEndOfStudyShifted <- list()

  for(shift in shiftedThetas){
    if(printProgress) print(paste0("theta shifted by: ",shift))

    o$inputs$mcmcData <- lapply(X = mcmcMonitoring, shiftTheta, shift = shift)
    o$inputs$outData  <- FALSE

    mcmcEndOfStudyShifted[[paste0("theta_",shift)]] <- do.call(sgpvAM, args=o$inputs)
  }

  # Append the original (unshifted) theta (if not already included in mcmcEndOfStudyShifted)
  if(! paste0("theta_",o$inputs$effectGeneration) %in% names(mcmcEndOfStudyShifted) ){
    mcmcEndOfStudyShifted[[paste0("theta_",o$inputs$effectGeneration)]] <- o
  }

  return(mcmcEndOfStudyShifted)

}

