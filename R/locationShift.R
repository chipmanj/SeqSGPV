#' @export
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

  # Make note of original value of theta
  effectOriginal <- o$inputs$effectGeneration


  for(shift in shiftedThetas){
    if(printProgress) print(paste0("theta shifted by: ",shift))


    # Use previously generated study design if provided in object o
    if(shift==0){
      mcmcEndOfStudyShifted[[paste0("theta_",effectOriginal)]]                     <- o
      mcmcEndOfStudyShifted[[paste0("theta_",effectOriginal)]][["mcmcMonitoring"]] <- NULL
      mcmcEndOfStudyShifted[[paste0("theta_",effectOriginal)]][["inputs"]][["mcmcData"]] <- NULL
      next
    }

    # New shifted study design
    oo                         <- o
    oo$inputs$mcmcData         <- lapply(X = mcmcMonitoring, shiftTheta, shift = shift)
    oo$inputs$effectGeneration <- effectOriginal + shift

    mcmcEndOfStudyShifted[[paste0("theta_",oo$inputs$effectGeneration)]] <- do.call(sgpvAM, args=oo$inputs)

    # If additional data was required for monitoring, keep for future use
    mcmcMonitoringNew <- mcmcEndOfStudyShifted[[paste0("theta_",oo$inputs$effectGeneration)]]$mcmcMonitoring
    mcmcShiftBack     <- lapply(X = mcmcMonitoringNew, shiftTheta, shift = -shift)
    mcmcMonitoring    <- lapply(X = mcmcShiftBack, keepColumns)

    # Don't keep generated data in output
    mcmcEndOfStudyShifted[[paste0("theta_",oo$inputs$effectGeneration)]]$mcmcMonitoring  <- NULL
    mcmcEndOfStudyShifted[[paste0("theta_",oo$inputs$effectGeneration)]]$inputs$mcmcData <- NULL

  }


  # Set class
  class(mcmcEndOfStudyShifted) <- append(class(mcmcEndOfStudyShifted),
                                         c("sgpvAM","sgpvAMlocationShift"))


  return(mcmcEndOfStudyShifted)

}
