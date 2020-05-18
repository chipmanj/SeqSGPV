#' @export
sgpvAMthetas <- function(o, thetas, printProgress=TRUE){

  # Only perform given data generated under fixed effect and of location-shift family
  if( is.function(o$inputs$effectGeneration) ){
    stop("Location shift only performed for data generated under fixed treatment effects.")
  }


  # Store results for each shifted effect
  mcmcEndOfStudyThetas <- list()

  # Make note of original value of theta
  effectOriginal <- o$inputs$effectGeneration

  # Save original effect
  mcmcEndOfStudyThetas[[paste0("theta_",effectOriginal)]]                     <- o
  mcmcEndOfStudyThetas[[paste0("theta_",effectOriginal)]][["mcmcMonitoring"]] <- NULL
  mcmcEndOfStudyThetas[[paste0("theta_",effectOriginal)]][["inputs"]][["mcmcData"]] <- NULL

  # Exclude the original effect if included in thetas
  thetas <- thetas[thetas != effectOriginal]


  if( any(grepl("rnorm|rt",deparse(o$inputs$dataGeneration))) ){
    ## A. If data were generated under normal or t-distributed, use previoulsy
    ##    generated data and shift the estimated mean and confidence intervals


    # Keep columns theta, n, y, trt, est, lo, up for each generated dataset
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


    # Obtain performance under shifted effect
    for(effect in thetas){

      if(printProgress) print(paste0("theta: ",effect))

      # Note the amoung of shift between current effect and original effect
      shift <- effect - effectOriginal

      # New shifted study design
      oo                         <- o
      oo$inputs$mcmcData         <- lapply(X = mcmcMonitoring, shiftTheta, shift = shift)
      oo$inputs$effectGeneration <- effect

      mcmcEndOfStudyThetas[[paste0("theta_",oo$inputs$effectGeneration)]] <- do.call(sgpvAM, args=oo$inputs)

      # If additional data was required for monitoring, keep for future use
      mcmcMonitoringNew <- mcmcEndOfStudyThetas[[paste0("theta_",oo$inputs$effectGeneration)]]$mcmcMonitoring
      mcmcShiftBack     <- lapply(X = mcmcMonitoringNew, shiftTheta, shift = -shift)
      mcmcMonitoring    <- lapply(X = mcmcShiftBack, keepColumns)

      # Don't keep generated data in output
      mcmcEndOfStudyThetas[[paste0("theta_",oo$inputs$effectGeneration)]]$mcmcMonitoring  <- NULL
      mcmcEndOfStudyThetas[[paste0("theta_",oo$inputs$effectGeneration)]]$inputs$mcmcData <- NULL

    }



  } else {
    ## B. Otherwise, generate new data for Thetas treatment effects


    # Drop saved data if any (because will generate new data)
    o$inputs$mcmcData <- NULL

    # Ensure outputs set to not keep output data
    o$inputs$outData <- FALSE

    # Store results for each theta of interest
    mcmcEndOfStudyThetas <- list()

    for(effect in thetas){

      if(printProgress) print(paste0("theta: ",effect))

      # New Thetas study design
      oo                         <- o
      oo$inputs$effectGeneration <- effect

      mcmcEndOfStudyThetas[[paste0("theta_",oo$inputs$effectGeneration)]] <- do.call(sgpvAM, args=oo$inputs)

    }

  }


  # Set class
  class(mcmcEndOfStudyThetas) <- append(c("sgpvAMthetas","sgpvAM"),class(mcmcEndOfStudyThetas))


  return(mcmcEndOfStudyThetas)

}
