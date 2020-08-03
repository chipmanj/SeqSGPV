#' @export
sgpvAMeffects <- function(o, effects, printProgress=TRUE,locationShift){

  # Only perform given data generated under fixed effect and of location-shift family
  if( is.function(o$inputs$effectGeneration) ){
    stop("Location shift only performed for data generated under fixed treatment effects.")
  }

  # Determine if inference can use location-shifts to speed up computations (rather than re-sampling for each new effect)
  # Identified through any of:
  # 1. Using the lmCI model fit assuming normality
  # 2. Setting locationShift == TRUE
  if( any(grepl("fastLmPure",deparse(o$inputs$modelFit))) ) {
    locationShift <- TRUE
  } else if (missing(locationShift)){
    locationShift <- FALSE
  } else if(!is.logical(locationShift)){
    stop("locationShift must be inputted as logical TRUE or FALSE.")
  }

  # Store results for each shifted effect
  mcmcEndOfStudyEffects <- list()

  # Make note of original value of effect
  effectOriginal <- o$inputs$effectGeneration

  # Save original effect
  mcmcEndOfStudyEffects[[paste0("effect_",effectOriginal)]]                     <- o
  mcmcEndOfStudyEffects[[paste0("effect_",effectOriginal)]][["mcmcMonitoring"]] <- NULL
  mcmcEndOfStudyEffects[[paste0("effect_",effectOriginal)]][["inputs"]][["mcmcData"]] <- NULL

  # Exclude the original effect if included in effects objects
  effects <- effects[effects != effectOriginal]


  if( locationShift ){

    # Keep columns theta0, effect1, n, y, trt, est, lo, up for each generated dataset
    keepColumns <- function(simData){
      simData <- simData[,c("theta0","effect1", "n", "y", "trt", "est", "lo", "up")]
      simData
    }
    mcmcMonitoring <- lapply(X=o[["mcmcMonitoring"]], keepColumns)


    # Function to location shift treatment effect
    shiftEffect <- function(simData, shift){

      simData[simData[,"trt"]==1,"y"]        <- simData[simData[,"trt"]==1,"y"]        + shift
      simData[,c("effect1","est","lo","up")] <- simData[,c("effect1","est","lo","up")] + shift
      simData

    }


    # Obtain performance under shifted effect
    for(effect in effects){

      if(printProgress) print(paste0("effect: ",effect))

      # Note the amoung of shift between current effect and original effect
      shift <- effect - effectOriginal

      # New shifted study design
      oo                         <- o
      oo$inputs$mcmcData         <- lapply(X = mcmcMonitoring, shiftEffect, shift = shift)
      oo$inputs$effectGeneration <- effect

      mcmcEndOfStudyEffects[[paste0("effect_",oo$inputs$effectGeneration)]] <- do.call(sgpvAM, args=oo$inputs)

      # If additional data was required for monitoring, keep for future use
      mcmcMonitoringNew <- mcmcEndOfStudyEffects[[paste0("effect_",oo$inputs$effectGeneration)]]$mcmcMonitoring
      mcmcShiftBack     <- lapply(X = mcmcMonitoringNew, shiftEffect, shift = -shift)
      mcmcMonitoring    <- lapply(X = mcmcShiftBack, keepColumns)

      # Don't keep generated data in output
      mcmcEndOfStudyEffects[[paste0("effect_",oo$inputs$effectGeneration)]]$mcmcMonitoring  <- NULL
      mcmcEndOfStudyEffects[[paste0("effect_",oo$inputs$effectGeneration)]]$inputs$mcmcData <- NULL

    }



  } else {


    # Drop saved data if any (because will generate new data)
    o$inputs$mcmcData <- NULL

    # Ensure outputs set to not keep output data
    o$inputs$outData <- FALSE

    for(effect in effects){

      if(printProgress) print(paste0("effect: ",effect))

      # New Effects study design
      oo                         <- o
      oo$inputs$effectGeneration <- effect

      mcmcEndOfStudyEffects[[paste0("effect_",oo$inputs$effectGeneration)]] <- do.call(sgpvAM, args=oo$inputs)

    }

  }


  # Set class
  class(mcmcEndOfStudyEffects) <- append(c("sgpvAMeffects","sgpvAM"),class(mcmcEndOfStudyEffects))


  return(mcmcEndOfStudyEffects)

}
