#' @title Function to assess operating characteristics under a range of shifted treatment effects
#'
#' @param o Object from SeqSGPV
#' @param shift Alternative effects to assess operating characteristics
#' @param printProgress TRUE to print progress
#'
#' @export
fixedDesignEffects <- function(o, shift, printProgress=TRUE){

  # Only perform given data generated under fixed effect and of location-shift family
  if( is.function(o$inputs$effectGeneration) ){
    stop("Location shift only performed for data generated under fixed treatment effects.")
  }

  if(length(o$inputs$allocation)==1){
    effectX <- "effect0"
  } else if(length(o$inputs$allocation)==2){
    effectX <- "effect1"
  }


  # Store results for each shifted effect
  mcmcEndOfStudyShifted <- list()

  # Make note of original value of parameter of interest
  effectOriginal <- o$inputs$effectGeneration



  if( any(grepl("rnorm|rt",deparse(o$inputs$dataGeneration))) ){
    ## A. If data were generated under normal or t-distributed, use previoulsy
    ##    generated data and shift the estimated mean and confidence intervals



    # Keep columns theta, n, y, trt, est, lo, up for each generated dataset
    keepColumns <- function(simData){
      simData <- simData[,c("theta0", effectX, "n", "y", "trt", "est", "lo", "up")]
      simData
    }
    mcmcMonitoring <- lapply(X=o[["mcmcMonitoring"]], keepColumns)


    # Function to location shift treatment effect (SE) if randomized or underlying central parameter if non-randomized
    shiftParam <- function(simData, SE){

      simData[simData[,"trt"]==1,"y"]      <- simData[simData[,"trt"]==1,"y"]      + SE
      simData[,c(effectX,"est","lo","up")] <- simData[,c(effectX,"est","lo","up")] + SE
      simData

    }


    # Obtain performance under shifted effect
    for(SE in shift){
      if(printProgress) print(paste0("effect: ",SE))


      # Use previously generated study design if provided in object o
      if(SE==effectOriginal){
        mcmcEndOfStudyShifted[[paste0(effectX,"_",effectOriginal)]]                     <- o
        mcmcEndOfStudyShifted[[paste0(effectX,"_",effectOriginal)]][["mcmcMonitoring"]] <- NULL
        mcmcEndOfStudyShifted[[paste0(effectX,"_",effectOriginal)]][["inputs"]][["mcmcData"]] <- NULL
        next
      }

      # New shifted study design
      oo                         <- o
      oo$inputs$mcmcData         <- lapply(X = mcmcMonitoring, shiftParam, SE = SE)
      oo$inputs$effectGeneration <- effectOriginal + SE

      mcmcEndOfStudyShifted[[paste0(effectX,"_",oo$inputs$effectGeneration)]] <- do.call(SeqSGPV, args=oo$inputs)

      # If additional data was required for monitoring, keep for future use
      mcmcMonitoringNew <- mcmcEndOfStudyShifted[[paste0(effectX,"_",oo$inputs$effectGeneration)]]$mcmcMonitoring
      mcmcShiftBack     <- lapply(X = mcmcMonitoringNew, shiftParam, SE = -SE)
      mcmcMonitoring    <- lapply(X = mcmcShiftBack, keepColumns)

      # Don't keep generated data in output
      mcmcEndOfStudyShifted[[paste0(effectX,"_",oo$inputs$effectGeneration)]]$mcmcMonitoring  <- NULL
      mcmcEndOfStudyShifted[[paste0(effectX,"_",oo$inputs$effectGeneration)]]$inputs$mcmcData <- NULL

    }



  } else {
    ## B. Otherwise, generate new data for shifted treatment effects


    # Drop saved data if any (because will generate new data)
    o$inputs$mcmcData <- NULL

    # Ensure outputs set to not keep output data
    o$inputs$outData <- FALSE

    # Store results for each theta of interest
    mcmcEndOfStudyShifted <- list()

    for(SE in shift){

      if(printProgress) print(paste0("effect: ",SE))

      # Use previously generated study design if provided in object o
      if(SE==effectOriginal){
        mcmcEndOfStudyShifted[[paste0(effectX,"_",effectOriginal)]]                     <- o
        mcmcEndOfStudyShifted[[paste0(effectX,"_",effectOriginal)]][["mcmcMonitoring"]] <- NULL
        mcmcEndOfStudyShifted[[paste0(effectX,"_",effectOriginal)]][["inputs"]][["mcmcData"]] <- NULL
        next
      }

      # New shifted study design
      oo                         <- o
      oo$inputs$effectGeneration <- SE

      mcmcEndOfStudyShifted[[paste0(effectX,"_",oo$inputs$effectGeneration)]] <- do.call(SeqSGPV, args=oo$inputs)

    }

  }


  # Set class
  class(mcmcEndOfStudyShifted) <- append(c("SeqSGPVeffects","SeqSGPV"),class(mcmcEndOfStudyShifted))


  return(mcmcEndOfStudyShifted)

}
