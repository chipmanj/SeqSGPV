#' @export
summary.sgpvAM <- function(am, alertK, waitTime, treatEffect, rd = 4){


  # Inputs stored in am object
  amInputs <- am$inputs


  # Select treatment effect to summarize and reduce am object
  if(is.element(el = "sgpvAMlocationShift",class(am))){

    if(missing(treatEffect)){
      cat(paste0(paste0("select treatment effect of interest from choices:\n",
                        paste0(unname(sapply(names(am),
                                             FUN = function(x) as.numeric(strsplit(x,split="_")[[1]][2]))),
                               collapse=", "))))
      te <- readline(prompt="input: ")
    } else {
      te <- treatEffect
    }

    # Reduce to specific out object
    o  <- am[[paste0("theta_",te)]]

  } else if(!is.function(am$inputs$effectGeneration)) {

    te  <- amInputs$effectGeneration
    o   <- am
  } else {
    return(NULL)
  }

  lag  <- amInputs$lagOutcomeN
  maxN <- amInputs$maxN


  # Select and reduce to wait time to summarize
  if(missing(waitTime)){
    cat(paste0(paste0("select wait time until applying monitoring rules:\n",
                      paste0(unname(sapply(names(o[["mcmcEndOfStudy"]]),
                                           FUN = function(x) as.numeric(strsplit(x,split="_")[[1]][2]))),
                             collapse=", "))))
    w <- readline(prompt="input: ")

  } else w <- waitTime

  o <- o[["mcmcEndOfStudy"]][[paste0("width_",w)]][["mcmcEndOfStudyAve"]]


  # Select and reduce to number of affirmation steps required for stopping
  if(missing(alertK)){
    cat(paste0(paste0("select required affirmation steps:\n",
                      paste0(o[,"alertK"], collapse=", "))))
    k <- readline(prompt="input: ")

  } else k <- alertK

  o <- o[o[,"alertK"]==k,]





  cat(paste0("Given: theta = ", te, ", wait time = ", w, ", and k = ",k))

  if(is.element("n",names(o))){
    cat("\n\nUnrestricted sample size (immediate outcomes):")
    cat(paste0("\nAverage sample size           = ", round(o["n"],rd)))
    cat(paste0("\nP( reject point null )        = ", o["rejPN"]))
    cat(paste0("\nP( conclude not ROPE effect ) = ", o["stopNotROPE"]))
    cat(paste0("\nP( conclude not ROME effect ) = ", o["stopNotROME"]))
    cat(paste0("\nP( conclude inconclusive )    = ", o["stopInconclusive"]))
    cat(paste0("\nCoverage                      = ", round(o["cover"],rd)))
    cat(paste0("\nBias                          = ", round(o["bias"],rd)))
  }


  if(is.element("lag.n",names(o))){
    cat(paste0("\n\nUnrestricted sample size (stopping and then observing ", lag, " lagged outcomes):"))
    cat(paste0("\nAverage Sample Size           = ", round(o["lag.n"],rd)))
    cat(paste0("\nP( reject point null )        = ", o["lag.rejPN"]))
    cat(paste0("\nP( conclude not ROPE effect ) = ", o["lag.stopNotROPE"]))
    cat(paste0("\nP( conclude not ROME effect ) = ", o["lag.stopNotROME"]))
    cat(paste0("\nP( conclude inconclusive )    = ", o["lag.stopInconclusive"]))
    cat(paste0("\nCoverage                      = ", round(o["lag.cover"],rd)))
    cat(paste0("\nBias                          = ", round(o["lag.bias"],rd)))
  }


  if(is.element("maxN.n",names(o))){
    cat(paste0("\n\nMaximum sample size of ",maxN," (immediate outcomes):"))
    cat(paste0("\nAverage Observed Sample Size  = ", round(o["maxN.n"],rd)))
    cat(paste0("\nP( reject point null )        = ", o["maxN.rejPN"]))
    cat(paste0("\nP( conclude not ROPE effect)  = ", o["maxN.stopNotROPE"]))
    cat(paste0("\nP( conclude not ROME effect ) = ", o["maxN.stopNotROME"]))
    cat(paste0("\nP( conclude inconclusive )    = ", o["maxN.stopInconclusive"]))
    cat(paste0("\nCoverage                      = ", round(o["maxN.cover"],rd)))
    cat(paste0("\nBias                          = ", round(o["maxN.bias"],rd)))
  }


  if(is.element("lagMaxN.n",names(o))){
    cat(paste0("\n\nMaximum sample size of ",maxN," (stopping and then observing ", lag, " lagged outcomes):"))
    cat(paste0("\nAverage Total Sample Size     = ", round(o["lagMaxN.n"],rd)))
    cat(paste0("\nP( reject point null )        = ", o["lagMaxN.rejPN"]))
    cat(paste0("\nP( conclude not ROPE effect ) = ", o["lagMaxN.stopNotROPE"]))
    cat(paste0("\nP( conclude not ROME effect ) = ", o["lagMaxN.stopNotROME"]))
    cat(paste0("\nP( conclude inconclusive )    = ", o["lagMaxN.stopInconclusive"]))
    cat(paste0("\nCoverage                      = ", round(o["lagMaxN.cover"],rd)))
    cat(paste0("\nBias                          = ", round(o["lagMaxN.bias"],rd)))
  }

}

