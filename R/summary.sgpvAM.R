#' @export
summary.sgpvAM <- function(am, alertK, waitTime, treatEffect, rd = 4){


  # Inputs stored in am object
  amInputs <- am$inputs


  # Select treatment effect to summarize and reduce am object
  if(is.element(el = "sgpvAMthetas",class(am))){

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





  cat(paste0("\nGiven: theta = ", te, ", wait time = ", w, ", and k = ",k))

  if(is.element("n",names(o))){
    cat("\n\no [Immediate outcomes] Unrestricted sample size")
    cat(paste0("\n  Average sample size           = ", round(o["n"],rd)))
    cat(paste0("\n  P( reject point null )        = ", round(o["rejPN"],rd)))
    cat(paste0("\n  P( conclude not ROPE effect ) = ", round(o["stopNotROPE"],rd)))
    cat(paste0("\n  P( conclude not ROME effect ) = ", round(o["stopNotROME"],rd)))
    cat(paste0("\n  P( conclude inconclusive )    = ", round(o["stopInconclusive"],rd)))
    cat(paste0("\n  Coverage                      = ", round(o["cover"],rd)))
    cat(paste0("\n  Bias                          = ", round(o["bias"],rd)))
  }


  if(is.element("maxN.n",names(o))){
    cat(paste0("\n\no [Immediate outcomes] Maximum sample size of ",maxN))
    cat(paste0("\n  Average Observed Sample Size  = ", round(o["maxN.n"],rd)))
    cat(paste0("\n  P( reject point null )        = ", round(o["maxN.rejPN"],rd)))
    cat(paste0("\n  P( conclude not ROPE effect)  = ", round(o["maxN.stopNotROPE"],rd)))
    cat(paste0("\n  P( conclude not ROME effect ) = ", round(o["maxN.stopNotROME"],rd)))
    cat(paste0("\n  P( conclude inconclusive )    = ", round(o["maxN.stopInconclusive"],rd)))
    cat(paste0("\n  Coverage                      = ", round(o["maxN.cover"],rd)))
    cat(paste0("\n  Bias                          = ", round(o["maxN.bias"],rd)))
  }


  if(is.element("lag.n",names(o))){

    cat(paste0("\n\no [Stopping and then observing ", lag, " lagged outcomes] Unrestricted sample size"))
    cat(paste0("\n  Average Sample Size                           = ", round(o["lag.n"],rd)))
    cat(paste0("\n  P( reject point null )                        = ", round(o["lag.rejPN"],rd)))
    cat(paste0("\n  P( conclude not ROPE effect )                 = ", round(o["lag.stopNotROPE"],rd)))
    cat(paste0("\n  P( conclude not ROME effect )                 = ", round(o["lag.stopNotROME"],rd)))
    cat(paste0("\n  P( conclude inconclusive )                    = ", round(o["lag.stopInconclusive"],rd)))
    cat(paste0("\n  P( conclusion changed with lagged outcomes )  = ", round(o["lag.stopInconsistent"],rd)))
    cat(paste0("\n  Coverage                                      = ", round(o["lag.cover"],rd)))
    cat(paste0("\n  Bias                                          = ", round(o["lag.bias"],rd)))
  }


  if(is.element("lagMaxN.n",names(o))){
    cat(paste0("\n\no [Stopping and then observing ", lag, " lagged outcomes] Maximum sample size of ",maxN))
    cat(paste0("\n  Average Total Sample Size                     = ", round(o["lagMaxN.n"],rd)))
    cat(paste0("\n  P( reject point null )                        = ", round(o["lagMaxN.rejPN"],rd)))
    cat(paste0("\n  P( conclude not ROPE effect )                 = ", round(o["lagMaxN.stopNotROPE"],rd)))
    cat(paste0("\n  P( conclude not ROME effect )                 = ", round(o["lagMaxN.stopNotROME"],rd)))
    cat(paste0("\n  P( conclude inconclusive )                    = ", round(o["lagMaxN.stopInconclusive"],rd)))
    cat(paste0("\n  P( conclusion changed with lagged outcomes )  = ", round(o["lagMaxN.stopInconsistent"],rd)))
    cat(paste0("\n  Coverage                                      = ", round(o["lagMaxN.cover"],rd)))
    cat(paste0("\n  Bias                                          = ", round(o["lagMaxN.bias"],rd)))
  }

}

