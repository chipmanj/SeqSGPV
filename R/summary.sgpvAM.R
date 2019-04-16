summary.sgpvAM <- function(am, alertK, waitTime, treatEffect, rd = 4){




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

    te <- am$inputs$effectGeneration
    o  <- am
  }


  # Select and reduce to wait time to summarize
  if(missing(waitTime)){
    cat(paste0(paste0("select wait time until confidence interval width of:\n",
                      paste0(unname(sapply(names(o[["mcmcEndOfStudy"]]),
                                           FUN = function(x) as.numeric(strsplit(x,split="_")[[1]][2]))),
                             collapse=", "))))
    w <- readline(prompt="input: ")

  } else w <- waitTime

  o <- o[["mcmcEndOfStudy"]][[paste0("width_",w)]][["mcmcEndOfStudyAve"]]


  # Select and reduce to number of affirmation steps required for stopping
  if(missing(alertK)){
    cat(paste0(paste0("select wait time until confidence interval width of:\n",
                      paste0(o[,"alertK"], collapse=", "))))
    k <- readline(prompt="input: ")

  } else k <- 50

  o <- o[o[,"alertK"]==k,]





  cat(paste0("Given: theta = ", te, ", wait CI width = ", w, ", and k = ",k))

  cat("\n\nUnrestricted Sample Size with immediate outcomes:")
  cat(paste0("\nP( reject point null )             = ", o["rejPN"]))
  cat(paste0("\nP( conclude not trivial effect )   = ", o["stopNotTrivial"]))
  cat(paste0("\nP( conclude not impactful effect ) = ", o["stopNotImpactful"]))
  cat(paste0("\nP( conclude inconclusive )         = ", o["stopInconclusive"]))
  cat(paste0("\nCoverage                           = ", round(o["cover"],rd)))
  cat(paste0("\nBias                               = ", round(o["bias"],rd)))


  if(length((grep("lag.",names(o))))>0){
    cat("\n\nUnrestricted Sample Size, stopping and then observing the lagged outcomes:")
    cat(paste0("\nP( reject point null )             = ", o["lag.rejPN"]))
    cat(paste0("\nP( conclude not trivial effect )   = ", o["lag.stopNotTrivial"]))
    cat(paste0("\nP( conclude not impactful effect ) = ", o["lag.stopNotImpactful"]))
    cat(paste0("\nP( conclude inconclusive )         = ", o["lag.stopInconclusive"]))
    cat(paste0("\nCoverage                           = ", round(o["lag.cover"],rd)))
    cat(paste0("\nBias                               = ", round(o["lag.bias"],rd)))
  }


  if(length((grep("maxN.",names(o))))>0){
    cat("\n\nMaximum Sample Size with immediate outcomes:")
    cat(paste0("\nP( reject point null )             = ", o["maxN.rejPN"]))
    cat(paste0("\nP( conclude not trivial effect )   = ", o["maxN.stopNotTrivial"]))
    cat(paste0("\nP( conclude not impactful effect ) = ", o["maxN.stopNotImpactful"]))
    cat(paste0("\nP( conclude inconclusive )         = ", o["maxN.stopInconclusive"]))
    cat(paste0("\nCoverage                           = ", round(o["maxN.cover"],rd)))
    cat(paste0("\nBias                               = ", round(o["maxN.bias"],rd)))
  }


  if(length((grep("lagMaxN.",names(o))))>0){
    cat("\n\nMaximum Sample Size, stopping and then observing the lagged outcomes:")
    cat(paste0("\nP( reject point null )             = ", o["lagMaxN.rejPN"]))
    cat(paste0("\nP( conclude not trivial effect )   = ", o["lagMaxN.stopNotTrivial"]))
    cat(paste0("\nP( conclude not impactful effect ) = ", o["lagMaxN.stopNotImpactful"]))
    cat(paste0("\nP( conclude inconclusive )         = ", o["lagMaxN.stopInconclusive"]))
    cat(paste0("\nCoverage                           = ", round(o["lagMaxN.cover"],rd)))
    cat(paste0("\nBias                               = ", round(o["lagMaxN.bias"],rd)))
  }

}

