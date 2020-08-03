#' @export
summary.sgpvAM <- function(am, affirmK, waitJ, treatEffect, rd = 4, unrestricted=TRUE,lag=TRUE,maxN=TRUE){


  # Inputs stored in am object
  if(is.element(el = "sgpvAMeffects",class(am))){
    amInputs <- am[[1]]$inputs
  } else {
    amInputs <- am$inputs
  }


  # Select treatment effect to summarize and reduce am object
  if(is.element(el = "sgpvAMeffects",class(am))){

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
    o  <- am[[paste0("effect_",te)]]

  } else if(!is.function(am$inputs$effectGeneration)) {

    te  <- amInputs$effectGeneration
    o   <- am
  } else {
    return(NULL)
  }


  # Select and reduce to wait time to summarize
  if(missing(waitJ)){
    cat(paste0(paste0("select wait time until applying monitoring rules:\n",
                      paste0(unname(sapply(names(o[["mcmcEndOfStudy"]]),
                                           FUN = function(x) as.numeric(strsplit(x,split="_")[[1]][2]))),
                             collapse=", "))))
    w <- readline(prompt="input: ")

  } else w <- waitJ

  o <- o[["mcmcEndOfStudy"]][[paste0("wait_",w)]][["mcmcEndOfStudyAve"]]


  # Select and reduce to number of affirmation steps required for stopping
  if(missing(affirmK)){
    cat(paste0(paste0("select required affirmation steps:\n",
                      paste0(o[,"affirmK"], collapse=", "))))
    k <- readline(prompt="input: ")

  } else k <- affirmK

  o <- o[o[,"affirmK"]==k,]






  if(is.na(amInputs$deltaL2) & is.na(amInputs$deltaL1)){
    H0label <- paste0("effect is less than or equal to ",am$inputs$effectPN)
  } else if(is.na(amInputs$deltaL2) & is.na(amInputs$deltaL1)){
    H0label <- paste0("effect is greater than or equal to ",am$inputs$effectPN)
  } else {
    H0label <- paste0("effect equals ",am$inputs$effectPN)
  }

  cat(paste0("\nGiven: effect = ", te, ", wait time = ", w, ", and k = ",k))
  cat(paste0("\nH0   : ",H0label))



  if(is.element("n",names(o)) & unrestricted==TRUE){
    cat("\n\no [Immediate outcomes] Unrestricted sample size")
    cat(paste0("\n  Average sample size              = ", round(o["n"],rd)))
    cat(paste0("\n  P( reject H0 )                   = ", round(o["rejH0"],rd)))
    cat(paste0("\n  P( conclude not ROPE effect )    = ", round(o["stopNotROPE"],rd)))
    cat(paste0("\n  P( conclude not ROME effect )    = ", round(o["stopNotROME"],rd)))
    cat(paste0("\n  P( conclude PRISM inconclusive ) = ", round(o["stopInconclusive"],rd)))
    cat(paste0("\n  Coverage                         = ", round(o["cover"],rd)))
    cat(paste0("\n  Bias                             = ", round(o["bias"],rd)))
  }


  if(is.element("maxN.n",names(o)) & maxN==TRUE){
    cat(paste0("\n\no [Immediate outcomes] Maximum sample size of ",amInputs$maxN))
    cat(paste0("\n  Average Observed Sample Size     = ", round(o["maxN.n"],rd)))
    cat(paste0("\n  P( reject H0 )                   = ", round(o["maxN.rejH0"],rd)))
    cat(paste0("\n  P( conclude not ROPE effect)     = ", round(o["maxN.stopNotROPE"],rd)))
    cat(paste0("\n  P( conclude not ROME effect )    = ", round(o["maxN.stopNotROME"],rd)))
    cat(paste0("\n  P( conclude PRISM inconclusive ) = ", round(o["maxN.stopInconclusive"],rd)))
    cat(paste0("\n  Coverage                         = ", round(o["maxN.cover"],rd)))
    cat(paste0("\n  Bias                             = ", round(o["maxN.bias"],rd)))
  }


  if(is.element("lag.n",names(o)) & lag==TRUE & unrestricted==TRUE){

    cat(paste0("\n\no [Stopping and then observing ", amInputs$lagN, " lagged outcomes] Unrestricted sample size"))
    cat(paste0("\n  Average Sample Size                                           = ", round(o["lag.n"],rd)))
    cat(paste0("\n  P( reject H0 )                                                = ", round(o["lag.rejH0"],rd)))
    cat(paste0("\n  P( conclude not ROPE effect )                                 = ", round(o["lag.stopNotROPE"],rd)))
    cat(paste0("\n  P( conclude not ROME effect )                                 = ", round(o["lag.stopNotROME"],rd)))
    cat(paste0("\n  P( conclude PRISM inconclusive )                              = ", round(o["lag.stopInconclusive"],rd)))
    cat(paste0("\n  P( conclusion on ROPE or ROME changed with lagged outcomes )  = ", round(o["lag.stopInconsistent"],rd)))
    cat(paste0("\n  P( newly reject H0 with lagged outcomes )                     = ", round(o["lag.stopRejH0_NY"],rd)))
    cat(paste0("\n  P( newly fail to reject H0 with lagged outcomes )             = ", round(o["lag.stopRejH0_YN"],rd)))
    cat(paste0("\n  Coverage                                                      = ", round(o["lag.cover"],rd)))
    cat(paste0("\n  Bias                                                          = ", round(o["lag.bias"],rd)))
  }


  if(is.element("lagMaxN.n",names(o)) & lag==TRUE & maxN==TRUE){
    cat(paste0("\n\no [Stopping and then observing ", amInputs$lagN, " lagged outcomes] Maximum sample size of ",amInputs$maxN))
    cat(paste0("\n  Average Total Sample Size                                     = ", round(o["lagMaxN.n"],rd)))
    cat(paste0("\n  P( reject H0 )                                                = ", round(o["lagMaxN.rejH0"],rd)))
    cat(paste0("\n  P( conclude not ROPE effect )                                 = ", round(o["lagMaxN.stopNotROPE"],rd)))
    cat(paste0("\n  P( conclude not ROME effect )                                 = ", round(o["lagMaxN.stopNotROME"],rd)))
    cat(paste0("\n  P( conclude PRISM inconclusive )                              = ", round(o["lagMaxN.stopInconclusive"],rd)))
    cat(paste0("\n  P( conclusion on ROPE or ROME changed with lagged outcomes )  = ", round(o["lagMaxN.stopInconsistent"],rd)))
    cat(paste0("\n  P( newly reject H0 with lagged outcomes )                     = ", round(o["lagMaxN.stopRejH0_NY"],rd)))
    cat(paste0("\n  P( newly fail to reject H0 with lagged outcomes )             = ", round(o["lagMaxN.stopRejH0_YN"],rd)))
    cat(paste0("\n  Coverage                                                      = ", round(o["lagMaxN.cover"],rd)))
    cat(paste0("\n  Bias                                                          = ", round(o["lagMaxN.bias"],rd)))
  }

  cat("\n")
}

