#' @title addStats
#'
#' @description Obtains bias, whether the interval rejects H0, coverage, sgpvROPE, and sgpvROME for a given interval.  addStats is called within SeqSGPV.
#'
#' @param o Single mcmc replicate
#' @param randomize TRUE if allocation is vector has length > 1
#' @param effectPN See SeqSGPV
#' @param null See SeqSGPV
#' @param deltaL2 See SeqSGPV
#' @param deltaL1 See SeqSGPV
#' @param deltaG1 See SeqSGPV
#' @param deltaG2 See SeqSGPV
#'
#' @export
addStats <- function(o, randomize, effectPN, null, deltaL2, deltaL1, deltaG1, deltaG2){

  # Add whether coverage and bias
  if(randomize==TRUE){
    cover <- as.numeric(o[,"lo"] < o[,"effect1"] & o[,"effect1"] < o[,"up"])
    bias  <- o[,"est"] - o[,"effect1"]
  } else{
    cover <- as.numeric(o[,"lo"] < (o[,"theta0"] + o[,"effect0"]) & (o[,"theta0"] + o[,"effect0"]) < o[,"up"])
    bias  <- o[,"est"] - (o[,"theta0"] + o[,"effect0"])
  }


  # Add whether reject H0 and obtain sgpv
  if(null=="two.sided"){
    # Two sided

    rejH0   <- as.numeric(o[,"lo"] < effectPN & o[,"up"] < effectPN |
                          o[,"lo"] > effectPN & o[,"up"] > effectPN)

    sgpvROPE  <-     sgpv::sgpvalue(est.lo = o[,"lo"], est.hi = o[,"up"], null.lo = deltaL1, null.hi = deltaG1)$p.delta
    sgpvROME  <- 1 - sgpv::sgpvalue(est.lo = o[,"lo"], est.hi = o[,"up"], null.lo = deltaL2, null.hi = deltaG2)$p.delta

  } else if(null=="greater"){
    # One sided: H0: effect >= null

    rejH0 <- as.numeric(o[,"lo"] < effectPN & o[,"up"] < effectPN)

    # suppress warning that at least one interval has infinite length
    defaultW <- getOption("warn")
    options(warn = -1)

    sgpvROPE  <-     sgpv::sgpvalue(est.lo = o[,"lo"], est.hi = o[,"up"], null.lo = deltaL1, null.hi =     Inf)$p.delta
    sgpvROME  <-     sgpv::sgpvalue(est.lo = o[,"lo"], est.hi = o[,"up"], null.lo =    -Inf, null.hi = deltaL2)$p.delta

    options(warn = defaultW)

  } else if(null=="less"){
    # One sided: H0: effect <= null

    rejH0 <- as.numeric(o[,"lo"] > effectPN & o[,"up"] > effectPN)

    # suppress warning that at least one interval has inviting length
    defaultW <- getOption("warn")
    options(warn = -1)

    sgpvROPE  <-  sgpv::sgpvalue(est.lo = o[,"lo"], est.hi = o[,"up"], null.lo =    -Inf, null.hi = deltaG1)$p.delta
    sgpvROME  <-  sgpv::sgpvalue(est.lo = o[,"lo"], est.hi = o[,"up"], null.lo = deltaG2, null.hi =     Inf)$p.delta

    options(warn = defaultW)
  }



  # Return appended statistics
  cbind(o, bias, rejH0, cover, sgpvROPE, sgpvROME)

}

