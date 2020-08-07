#' @export

# addStats
#
# adds statistics following each new observation
# - reject point null
# - bias
# - coverage
# - sgpv Non Trivial and Futility
#
# Inputs
# o        : matrix with at least columns: theta, est, lo (lower CI bound), and up (upper CI bound)
# effectPN : specified effect point null
# deltaL2  : lower delta furthest from point null
# deltaL1  : lower delta closest    to point null
# deltaG1  : upper delta furthest from point null
# deltaG2  : upper delta closest    to point null
#
# Returns
# o with new columns appended
# - rejectPN (reject point null)
# - bias
# - cover (does CI cover theta)
# - sgpvTrivial
# - sgpvImpactful



addStats <- function(o, effectPN, deltaL2, deltaL1, deltaG1, deltaG2){

  # Add whether coverage and bias
  cover <- as.numeric(o[,"lo"] < o[,"effect1"] & o[,"effect1"] < o[,"up"])
  bias  <- o[,"est"] - o[,"effect1"]


  # Add whether reject H0 and obtain sgpv
  if(!anyNA(c(deltaL2, deltaL1, deltaG1, deltaG2))){
    # Two sided

    rejH0   <- as.numeric(o[,"lo"] < effectPN & o[,"up"] < effectPN |
                          o[,"lo"] > effectPN & o[,"up"] > effectPN)

    sgpvROPE  <-     sgpv::sgpvalue(est.lo = o[,"lo"], est.hi = o[,"up"], null.lo = deltaL1, null.hi = deltaG1)$p.delta
    sgpvROME  <- 1 - sgpv::sgpvalue(est.lo = o[,"lo"], est.hi = o[,"up"], null.lo = deltaL2, null.hi = deltaG2)$p.delta

  } else if(!anyNA(c(deltaL2, deltaL1))){
    # One sided: H0: effect >= null

    rejH0 <- as.numeric(o[,"lo"] < effectPN & o[,"up"] < effectPN)

    # suppress warning that at least one interval has inviting length
    defaultW <- getOption("warn")
    options(warn = -1)

    sgpvROPE  <-     sgpv::sgpvalue(est.lo = o[,"lo"], est.hi = o[,"up"], null.lo = deltaL1, null.hi =     Inf)$p.delta
    sgpvROME  <-     sgpv::sgpvalue(est.lo = o[,"lo"], est.hi = o[,"up"], null.lo =    -Inf, null.hi = deltaL2)$p.delta

    options(warn = defaultW)

  } else if(!anyNA(c(deltaG1, deltaG2))){
    # One sided: H0: effect <= null

    rejH0 <- as.numeric(o[,"lo"] > effectPN & o[,"up"] > effectPN)

    # suppress warning that at least one interval has inviting length
    defaultW <- getOption("warn")
    options(warn = -1)

    sgpvROPE  <-     sgpv::sgpvalue(est.lo = o[,"lo"], est.hi = o[,"up"], null.lo =    -Inf, null.hi = deltaG1)$p.delta
    sgpvROME  <-     sgpv::sgpvalue(est.lo = o[,"lo"], est.hi = o[,"up"], null.lo = deltaG2, null.hi =     Inf)$p.delta

    options(warn = defaultW)
  } else{

    stop("A one sided study requires both deltas to be strictly greater or lower than point null")
  }


  # Return appended statistics
  cbind(o, bias, rejH0, cover, sgpvROPE, sgpvROME)

}

