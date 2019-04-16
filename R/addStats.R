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
# pointNull: specified point null
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



addStats <- function(o, pointNull, deltaL2, deltaL1, deltaG1, deltaG2){

  # Add whether reject point null, coverage, and bias
  rejPN <- as.numeric(o[,"lo"] < pointNull & o[,"up"] < pointNull |
                      o[,"lo"] > pointNull & o[,"up"] > pointNull)
  cover <- as.numeric(o[,"lo"] < o[,"theta"] & o[,"theta"] < o[,"up"])
  bias  <- o[,"est"] - o[,"theta"]


  # Obtain sgpv
  if(!anyNA(c(deltaL2, deltaL1, deltaG1, deltaG2))){

    # Two sided
    sgpvTrivial   <-     sgpv::sgpvalue(est.lo = o[,"lo"], est.hi = o[,"up"], null.lo = deltaL1, null.hi = deltaG1)$p.delta
    sgpvImpactful <- 1 - sgpv::sgpvalue(est.lo = o[,"lo"], est.hi = o[,"up"], null.lo = deltaL2, null.hi = deltaG2)$p.delta
  } else if(!anyNA(c(deltaL2, deltaL1))){

    # One sided: efficacy when less than null
    sgpvTrivial   <-     sgpv::sgpvalue(est.lo = o[,"lo"], est.hi = o[,"up"], null.lo =  deltaL1, null.hi =   10^10)$p.delta
    sgpvImpactful <-     sgpv::sgpvalue(est.lo = o[,"lo"], est.hi = o[,"up"], null.lo =   -10^10, null.hi = deltaL2)$p.delta
  } else if(!anyNA(c(deltaG1, deltaG2))){

    # One sided: efficacy when greater than null
    sgpvTrivial   <-     sgpv::sgpvalue(est.lo = o[,"lo"], est.hi = o[,"up"], null.lo =  -10^10, null.hi = deltaG1)$p.delta
    sgpvImpactful <-     sgpv::sgpvalue(est.lo = o[,"lo"], est.hi = o[,"up"], null.lo = deltaG2, null.hi =   10^10)$p.delta
  } else{

    stop("A one sided study requires both deltas to be strictly greater or lower than point null")
  }


  # Return appended statistics
  cbind(o, bias, rejPN, cover, sgpvTrivial, sgpvImpactful)

}

