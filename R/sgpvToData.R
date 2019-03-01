# sgpvToData.R
# J Chipman
#
# Adds sgpv to data which includes interval endpoints: "lo" and "up"
#
# data:       matrix from sgpvAMdata.  Or a similar matrix with columns:
#             n, y, lo, hi, z
# deltaL1:    The delta less than and closest to the point null.
# deltaL2:    The delta less than and furthest from the point null.
# deltaG1:    The delta greater than and closest to the point null.
# deltaG2:    The delta greater than and furthest from the point null.
#
# Returns original dataframe with two new columns: sgpvNonTrivial and sgpvFutility


sgpvToData <- function(data, deltaL1=NA, deltaL2=NA,deltaG1=NA,deltaG2=NA){


  # 1 Obtain sgpv
  if(!anyNA(c(deltaL2, deltaL1, deltaG1, deltaG2))){

    # Two sided
    sgpvNonTrivial <- 1 - sgpv::sgpvalue(est.lo = data[,"lo"], est.hi = data[,"up"], null.lo = deltaL1, null.hi = deltaG1)$p.delta
    sgpvFutility   <-     sgpv::sgpvalue(est.lo = data[,"lo"], est.hi = data[,"up"], null.lo = deltaL2, null.hi = deltaG2)$p.delta
  } else if(!anyNA(c(deltaL2, deltaL1))){

    # One sided: efficacy when less than null
    sgpvNonTrivial <- 1 - sgpv::sgpvalue(est.lo = data[,"lo"], est.hi = data[,"up"], null.lo =  -10^10, null.hi = deltaL1)$p.delta
    sgpvFutility   <-     sgpv::sgpvalue(est.lo = data[,"lo"], est.hi = data[,"up"], null.lo = deltaL2, null.hi = 10^10)$p.delta
  } else if(!anyNA(c(deltaG1, deltaG2))){

    # One sided: efficacy when less than null
    sgpvNonTrivial <- 1 - sgpv::sgpvalue(est.lo = data[,"lo"], est.hi = data[,"up"], null.lo = deltaG1, null.hi = 10^10)$p.delta
    sgpvFutility   <-     sgpv::sgpvalue(est.lo = data[,"lo"], est.hi = data[,"up"], null.lo =  -10^10, null.hi = deltaG2)$p.delta
  } else{

    stop("A one sided study requires both deltas to be strictly greater or lower than point null")
  }

  data <- cbind(data,sgpvNonTrivial,sgpvFutility)


}
