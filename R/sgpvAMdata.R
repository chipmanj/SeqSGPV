# sgpvAMdata.R
# J Chipman
#
# Generates data, estimates fully sequential confidence intervals,
#  and calculates sgpv.
#
# Currently fits t.test and logistic model
#
# dataGeneration: function for generating data (such as rnorm and rbinom)
# dataGenArgs   : arguments to be included when generating data (such as
#                 n, mean, sd for rnorm and n, size, and prob for rbinom)
# effectGeneration: either a single numeric for a exploring a fixed effect
#                   or a function for generating a single treatment effect
# effectGenArgs:    arguments to be included when generating a treatment effect
# modelFit:         a user provided model fit for estimating the treatment effect
#                   must return a 2 length vector of interval lower and upper bound.
#                   If missing, will default to t.test for normal data and glm for
#                   binomial data.
# modelFitArgs:     Arguments for modelFit
# deltaL1:          The delta less than and closest to the point null.
# deltaL2:          The delta less than and furthest from the point null.
# deltaG1:          The delta greater than and closest to the point null.
# deltaG2:          The delta greater than and furthest from the point null.
# monitoringIntervalLevel: The traditional alpha in the (1-alpha) monitoring interval.
#                   These intervals are purely for monitoring, and we do not report the
#                   frequency properties of monitoring intervals.
# existingData:     If previous generation of data were insufficient to follow to completion
#                   they can be supplied such that current call of sgpvAMdata appends to
#                   previous data.


sgpvAMdataSingle <- function(dataGeneration,   dataGenArgs,
                             effectGeneration, effectGenArgs,
                             modelFit,         modelFitArgs,
                             deltaL2, deltaL1, deltaG1, deltaG2,
                             monitoringIntervalLevel,
                             existingData=NULL, ...){


  # 0 Set arguments for dataGeneration and effectGeneration
  dataType <- deparse(substitute(dataGeneration))


  # 1 Generate data under null of no difference
  y <- do.call(dataGeneration,dataGenArgs)



  # 2 Generate treatment effect if not provided
  if( is.function(effectGeneration) ){
         z <- do.call(effectGeneration,effectGenArgs)
  } else z <- effectGeneration



  # 3 Randomly treat half with effect
  # Use block 2 randomization for sims but not in practice
  # Enough enough treatments if generating an odd number of observations
  trt <- c(replicate(n=ceiling(length(y)/2), sample(c(0,1))))[1:dataGenArgs$n]


  if(dataType=="rnorm"){
    y[trt==1] <- y[trt==1] + z
  } else if(dataType=="rbinom"){
    oddsNull  <- dataGenArgs[["prob"]] / (1 - dataGenArgs[["prob"]])
    y[trt==1] <- rbinom(n = length(trt)/2,size = 1,prob = z * oddsNull / (1 + z * oddsNull))
  }


  # 3.5 Add to previous data if any (ie previous data was insufficient)
  if(! is.null(existingData) ) {
    y   <- c(existingData[,"y"], y)
    trt <- c(existingData[,"trt"], trt)
  }



  # 4 Obtain 1-alpha/2 monitoring confidence intervals
  # Wait until at least two observations in each group
  fullySequentialCIs <- function(look, miLevel){

    # For first for looks provide essentially infinite intervals
    if(look < 4){ return( c(NA, -10^10, 10^10, NA, NA) ) }

    if(dataType=="rnorm"){
      f     <- lm(y[1:look] ~ trt[1:look])
      coefs <- summary(f)$coefficients
      est   <- coefs[2,"Estimate"]
      ci    <- est + c(-1,1) * qt(1-miLevel/2, df = f$df.residual) * coefs[2,"Std. Error"]
      rejPN <- as.numeric(ci[1] < 0 & ci[2] < 0 | ci[1] > 0 & ci[2] > 0)

    } else if(dataType=="rbinom"){
      f     <- glm(y[1:look] ~ trt[1:look], family = binomial)
      coefs <- summary(f)$coefficients
      est   <- exp( coefs[2,"Estimate"] )
      ci    <- exp( coefs[2,"Estimate"] + c(-1,1) * qnorm(1-miLevel/2) * coefs[2,"Std. Error"] )
      # Infinite CI bounds may occur with binomial data and insufficient
      #  data to estimate an effect and error.  Set infinite bounds to 10^10
      ci[is.infinite(ci)] <- 10^10
      rejPN <- as.numeric(ci[1] < 1 & ci[2] < 1 | ci[1] > 1 & ci[2] > 1)
    }

    cover <- as.numeric(ci[1] < z & z < ci[2])


    return(c(est,ci,rejPN,cover))

  }

  if(! is.null(existingData) ) {
    eci <- rbind(existingData[,c("est","lo","hi","rejPN","cover")],
                t(sapply( (length(y)-dataGenArgs[["n"]] + 1):length(y),
                          fullySequentialCIs, miLevel = monitoringIntervalLevel)))
  } else {
    eci <- t(sapply(1:length(y),fullySequentialCIs, miLevel = monitoringIntervalLevel))
  }
  colnames(eci) <- c("est","lo","hi","rejPN","cover")



  # 5 Obtain sgpv
  if(!anyNA(c(deltaL2, deltaL1, deltaG1, deltaG2))){

    # Two sided
    sgpvNonTrivial <- 1 - sgpv::sgpvalue(est.lo = eci[,"lo"], est.hi = eci[,"hi"], null.lo = deltaL1, null.hi = deltaG1)$p.delta
    sgpvFutility   <-     sgpv::sgpvalue(est.lo = eci[,"lo"], est.hi = eci[,"hi"], null.lo = deltaL2, null.hi = deltaG2)$p.delta
  } else if(!anyNA(c(deltaL2, deltaL1))){

    # One sided: efficacy when less than null
    sgpvNonTrivial <- 1 - sgpv::sgpvalue(est.lo = eci[,"lo"], est.hi = eci[,"hi"], null.lo =  -10^10, null.hi = deltaL1)$p.delta
    sgpvFutility   <-     sgpv::sgpvalue(est.lo = eci[,"lo"], est.hi = eci[,"hi"], null.lo = deltaL2, null.hi = 10^10)$p.delta
  } else if(!anyNA(c(deltaG1, deltaG2))){

    # One sided: efficacy when less than null
    sgpvNonTrivial <- 1 - sgpv::sgpvalue(est.lo = eci[,"lo"], est.hi = eci[,"hi"], null.lo = deltaG1, null.hi = 10^10)$p.delta
    sgpvFutility   <-     sgpv::sgpvalue(est.lo = eci[,"lo"], est.hi = eci[,"hi"], null.lo =  -10^10, null.hi = deltaG2)$p.delta
  } else{

    stop("A one sided study requires both deltas to be strictly greater or lower than point null")
  }




  # 6 Return matrix of data, confidence interval, estimate, errors, and sgpvs
  cbind(n = 1:length(y), y, trt, eci, bias=eci[,"est"] - z, sgpvNonTrivial, sgpvFutility, z)

}


sgpvAMdata <- function(nreps, ...){
  mcmcMonitoring <- plyr::rlply(.n = nreps, .expr = { sgpvAMdataSingle( ... ) })
}


# library(doParallel)


sgpvAMdataPP <- function(nreps, ...){

  cl <- makeCluster(detectCores())
  registerDoParallel(cl)

  mcmcMonitoring <- foreach(i = 1:nreps) %dopar% {
    sgpvAMdataSingle( ... )
  }

  # i       <- seq_len(nreps)
  # fe_call <- as.call(c(list(quote(foreach::foreach), i = i)))
  # fe      <- eval(fe_call)
  #
  # result <- foreach::`%dopar%`(fe, sgpvAMdataSingle( ... ))
}


# a <- sgpvAMdataPP(nreps=10, rnorm, dataGenArgs = list(n=70), effectGeneration = 0.5,
#                   deltaL2 = -0.5, deltaL1 = -0.15, deltaG1 = 0.15, deltaG2 = 0.5,
#                   monitoringIntervalLevel = 0.05)


# else {
#   result <- loop_apply(n, do.ply)
# }
# cl <- makeCluster(detectCores())
# registerDoParallel(cl)
#
# mcmcMonitoring <- foreach(i = 1:nreps) %dopar% {
#   sgpvAMdataSingle( ... )
# }
#
# stopCluster(cl)
# mcmcMonitoring
#
# }


# Examples
mcmcMonotiringFixed <- sgpvAMdata(rnorm, dataGenArgs = list(n=70), effectGeneration = 0.5,
                                  deltaL2 = -0.5, deltaL1 = -0.15, deltaG1 = 0.15, deltaG2 = 0.5,
                                  monitoringIntervalLevel = 0.05)
# head(mcmcMonotiringFixed)
#
# mcmcMonotiringNorm <- sgpvAMdata(rnorm, dataGenArgs = list(n=700),
#                                  effectGeneration = rnorm, effectGenArgs = list(n=1),
#                                   deltaL2 = -0.5, deltaL1 = -0.15, deltaG1 = 0.15, deltaG2 = 0.5,
#                                   monitoringIntervalLevel = 0.05)
# head(mcmcMonotiringNorm)


# development
# dataGeneration <- rnorm
# dataGenArgs <- list(n=700)
# effectGeneration <- 0
# deltaL2 <- -0.5
# deltaL1 <- -0.15
# deltaG1 <-  0.15
# deltaG2 <-  0.5
# monitoringIntervalLevel <- 0.05
