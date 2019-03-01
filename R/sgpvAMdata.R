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
                             existingData=NULL){


  # 0 Set arguments for dataGeneration and effectGeneration
  dataType <- class(modelFit)


  # 1 Generate data under null of no difference
  y <- do.call(dataGeneration,dataGenArgs)



  # 2 Generate treatment effect if not provided
  if( is.function(effectGeneration) ){
         z <- do.call(effectGeneration,effectGenArgs)
  } else z <- effectGeneration



  # 3 Randomly treat half with effect
  # Use block 2 randomization for sims but not in practice
  # Ensure enough treatments if generating an odd number of observations
  # Indicate the point null
  trt <- c(replicate(n=ceiling(length(y)/2), sample(c(0,1))))[1:dataGenArgs$n]


  if(dataType=="normal"){
    y[trt==1] <- y[trt==1] + z
    pointNull <- 0
  } else if(dataType=="binomial"){
    oddsNull  <- dataGenArgs[["prob"]] / (1 - dataGenArgs[["prob"]])
    y[trt==1] <- rbinom(n = length(trt)/2,size = 1,prob = z * oddsNull / (1 + z * oddsNull))
    pointNull <- 1
  }


  # 3.5 Add to previous data if any (ie previous data was insufficient)
  if(! is.null(existingData) ) {
    y   <- c(existingData[,"y"], y)
    trt <- c(existingData[,"trt"], trt)
  }


  # 4 Obtain (or add to) 1-alpha/2 monitoring confidence intervals
  #   Wait until at least two observations in each group
  if(! is.null(existingData) ) {
    eci <- rbind(existingData[,c("est","lo","up")],
                t(sapply( (length(y)-dataGenArgs[["n"]] + 1):length(y),
                          modelFit, y=y, trt=trt, miLevel = monitoringIntervalLevel, ... )))
  } else {
    eci <- rbind(matrix(rep(c(NA,-10^10,10^10),4),byrow = TRUE,nrow=4),
                 t(sapply(5:length(y), modelFit, y=y, trt=trt, miLevel = monitoringIntervalLevel, ... )))
    colnames(eci) <- c("est", "lo", "up")
  }


  # Add whether reject point null, coverage, and bias
  rejPN <- as.numeric(eci[,"lo"] < pointNull & eci[,"up"] < pointNull |
                      eci[,"lo"] > pointNull & eci[,"up"] > pointNull)
  cover <- as.numeric(eci[,"lo"] < z & z < eci[,"up"])
  bias  <- eci[,"est"] - z
  eci   <- cbind(eci,bias,rejPN,cover)




  # 5 Obtain sgpv
  if(!anyNA(c(deltaL2, deltaL1, deltaG1, deltaG2))){

    # Two sided
    sgpvNonTrivial <- 1 - sgpv::sgpvalue(est.lo = eci[,"lo"], est.hi = eci[,"up"], null.lo = deltaL1, null.hi = deltaG1)$p.delta
    sgpvFutility   <-     sgpv::sgpvalue(est.lo = eci[,"lo"], est.hi = eci[,"up"], null.lo = deltaL2, null.hi = deltaG2)$p.delta
  } else if(!anyNA(c(deltaL2, deltaL1))){

    # One sided: efficacy when less than null
    sgpvNonTrivial <- 1 - sgpv::sgpvalue(est.lo = eci[,"lo"], est.hi = eci[,"up"], null.lo =  -10^10, null.hi = deltaL1)$p.delta
    sgpvFutility   <-     sgpv::sgpvalue(est.lo = eci[,"lo"], est.hi = eci[,"up"], null.lo = deltaL2, null.hi = 10^10)$p.delta
  } else if(!anyNA(c(deltaG1, deltaG2))){

    # One sided: efficacy when less than null
    sgpvNonTrivial <- 1 - sgpv::sgpvalue(est.lo = eci[,"lo"], est.hi = eci[,"up"], null.lo = deltaG1, null.hi = 10^10)$p.delta
    sgpvFutility   <-     sgpv::sgpvalue(est.lo = eci[,"lo"], est.hi = eci[,"up"], null.lo =  -10^10, null.hi = deltaG2)$p.delta
  } else{

    stop("A one sided study requires both deltas to be strictly greater or lower than point null")
  }




  # 6 Return matrix of data, confidence interval, estimate, errors, and sgpvs
  cbind(n = 1:length(y), y, trt, eci, sgpvNonTrivial, sgpvFutility, z)

}


sgpvAMdata <- function(nreps, fork=TRUE, socket = TRUE, cores = detectCores(), ...){

  if(fork==TRUE){
    # Only works on POSIX systems (Mac, Linux, Unix, BSD) and not Windows.
    mcmcMonitoring <- parallel::mclapply(1:nreps, sgpvAMdataSingle, ... , mc.cores = cores)
  } else if(socket==TRUE){
    # Works on Mac and Windows; only slightly slower
    cl             <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    mcmcMonitoring <- parallel::parLapply(cl, 1:nreps, sgpvAMdataSingle, ...)
  } else {
    mcmcMonitoring <- plyr::rlply(.n = nreps, .expr = { sgpvAMdataSingle( ... ) })
  }

  return(mcmcMonitoring)
}




# system.time(a <- sgpvAMdata(nreps=1000, dataGeneration = rnorm, dataGenArgs = list(n=1600), effectGeneration = 0.5,
#                             modelFit = lmCI,
#                             deltaL2 = -0.5, deltaL1 = -0.15, deltaG1 = 0.15, deltaG2 = 0.5,
#                             monitoringIntervalLevel = 0.05))
#
#
# a <- sgpvAMdataSingle(dataGeneration = rnorm, dataGenArgs = list(n=120), effectGeneration = 0.5,
#                       modelFit = lmCI,
#                       deltaL2 = -0.5, deltaL1 = -0.15, deltaG1 = 0.15, deltaG2 = 0.5,
#                       monitoringIntervalLevel = 0.05)


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
# mcmcMonotiringFixed <- sgpvAMdata(rnorm, dataGenArgs = list(n=70), effectGeneration = 0.5,
#                                   deltaL2 = -0.5, deltaL1 = -0.15, deltaG1 = 0.15, deltaG2 = 0.5,
#                                   monitoringIntervalLevel = 0.05)
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
