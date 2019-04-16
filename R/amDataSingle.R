#' amDataSingle
#'
#' Generate a single adaptive monitoring simulation
#'
#' @export
amDataSingle <- function(dataGeneration,   dataGenArgs,
                         effectGeneration, effectGenArgs,
                         modelFit,         modelFitArgs,
                         monitoringIntervalLevel,
                         existingData=NULL){


  # 0 Set arguments for dataGeneration and effectGeneration
  dataType <- class(modelFit)


  # 1 Generate data under null of no difference
  y <- do.call(dataGeneration,dataGenArgs)



  # 2 Generate treatment effect if not provided
  if( is.function(effectGeneration) ){
    theta <- do.call(effectGeneration,effectGenArgs)
  } else theta <- effectGeneration



  # 3 Randomly treat half with effect
  # Use block 2 randomization for sims but not in practice
  # Ensure enough treatments if generating an odd number of observations
  # Indicate the point null
  trt <- c(replicate(n=ceiling(length(y)/2), sample(c(0,1))))[1:dataGenArgs$n]


  if(dataType=="normal"){
    y[trt==1] <- y[trt==1] + theta
    pointNull <- 0
  } else if(dataType=="binomial"){
    oddsNull  <- dataGenArgs[["prob"]] / (1 - dataGenArgs[["prob"]])
    y[trt==1] <- rbinom(n = length(trt)/2,size = 1,prob = theta * oddsNull / (1 + theta * oddsNull))
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
                           modelFit, y=y, trt=trt, miLevel = monitoringIntervalLevel )))
  } else {
    eci <- rbind(matrix(rep(c(NA,-10^10,10^10),4),byrow = TRUE,nrow=4),
                 t(sapply(5:length(y), modelFit, y=y, trt=trt, miLevel = monitoringIntervalLevel )))
    colnames(eci) <- c("est", "lo", "up")
  }




  # 5 Return matrix of data, confidence interval, estimate, errors, and sgpvs
  cbind(theta, n = 1:length(y), y, trt, eci)

}
