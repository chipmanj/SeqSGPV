#' amDataSingle
#'
#' Generate a single adaptive monitoring simulation
#'
#' @export
amDataSingle <- function(dataGeneration,   dataGenArgs,
                         effectGeneration, effectGenArgs,
                         effectScale="identity",
                         randomize=FALSE,
                         modelFit,         modelFitArgs,
                         monitoringIntervalLevel,
                         existingData=NULL){


  # 0 Set arguments for dataGeneration and effectGeneration
  if(any(grepl("rbinom",deparse(dataGeneration)))){
    dataType <- "dichotomous"
  } else {
    dataType <- NA
  }


  # 1 Generate data under null of no difference
  if(!is.na(dataType) & dataType == "dichotomous"){
    rUniform <- runif(dataGenArgs$n)
    y        <- as.numeric(rUniform <= dataGenArgs$prob)
  } else {
    rUniform <- NA
    y <- do.call(dataGeneration,dataGenArgs)
  }



  # 2 Generate treatment effect if not provided
  if( is.function(effectGeneration) ){
         theta <- do.call(effectGeneration,effectGenArgs)
  } else theta <- effectGeneration


  # 2.1 Transform treatment effect to identity scale
  #     Example, if effect is on the odds ratio scale, transform to probability
  if( toupper(effectScale) == "IDENTITY" ){
    thetaIdentity <- theta
  } else if( dataType=="dichotomous"){
    # transform to probabilities
    if(toupper(effectScale) %in% c("ODDSRATIO","OR")){
      dataOR <- dataGenArgs$prob / (1 - dataGenArgs$prob)
      thetaIdentity <- (dataOR * theta) / (1 + dataOR * theta) - dataGenArgs$prob
    } else if (toupper(effectScale) %in% c("ODDS")){
      thetaIdentity <- theta / (1 + theta)
    }
  } else if(toupper(effectScale) %in% c("LOG")){
      thetaIdentity <- exp(theta)
  }




  # 3 Set treatment variable
  # If randomize == FALSE, then all patients will have effect shifted by theta
  # Else randomly treat half with effect who will have effect shifted by theta
  # Use block 2 randomization for sims but not in practice
  # Ensure enough treatments if generating an odd number of observations

  trt <- rep(1,length(y))
  if(randomize) trt[seq(2,length(y),2)] <- 0

  if(!is.na(dataType) & dataType=="dichotomous"){
    y[trt==1] <- as.numeric( as.numeric(rUniform[trt==1] <= dataGenArgs$prob + thetaIdentity) )
  } else {
    y[trt==1] <- y[trt==1] + thetaIdentity
  }


  # 3.5 Add to previous data if any (ie previous data was insufficient)
  if(! is.null(existingData) ) {
    y   <- c(existingData[,"y"], y)
    trt <- c(existingData[,"trt"], trt)
  }


  # 3.75 Create design matrix
  if(randomize==FALSE){
    XD <- as.matrix(trt,ncol=1)
  } else {
    XD <- as.matrix(cbind(1,trt),ncol=2)
  }

  # 4 Obtain (or add to) estimate and 1-alpha/2 monitoring confidence intervals
  #   Wait until at least two observations in each group
  if(! is.null(existingData) ) {
    eci <- rbind(existingData[,c("est","lo","up")],
                 t(sapply( (length(y)-dataGenArgs[["n"]] + 1):length(y),
                           modelFit, y=y, XD=XD, miLevel = monitoringIntervalLevel )))
  } else {
    eci <- rbind(matrix(rep(c(NA,-10^10,10^10),4),byrow = TRUE,nrow=4),
                 t(sapply(5:length(y), modelFit, y=y, XD=XD, miLevel = monitoringIntervalLevel )))
    colnames(eci) <- c("est", "lo", "up")
  }




  # 5 Return matrix of data, confidence interval, estimate, errors, and sgpvs
  cbind(theta, n = 1:length(y), y, rUniform, trt, eci)

}
