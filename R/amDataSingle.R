#' @title amDataSingle
#'
#' @description Generate a single adaptive monitoring simulation.  amDataSingle is called within SeqSGPV.
#'
#' @param dataGeneration See SeqSGPV
#' @param dataGenArgs See SeqSGPV
#' @param effectGeneration See SeqSGPV
#' @param effectGenArgs See SeqSGPV
#' @param effectScale See SeqSGPV
#' @param allocation See SeqSGPV
#' @param randomize TRUE if length(allocation) > 1
#' @param modelFit See SeqSGPV
#' @param modelFitArgs See SeqSGPV
#' @param existingData Previously provided data
#'
#' @export
amDataSingle <- function(dataGeneration,   dataGenArgs,
                         effectGeneration, effectGenArgs,
                         effectScale,
                         allocation,
                         randomize,
                         modelFit,
                         modelFitArgs,
                         existingData){


  # 0.1 if data type is dichotomous, use uniform data with cutoffs for dataGeneration
  #     Note: uniform data makes it easier to take the same data and impose a treatment effect (such as on odds ratio scale)
  # 0.2 set theta0 from data generation arguments
  if(any(grepl("rbinom",deparse(dataGeneration)))){
    dichotomous <- TRUE
    theta0      <- dataGenArgs$prob

  } else if (any(grepl("rnorm",deparse(dataGeneration)))){
    dichotomous <- FALSE
      if(is.element("mean",names(dataGenArgs))) {
        theta0 <- dataGenArgs$mean
      } else {
        theta0 <- 0
      }

  } else if( any(grepl("runif",deparse(dataGeneration))) ) {
    dichotomous <- FALSE
    theta0      <- mean(c(dataGenArgs$min,dataGenArgs$max))

  } else {
    dichotomous <- NA
    theta0      <- NA

  }


  # 1 Generate data under null of no difference
  #  -- If binomiail generate random uniform for y == 0 and y == 1
  #  -- To be used when shifting treatment effects
  y <- do.call(dataGeneration,dataGenArgs)
  rUniform <- NA
  if(!is.na(dichotomous) & dichotomous == TRUE){
    rUniform[y==1] <- runif(sum(y==1),min = 0,max=dataGenArgs$prob)
    rUniform[y==0] <- runif(sum(y==0),min = dataGenArgs$prob,max=1)
  }



  # 2 Generate treatment effect if not provided
  if( is.function(effectGeneration) ){
         effect <- do.call(effectGeneration,effectGenArgs)
  } else effect <- effectGeneration


  # 2.1 Transform treatment effect to identity scale
  #     Example, if effect is on the odds ratio scale, transform to probability
  if( toupper(effectScale) == "IDENTITY" ){

    effectIdentity <- effect

  } else if (toupper(effectScale) %in% c("ODDSRATIO","OR")){

    dataOR         <-  dataGenArgs$prob / (1 - dataGenArgs$prob)
    effectIdentity <- (dataOR * effect) / (1 + dataOR * effect) - dataGenArgs$prob

  } else if (toupper(effectScale) %in% c("ODDS")){

    effectIdentity <- effect / (1 + effect)

  } else if (toupper(effectScale) %in% c("LOG")){

    effectIdentity <- exp(effect)

  }




  # 2.2 Set treatment variable and apply effect
  # If length(allocation) == 0, then all patients will have effect shifted by theta
  # Else randomly treat using block allocation with effect who will have effect shifted by theta
  # Ensure enough treatments if generating an odd number of observations
  if(randomize==FALSE){
    trt <- rep(1, length(y))
  } else {
    trt <- rep(0:1, times = length(y) / sum(allocation))
  }

  if(!is.na(dichotomous) & dichotomous==TRUE){
    y[trt==1] <- as.numeric( as.numeric(rUniform[trt==1] <= dataGenArgs$prob + effectIdentity) )
  } else {
    y[trt==1] <- y[trt==1] + effectIdentity
  }


  # 2.3 Record effect (thetas) for output
  # theta0  : Expectation of theta used in generation data (set in SeqSGPV)
  # effectX : Arm X Effect (on scale specified by effectScale; ex. could be odds ratio scale)
  if(randomize==FALSE){
    effect0 <- effect
  } else if(randomize==TRUE){
    effect1 <- effect
  }




  # 3 Add to previous data if any (ie previous data was insufficient)
  if(! is.null(existingData) ) {
    y   <- c(existingData[,"y"], y)
    trt <- c(existingData[,"trt"], trt)
    if(!is.na(dichotomous) & dichotomous == TRUE){
      rUniform <- c(existingData[,"rUniform"], rUniform)
    }

  }


  # 3.5 Create design matrix
  if(randomize==FALSE){
    XD <- as.matrix(trt,ncol=1)
  } else {
    XD <- as.matrix(cbind(1,trt),ncol=2)
  }

  # 4 Obtain (or add to) estimate and 1-alpha/2 monitoring confidence intervals
  #   Wait until at least two observations in each group

  # Common inputs for modelFit function
  modelFitArgs$y       <- y
  modelFitArgs$XD      <- XD



  if(! is.null(existingData) ) {
    eci <- rbind(existingData[,c("est","lo","up")],
                 t(sapply( (length(y)-dataGenArgs[["n"]] + 1):length(y),
                           modelFit, modelFitArgs=modelFitArgs )))
  } else {
    # First assess after 2 observations in each group
    if(randomize==FALSE){
      firstAssess <- 2
    } else {
      firstAssess <- max(which(trt==0)[2], which(trt==1)[2])
    }
    eci <- rbind(matrix(rep(c(NA,-Inf,Inf),firstAssess),byrow = TRUE,nrow=firstAssess),
                 t(sapply((firstAssess+1):length(y), modelFit, modelFitArgs=modelFitArgs )))
    colnames(eci) <- c("est", "lo", "up")
  }




  # 5 Return matrix of data, confidence interval, estimate, errors, and sgpvs
  if(randomize==FALSE){
    cbind(theta0, effect0, n = 1:length(y), y, rUniform, trt, eci)
  } else if (randomize==TRUE){
    cbind(theta0, effect1, n = 1:length(y), y, rUniform, trt, eci)
  }

}
