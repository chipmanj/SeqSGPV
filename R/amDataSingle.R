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
                         miLevel,
                         existingData=NULL){


  # 0.1 if data type is dichotomous, use uniform data with cutoffs for dataGeneration
  #     Note: uniform data makes it easier to take the same data and impose a treatment effect (such as on odds ratio scale)
  # 0.2 set theta0 from data generation arguments
  if(any(grepl("rbinom",deparse(dataGeneration)))){
    dichotomous <- TRUE
    theta0      <- dataGenArgs$probs

  } else if (any(grepl("rnorm",deparse(dataGeneration)))){
    dichotomous <- FALSE
      if(is.element("mean",dataGenArgs)) {
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
  if(!is.na(dichotomous) & dichotomous == TRUE){
    rUniform <- runif(dataGenArgs$n)
    y        <- as.numeric(rUniform <= dataGenArgs$prob)
  } else {
    rUniform <- NA
    y        <- do.call(dataGeneration,dataGenArgs)
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
  # If randomize == FALSE, then all patients will have effect shifted by theta
  # Else randomly treat half with effect who will have effect shifted by theta
  # Use block 2 randomization for sims but not in practice
  # Ensure enough treatments if generating an odd number of observations

  trt <- rep(1,length(y))
  if(randomize) trt[seq(2,length(y),2)] <- 0

  if(!is.na(dichotomous) & dichotomous==TRUE){
    y[trt==1] <- as.numeric( as.numeric(rUniform[trt==1] <= dataGenArgs$prob + effectIdentity) )
  } else {
    y[trt==1] <- y[trt==1] + effectIdentity
  }


  # 2.3 Record effect (thetas) for output
  # theta0  : Expectation of theta used in generation data (set in sgpvAM)
  # effectX : Effect (on scale specified by effectScale; ex. could be odds ratio scale)
  effect1   <- effect

  # 3 Add to previous data if any (ie previous data was insufficient)
  if(! is.null(existingData) ) {
    y   <- c(existingData[,"y"], y)
    trt <- c(existingData[,"trt"], trt)
  }


  # 3.5 Create design matrix
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
                           modelFit, y=y, XD=XD, miLevel = miLevel )))
  } else {
    eci <- rbind(matrix(rep(c(NA,-10^10,10^10),4),byrow = TRUE,nrow=4),
                 t(sapply(5:length(y), modelFit, y=y, XD=XD, miLevel = miLevel )))
    colnames(eci) <- c("est", "lo", "up")
  }




  # 5 Return matrix of data, confidence interval, estimate, errors, and sgpvs
  cbind(theta0, effect1, n = 1:length(y), y, rUniform, trt, eci)

}
