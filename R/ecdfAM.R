#'@export
ecdfAM <- function(am, stat, sizeRestrictions,
                        waitTime = NULL, alertK = NULL, treatEffect = NULL,
                        xlim,             ylim,          maxVary     = 10 ,
                        doPlot = TRUE){


  # if(missing(xlim)) { stop("Specify xlim") }
  # if(length(treatEffect)!=1){
  #   stop("Must specify a treatment effect")
  # }

  if((length(alertK)>1 & length(waitTime)>1) | length(c(alertK,waitTime))==0){
    stop("At least one of alertK and waitWidth must be 1")
  }

  # For naming convention with ECDF elements
  if(!is.element(stat,c("n","bias"))){
    stop("Stat must be one of 'n' or 'bias'")
  }
  if(stat=="n")    stat <- "Size"
  if(stat=="bias") bias <- "Bias"



  # Any restrictions
  if (!missing(sizeRestrictions) & is.element(sizeRestrictions,c("maxN","lag","lagMaxN"))){

    # For naming convention with ECDF elements
    if(sizeRestrictions=="lagMaxN") sizeRestrictions <- "LagMaxN"

    stat    <- paste0(stat,sizeRestrictions)
  } else stop("If provided, sizeRestrictions parameter must be one of: maxN, lag, lagMaxN")


  # Select treatment effect to summarize and reduce am object
  if(is.element(el = "sgpvAMlocationShift",class(am))){

    if(length(treatEffect)!=1){
      cat(paste0(paste0("select treatment effect of interest from choices:\n",
                        paste0(unname(sapply(names(am),
                                             FUN = function(x) as.numeric(strsplit(x,split="_")[[1]][2]))),
                               collapse=", "))))
      treatEffect <- readline(prompt="input: ")
    }

    o <- am[[paste0("theta_",treatEffect)]]

  } else {

    o <- am

  }

  amInputs <- o$inputs

  if(!is.function(o$inputs$effectGeneration)){
    mainTE <- paste0("theta = ", o$inputs$effectGeneration," ")
  }


  # If waitWidth not provided, get from inputs
  # Limit number of waitWidths if specified by maxInputs
  # Set colors based on number of alertK explored
  # if(is.null(waitWidth)){
  #
  #   waitWidth <- o$inputs$waitWidths
  #
  #   if(!is.na(maxVary)){
  #     s         <- ceiling(length(waitWidth)/maxVary)
  #     waitWidth <- waitWidth[unique(c(1,seq(s,length(waitWidth), by=s)))]
  #   }
  #
  #   cols      <- brewer.pal(length(waitWidth), name="Spectral")
  #
  #   mainWW <- NULL
  #
  # } else if(length(waitWidth) == 1){
  #
  #   mainWW <- paste0("wait width = ", waitWidth," ")
  #
  # }


  # Get wait width and alert k parameters
  if(is.null(waitTime)){
    waitTime <- amInputs$waitTime
  }
  if(length(waitTime)>10) stop("Please select at most 10 wait times to investigate")

  if(length(waitTime) == 1){

    mainWait <- paste0("wait time = ", waitTime," ")

  } else mainWait <- NULL



  if(is.null(alertK)){
    alertK    <- seq(0, amInputs$maxAlertSteps, by = amInputs$kSteps)
  }
  if(length(alertK)>11) stop("Please select at most 11 required affirmation steps to investigate")


  if(length(alertK) == 1){

    mainK <- paste0("alert k = ", alertK)

  } else mainK <- NULL



  # if(length(waitTime) == 1 & length(alertK) >  1){
  #   mainGiven <- paste0("Wait time = ", waitTime)
  # } else if(length(waitTime) >  1 & length(alertK) == 1){
  #   mainGiven <- paste0("Required affirmation steps = ", alertK)
  # } else {
  #   mainGiven <- paste0("Wait time = ", waitTime,
  #                       "; Required affirmation steps = ", alertK)
  # }




  # If alertK not provided, get from inputs
  # Limit number of alertK if specified by maxInputs
  # Set colors based on number of alertK explored
  # if(is.null(alertK)){
  #
  #   alertK    <- seq(0, o$inputs$maxAlertSteps, by = o$inputs$lookSteps)
  #
  #   if(!is.na(maxVary)){
  #
  #     s      <- ceiling(length(alertK)/maxVary)
  #     alertK <- alertK[unique(c(1,seq(s,length(alertK), by=s)))]
  #
  #   }
  #
  #   cols      <- brewer.pal(length(alertK), name="Spectral")
  #
  #   mainK <- NULL
  #
  # } else if(length(alertK) == 1){
  #
  #   mainK <- paste0("alert k = ", alertK)
  #
  # }




  # Set colors
  nVary <- max(c(length(waitTime),length(alertK)))
  if(nVary >= 3) {
    cols <- brewer.pal(nVary, name="Spectral")
  } else {
    cols <- 1:nVary
  }






  if(doPlot==TRUE){

    # One of the ECDF functions to be used for creating the plot template
    tempECDF <- o[["mcmcEndOfStudy"]][[paste0("width_",waitTime[1])]][["mcmcECDFs"]][[paste0("mcmcEndOfStudyEcdf",stat)]][[1]]

    if(missing(xlim)) {

      xlim <- c(0,ceiling(quantile(tempECDF,1)*1.2))

    }


    # If more than one study design to plot, place key on right of plot
    if(length(c(alertK,waitTime))>1){
      par(mar=c(5,4,4,6)+.1)
    }


    # Blank plot canvas
    plot(tempECDF, col="white", las = 1, xlim=xlim,
         main=paste0("ECDF of ", stat,"\n",paste0(mainTE,mainWait,mainK,collapse = ", ")))
    colIter <- 1

    for(w in waitTime){
      oo <- o[["mcmcEndOfStudy"]][[paste0("width_",waitTime[1])]][["mcmcECDFs"]][[paste0("mcmcEndOfStudyEcdf",stat)]]

      for (k in alertK){

        plot(oo[[paste0("alertK_",k)]], col=cols[colIter], cex=0.5, add=TRUE)

        colIter <- colIter + 1
      }

      colIter <- colIter + 1

    }

    # If more than one study design to plot, place key on right of plot
    if(length(c(waitTime,alertK))>1){
      if(length(waitTime)>1){

        legend(x = max(xlim) + abs(diff(xlim)) * .35,
               y = 1,
               legend=waitTime,
               xjust = 1,
               col=cols, pch=19, bty="n",xpd=TRUE,title="Wait Time")

      } else {
        legend(x = max(xlim) + abs(diff(xlim)) * .35,
               y = 1,
               legend=alertK,
               xjust = 1,
               col=cols, pch=19, bty="n",xpd=TRUE,title="Required\nAffirmation\nSteps")      }
    }



  }


  # If one K and one wait width return ecdf function
  if(length(alertK)==1 & length(waitTime)==1){
    return(oo[[paste0("alertK_",alertK)]])
  }

}



