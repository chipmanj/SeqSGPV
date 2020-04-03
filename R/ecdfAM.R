#'@export
ecdfAM <- function(am, stat, sizeRestrictions,
                        waitTime = NULL, alertK = NULL, treatEffect = NULL,
                        xlim,             ylim,         doPlot      = TRUE, ...){


  # 1. Select treatment effect to summarize and reduce am object if needed
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


  # 2. Inputs from main object
  amInputs <- o$inputs




  # Get wait width and alert k parameters
  if(is.null(waitTime)){
    waitTime <- amInputs$waitTimes
  }
  if(length(waitTime)>10) stop("Please select at most 10 wait times to investigate")


  if(is.null(alertK)){
    alertK    <- seq(0, amInputs$maxAlertSteps, by = amInputs$kSteps)
  }
  if(length(alertK)>11) stop("Please select at most 11 required affirmation steps to investigate")


  if(length(waitTime)>1 & length(alertK) > 1){
    stop("Must fix at least waitTime or alertK to one value.")
  }

  # 3. Translate stat to naming convention used with ECDF elements
  if(!is.element(stat,c("n","bias"))){
    stop("Stat must be one of 'n' or 'bias'")
  }
  if(stat=="n")    stat <- "Size"
  if(stat=="bias") stat <- "Bias"

  # Any maximum sample size restrictions and/or lag times
  if (!missing(sizeRestrictions)){
    if(is.element(sizeRestrictions,c("maxN","lag","lagMaxN"))){
      stat    <- paste0(stat,toupper(substr(sizeRestrictions,1,1)),substring(sizeRestrictions,2))
    } else stop("If provided, sizeRestrictions parameter must be one of: maxN, lag, lagMaxN")
  }



  # 4.1 Get the wait times and alert ks used in figure
  #     Also get the figure titles associated with wiat time and alert k

  # Get wait width and title for wait width
  if(is.null(waitTime)){
    waitTime <- amInputs$waitTime
  }
  if(length(waitTime)>10) stop("Please select at most 10 wait times to investigate")

  if(length(waitTime) == 1){

    mainWait <- paste0("wait time = ", waitTime)

  } else mainWait <- NULL


  # Get alert K and title for alert K
  if(is.null(alertK)){
    alertK    <- seq(0, amInputs$maxAlertSteps, by = amInputs$kSteps)
  }
  if(length(alertK)>11) stop("Please select at most 11 required affirmation steps to investigate")

  if(length(alertK) == 1){

    mainK <- paste0("alert k = ", alertK)

  } else mainK <- NULL



  # 4.2 Title for treatment effect
  if(!is.function(amInputs$effectGeneration)){
    mainTE <- paste0("theta = ", amInputs$effectGeneration)
  }

  # 4.3 Title if there is a lag time
  if(grepl("LAG",toupper(stat))){
    mainLag      <- paste0("lag time = ", amInputs$lagOutcomeN)
  } else mainLag <- NULL

  # 4.3 Title if there is a lag time
  if(grepl("MAX",toupper(stat))){
    mainMax      <- paste0("max n = ", amInputs$maxN)
  } else mainMax <- NULL


  # 5. Colors for figure
  nVary <- max(c(length(waitTime),length(alertK)))
  if(nVary >= 3) {
    cols <- brewer.pal(nVary, name="Spectral")
  } else {
    cols <- 1:nVary
  }


  # 6. Plot
  if(doPlot==TRUE){

    # One of the ECDF functions to be used for creating the plot template
    tempECDF <- o[["mcmcEndOfStudy"]][[paste0("width_",waitTime[1])]][["mcmcECDFs"]][[paste0("mcmcEndOfStudyEcdf",stat)]][[1]]


    # Xlim and xlab
    if(grepl("Size",stat)){
      if(missing(xlim)) xlim <- c(0,ceiling(quantile(tempECDF,1)*1.2))
      xlab <- "sample size"
    } else if (grepl("Bias",stat)){
      if(missing(xlim)) xlim <- c(floor(quantile(tempECDF,0)*1.2),ceiling(quantile(tempECDF,1)*1.2))
      xlab <- "estimate - theta"
    }


    # If more than one study design to plot, place key on right of plot
    if(length(c(alertK,waitTime))>1){
      par(mar=c(5,4,4,6)+.1)
    }


    # Blank plot canvas
    plot(tempECDF, col="white", las = 1, xlim=xlim, xlab=xlab,
         main=paste0("ECDF of ", xlab,"\n",paste0(c(mainTE,mainWait,mainK,mainLag,mainMax),collapse = ", ")),
         ...)
    colIter <- 1


    if(length(waitTime)>1){

      for(w in waitTime){
        oo <- o[["mcmcEndOfStudy"]][[paste0("width_",w)]][["mcmcECDFs"]][[paste0("mcmcEndOfStudyEcdf",stat)]]
        plot(oo[[paste0("alertK_",alertK)]], col=cols[colIter], cex=0.5, add=TRUE)
        colIter <- colIter + 1

      }

    } else if(length(alertK)>1){

      oo <- o[["mcmcEndOfStudy"]][[paste0("width_",waitTime)]][["mcmcECDFs"]][[paste0("mcmcEndOfStudyEcdf",stat)]]

      for(k in alertK){

        plot(oo[[paste0("alertK_",k)]], col=cols[colIter], cex=0.5, add=TRUE)
        colIter <- colIter + 1

      }

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
    ECDF <- o[["mcmcEndOfStudy"]][[paste0("width_",waitTime)]][["mcmcECDFs"]][[paste0("mcmcEndOfStudyEcdf",stat)]][[paste0("alertK_",alertK)]]
    return(ECDF)
  }


  # Reset margins
  par(mar=c(5,4,4,2)+.1)

}



