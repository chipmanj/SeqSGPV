#'@export
ecdf.sgpvAM <- function(am,        stat,
                        waitWidth = NULL, alertK = NULL, treatEffect = NULL,
                        xlim,             ylim,          maxVary     = 10 ,
                        doPlot = TRUE){



  if(length(treatEffect)!=1){
    stop("Must specify a treatment effect")
  }

  if((length(alertK)>1 & length(waitWidth)>1) | length(c(alertK,waitWidth))==0){
    stop("At least one of alertK and waitWidth must be 1")
  }

  if(!is.element(stat,c("Size","Bias"))){
    stop("Stat must be one of 'Size' or 'Bias'")
  }


  # Select treatment effect to summarize and reduce am object
  if(is.element(el = "sgpvAMlocationShift",class(am))){

    o <- am[[paste0("theta_",treatEffect)]]

  } else {

    o <- am

  }

  if(!is.function(o$inputs$effectGeneration)){
    mainTE <- paste0("theta = ", o$inputs$effectGeneration," ")
  }


  # If waitWidth not provided, get from inputs
  # Limit number of waitWidths if specified by maxInputs
  # Set colors based on number of alertK explored
  if(is.null(waitWidth)){

    waitWidth <- o$inputs$waitWidths

    if(!is.na(maxVary)){
      s         <- ceiling(length(waitWidth)/maxVary)
      waitWidth <- waitWidth[unique(c(1,seq(s,length(waitWidth), by=s)))]
    }

    cols      <- brewer.pal(length(waitWidth), name="Spectral")

    mainWW <- NULL

  } else if(length(waitWidth) == 1){

    mainWW <- paste0("wait width = ", waitWidth," ")

  }




  # If alertK not provided, get from inputs
  # Limit number of alertK if specified by maxInputs
  # Set colors based on number of alertK explored
  if(is.null(alertK)){

    alertK    <- seq(0, o$inputs$maxAlertSteps, by = o$inputs$lookSteps)

    if(!is.na(maxVary)){

      s      <- ceiling(length(alertK)/maxVary)
      alertK <- alertK[unique(c(1,seq(s,length(alertK), by=s)))]

    }

    cols      <- brewer.pal(length(alertK), name="Spectral")

    mainK <- NULL

  } else if(length(alertK) == 1){

    mainK <- paste0("alert k = ", alertK)

  }


  if(doPlot==TRUE){

    if(missing(xlim)) { stop("Specify xlim") }


    # If more than one study design to plot, place key on right of plot
    if(length(c(alertK,waitWidth))>1){
      par(mar=c(5,4,4,6)+.1)
    }


    # Blank plot canvas
    plot(o[["mcmcEndOfStudy"]][[paste0("width_",w)]][[paste0("mcmcEndOfStudyEcdf",stat)]][[1]],
         col="white", las = 1, xlim=xlim,
         main=paste0("ECDF of ", stat,"\n",paste0(mainTE,mainWW,mainK,collapse = ", ")))
    colIter <- 1

    for(w in waitWidth){
      oo <- o[["mcmcEndOfStudy"]][[paste0("width_",w)]][[paste0("mcmcEndOfStudyEcdf",stat)]]

      for (k in alertK){

        plot(oo[[paste0("alertK_",k)]], col=cols[colIter], cex=0.5, add=TRUE)

      }

      colIter <- colIter + 1

    }

    # If more than one study design to plot, place key on right of plot
    if(length(c(waitWidth,alertK))>1){
      if(length(waitWidth)>1){
        legend("topright", inset=c(-.2, .05),legend="Wait Width", bty="n",xpd=TRUE)
        legend("topright", inset=c(-.2, .15),legend=waitWidth,col=cols, pch=19, bty="n",xpd=TRUE)
      } else {
        legend("topright", inset=c(-.2, .05),legend="Alert K", bty="n",xpd=TRUE)
        legend("topright", inset=c(-.2, .15),legend=alertK,col=cols, pch=19, bty="n",xpd=TRUE)
      }
    }



  }


  # If one K and one wait width return ecdf function
  if(length(alertK)==1 & length(waitWidth)==1){
    return(oo[[paste0("alertK_",alertK)]])
  }

}



