#' @export
plot.sgpvAMlocationShift <- function( amShifted,        stat,
                                      waitWidth = NULL, alertK = NULL, treatEffect = NULL,
                                      xlim,             ylim,          maxVary     = 10 ,
                                      sizeRestrictions ){


  if(missing(xlim)) { stop("Specify xlim") }
  if(missing(ylim)) { stop("Specify ylim") }


  # Set parameters to explore
  if(length(waitWidth)!=1 & length(alertK)!= 1 & length(treatEffect) != 1 &
     length(c(waitWidth, alertK))      != 2 &
     length(c(waitWidth, treatEffect)) != 2 &
     length(c(alertK,    treatEffect)) != 2   ){
    stop("Exactly one of waitWidth, alertK, and treatEffect must be a fixed value.")
  }


  # If waitWidth not provided, get from inputs
  # Limit number of waitWidths if specified by maxInputs
  # Set colors based on number of alertK explored
  if(is.null(waitWidth)){
    waitWidth <- amShifted[[1]]$inputs$waitWidths
    if(!is.na(maxVary)){
      s         <- ceiling(length(waitWidth)/maxVary)
      waitWidth <- waitWidth[unique(c(1,seq(s,length(waitWidth), by=s)))]
    }
    cols      <- brewer.pal(length(waitWidth), name="Spectral")
  } else if(length(waitWidth) == 1){
    mainGiven <- paste0(" | wait width = ", waitWidth)
  }

  # If alertK not provided, get from inputs
  # Limit number of alertK if specified by maxInputs
  # Set colors based on number of alertK explored
  if(is.null(alertK)){
    alertK    <- seq(0, amShifted[[1]]$inputs$maxAlertSteps, by = amShifted[[1]]$inputs$kSteps)
    if(!is.na(maxVary)){
      s      <- round(length(alertK)/maxVary)
      alertK <- alertK[unique(c(1,seq(s,length(alertK), by=s)))]
    }
    cols      <- brewer.pal(length(alertK), name="Spectral")
  } else if(length(alertK) == 1){
    mainGiven <- paste0(" | alert k = ", alertK)
  }

  # If alertK not provided, get from names of amShifted
  # Number of treatment effects not limited
  if(is.null(treatEffect)){
    treatEffect <- unname(sapply(names(amShifted), FUN = function(x) as.numeric(strsplit(x,split="_")[[1]][2])))
  } else if(length(treatEffect) == 1){
    mainGiven <- paste0(" | treatment effect = ", treatEffect)
  }


  # plot parameters (if not specified)
  if(stat=="rejPN"){
    main <- "P( reject point null )"
  } else if (stat=="cover"){
    main <- "CI Coverage"
  } else if (stat=="n"){
    main <- "Average Sample Size"
  } else if (stat=="mse"){
    main <- "MSE"
  } else if (stat=="bias"){
    main <- "Bias"
  } else if(stat=="stopInconclusive"){
    main <- "P( inconclusive )"
  } else if(stat=="stopNotTrivial"){
    main <- "P( not trivial effect )"
  } else if(stat=="stopNotImpactful"){
    main <- "P( not impactful effect )"
  }


  # Any restrictions
  if (missing(sizeRestrictions)) {
    mainSub <- "Unrestricted Sample Size"
  } else if (is.element(sizeRestrictions,c("maxN","lag","lagMaxN"))){
    stat    <- paste0(sizeRestrictions,".",stat)
    if(sizeRestrictions=="maxN")    mainSub <- "Maximum sample size specified" else
      if(sizeRestrictions == "lag") mainSub <- "Lag outcome specified" else
        if(sizeRestrictions == "lagMaxN") mainSub <- "Maximum observed sample size before lag outcomes"
  } else stop("restrictions parameter must be any of NULL, maxN, lag, lagMaxN")




  # Blank plot
  par(mar=c(5,4,4,6)+.1)
  plot(x=0,y=0,xlim=xlim,ylim=ylim,type="n",
       las=1,xlab="",ylab="",
       main=paste0(main,mainGiven,"\n",mainSub))


  # Add line for point null
  abline(v=amShifted[[1]]$inputs$pointNull,lty=2)
  abline(v=c(amShifted[[1]]$inputs$deltaL2,
             amShifted[[1]]$inputs$deltaL1,
             amShifted[[1]]$inputs$deltaG1,
             amShifted[[1]]$inputs$deltaG2),lwd=2)


  # Collect average effects across specified parameters
  toPlot           <- matrix(NA,nrow=length(treatEffect)*length(waitWidth)*length(alertK),ncol=4)
  colnames(toPlot) <- c("te","w","k","y")
  iter   <- 1

  for(te in treatEffect){

    for(w in waitWidth){
      o     <- amShifted[[paste0("theta_",te)]][["mcmcEndOfStudy"]][[paste0("width_",w)]][["mcmcEndOfStudyAve"]]

      for (k in alertK){

        toPlot[iter,] <- c(te,w,k,o[o[,"alertK"]==k,stat])
        iter <- iter + 1

      }

    }

  }
  toPlot <- toPlot[order(toPlot[,"te"], toPlot[,"w"], toPlot[,"k"]),]

  # Plot stats across varying parameter
  colIter <- 1
  if(length(alertK)==1){
    for(w in waitWidth){
      lines( x=toPlot[toPlot[,"w"]==w,"te"], toPlot[toPlot[,"w"]==w,"y"], col=cols[colIter])
      points(x=toPlot[toPlot[,"w"]==w,"te"], toPlot[toPlot[,"w"]==w,"y"], col=cols[colIter])
      colIter <- colIter + 1
    }

    legend("topright", inset=c(-.2, .05),legend="Wait Width", bty="n",xpd=TRUE)
    legend("topright", inset=c(-.2, .15),legend=waitWidth,col=cols, pch=19, bty="n",xpd=TRUE)

  } else if(length(waitWidth)==1){
    for(k in alertK){
      lines( x=toPlot[toPlot[,"k"]==k,"te"], toPlot[toPlot[,"k"]==k,"y"], col=cols[colIter])
      points(x=toPlot[toPlot[,"k"]==k,"te"], toPlot[toPlot[,"k"]==k,"y"], col=cols[colIter])
      colIter <- colIter + 1
    }

    legend("topright", inset=c(-.2, .05),legend="Alert K", bty="n",xpd=TRUE)
    legend("topright", inset=c(-.2, .15),legend=alertK,col=cols, pch=19, bty="n",xpd=TRUE)

  }



  # Reset margins
  par(mar=c(5,4,4,2)+.1)

}
