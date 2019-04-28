#' @export
plot.sgpvAMlocationShift <- function( amShifted,        stat,
                                      waitWidth = NULL, alertK = NULL,
                                      xlim,             ylim,
                                      sizeRestrictions ){

  amInputs <- amShifted[[1]]$inputs


  if(missing(xlim)) { xlim <- readline(prompt="xlim: ") }
  if(missing(ylim)) { ylim <- readline(prompt="ylim: ") }


  # Get wait width and alert k parameters
  if(is.null(waitWidth)){
    waitWidth <- amInputs$waitWidths
  }
  if(length(waitWidth)>10) stop("Please select at most 10 wait times to investigate")


  if(is.null(alertK)){
    alertK    <- seq(0, amInputs$maxAlertSteps, by = amInputs$kSteps)
  }
  if(length(alertK)>11) stop("Please select at most 11 required affirmation steps to investigate")


  if(length(waitWidth)>1 & length(alertK) > 1){
    stop("At least one of waitWidth and alertK must be a singular.")
  } else if(length(waitWidth) == 1 & length(alertK) >  1){
    mainGiven <- paste0("Effects within bound periphery wait for CI width = ", waitWidth)
  } else if(length(waitWidth) >  1 & length(alertK) == 1){
    mainGiven <- paste0("Required affirmation steps = ", alertK)
  } else {
    mainGiven <- paste0("Effects within bound periphery wait for CI width = ", waitWidth,
                        "; Required affirmation steps = ", alertK)
  }

  # Set colors
  nVary <- max(c(length(waitWidth),length(alertK)))
  if(nVary >= 3) {
    cols <- brewer.pal(nVary, name="Spectral")
  } else {
    cols <- 1:nVary
  }


  # If alertK not provided, get from names of amShifted
  # Number of treatment effects not limited
  treatEffect <- unname(sapply(names(amShifted), FUN = function(x) as.numeric(strsplit(x,split="_")[[1]][2])))


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
      if(sizeRestrictions == "lag") mainSub <- paste0(amInputs$lagOutcomeN," Lag outcomes") else
        if(sizeRestrictions == "lagMaxN") mainSub <- paste0("Maximum observed sample size before ", amInputs$lagOutcomeN," lag outcomes")
  } else stop("restrictions parameter must be any of NULL, maxN, lag, lagMaxN")




  # Blank plot
  par(mar=c(5,4,5,8)+.1)
  plot(x=0,y=0,xlim=xlim,ylim=ylim,type="n",
       las=1,xlab="",ylab="",
       main=paste0(main,"\n",mainGiven,"\n",mainSub))


  # Add line for point null
  abline(v=amInputs$pointNull,lty=2)
  abline(v=c(amInputs$deltaL2,
             amInputs$deltaL1,
             amInputs$deltaG1,
             amInputs$deltaG2),lwd=2)


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

    legend("topright", inset=c(-.35, .05),legend="Effects within\nbound periphery\nwait for CI width", bty="n",xpd=TRUE)
    legend("topright", inset=c(-.25, .35),legend=waitWidth,col=cols, pch=19, bty="n",xpd=TRUE)

  } else if(length(waitWidth)==1){
    for(k in alertK){
      lines( x=toPlot[toPlot[,"k"]==k,"te"], toPlot[toPlot[,"k"]==k,"y"], col=cols[colIter])
      points(x=toPlot[toPlot[,"k"]==k,"te"], toPlot[toPlot[,"k"]==k,"y"], col=cols[colIter])
      colIter <- colIter + 1
    }

    legend("topright", inset=c(-.35, .05),legend="Required\nAffirmation\nSteps", bty="n",xpd=TRUE)
    legend("topright", inset=c(-.25, .275),legend=alertK,col=cols, pch=19, bty="n",xpd=TRUE)

  }



  # Reset margins
  par(mar=c(5,4,4,2)+.1)

}
