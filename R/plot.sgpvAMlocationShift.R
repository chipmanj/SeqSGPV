#' @export
plot.sgpvAMlocationShift <- function( amShifted,        stat,
                                      waitTime = NULL,  alertK = NULL,
                                      xlim,             ylim,
                                      pts = FALSE,
                                      sizeRestrictions ){

  amInputs <- amShifted[[1]]$inputs


  if(missing(xlim)) { xlim <- readline(prompt="xlim: ") }
  if(missing(ylim)) { ylim <- readline(prompt="ylim: ") }


  # Get wait width and alert k parameters
  if(is.null(waitTime)){
    waitTime <- amInputs$waitTimes
  }
  if(length(waitTime)>10) stop("Please select at most 10 wait times to investigate")


  if(is.null(alertK)){
    alertK    <- seq(0, amInputs$maxAlertSteps, by = amInputs$kSteps)
  }
  if(length(alertK)>11) stop("Please select at most 11 required affirmation steps to investigate")

  # Effects within bound periphery wait for CI width
  mainWaitTime <- "Wait until n = "
  # if(amInputs$waitEmpirical==TRUE){
  #   mainME <- "Wait until margin of error < "
  # } else {
  #   mainME <- "Wait for expected CI width = "
  # }
  if(length(waitTime)>1 & length(alertK) > 1){
    stop("Must fix at least waitTime or alertK to one value.")
  } else if(length(waitTime) == 1 & length(alertK) >  1){
    mainGiven <- paste0(mainWaitTime, waitTime)
  } else if(length(waitTime) >  1 & length(alertK) == 1){
    mainGiven <- paste0("Required affirmation steps = ", alertK)
  } else {
    mainGiven <- paste0(mainWaitTime, waitTime,
                        "; Required affirmation steps = ", alertK)
  }

  # Set colors
  nVary <- max(c(length(waitTime),length(alertK)))
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
    main <- "Interval Coverage"
  } else if (stat=="n"){
    main <- "Average Sample Size"
  } else if (stat=="mse"){
    main <- "MSE"
  } else if (stat=="bias"){
    main <- "Bias"
  } else if(stat=="stopInconclusive"){
    main <- "P( inconclusive )"
  } else if(stat=="stopNotROPE"){
    main <- "P( not ROPE )"
  } else if(stat=="stopNotROME"){
    main <- "P( not ROME )"
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
  par(mar=c(5,4,5,6)+.1)
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
  toPlot           <- matrix(NA,nrow=length(treatEffect)*length(waitTime)*length(alertK),ncol=4)
  colnames(toPlot) <- c("te","w","k","y")
  iter   <- 1

  for(te in treatEffect){

    for(w in waitTime){
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
    for(w in waitTime){
      lines( x=toPlot[toPlot[,"w"]==w,"te"], toPlot[toPlot[,"w"]==w,"y"], col=cols[colIter])
      if(pts) points(x=toPlot[toPlot[,"w"]==w,"te"], toPlot[toPlot[,"w"]==w,"y"], col=cols[colIter])
      colIter <- colIter + 1
    }

    legend(x = max(xlim) + abs(diff(xlim)) * .35,
           y = max(ylim) * abs(diff(xlim)) * 0.9,
           legend=waitTime,
           xjust = 1,
           col=cols, pch=19, bty="n",xpd=TRUE, title = "Wait Time")

  } else if(length(waitTime)==1){
    for(k in alertK){
      lines( x=toPlot[toPlot[,"k"]==k,"te"], toPlot[toPlot[,"k"]==k,"y"], col=cols[colIter])
      if(pts) points(x=toPlot[toPlot[,"k"]==k,"te"], toPlot[toPlot[,"k"]==k,"y"], col=cols[colIter])
      colIter <- colIter + 1
    }

    legend(x = max(xlim) + abs(diff(xlim)) * .35,
           y = max(ylim) * abs(diff(xlim)) * 0.9,
           legend=alertK,
           xjust = 1,
           col=cols, pch=19, bty="n",xpd=TRUE,title="Required\nAffirmation\nSteps")

  }



  # Reset margins
  par(mar=c(5,4,4,2)+.1)

}
