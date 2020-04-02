#' @export
plot.sgpvAMlocationShift <- function( amShifted,        stat,
                                      waitTime = NULL,  alertK = NULL,
                                      xlim,             ylim,
                                      xlab = "", ylab = "",
                                      pts = FALSE,
                                      addLegend = TRUE, addMain = TRUE,
                                      addRegions = TRUE, addRegionLines = TRUE,
                                      sizeRestrictions , setMargins, ...){

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


  # plot parameters
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
        if(sizeRestrictions == "lagMaxN") mainSub <- paste0("Maximum sample size after ", amInputs$lagOutcomeN," lag outcomes")
  } else stop("If provided, sizeRestrictions parameter must be one of: maxN, lag, lagMaxN")


  if(addMain) {
    figMain <- paste0(main,"\n",mainGiven,"\n",mainSub)
  } else {
    figMain <- NULL
  }


  if(!missing(setMargins)){
    par(mar=setMargins)
  } else {
    margins <- par("mar")
    if (addLegend){
     margins[4] <- margins[4] + 4
    } else if (addMain){
      margins[3] <- margins[3] + 1
    }
    par(mar=margins)
  }

  # Blank plot
  plot(x=0,y=0,xlim=xlim,ylim=ylim,type="n",
       las=1, xlab=xlab, ylab=ylab,
       main=figMain, ... )

  if(addRegions){

    # Colors for figures
    cbPalette <- c("#d0d7e8","#f7f5f3", "#e8dbd0")
    cbAlpha   <- paste0(cbPalette,sep = "40")
    colsZone  <- cbPalette

    ptCols <- c("#ff0000","#0080ff","#bc03ff")
    ptSymb <- c(19,17,15)

    # Shade clinically meaningful zones
    xLimBuff1 <- min(xlim) - abs(diff(xlim))
    xLimBuff2 <- max(xlim) + abs(diff(xlim))
    yLimBuff1 <- min(ylim) - abs(diff(ylim))
    yLimBuff2 <- max(ylim) + abs(diff(ylim))

    dL2 <- amInputs$deltaL2
    dL1 <- amInputs$deltaL1
    dG1 <- amInputs$deltaG1
    dG2 <- amInputs$deltaG2

    # ROME
    polygon(x=c(xLimBuff1, xLimBuff1, dL2, dL2), y=c(yLimBuff1,yLimBuff2,yLimBuff2,yLimBuff1),col=cbAlpha[3], xpd=FALSE, border=NA)
    polygon(x=c(xLimBuff2, xLimBuff2, dG2, dG2), y=c(yLimBuff1,yLimBuff2,yLimBuff2,yLimBuff1),col=cbAlpha[3], xpd=FALSE, border=NA)

    # Grey Zone
    polygon(x=c(dL2,dL2,dL1,dL1), y=c(yLimBuff1,yLimBuff2,yLimBuff2,yLimBuff1),col=cbAlpha[2], xpd=FALSE, border=NA)
    polygon(x=c(dG1,dG1,dG2,dG2), y=c(yLimBuff1,yLimBuff2,yLimBuff2,yLimBuff1),col=cbAlpha[2], xpd=FALSE, border=NA)

    # ROPE / ROWPE
    # ROWPE: For figure plotting purposes only, set dL2 (or DG2) from NA to a point beyond plot region
    if(is.na(dL2) & is.na(dL1)) dL2 <- xLimBuff1
    if(is.na(dG2) & is.na(dG1)) dG2 <- xLimBuff2

    polygon(x=c(dL1,dL1, dG1, dG1), y=c(yLimBuff1,yLimBuff2,yLimBuff2,yLimBuff1),col=colsZone[1], border=NA)


  }


  if(addRegionLines){
    # Add line for point null
    abline(v=amInputs$pointNull,lty=2)
    abline(v=c(amInputs$deltaL2,
               amInputs$deltaL1,
               amInputs$deltaG1,
               amInputs$deltaG2),lwd=2)
  }

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
      lines( x=toPlot[toPlot[,"w"]==w,"te"], toPlot[toPlot[,"w"]==w,"y"], col=cols[colIter], ... )
      if(pts) points(x=toPlot[toPlot[,"w"]==w,"te"], toPlot[toPlot[,"w"]==w,"y"], col=cols[colIter], ... )
      colIter <- colIter + 1
    }

    if(addLegend){
      legend(x = max(xlim) + abs(diff(xlim)) * .35,
           y = max(ylim),
           legend=waitTime,
           xjust = 1,
           col=cols, pch=19, bty="n",xpd=TRUE, title = "Wait Time")
    }

  } else if(length(waitTime)==1){
    for(k in alertK){
      lines( x=toPlot[toPlot[,"k"]==k,"te"], toPlot[toPlot[,"k"]==k,"y"], col=cols[colIter], ... )
      if(pts) points(x=toPlot[toPlot[,"k"]==k,"te"], toPlot[toPlot[,"k"]==k,"y"], col=cols[colIter], ... )
      colIter <- colIter + 1
    }

    if(addLegend){
      legend(x = max(xlim) + abs(diff(xlim)) * .35,
           y = max(ylim),
           legend=alertK,
           xjust = 1,
           col=cols, pch=19, bty="n",xpd=TRUE,title="Required\nAffirmation\nSteps")
    }

  }



  # Reset margins
  par(mar=c(5,4,4,2)+.1)

}
