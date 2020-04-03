#' @export
getStats <- function(am, stat, k, wNs){
  # Function to obtain statistics for plotting

  if(is.element("sgpvAMlocationShift",class(am))){
    return(invisible())
  }

  s <- vector()
  for (w in wNs){
    o <- am[["mcmcEndOfStudy"]][[paste0("width_",w)]]$mcmcEndOfStudyAve
    s <- c(s,o[o[,"alertK"]==k,stat])
  }
  s

}

#' @export
plot.sgpvAM <- function(am, stat, xlim, ylim, addMain = TRUE, sizeRestrictions, setMargins, ...){

  # Wait times and k steps
  amInputs <- am$inputs
  wNs <- amInputs$waitTimes
  ks  <- seq(0,amInputs$maxAlertSteps,by=amInputs$kSteps)

  if(missing(xlim)) xlim <- c(min(wNs),max(wNs))

  if(missing(ylim)) {
    if(stat%in%c("rejPN","stopNotROPE")){
      ylim <- c(0,0.10)
    } else if(stat%in%c("stopNotROME")){
      ylim <- c(0.75,1)
    } else {
      stop("Please enter ylim")
    }
  }


  # Set colors
  nVary <- length(ks)
  if(nVary >= 3) {
    cols <- brewer.pal(nVary, name="Spectral")
  } else {
    cols <- 1:nVary
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
  } else if(stat=="stopInconsistent"){
    main <- "P( Conclusion changed after lagged outcomes )"
  }



  # Any restrictions
  if (missing(sizeRestrictions)) {
    mainLag <- NULL
    mainMax <- NULL
  } else if (is.element(sizeRestrictions,c("maxN","lag","lagMaxN"))){
    stat    <- paste0(sizeRestrictions,".",stat)

    # Title if there is a lag time
    if(grepl("LAG",toupper(stat))){
      mainLag      <- paste0("lag time = ", amInputs$lagOutcomeN)
    } else mainLag <- NULL

    # Title if there is a lag time
    if(grepl("MAX",toupper(stat))){
      mainMax      <- paste0("max n = ", amInputs$maxN)
    } else mainMax <- NULL

  } else stop("restrictions parameter must be any of NULL, maxN, lag, lagMaxN")


  if(addMain) {
    figMain <- paste0(main,"\n",paste0(c(mainLag,mainMax),collapse = ", "))
  } else {
    figMain <- NULL
  }


  if(!missing(setMargins)){
    par(mar=setMargins)
  } else {

    margins <- par("mar")
    margins[4] <- margins[4] + 4

    if (addMain){ margins[3] <- margins[3] + 1 }

    par(mar=margins)
  }



  # Blank plot canvas
  plot(x = wNs, y = getStats(am, stat, k=0, wNs = wNs),
       xlim=xlim, ylim=ylim,
       las = 1, type="l", col = "white",
       ylab="", xlab = "Wait time", main=figMain, ...)

  # for k = 10 to 70 by 10
  for (i in 1:length(ks)){
    rpn <- getStats(am, stat, k=ks[i], wNs = wNs)
    lines(wNs, rpn, col=cols[i])
  }

  legend(x = max(xlim) + abs(diff(xlim)) * .35,
         y = max(ylim),
         legend=ks,
         xjust = 1,
         col=cols, pch=19, bty="n",xpd=TRUE,title="Required\nAffirmation\nSteps")



  # Reset margins
  par(mar=c(5,4,4,2)+.1)

}
