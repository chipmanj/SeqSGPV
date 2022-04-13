#' @title Plot the design-based operating characteristics conditioned across a range of effects.
#'
#' @description Condition all but two monitoring frequency parameters (wait, steps, affirm, lag, N) must have a single value.
#'
#' @param am SeqSGPV object
#' @param stat Any of: rejH0, cover, bias, stopInconclusive, stopNotROPE, stopNotROME, stopInconsistent.  Use prefix lag.x when some outcomes were lagged/delayed.
#' @param wait Vector of possible wait times (W) before monitoring.
#' @param steps Vector of the number of observations (S) in between monitoring assessments.
#' @param affirm Vector of the number of observations required for affirming a stopping rule (A) in between raising an alert of stopping.
#' @param lag Vector of the number of delayed outcomes between enrolling a single observation and observing its outcome.
#' @param N Vector of maximum sample size. Can be set as Inf (indefinite).
#' @param xlim Optional limits for x-axis.
#' @param ylim Optional limits for y-axis.
#' @param addRegions Defaults to TRUE to add PRISM boundaries to plot.
#' @param addRegionLines Defaults to TRUE to add PRISM boundary lines to plot.
#' @param addMain Provides default title for figure. Defaults to TRUE.
#' @param addLegend Provides default legend for figure. Defaults to TRUE.
#' @param ablineH Add any set of horizontal reference lines.
#' @param ablineV Add any set of vertical reference lines.
#' @param setMargins Optional figure margins
#' @param ylab ylab
#' @param ... additional plot parameters
#'
#' @export
plot.SeqSGPVeffects <- function( am, stat, wait, steps, affirm, lag, N, xlim, ylim, log,
                                 addRegions=TRUE, addRegionLines=TRUE, addMain = TRUE,
                                 addLegend=TRUE, ablineH=NULL,ablineV=NULL, setMargins, xlab,ylab="", ...){



  amInputs <- am[[1]]$inputs

  # 1 or 2 arm trial
  if(length(amInputs$allocation)==1){
    effectX <- "effect0"
    if(amInputs$effectScale=="identity" & missing(xlab)){
      xlab <- "Theta"
    } else if(missing(xlab)){
      xlab <- "Effect"
    }

  } else {
    effectX <- "effect1"

    if(missing(xlab)){
      xlab <- "Effect"
    }

  }


  # Wait time
  if(missing(wait)){

    w <- sort(unique(amInputs$wait))

  } else w <- wait

  # Monitoring steps
  if(missing(steps)){

    s <- sort(unique(amInputs$steps))

  } else s <- steps

  # Affirmation steps
  if(missing(affirm)){

    a <- sort(unique(amInputs$affirm))

  } else a <- affirm

  # Maximum sample size
  if(missing(N)){

    n <- sort(unique(amInputs$N))

  } else n <- N

  # Number of lag (delayed) outcomes
  if(missing(lag)){

    l <- sort(unique(amInputs$lag))

  } else l <- lag

  # List of monitoring frequency and lag outcomes
  mfl         <- list(wait=w, steps=s, affirm=a, N=n, lag=l)
  mfl_lengths <- sapply(mfl,length)

  # effects
  effects <- paste0(unname(sapply(names(am),FUN = function(x) as.numeric(strsplit(x,split="_")[[1]][2]))))


  # Ensure all, or all but one, parameter(s) held constant
  if( sum(mfl_lengths==1) < 4){

    print(mfl)
    stop("Set all, or all but one, monitoring frequency or lag (delayed) settings to a single value.")

  } else {

    whichVary   <- which(mfl_lengths>1)
    nVariations <- mfl_lengths[whichVary]

    if(length(whichVary)==0){
      nVariations <- 1
    }

    y <- matrix(nrow=length(effects),ncol=nVariations)



    for(te in 1:length(effects)){

      o <- am[[paste0(effectX,"_",effects[te])]][["mcmcOC"]]

      y[te,] <- o[o[,"wait"] %in% mfl[["wait"]] &
                o[,"steps"]  %in% mfl[["steps"]] &
                o[,"affirm"] %in% mfl[["affirm"]] &
                o[,"N"]      %in% mfl[["N"]] &
                o[,"lag"]    %in% mfl[["lag"]], stat]

    }

  }


  # Main header
  if(grepl(pattern = "rejH0",x = stat)){
    main <- "P( reject H0 )"
  } else if (grepl(pattern="cover", x = stat)){
    main <- "Interval Coverage"
  } else if (stat=="n" | stat=="lag.n"){
    main <- "Average Sample Size"
  } else if (grepl(pattern="mse",x = stat)){
    main <- "MSE"
  } else if (grepl(pattern="bias",x = stat)){
    main <- "Bias (standardized by SD)"
  } else if(grepl(pattern="stopInconclusive", x = stat)){
    main <- "P( inconclusive )"
  } else if(grepl(pattern="stopNotROPE",x = stat)){
    main <- "P( not ROPE )"
  } else if(grepl(pattern="stopNotROME", x= stat)){
    main <- "P( not ROME )"
  } else if(grepl(pattern="stopInconsistent", x = stat)){
    main <- "P( Conclusion changed after lagged outcomes )"
  } else {
    main <- ""
  }

  if(addMain) {
    figMain <- paste(main,"\n",
                     paste0(toupper(substr(names(unlist(mfl[which(mfl_lengths==1)])),start=1,stop = 1)),
                            " = ",
                            unlist(mfl[which(mfl_lengths==1)]), collapse=", "))
  } else {
    figMain <- NULL
  }


  # Set colors
  if(nVariations >= 3) {
    cols <- RColorBrewer::brewer.pal(nVariations, name="Set1")
  } else {
    cols <- 1:nVariations
  }


  # Plot boundaries
  if(length(amInputs$allocation) == 1){

    # If baseline effect is scalar, then it is the same for all simulations
    # -- Use last mcmcOC object (o) to obtain scalar
    if(length(amInputs$effectGenArgs)==0){
      theta0 <- o[1,"theta0"]
    }
    effects <- as.numeric(effects) + theta0

  }


  if(missing(xlim)) xlim <- range(as.numeric(effects))
  if(missing(ylim)) ylim <- range(y)

  if(length(unique(ylim))==1){
    ylim[2] <- ylim[1] + 1.25 * ylim[1]
  }


  if(!missing(setMargins)){
    par(mar=setMargins)
  } else {
    margins <- par("mar")
    if (addLegend & nVariations>1){
     margins[4] <- margins[4] + 4
    } else if (addMain){
      margins[3] <- margins[3] + 1
    }
    par(mar=margins)
  }


  # Blank plot
  if(!missing(log)){
    plot(x=1,y=0,xlim=xlim,ylim=ylim,type="n",
         las=1, xlab=xlab, ylab=ylab,
         main=figMain, log=log, ... )
  } else {
    plot(x=0,y=0,xlim=xlim,ylim=ylim,type="n",
         las=1, xlab=xlab, ylab=ylab,
         main=figMain, ... )
  }

  if(addRegions){

    # Colors for figures
    cbPalette <- c("#d0d7e8","#f7f5f3", "#e8dbd0")
    cbAlpha   <- paste0(cbPalette,sep = "40")

    # Shade clinically meaningful zones
    xLimBuff1 <- min(xlim) - abs(diff(xlim))
    xLimBuff2 <- max(xlim) + abs(diff(xlim))
    yLimBuff1 <- min(ylim) - abs(diff(ylim))
    yLimBuff2 <- max(ylim) + abs(diff(ylim))

    # When plotting on the log scale, set lower plotting region bound to just greater than zero
    if(!missing(log)){
      if(log=="x"){
        if(xLimBuff1<0) xLimBuff1 <- 0.001
        if(xLimBuff2<0) xLimBuff2 <- 0.001
      } else if (log=="y"){
        if(yLimBuff1<0) yLimBuff1 <- 0.001
        if(yLimBuff2<0) yLimBuff2 <- 0.001
      }
    }

    dL2 <- amInputs[["PRISM"]]["deltaL2"]
    dL1 <- amInputs[["PRISM"]]["deltaL1"]
    dG1 <- amInputs[["PRISM"]]["deltaG1"]
    dG2 <- amInputs[["PRISM"]]["deltaG2"]

    # ROME
    polygon(x=c(xLimBuff1, xLimBuff1, dL2, dL2), y=c(yLimBuff1,yLimBuff2,yLimBuff2,yLimBuff1),col=cbAlpha[3], xpd=FALSE, border=NA)
    polygon(x=c(xLimBuff2, xLimBuff2, dG2, dG2), y=c(yLimBuff1,yLimBuff2,yLimBuff2,yLimBuff1),col=cbAlpha[3], xpd=FALSE, border=NA)

    # Grey Zone
    polygon(x=c(dL2,dL2,dL1,dL1), y=c(yLimBuff1,yLimBuff2,yLimBuff2,yLimBuff1),col=cbAlpha[2], xpd=FALSE, border=NA)
    polygon(x=c(dG1,dG1,dG2,dG2), y=c(yLimBuff1,yLimBuff2,yLimBuff2,yLimBuff1),col=cbAlpha[2], xpd=FALSE, border=NA)

    # ROPE / ROWPE
    # ROWPE: For figure plotting purposes only, set dL2 (or DG2) from NA to a point beyond plot region
    if(is.na(dL2) & is.na(dL1)) dL1 <- xLimBuff1
    if(is.na(dG2) & is.na(dG1)) dG1 <- xLimBuff2

    polygon(x=c(dL1,dL1, dG1, dG1), y=c(yLimBuff1,yLimBuff2,yLimBuff2,yLimBuff1),col=cbPalette[1], border=NA)


  }


  if(addRegionLines){
    # Add line for point null
    abline(v=amInputs$effectPN,lty=2)
    abline(v=c(dL2, dL1, dG1, dG2),lwd=2)
  }



  # Plot stat across effects
  for (i in 1:ncol(y)){
    lines( x=effects, y[,i], col=cols[i], ... )
  }

  if(nVariations>1){
    legend(x = max(xlim) + abs(diff(xlim)) * .35,
           y = max(ylim),
           legend=mfl[[whichVary]],
           xjust = 1,
           col=cols, pch=19, bty="n",xpd=TRUE, title = names(mfl)[whichVary])
  }

  # add reference lines
  abline(h=ablineH)
  abline(v=ablineV)


  # Reset margins
  par(mar=c(5,4,4,2)+.1)

}
