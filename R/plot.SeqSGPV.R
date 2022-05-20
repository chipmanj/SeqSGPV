#' @title Plot the design-based operating characteristics conditioned on a single effect.
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
#' @param addMain Provides default title for figure. Defaults to TRUE.
#' @param addLegend Provides default legend for figure. Defaults to TRUE.
#' @param ablineH Add any set of horizontal reference lines.
#' @param ablineV Add any set of vertical reference lines.
#' @param setMargins Optional figure margins
#' @param ylab Label for y-axis
#' @param ... additional plot parameters
#'
#' @export
plot.SeqSGPV <- function( am, stat, wait, steps, affirm, lag, N, xlim, ylim, log,
                          addMain = TRUE, addLegend=TRUE, ablineH=NULL,ablineV=NULL,
                          setMargins, ylab="", ...){

  if(length(am$inputs$effectGenArgs)!=0){

    stop("plot function is currently only for when conditioning on a fixed parameter value")

  }

  # Inputs stored in am object
  if(is.element(el = "SeqSGPVeffects",class(am))){
    amInputs <- am[[1]]$inputs
  } else {
    amInputs <- am$inputs
  }


  # 1 or 2 arm trial
  if(length(amInputs$allocation)==1){
    effectX <- "effect0"
  } else {
    effectX <- "effect1"
  }


  # Select treatment effect to summarize and reduce am object
  if(is.element(el = "SeqSGPVeffects",class(am))){

    if(missing(effect)){
      cat(paste0(paste0("select treatment effect of interest from choices:\n",
                        paste0(unname(sapply(names(am),
                                             FUN = function(x) as.numeric(strsplit(x,split="_")[[1]][2]))),
                               collapse=", "))))
      e <- readline(prompt="input: ")
    } else {
      e <- effect
    }

    # Reduce to specific out object
    o  <- am[[paste0(effectX,"_",e)]][["mcmcOC"]]

  } else if(!is.function(am$inputs$effectGeneration)) {

    e  <- amInputs$effectGeneration
    o  <- am[["mcmcOC"]]

  } else {
    return(NULL)
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

  # # effects
  # effects <- paste0(unname(sapply(names(am),FUN = function(x) as.numeric(strsplit(x,split="_")[[1]][2]))))


  # Ensure all but one or two parameter(s) held constant
  if( sum(mfl_lengths>1) > 2){

    print(mfl)
    stop("Set all but one or two monitoring frequency or lag (delayed) settings to a single value.")

  } else {

    namesVary <- names(mfl_lengths)[which(mfl_lengths>1)]
    whichVary   <- mfl_lengths[namesVary]

    # grouping variable for legend and x-axis variable
    if(length(namesVary)>1){

      if(whichVary[1] == whichVary[2]){

        glab  <- names(whichVary[1])
        xlab  <- names(whichVary[2])

      } else {

        glab <- names(which.min(whichVary))
        xlab <- names(which.max(whichVary))

      }

      # if(glab=="lag (delayed) outcomes") glab <- "lag"

      gVary       <- mfl[[glab]]
      gVariations <- length(gVary)

    } else {
      gVariations <- 1
      xlab        <- namesVary
    }
    xVary       <- mfl[[xlab]]


    # collect summary outcomes
    y <- matrix(nrow=gVariations,ncol=length(xVary))

    for(g in 1:gVariations){

      if(length(namesVary)>1){
        oo <- o[o[,glab]==gVary[g], ]
      } else {
        oo <- o
      }

      y[g,] <- oo[oo[,"wait"]     %in% w &
                    oo[,"steps"]  %in% s &
                    oo[,"affirm"] %in% a &
                    oo[,"N"]      %in% n &
                    oo[,"lag"]    %in% l, stat]

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
    main <- "Bias"
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




    # 1 or 2 arm trial
    if(length(amInputs$allocation)==1){

      # 1 arm trial with effect on identity scale, report theta
      if(amInputs$effectScale=="identity"){
        param      <- "Theta"
        fixedParam <- e + o[1,"theta0"]
      } else {
        param      <- "Effect"
        fixedParam <- e
      }

    } else {
      param <- "Effect"
      fixedParam <- e
    }


    figMain <- paste0(main,"\n",param," = ", fixedParam, ", ",
                     paste0(toupper(substr(names(unlist(mfl[which(mfl_lengths==1)])),start=1,stop = 1)),
                            " = ",
                            unlist(mfl[which(mfl_lengths==1)]), collapse=", "))
  } else {
    figMain <- NULL
  }


  # Set colors
  if(gVariations >= 3) {
    cols <- RColorBrewer::brewer.pal(gVariations, name="Set1")
  } else {
    cols <- 1:gVariations
  }


  # Plot boundaries
  #
  # if(length(amInputs$allocations)==1){
  #
  #   # If baseline effect is scalar, then it is the same for all simulations
  #   # -- Use last mcmcOC object (o) to obtain scalar
  #   if(length(amInputs$effectGenArgs)==0){
  #     theta0 <- o[1,"theta0"]
  #   }
  #   effects <- as.numeric(effects) + theta0
  #
  # }


  if(missing(xlim)) {
    if(any(is.infinite(xVary))){
      xlim <- c(1,length(xVary))
    } else xlim <- range(xVary)
  }
  if(missing(ylim)) ylim <- range(y)

  # Add to ylim if no range in y
  if(length(unique(ylim))==1){
    ylim[2] <- ylim[1] + 1.25 * ylim[1]
  }


  if(!missing(setMargins)){
    par(mar=setMargins)
  } else {
    margins <- par("mar")
    if (addLegend & gVariations>1){
      margins[4] <- margins[4] + 4
    } else if (addMain){
      margins[3] <- margins[3] + 1
    }
    par(mar=margins)
  }


  # Blank plot
  if(!missing(log)){
    plot(x=1,y=0,xlim=xlim,ylim=ylim,type="n",
         las=1, xlab=xlab, ylab=ylab,xaxt="n",
         main=figMain, log=log, ... )
  } else {
    plot(x=0,y=0,xlim=xlim,ylim=ylim,type="n",
         las=1, xlab=xlab, ylab=ylab,xaxt="n",
         main=figMain, ... )
  }

  if(any(is.infinite(xVary))){
    axis(1, at = 1:length(xVary), labels = xVary)
  } else {
    axis(1, at = xVary, labels = xVary)
  }


  # Plot stat across effects
  for (i in 1:nrow(y)){
    if(any(is.infinite(xVary))){
      lines( x=1:length(xVary), y[i,], col=cols[i], ... )
    }
    lines( x=xVary, y[i,], col=cols[i], ... )

  }

  if(gVariations>1){
    legend(x = max(xlim) + abs(diff(xlim)) * .35,
           y = max(ylim),
           legend=gVary,
           xjust = 1,
           col=cols, pch=19, bty="n",xpd=TRUE, title = glab)
  }

  # add reference lines
  abline(h=ablineH)
  abline(v=ablineV)


  # Reset margins
  par(mar=c(5,4,4,2)+.1)

}
