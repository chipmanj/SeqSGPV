#' S@title summary.SeqSGPV
#'
#' @description Summary of the design-based operating characteristics conditioned on a single effect and a single combination of monitoring frequency parameters (wait, steps, affirm, lag, N).
#'
#' @param am SeqSGPV object.
#' @param wait Vector of possible wait times (W) before monitoring.
#' @param steps Vector of the number of observations (S) in between monitoring assessments.
#' @param affirm Vector of the number of observations required for affirming a stopping rule (A) in between raising an alert of stopping.
#' @param lag Vector of the number of delayed outcomes between enrolling a single observation and observing its outcome.
#' @param N Vector of maximum sample size. Can be set as Inf (indefinite).
#' @param effect Effect of interest.
#' @param rd Number of rounding decimals.
#'
#' @export
summary.SeqSGPV <- function(am, wait, steps, affirm, N, lag, effect, rd = 4){

  # Inputs stored in am object
  if(is.element(el = "SeqSGPVeffects",class(am))){
    amInputs <- am[[1]]$inputs
  } else {
    amInputs <- am$inputs
  }

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
    o   <- am[["mcmcOC"]]
  } else {
    return(NULL)
  }


  # Select and reduce to wait time to summarize
  ws <- sort(unique(amInputs$wait))

  if(missing(wait)){

    if(length(ws)>1){

      cat(paste0("select wait time until applying monitoring rules: ", paste0(ws,collapse=", ")))
      w <- readline(prompt="input: ")

    } else w <- ws

  } else if (!is.element(wait,ws)) {

    cat(paste0("select wait time until applying monitoring rules: ", paste0(ws,collapse=", ")))
    w <- readline(prompt="input: ")

  } else w <- wait

  # o <- o[o[,"wait"]==w,]


  # Select and reduce to wait time to summarize
  ss <- sort(unique(amInputs$steps))

  if(missing(steps)){

    if(length(ss)>1){

      cat(paste0("select monitoring steps (number of observations between SGPV assessments): ", paste0(ss,collapse=", ")))
      s <- readline(prompt="input: ")

    } else s <- ss

  } else if (!is.element(steps,ss)) {

    cat(paste0("select monitoring steps (number of observations between SGPV assessments): ", paste0(ss,collapse=", ")))
    s <- readline(prompt="input: ")

  } else {

    s <- steps

  }

  # o <- o[o[,"steps"]==s,]



  # Select and reduce to number of affirmation steps required for stopping
  as <- sort(unique(amInputs$affirm))

  if(missing(affirm)){

    if(length(as)>1){

      cat(paste0("select required affirmation steps: ", paste0(as,collapse=", ")))
      a <- readline(prompt="input: ")

    } else a <- as

  } else if (!is.element(affirm,as)) {

    cat(paste0("select required affirmation steps: ", paste0(as,collapse=", ")))
    a <- readline(prompt="input: ")

  } else {

    a <- affirm

  }

  # o <- o[o[,"affirm"]==a,]


  # Select and reduce to number of lag (delayed) outcomes
  ls <- sort(unique(amInputs$lag))

  if(missing(lag)){

    if(length(ls)>1){

      cat(paste0("select number of lag (delayed) outcomes: ", paste0(ls,collapse=", ")))
      l <- readline(prompt="input: ")

    } else l <- ls

  } else if (!is.element(lag,ls)) {

    cat(paste0("select number of lag (delayed) outcomes: ", paste0(ls,collapse=", ")))
    l <- readline(prompt="input: ")

  } else l <- lag

  # o <- o[o[,"lag"]==l,]


  # Select and reduce to maximum sample size
  ns <- sort(unique(amInputs$N))

  if(missing(N)){

    if(length(ns)>1){

      cat(paste0("select maximum sample size: ", paste0(ns,collapse=", ")))
      n <- readline(prompt="input: ")

    } else n <- ns

  } else if (!is.element(N,ns)) {

    cat(paste0("select maximum sample size: ", paste0(ns,collapse=", ")))
    n <- readline(prompt="input: ")

  } else n <- N



  # List of monitoring frequency and lag outcomes
  mfl         <- list(wait=w, steps=s, affirm=a, N=n, "lag (delayed) outcomes"=l)
  mfl_lengths <- sapply(mfl,length)

  # Ensure all but one parameter held constant
  if( sum(mfl_lengths==1) != 5){

    print(mfl)
    stop("Set all all monitoring frequency and lag (delayed) settings to a single value.")

  }


  o <-o[o[,"wait"]   %in% mfl[["wait"]] &
        o[,"steps"]  %in% mfl[["steps"]] &
        o[,"affirm"] %in% mfl[["affirm"]] &
        o[,"N"]      %in% mfl[["N"]] &
        o[,"lag"]    %in% mfl[["lag (delayed) outcomes"]],]


  # Differences if single vs mx arm
  if(length(amInputs$allocation) == 1){

    # Plot boundaries
    # If baseline effect is scalar, then it is the same for all simulations
    # -- Use last mcmcOC object (o) to obtain scalar
    if(length(amInputs$effectGenArgs)==0){
      theta0 <- o["theta0"]
    }

    # 1 arm trial with effect on identity scale, report theta
    if(amInputs$effectScale=="identity"){
      param      <- "theta"
      fixedParam <- e + theta0
    } else {
      param      <- "effect"
      fixedParam <- e
    }

  } else {

    # 2 arm trial, report effect
    fixedParam <- e
    param <- "effect"

  }


  PRISM <- amInputs$PRISM

  if(is.na(PRISM$deltaL2) & is.na(PRISM$deltaL1)){
    H0label <- paste0(param," is less than or equal to ",amInputs$effectPN)
  } else if(is.na(PRISM$deltaL2) & is.na(PRISM$deltaL1)){
    H0label <- paste0(param," is greater than or equal to ",amInputs$effectPN)
  } else {
    H0label <- paste0(param," equals ",amInputs$effectPN)
  }

  cat(paste0("\nGiven: ",param," = ", fixedParam, ", W = ", w, " S = ", s ,", A = ",a, " and N = ", n, ", with ", l, " lag (delayed) outcomes"))
  cat(paste0("\nH0   : ",H0label))



  if(l > 0){

    cat(paste0("\n  Average Sample Size                                           = ", round(o["lag.n"],rd)))
    cat(paste0("\n  P( reject H0 )                                                = ", round(o["lag.rejH0"],rd)))
    cat(paste0("\n  P( conclude not ROPE effect )                                 = ", round(o["lag.stopNotROPE"],rd)))
    cat(paste0("\n  P( conclude not ROME effect )                                 = ", round(o["lag.stopNotROME"],rd)))
    cat(paste0("\n  P( conclude PRISM inconclusive )                              = ", round(o["lag.stopInconclusive"],rd)))
    cat(paste0("\n  P( conclusion on ROPE or ROME changed with lagged outcomes )  = ", round(o["lag.stopInconsistent"],rd)))
    cat(paste0("\n  P( newly reject H0 with lagged outcomes )                     = ", round(o["lag.stopRejH0_NY"],rd)))
    cat(paste0("\n  P( newly fail to reject H0 with lagged outcomes )             = ", round(o["lag.stopRejH0_YN"],rd)))
    cat(paste0("\n  Coverage                                                      = ", round(o["lag.cover"],rd)))
    cat(paste0("\n  Bias                                                          = ", round(o["lag.bias"],rd)))

  } else {

    cat(paste0("\n  Average sample size              = ", round(o["n"],rd)))
    cat(paste0("\n  P( reject H0 )                   = ", round(o["rejH0"],rd)))
    cat(paste0("\n  P( conclude not ROPE effect )    = ", round(o["stopNotROPE"],rd)))
    cat(paste0("\n  P( conclude not ROME effect )    = ", round(o["stopNotROME"],rd)))
    cat(paste0("\n  P( conclude PRISM inconclusive ) = ", round(o["stopInconclusive"],rd)))
    cat(paste0("\n  Coverage                         = ", round(o["cover"],rd)))
    cat(paste0("\n  Bias                             = ", round(o["bias"],rd)))

  }

  cat("\n\n")
}

