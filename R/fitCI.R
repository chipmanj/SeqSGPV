#' @title lmCI
#'
#' @param look Current number of observations accrued
#' @param modelFitArgs Arguments specific to the interval function.  lmCI and lrCI require y (outcomes), XD (design matrix), and miLevel (monitoring interval level)  which are set in SeqSGPV. binomCI requires x (successes), n (observations), and any other inputs used in binom::binconf.
#'
#' @return est, lower bound, and upper bound
#' @export
lmCI <- function(look, modelFitArgs){

  # model fit specification
  y       <- modelFitArgs$y
  XD      <- modelFitArgs$XD
  miLevel <- modelFitArgs$miLevel

  f <- RcppEigen::fastLmPure(X = as.matrix(XD[1:look,]), y = y[1:look])

  betaCol <- ncol(XD)
  eci     <- c(f$coefficients[betaCol],
               f$coefficients[betaCol] + c(-1,1) * qt(1-(1-miLevel)/2, df = f$df.residual) * f$se[betaCol])

  return(eci)

}

#' @title lrCI
#'
#' @param look Current number of observations accrued
#' @param modelFitArgs Arguments specific to the interval function.  lmCI and lrCI require y (outcomes), XD (design matrix), and miLevel (monitoring interval level)  which are set in SeqSGPV. binomCI requires x (successes), n (observations), and any other inputs used in binom::binconf.
#'
#' @return est, lower bound, and upper bound
#' @export
lrCI <- function(look, modelFitArgs){

  # model fit specification
  y       <- modelFitArgs$y
  XD      <- modelFitArgs$XD
  miLevel <- modelFitArgs$miLevel

  # logistic regression CIs
  f       <- fastglm::fastglmPure(x = as.matrix(XD[1:look,]), y=y[1:look],family = binomial())
  betaCol <- ncol(XD)
  eci     <- exp(c(f$coefficients[betaCol],
                   f$coefficients[betaCol] + c(-1,1) * qnorm(1-(1-miLevel)/2) * f$se[betaCol]))
  # Infinite upper CI bound may occur with binomial data and insufficient
  #  data to estimate an effect and error.  Set infinite bounds to 10^10
  eci[is.infinite(eci)] <- 10^10
  return(eci)

}


#' @title binomCI
#'
#' @param look Current number of observations accrued
#' @param modelFitArgs Arguments specific to the interval function.  lmCI and lrCI require y (outcomes), XD (design matrix), and miLevel (monitoring interval level)  which are set in SeqSGPV. binomCI requires x (successes), n (observations), and any other inputs used in binom::binconf.
#'
#' @return est, lower bound, and upper bound
#' @export
binomCI <- function(look, modelFitArgs){

  # model fit specification
  modelFitArgs$x <- sum(modelFitArgs$y[1:look])
  modelFitArgs$n <- look

  # obtain estimates
  ests <- do.call(binom::binom.confint, modelFitArgs)

  # Return mean and CI
  eci <- unlist(c(ests["mean"],ests["lower"],ests["upper"]))
  return(eci)

}
