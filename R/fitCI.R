##' @name fitCI
##' @rdname fitCI
##'
##' @title Functions for obtaining Confidence Intervals under a linear model and under logistic regression
##'
##' @param y outcomes
##' @param trt treatment assignments for each patient
##' @param look current number of patients accrued
##' @param miLevel alpha in traditional (1-alpha) confidence interval used in monitoring intervals
##'
##' @return est, lowwer bound, and upper bound

##' @rdname fitCI
##' @export
lmCI <- function(y, trt, look, miLevel){
  # ols CIs
  f     <- lm(y[1:look] ~ trt[1:look])
  coefs <- summary(f)$coefficients
  est   <- coefs[2,"Estimate"]
  eci   <- c(est, est + c(-1,1) * qt(1-miLevel/2, df = f$df.residual) * coefs[2,"Std. Error"])
  return(eci)
}
class(lmCI) <- "normal"


##' @rdname fitCI
##' @export
lrCI <- function(y, trt, look, miLevel){
  # logistic regression CIs
  f     <- glm(y[1:look] ~ trt[1:look], family = binomial)
  coefs <- summary(f)$coefficients
  est   <- exp( coefs[2,"Estimate"] )
  eci   <- c(est,exp( coefs[2,"Estimate"] + c(-1,1) * qnorm(1-miLevel/2) * coefs[2,"Std. Error"] ))
  # Infinite CI bounds may occur with binomial data and insufficient
  #  data to estimate an effect and error.  Set infinite bounds to 10^10
  eci[is.infinite(eci)] <- 10^10
  return(ci)

}
class(lrCI) <- "binomial"
