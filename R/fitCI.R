##' @name fitCI
##' @rdname fitCI
##'
##' @title Functions for obtaining Confidence Intervals under a linear model and under logistic regression
##'
##' @param y outcomes
##' @param X X matrix with last column corresponding to coefficient of interest
##' @param look current number of patients accrued
##' @param miLevel alpha in traditional (1-alpha) confidence interval used in monitoring intervals
##'
##' @return est, lowwer bound, and upper bound

##' @rdname fitCI
##' @export
lmCI <- function(y, X, look, miLevel){
  # ols CIs
  f <- RcppEigen::fastLmPure(X = as.matrix(X[1:look,]), y = y[1:look])
  eci <- c(f$coefficients[2],
           f$coefficients[2] + c(-1,1) * qt(1-miLevel/2, df = f$df.residual) * f$se[2])
  return(eci)
}
class(lmCI) <- "normal"


##' @rdname fitCI
##' @export
lrCI <- function(y, trt, look, miLevel){
  # logistic regression CIs
  f <- fastglm::fastglmPure(x = as.matrix(X[1:look,]), y=y[1:look],family = binomial())
  eci <- exp(c(f$coefficients[2],
               f$coefficients[2] + c(-1,1) * qnorm(1-miLevel/2) * f$se[2]))
  # Infinite upper CI bound may occur with binomial data and insufficient
  #  data to estimate an effect and error.  Set infinite bounds to 10^10
  eci[is.infinite(eci)] <- 10^10
  return(eci)

}
class(lrCI) <- "binomial"
