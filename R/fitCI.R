# fitCI.R
# J Chipman
#
# Model fits for ordinary least squares and logistic regression

# ols CIs
lmCI <- function(y, trt, look, miLevel){
  f     <- lm(y[1:look] ~ trt[1:look])
  coefs <- summary(f)$coefficients
  est   <- coefs[2,"Estimate"]
  eci   <- c(est, est + c(-1,1) * qt(1-miLevel/2, df = f$df.residual) * coefs[2,"Std. Error"])
  return(eci)
}
class(lmCI) <- "normal"

# logistic regression CIs
lrCI <- function(y, trt, look, miLevel){

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
