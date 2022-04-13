##' @title ecdfDataReduction.R
##'
##' @description Order and reduce x observations before obtaining ecdf.  Helps save storage space when saving ecdf function to object and when drawing upon many observations.  ecdfDataReduction is called within SeqSGPV.
##'
##' @param x observations
##' @param n desired number of reduced, ordered observations
##'
##' @return ecdf after reducing (if applicable) the number of observations
##'
##' @export
ecdfDataReduction <- function(x,n=200){
  # Data reduction not necessary if number of x observations < n
  if(length(x) < n) n <- length(x)

  # Sort and reduce x to n observations
  xN <- sort(x)[seq(1,length(x),length.out = n)]

  # Return ECDF
  return(ecdf(xN))
}
