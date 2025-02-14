#' @noRd
Laplace_cdf = function(x){ #standard Laplace cdf function for HT model
  u = c()
  u[x<0] = exp(x[x<0])/2
  u[x>=0] = 1-exp(-x[x>=0])/2
  return(u)
}