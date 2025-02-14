#' @noRd
Laplace_inverse = function(u){ #inverse function of standard Laplace cdf for HT model
  x = c()
  x[u<=0.5] = log(2*u[u<=0.5])
  x[u>0.5] = -log(2*(1-u[u>0.5]))
  return(x)
}