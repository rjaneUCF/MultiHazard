#' @noRd
heff_tawn_alphas <- function(data, q){
  if(any(q < 0) | any(q > 1)){
    stop("Marginal quantile needs to be in [0, 1].")
  }
  u <- sapply(1:dim(data)[2], function(i) quantile(data[, i], probs = q[i]))
  excdata <- sapply(1:dim(data)[2], function(i) data[data[, i] > u[i], ], simplify = F)
  par <- rep(1/2, 4)
  Yopt <- optim(fn = HeffTawnNegLL, X = excdata[[2]][, 2], Y = excdata[[2]][, 1], par = par, control = list(maxit = 100000)) 
  Ypar <- Yopt$par
  Xopt <- optim(fn = HeffTawnNegLL, X = excdata[[1]][, 1], Y = excdata[[1]][, 2], par = par, control = list(maxit = 100000))
  Xpar <- Xopt$par
  return(c(Ypar[1], Xpar[1]))
}
