#' @noRd
transFun.HT04 <- function(x, mod) {
  
  #Empirical distribution function
  ox <- order(x)
  names(x) <- 1:length(x)
  x <- sort(x)
  run <- rle(x)
  p <- cumsum(run$lengths)/length(x)  #divisor
  p <- rep(p, run$lengths)
  p <- p[order(as.character(names(x)))]
  x <- x[order(as.character(names(x)))]
  Femp <- p
  
  #dpareto
  th <-mod$threshold
  sigma <- exp(mod$coefficients[1])
  xi <- mod$coefficients[2]
  Para <- (1 + xi * (x - th)/sigma)^(-1/xi)
  Para <- 1 - mean(x > th) * Para
  
  #Pareto distribution above threshold and emprical distribution below
  res <- ifelse(x <= th, Femp, Para)
  
  #Ensurig results are in the order of x
  res[ox] <- as.numeric(sort(res))
  res
}