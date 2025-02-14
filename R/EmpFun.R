#' @noRd
EmpFun <- function(x, r, mod) {
  
  X<-x
  #Empirical distribution function
  F <-ecdf(x)
  Femp<-F(r)
  
  #Empirical distribution function
  ox <- order(x)
  names(x) <- 1:length(x)
  x <- sort(x)
  run <- rle(x)
  p <- cumsum(run$lengths)/length(x)  #divisor
  p <- rep(p, run$lengths)
  p <- p[order(as.character(names(x)))]
  x <- x[order(as.character(names(x)))]
  #Femp <- p
  
  #dpareto
  th <-mod$threshold
  sigma <- exp(mod$coefficients[1])
  xi <- mod$coefficients[2]
  #Para <- (1 + xi * (r - th)/sigma)^(-1/xi)
  #Para <- 1 - (1-th) * Para
  Para <- 1-(1-F(th)) * pgpd(r,u=th,sigma=exp(mod$coefficients[1]),xi=mod$coefficients[2],lower.tail = FALSE)
  
  #Pareto distribution above threshold and emprical distribution below
  res <- ifelse(Femp <= th, Femp, Para)
  
  #Ensuring results are in the order of x
  res <- as.numeric(res)
  res
}
