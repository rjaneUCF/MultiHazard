#' @noRd
HeffTawnNegLL <- function(X, Y, par){ 
  alpha <- par[1]
  beta <- par[2]
  sig <- par[3]
  mu <- par[4]
  if(alpha < 0 || alpha > 1 || beta > 1 || beta < 0 || sig <= 0){
    return(1e10)
  }
  negloglik <- -sum(dnorm(Y, alpha*X + mu*((X)^beta), sig*((X)^beta), log = T))
  if(is.finite(negloglik)){
    return(negloglik)
  }
  else{
    return(1e10)
  }
}
