#' @noRd
PosGumb.Laplace.negProfileLogLik <- function(yex, ydep, a, b, constrain, v, aLow) { 
  Z <- (ydep - yex * a) / (yex^b)
  
  m <- mean(Z)
  s <- sd(Z)
  
  res <- PosGumb.Laplace.negloglik(yex,ydep,a,b,m=m,s=s,constrain,v,aLow=aLow)
  res <- list(profLik=res,m=m, s=s)
  res
}