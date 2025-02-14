#' @noRd
PosGumb.Laplace.negloglik <- function(yex, ydep, a, b, m, s, constrain, v, aLow) {
  BigNumber <- 10^40
  WeeNumber <- 10^(-10)
  
  if(a < aLow[1] | s < WeeNumber | a > 1-WeeNumber  | b > 1-WeeNumber) {
    res <- BigNumber
  } else {
    mu <- a * yex + m * yex^b
    sig <- s * yex^b
    
    res <- sum(0.5 * log(2*pi) + log(sig) + 0.5 * ((ydep - mu)/sig)^2)
    
    if (is.infinite(res)){
      if (res < 0){ 
        res <- -BigNumber 
      } else {
        res <- BigNumber
      }
      warning("Infinite value of Q in mexDependence")
    } else if (constrain){
      #v <- v * max(yex)
      zpos <- range(ydep - yex) # q0 & q1
      z <- range((ydep - yex * a) / (yex^b)) # q0 & q1 
      zneg <- range(ydep + yex) # q0 & q1
      
      if (!ConstraintsAreSatisfied(a,b,z,zpos,zneg,v)){
        res <- BigNumber
      }
    } 
  }    
  res
}