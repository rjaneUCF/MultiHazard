#' @noRd
bbp <- function(w, k, a, b){
  if(a >= b){
    stop("[0, a] need to be disjoint from [b, 1].")
  } 
  l <- sum(w >= a & w <= b)
  v <- seq(a, b, length.out = l)
  vnew <- (v - a)/(b - a)
  basis <- array(0, dim = c(l, k+1))
  for(i in 0:k){
    basis[, i + 1] <- (choose(k, i)) * vnew^i * (1 - vnew)^(k - i)
  }
  return(list("basis" = basis, "angles" = v)) 
}