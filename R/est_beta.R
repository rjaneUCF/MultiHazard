#' @noRd
est_beta <- function(par, basis, t, len_vec, w, lam_end){ 
  if(any(lam_end > 1)){
    stop("Value of the ADF at the endpoints not supported.")
  }
  beta <- c(lam_end[1], exp(par), lam_end[2])
  lam <- basis %*% beta
  if(sum(is.infinite(beta)) > 0 | sum(is.na(beta)) > 0){
    return(1e10)
  }else{
    lam2 <- c() 
    for(i in 1:length(len_vec)){
      lam2[(length(lam2) + 1):(length(lam2) + len_vec[i])] <- rep(lam[i], len_vec[i])
    }
    loglike <- sum(log(lam2)) - sum(lam2 * t) 
    return(-loglike)
  }
}