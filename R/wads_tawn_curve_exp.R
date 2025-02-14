#' @noRd
wads_tawn_curve_exp = function(data_exp,prob,q){ #function for estimating return curve for data on standard exponential margins using HT model
  #prob is curve survival probability (small)
  #qalphas is quantile level for estimating heffernan tawn parameters
  
  w <- seq(0, 1, by = 0.01)
  
  alphas <- heff_tawn_alphas(data = data_exp, q = rep(q,2))
  
  a <- alphas[1]/(1 + alphas[1])
  b <- 1/(1 + alphas[2])
  interval <- c(a, b)
  indx <- w < a | w > b
  
  lambda_cl <- c()
  if(sum(!indx) < 2){
    lambda_cl <- pmax(w, 1 - w)
  }
  else{
    lambda_cl[indx] <- pmax(w, 1 - w)[indx]
    basis <- bbp(w = w, k = 7, a = a, b = b)$basis
    lam_end <- c(max(a, 1 - a), max(b, 1 - b))
    betacl <- minfunction_mle(w = w, data = data_exp, a = a, b = b, lam_end = lam_end, k = 7, q_minproj = q, tol = 0.0001, par_init = rep(0, 6))
    lambda_cl[!indx] <- basis %*% betacl
    lambda_cl <- properties(w = w, lambda = as.vector(lambda_cl))
  }
  
  nw = length(w)
  
  xp = qexp(1 - prob)
  
  thresh <- sapply(w, function(i) minproj_lambda(data_exp, i, q_minproj = q)$thresh)
  r <- sapply(1:nw, function(i) thresh[i] - log(prob/(1 - q))/lambda_cl[i])
  x <- sapply(1:nw, function(i) r[i] * w[i])
  y <- sapply(1:nw, function(i) r[i] * (1 - w[i]))
  for(i in 1:length(x)){
    if(x[i] > xp){
      x[i] <- xp
    }
    if(x[i] < 0){
      x[i] <- 0
    }
    if(y[i] > xp){
      y[i] <- xp
    }
    if(y[i] < 0){
      y[i] <- 0
    }
  }
  x[1] <- 0
  y[1] <- xp
  x[nw] <- xp
  y[nw] <- 0
  for(i in length(w[w < 0.5]):1){
    if(x[i] > x[i + 1]){
      x[i] <- x[i + 1]
    }
    if(y[i] < y[i + 1]){
      y[i] <- y[i + 1]
    }
  }
  for(i in length(w[w > 0.5]):(length(w))){
    if(x[i] < x[i - 1]){
      x[i] <- x[i - 1]
    }
    if(y[i] > y[i - 1]){
      y[i] <- y[i - 1]
    }
  }
  
  return(cbind(x,y))#returns matrix containing curve estimate
}