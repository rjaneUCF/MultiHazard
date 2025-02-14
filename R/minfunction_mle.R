#' @noRd
minfunction_mle <- function(w, data, a, b, lam_end, k, q_minproj, tol, par_init){
  if(tol < 0){
    stop("Convergence tolerance needs to be positive.")
  }
  polynomials <- bbp(w = w, k = k, a = a, b = b)
  basis <- polynomials$basis
  angles <- polynomials$angles
  min_proj <- sapply(angles, function(i) minproj_lambda(data = data, w = i, q_minproj = q_minproj))
  t <- c()
  len_vec <- c()
  for(i in 1:length(angles)){
    aux <- min_proj[1, ][[i]][min_proj[1, ][[i]] > min_proj[2, ][[i]]]
    len_vec[i] <- length(aux)
    t[(length(t) + 1):(length(t) + len_vec[i])] <- aux - min_proj[2, ][[i]]
  }
  results <- tryCatch(optim(par = par_init, fn = est_beta, basis = basis, t = t, len_vec = len_vec, w = angles, lam_end = lam_end, method = "BFGS", control = list(maxit = 100000)),
                      error = function(e){1})
  if(is.list(results)){optim_output <- results}
  else{
    optim_output <- optim(par = par_init, fn = est_beta, basis = basis, t = t, len_vec = len_vec, w = angles, lam_end = lam_end, control = list(maxit = 100000))
  }
  results <- tryCatch(optim(par = optim_output$par, fn = est_beta, basis = basis, t = t, len_vec = len_vec, w = angles, lam_end = lam_end, control = list(maxit = 100000)),
                      error = function(e){1})
  if(is.list(results)){optim_output2 <- results}
  else{
    optim_output2 <- optim(par = optim_output$par, fn = est_beta, basis = basis, t = t, len_vec = len_vec, w = angles, lam_end = lam_end, control = list(maxit = 100000))
  }
  
  while(abs(optim_output2$val - optim_output$val) >= tol){
    optim_output <- optim_output2
    results <- tryCatch(optim(par = optim_output$par, fn = est_beta, basis = basis, t = t, len_vec = len_vec, w = angles, lam_end = lam_end, control = list(maxit = 100000)),
                        error = function(e){1})
    if(is.list(results)){optim_output2 <- results}
    else{
      optim_output2 <- optim(par = optim_output$par, fn = est_beta, basis = basis, t = t, len_vec = len_vec, w = angles, lam_end = lam_end, control = list(maxit = 100000))
    }
  }
  return(c(lam_end[1], exp(optim_output2$par), lam_end[2]))
}