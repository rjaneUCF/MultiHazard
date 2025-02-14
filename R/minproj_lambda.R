#' @noRd
minproj_lambda <- function(data, w, q_minproj){
  if(q_minproj < 0 | q_minproj > 1){
    stop("Marginal quantile needs to be in [0, 1].")
  }
  t <- pmin(data[, 1]/w, data[, 2]/(1 - w))
  u <- quantile(t, q_minproj)
  lambda <- 1/(mean(t[t > u] - u))
  return(list("minproj" = t, "thresh" = u, "lambdahill" = lambda))
}