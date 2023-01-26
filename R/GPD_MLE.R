#' @noRd
GPD_MLE<-function(Data){
  mod = evm(na.omit(Data), th = quantile(Data, 0))
  xi = mod[[1]][2]
  sigma = exp(mod[[1]][1])
  u = as.numeric(quantile(Data,0))
  res<-list(xi=xi,sigma=sigma,u=u)
  return(res)
}