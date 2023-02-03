#' @noRd
inverse_pit_gpd = function(u,Data,Data_Declust,q){ #function for moving marginal components of curve estimates on uniform margins back to original
  thres = quantile(Data,q)
  mod = evm(Data_Delcust,th = thres,show=FALSE)
  x = c()
  x[u>q] = qgpd((u[u>q]-q)/(1-q),mod$threshold, sigma = exp(mod$par[1]), xi = mod$par[2])
  x[u<=q] = quantile(Data,u[u<=q])
  return(nvec)
}
