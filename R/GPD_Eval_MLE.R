#' @noRd
GPD_Eval_MLE<-function(Data,RP,Years){
  # parameters
  mod = evm(na.omit(Data), th = quantile(Data, 0))
  xi = mod[[1]][2]
  sigma = exp(mod[[1]][1])
  u = as.numeric(quantile(Data,0))
  # MRLP
  MRLP       = mean(Data - u)
  # Modifiewd shape parameter
  mod_sigma   = sigma - xi*u
  # Number of events per year
  rate         = length(Data)/Years
  # Return period
  PRI        = 1-1/(RP*rate)
  RP_Est        = qgpd(PRI,xi=xi,sigma=sigma,u=u)
  # Output
  res        = list(xi=xi,sigma=sigma,u=u,MRLP=MRLP,mod_sigma=mod_sigma,rate=rate,RP_Est=RP_Est)
  return(res)
}
