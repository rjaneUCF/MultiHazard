#' @noRd
GPD_MLE_Boot<-function(Data,RP,Years){
  mod        = GPD_MLE(Data)
  xi         = mod$K
  sigma      = mod$SIG
  u          = mod$u
  rate       = length(Data)/Years
  PRI        = 1-1/(RP*rate)
  RP_Est     = qgpd(PRI,xi=K,sigma=SIG,u=u)
  res        = list(xi=xi,sigma=sigma,u=u,rate=rate,RP_Est=RP_Est)
  return(res)
}
