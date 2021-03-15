#' @noRd
AR2<-function(Par,Data){
  N     = length(Data)
  CDF      = sort(pgpd(q=Data,sigma=as.numeric(Par[2]),xi=as.numeric(Par[1]),u=as.numeric(Par[3])))-0.0001
  Est = N/2 + -2*sum(CDF)-sum((2-(2*(1:N)-1)/N)*log(1-CDF))
  return(Est)
}