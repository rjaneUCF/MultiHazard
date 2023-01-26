#' @noRd
local_minimum<-function(Data){
  max.element<-rep(NA,length(Data))
  for(j in 4:(length(Data)-4)){
    max.element[j]<-ifelse(Data[j-3]>Data[j] & Data[j-2]>Data[j] & Data[j-1]>Data[j] & Data[j+1]>Data[j] & Data[j+2]>Data[j] & Data[j+3]>Data[j],j,NA)
  }
  return(max.element)
}