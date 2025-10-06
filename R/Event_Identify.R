#' Event_Identify
#' @noRd
Event_Identify<-function(Data,Threshold,SeparationPeriod=3){  #Threshold is a quantile
  x.exceed<-order(Data,decreasing = TRUE)[1:length(which(Data>=Threshold))]
  x.exceed.position<-numeric(length(x.exceed))
  for(i in 1:length(x.exceed)){
    x.exceed.position[i]<-ifelse(max(Data[(x.exceed[i]+1):min((x.exceed[i]+SeparationPeriod),length(Data))])>=Threshold,0,order(Data,decreasing = TRUE)[i])
  }
  x.exceed.position<-x.exceed.position[which(x.exceed.position>0)]
  return(x.exceed.position)
}
