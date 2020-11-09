#' @noRd
Event_Start<-function(Data,Threshold,Events,Event.Max){
  Events<-c(Events,1)
  x.exceed.lower.bound<-numeric(length(Events))
  for(i in 1:(length(x.exceed.lower.bound))){
    x.exceed.lower.bound[i]<-ifelse(Events[i]==1,1,max(Events[(which(Events<Events[i]))]))
  }

  x.exceed.lower.bound<-x.exceed.lower.bound[-length(x.exceed.lower.bound)]
  #return(x.exceed.lower.bound)

  x.exceed.start.position<-numeric(length(x.exceed.lower.bound))
  for(i in 1:(length(x.exceed.lower.bound))){
    x.exceed.start.position[i]<-ifelse(Events[i]==1,1,x.exceed.lower.bound[i]+min(which(Data[(x.exceed.lower.bound[i]+1):Events[i]]>=Thres))))
  }
  #x.exceed.start.position<-ifelse(x.exceed.start.position==0,min(Event.Max),x.exceed.start.position)
  return(x.exceed.start.position)
}
