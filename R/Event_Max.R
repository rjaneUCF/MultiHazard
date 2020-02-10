#' @noRd
Event_Max<-function(Data,Events){
  Events<-c(Events,1)
  x.exceed.max.position<-numeric(length(Events)-1)
  for(i in 1:(length(Events)-1)){
    x.exceed.lower.bound<-ifelse(Events[i]==1,1,max(Events[(which(Events<Events[i]))]))
    x.exceed.value<-max(Data[(x.exceed.lower.bound+1):Events[i]])
    x.exceed.max.position[i]<-(x.exceed.lower.bound)+which(Data[(x.exceed.lower.bound+1):Events[i]]==x.exceed.value)
  }
  return(x.exceed.max.position)
}