#' Declusters a time series
#'
#' Identify cluster maxima above a threshold, using the runs method of Smith and Weissman (1994).
#'
#' @param Data Numeric vector of the time series.
#' @param u Numeric vector of length one specifying the declustering threshold; as a quantile \code{[0,1]} of \code{Data} vector. Default is \code{0.95}.
#' @param SepCrit Integer; specifying the separation criterion under which events are declustered. Default is \code{3} corresponding to a storm window of three days in the case of daily data.
#' @param mu (average) occurrence frequency of events in \code{Data}. Numeric vector of length one. Default is \code{365.25}, daily data.
#' @return List comprising the \code{Threshold} above which cluster maxima are identified, rate of cluster maxima \code{Rate}, a vector containing the original time series \code{Detrended} and the \code{Declustered} series.
#' @seealso \code{\link{Detrend}}
#' @export
#' @examples
#' Decluster(data=S20_T_MAX_Daily_Completed_Detrend$Detrend)
Decluster<-function(Data,u=0.95,SepCrit=3,mu=365.25){

  Thres<-quantile(na.omit(Data),u)

  if(length(which(is.na(Data)==T))>0){
  z<-which(is.na(Data)==T)
  Data[z]<-min(Data,na.rm=T)-1000
  }

  Events<-Event_Identify(Data=Data,Threshold=Thres,SeparationPeriod=SepCrit)
  Events.Max<-Event_Max(Data=Data,Events=Events)
  Events.Start<-Event_Start(Data=Data,Threshold=Thres,Events=Events,Event.Max=Events.Max)
  print(Events)
  print(Events.Start)

  Threshold<-Thres
  Rate<-length(Events)/(length(Data)/mu)

  #Declustered data as a vector
  data_Detrend_Declustered<-Data
  for(i in 1:length(Events)){
    data_Detrend_Declustered[Events.Start[i]:Events[i]]<-NA
  }
  data_Detrend_Declustered[Events.Max]<-Data[Events.Max]
  data_Detrend_Declustered[z]<-NA
  Data[z]<-NA

  res<-list("Threshold" = Threshold, "Rate" = Rate, "EventsMax" = Events.Max, "Detrended" = Data, "Declustered" = data_Detrend_Declustered)
  return(res)
}

