#' Declusters a time series
#'
#' Identify cluster maxima above a threshold, using the runs method of Smith and Weissman (1994).
#'
#' @param Data Numeric vector of the time series.
#' @param u Numeric vector of length one specifying the declustering threshold; as a quantile \code{[0,1]} of \code{Data} vector. Default is \code{0.95}.
#' @param SepCrit Integer; specifying the separation criterian under which events are declustered. Default is \code{3} corresponding to a storm window of three days in the case of daily data.
#' @param mu (average) Number of events per year. Numeric vector of length one. Default is \code{365.25}, daily data.
#' @return List comprising the \code{Threshold} above which cluster maxima are identifed, average number of declustered excesses per year \code{EventsPerYear}, a vector containing the origional time series \code{Detrended} and the \code{Declustered} series.
#' @seealso \code{\link{Detrend}}
#' @export
#' @examples
#' Decluster(data=S28_T_MAX_Daily_Completed_Detrend$Detrend)
Decluster<-function(Data,u=0.95,SepCrit=3,mu=365.25){

  Events<-Event_Identify(Data=Data,Threshold=u,SeparationPeriod = SepCrit)
  Events.Max<-Event_Max(Data=Data,Events=Events)
  Events.Start<-Event_Start(Data=Data,Threshold=u,Events=Events,Event.Max=Events.Max)

  Threshold<-as.numeric(quantile(Data,u))
  EPY<-length(Events)/(length(Data)/mu)

  #Declustered data as a vector
  data_Detrend_Declustered<-Data
  for(i in 1:length(Events)){
    data_Detrend_Declustered[Events.Start[i]:Events[i]]<-NA
  }
  data_Detrend_Declustered[Events.Max]<-Data[Events.Max]

  res<-list("Threshold" = Threshold, "EventsPerYear" = EPY, "EventsMax" = Events.Max, "Detrended" = Data, "Declustered" = data_Detrend_Declustered)
  return(res)
}

