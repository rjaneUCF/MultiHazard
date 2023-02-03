#' Ocean-side Water Level Intensity
#'
#' Calculates the "intensity" of extreme water levels, as defined in Wahl et al. (2011).
#'
#' @param Data A data frame with co-occurring rainfall and O-sWL time series in two columns labeled \code{"Rainfall"} and \code{"OsWL"}, respectively.
#' @param Cluster_Max Numeric vector containing indexes of peaks in the O-sWL column of \code{Data}. If analyzing a sample conditioned on O-sWL derived using \code{Con_Sample_2D()} set equal to the \code{$xcon} output.
#' @param Base_Line Numeric vector of length one, specifying water level about which to calculate the intensity. Default is the mean O-sWL.
#' @param Rainfall_Interval Numeric vector of length one, specifying length of time before and after a peak over which to sum rainfall totals. Total window width is 2*\code{Rainfall_Interval}+1. Default is \code{24}.
#' @return A data frame with the following columns:  \itemize{
#' \item \code{Pre.High} Index of the OsWL column of \code{Data} containing the preceding high water level.
#' \item \code{Fol.High} Index of the OsWL column of \code{Data} containing the following high water level.
#' \item \code{Pre.Low} Index of the OsWL column of \code{Data} containing the preceding  low water level.
#' \item \code{Fol.Low} Index of the OsWL column of \code{Data} containing the following low water level.
#' \item \code{Intensity} Intensity of the O-sWL.
#' \item \code{V} Total rainfall volume within \code{Rainfall_Interval} before and after the peak.
#' }
#' @seealso \code{\link{Decluster}} \code{\link{WL_Curve}}
#' @export
#' @examples
#' #Decluster O-sWL series at S-13 using a runs method
#' S13.OsWL.Declust = Decluster(Data=S13.Detrend.df$OsWL,
#'                             SepCrit=24*7, u=0.99667)
#' #Calculate O-sWL of the identified cluster maximum
#' intensity = OsWL_Intensity(Data=S13.Detrend.df,Cluster_Max=S13.OsWL.Declust$EventsMax)
#' #Plot O-sWL series identifying cluster maximum (in red) and print "intensity" above each maximum
#' plot(as.Date(S13.Detrend.df$Date_Time),
#'      S13.Detrend.df$OsWL)
#' points(as.Date(S13.Detrend.df$Date_Time[S13.OsWL.Declust$EventsMax]),
#'        S13.Detrend.df$OsWL[S13.OsWL.Declust$EventsMax],pch=16,col=2)
#' text(as.Date(S13.Detrend.df$Date_Time[S13.OsWL.Declust$EventsMax]),
#'      S13.Detrend.df$OsWL[S13.OsWL.Declust$EventsMax]+0.2,
#'      round(intensity$Intensity,0),cex=0.5)
OsWL_Intensity<-function(Data,Cluster_Max,Base_Line=mean(Data$OsWL,na.rm=T),Rainfall_Interval=24){

 #Find local maximum in the O-sWL series
 x.max = local_maximum(Data$OsWL)
 #Find local minimum in the O-sWL series
 x.min = local_minimum(Data$OsWL)

 #Result vectors
 pre.low = numeric(length(Cluster_Max))
 fol.low = numeric(length(Cluster_Max))
 pre.high = numeric(length(Cluster_Max))
 fol.high = numeric(length(Cluster_Max))
 intensity = numeric(length(Cluster_Max))
 vol = numeric(length(Cluster_Max))

 #Loop repeated for each of the cluster maximum
 for(i in 1:length(Cluster_Max)){

  #Preceding high water level
  pre.high[i] = max(x.max[x.max<Cluster_Max[i]],na.rm=T)

  #Following high water level
  fol.high[i] = min(x.max[x.max>Cluster_Max[i]],na.rm=T)

  #Preceding low water level
  pre.low[i] = max(x.min[x.min<pre.high[i]],na.rm=T)

  #Following low water level
  fol.low[i] = min(x.min[x.min>fol.high[i]],na.rm=T)

  #Calculate surge "intensity"
  intensity[i] = sum(Data$OsWL[(pre.low[i]):(fol.low[i])][
    which(Data$OsWL[(pre.low[i]):(fol.low[i])]>
            mean(Data$OsWL,na.rm=T))]-
      mean(Data$OsWL,na.rm=T))

  #Calculate rainfall volume within +/- Rainfall_Interval hrs about water level peak
  vol[i] = sum(Data$Rainfall[(Cluster_Max[i]-Rainfall_Interval):
                             (Cluster_Max[i]+Rainfall_Interval)],na.rm = T)
 }

 #Put results in a data frame
 res = data.frame(pre.high,fol.high,pre.low,fol.low,intensity,vol)
 colnames(res) = c("Pre.High","Fol.High","Pre.Low","Fol.Low","Intensity","V")

 #Output data frame containing the results
 return(res)
}
