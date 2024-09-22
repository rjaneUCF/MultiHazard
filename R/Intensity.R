#' Intensity
#'
#' Calculates the "intensity" of extreme water levels, as defined in Wahl et al. (2011).
#'
#' @param Data A data frame containing the water level time series. First column may be a \code{"Date"} object.
#' @param Cluster_Max Numeric vector containing indexes of the water level peaks in \code{Data}. If analyzing a sample conditioned on water level derived using \code{Con_Sample_2D()} set equal to the \code{$xcon} output.
#' @param Base_Line Vector of length one, specifying water level about which to calculate the intensity. Default is \code{"Mean"} where the mean of the entire time series is used as the baseline water level above which intensity is calculated.
#' @return A data frame with the following columns:  \itemize{
#' \item \code{Pre.High} Index of the water level column of \code{Data} containing the preceding high water level.
#' \item \code{Fol.High} Index of the water level column of \code{Data} containing the following high water level.
#' \item \code{Pre.Low} Index of the water level column of \code{Data} containing the preceding  low water level.
#' \item \code{Fol.Low} Index of the water level column of \code{Data} containing the following low water level.
#' \item \code{Intensity} Intensity of the extreme water level event.
#' }
#' @seealso \code{\link{Decluster}} \code{\link{WL_Curve}}
#' @export
#' @examples
#' #Decluster O-sWL series at S-13 using a runs method
#' S13.OsWL.Declust = Decluster(Data=S13.Detrend.df$OsWL,
#'                             SepCrit=24*7, u=0.99667)
#' #Calculate O-sWL of the identified cluster maximum
#' intensity = Intensity(Data=S13.Detrend.df,Cluster_Max=S13.OsWL.Declust$EventsMax)
#' #Plot O-sWL series identifying cluster maximum (in red) and print "intensity" above each maximum
#' plot(as.Date(S13.Detrend.df$Date_Time),
#'      S13.Detrend.df$OsWL)
#' points(as.Date(S13.Detrend.df$Date_Time[S13.OsWL.Declust$EventsMax]),
#'        S13.Detrend.df$OsWL[S13.OsWL.Declust$EventsMax],pch=16,col=2)
#' text(as.Date(S13.Detrend.df$Date_Time[S13.OsWL.Declust$EventsMax]),
#'      S13.Detrend.df$OsWL[S13.OsWL.Declust$EventsMax]+0.2,
#'      round(intensity$Intensity,0),cex=0.5)
Intensity<-function(Data,Cluster_Max,Base_Line="Mean"){

  if(class(Data[,1])[1]=="Date" | class(Data[,1])[1]=="factor" | class(Data[,1])[1]=="POSIXct" | class(Data[,1])[1]=="character"){
    Data<-Data[,-1]
  }

 #Calculating Base Line from which intensity is calculated
 Base_Line = ifelse(Base_Line == "Mean", mean(Data,na.rm=T),Base_Line)

 #Find local maximum in the time series
 x.max = local_maximum(Data)
 #Find local minimum in the time series
 x.min = local_minimum(Data)

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
  intensity[i] = sum(Data[(pre.low[i]):(fol.low[i])][
    which(Data[(pre.low[i]):(fol.low[i])]>Base_Line)]- Base_Line)
 }

 #Put results in a data frame
 res = data.frame(pre.high,fol.high,pre.low,fol.low,intensity)
 colnames(res) = c("Pre.High","Fol.High","Pre.Low","Fol.Low","Intensity")

 #Output data frame containing the results
 return(res)
}
