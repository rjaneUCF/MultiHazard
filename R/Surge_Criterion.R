#' Surge identification criterion
#'
#' Classify extreme water levels as either tidally dominated or surge driven.
#'
#' @param Data A data frame with co-occurring rainfall and O-sWL time series in two columns labeled \code{"Rainfall"} and \code{"OsWL"}, respectively.
#' @param Cluster_Max Numeric vector containing indexes of peaks in the O-sWL column of \code{Data}. If analyzing a sample conditioned on O-sWL derived using \code{Con_Sample_2D()} set equal to the \code{$xcon} output.
#' @param Criterion_Number Numeric vector of length one, specifying which of the five criterion detailed in the report to adopt. If a user-defined criterion is adopted set to \code{NA} which is the default.
#' @param Surge_Thres Numeric vector of length one, specifying the minimum elevation difference between a peak and prior maximum or minimum for the peak to be classified as surge driven. Default is \code{0.25}.
#' @param Rainfall_Thres Numeric vector of length one, specifying minimum rainfall within a +/- \code{Rainfall_Interval} period of a peak for the peak to be classified as surge driven. Default is \code{NA}.
#' @param Pre_Sur Numeric vector of length one, specifying, minimum length of time allowed between preceding maximum or minimum and the peak. Default is \code{7}.
#' @param MaxMin Character vector of length one, specifying whether elevation difference refers to the preceding minimum (\code{"Min"}) or maximum (\code{"Max"}). Default is \code{"Max"}.
#' @param Rainfall_Interval Numeric vector of length one, specifying length of time before and after a peak over which to sum rainfall totals. Total window width is 2*\code{Rainfall_Interval}+1. Default is \code{NA}.
#' @return A vector with each cluster maximum classified as either \code{Tide} or \code{Surge} driven.
#' @seealso \code{\link{Con_Sampling_2D}}
#' @export
#' @examples
#' #Decluster O-sWL series at S-13 using a runs method
#' S13.OsWL.Declust = Decluster(Data=S13.Detrend.df$OsWL,
#'                              SepCrit=24*7, u=0.99667)
#' #Classify peak water levels as either surge or tidally driven
#' surge_class = Surge_Criterion(Data = S13.Detrend.df,
#'                              Cluster_Max = S13.OsWL.Declust$EventsMax,
#'                              Criterion_Number = 5)
#' #Plot O-sWL time series with peaks the color of peaks representing classification
#'  S13.Detrend.df$Date_Time =  as.POSIXct(S13.Detrend.df$Date_Time)
#' plot(S13.Detrend.df$Date_Time,S13.Detrend.df$OsWL)
#' points(S13.Detrend.df$Date_Time[S13.OsWL.Declust$EventsMax],
#'       S13.Detrend.df$OsWL[S13.OsWL.Declust$EventsMax],
#'       col=ifelse(surge_class=="Tide","Blue","Red"),pch=16)
#' legend("topleft",c("Tide","Surge"),pch=16,col=c("Blue","Red"))
#' #Example of a custom surge criterion. Peak is classified as tidal if
#' #Elevation difference between peak and preceding minimum at least 7 hrs before is less than 0.25.
#' #Total rainfall from 72 hours before and to 72 hrs after the peak is less than 2 Inches
#' surge_class = Surge_Criterion(Data = S13.Detrend.df,
#'                               Cluster_Max = S13.OsWL.Declust$EventsMax,
#'                               Surge_Thres=2.5,Rainfall_Thres=2,Pre_Sur=7,
#'                               MaxMin="Min",Rainfall_Interval=72)
Surge_Criterion<-function(Data,Cluster_Max,Criterion_Number="NA",Surge_Thres=0.25,Rainfall_Thres=NA,Pre_Sur=7,MaxMin="Max",Rainfall_Interval=NA){

 #Vectors for storing the results
 surge_est_max<-rep(NA,length(Cluster_Max))
 surge_est_min<-rep(NA,length(Cluster_Max))
 surge_est<-rep(NA,length(Cluster_Max))
 Total_Rainfall<-rep(NA,length(Cluster_Max))

 #Minimum length of time allowed between preceding maximum or minimum and the peak O-sWL.
 Pre_Sur = ifelse(Criterion_Number == 2 | Criterion_Number == 3, 12,
              ifelse(Criterion_Number == 4 | Criterion_Number == 5, 7, Pre_Sur))

 #Time interval over which total rainfall volumes are calculated for relevant pre-defined criterion
 Rainfall_Interval = ifelse(Criterion_Number == 5, 72, Rainfall_Interval)

 #Rainfall threshold for an O-sWL peak to be classified as a surge event for relevant pre-defined criterion
 Rainfall_Thres = ifelse(Criterion_Number == 5, 2.5, Rainfall_Thres)

 #Identify local maximum in the O-sWL series
 local_max = local_maximum(Data$OsWL)

 #Identify local minimum in the O-sWL series
 local_min = local_minimum(Data$OsWL)

 for(i in 1:length(Cluster_Max)){

  #Difference in elevation between peak and previous peak that occurred at least Pre_Sur observations prior
  surge_est_max[i] = Data$OsWL[Cluster_Max[i]]-Data$OsWL[max(which(local_max<Cluster_Max[i]-Pre_Sur))]

  #Difference in elevation between peak and previous low water level that occurred at least Pre_Sur observations prior
  surge_est_min[i] = Data$OsWL[Cluster_Max[i]]-Data$OsWL[max(which(local_min<Cluster_Max[i]-Pre_Sur))]
 }

 #Keep set of elevation differences specified by MaxMin input
  surge_est = ifelse(MaxMin=="Max",1,0)*surge_est_max+ifelse(MaxMin=="Min",1,0)*surge_est_min

 #Calculate rainfall about the time of peak water level if instructed to by set of inputs
 if(is.na(Rainfall_Interval)==F){
  for(i in 1:length(Cluster_Max)){

   #Total rainfall volume within 3-days of the peak
   Total_Rainfall[i]<-sum(Data$Rainfall[(Cluster_Max[i]-Rainfall_Interval):(Cluster_Max[i]+Rainfall_Interval)],na.rm=T)

  }

  #Apply criterion specified in input when it involves rainfall
  res = ifelse(surge_est<Surge_Thres & Total_Rainfall<Rainfall_Thres, "Tide", "Surge")

  } else{

  #Apply criterion specified in input when it doesn't involve rainfall
  res = ifelse(surge_est<Surge_Thres, "Tide", "Surge")

  }

 #Overwriting res with pre-defined surge criterion if Criterion_Number is given in input
 if(Criterion_Number==1){ res = ifelse(surge_est<1.5, "Tide", "Surge") }
 if(Criterion_Number==2){ res = ifelse(surge_est<1.5, "Tide", "Surge") }
 if(Criterion_Number==3){ res = ifelse(surge_est<2, "Tide", "Surge") }
 if(Criterion_Number==4){ res = ifelse(surge_est<0.25, "Tide", "Surge") }
 if(Criterion_Number==5){ res = ifelse(surge_est<0.25 & Total_Rainfall<2.5, "Tide", "Surge") }

 #Output results
 return(res)
}
