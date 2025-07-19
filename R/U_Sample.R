#' Implements (unconditional) bootstrap procedure in Serinaldi and Kilsby (2013)
#'
#' Implements the unconditional bootstrap procedure i.e. peak is not conditioned on duration outlined in Serinaldi and Kilsby (2013) to generate non-peak rainfall totals for a simulated peak. The function also calculates hyetograph properties including net characteristics.
#'
#' @param Data Vector of the rainfall time series.
#' @param Cluster_Max Vector of the index of \code{Data} containing the cluster maximum. If declustering is carried out using \code{Decluster_SW()} set equal to \code{$EventsMax} output.
#' @param D Numeric vector of the duration of the cluster maximum events.
#' @param Xp Numeric vector of simulated peaks. To implement the method exactly as in Serinaldi and Kilsby (2013), set equal to a sample (taken with replacement) of the observed cluster maximum (peaks).
#' @param Start Numeric vector of the index of \code{Data} where each cluster maximum event begins.
#' @param End Numeric vector of the index of \code{Data} where each cluster maximum event ends.
#' @return A data frame with the following columns:  \itemize{
#' \item \code{Xp} Simulated event peaks i.e. input \code{Xp}.
#' \item \code{D} Duration sampled from the duration vector \code{D} for each simulated event.
#' \item \code{Samp} Index of the cluster maximum event, sampled conditionally on \code{D}, that provides non-peak rainfall depths.
#' \item \code{V} Volume of simulated events.
#' \item \code{Vn} Net volume of simulated events.
#' \item \code{I} Intensity of simulated events.
#' \item \code{In} Net intensity of simulated events.
#' \item \code{Start} Index of \code{Data} where the sampled (\code{Samp}) event begins.
#' \item \code{End} Index of \code{Data} where the sampled (\code{Samp}) event ends.
#' }
#' @seealso \code{\link{Decluster}} \code{\link{Time_Series_Plot}} \code{\link{WL_Curve}}
#' @export
#' @examples
#' #First decluster the rainfall series to find the 500 events
#' #with the highest peaks
#' S13.Rainfall.Declust = Decluster(Data=S13.Detrend.df$Rainfall,
#'                                  SepCrit=24*3, u=0.99667)
#' #Set very small rainfall measurements to zero.
#' #Assumed to be the result of uncertainty in measuring equipment.
#' S13_Rainfall$Rainfall[which(S13_Rainfall$Rainfall<0.01)] = 0
#' #Find NAs in rainfall series
#' z = which(is.na(S13_Rainfall$Rainfall)==T)
#' #Temporarily set NAs to zero
#' S13_Rainfall$Rainfall[z] = 0
#' #Find times where there is 6-hours of no rainfall
#' no.rain = rep(NA,length(S13_Rainfall$Rainfall))
#' for(i in 6:length(S13_Rainfall$Rainfall)){
#'   no.rain[i] = ifelse(sum(S13_Rainfall$Rainfall[(i-5):i])==0,i,NA)
#' }
#' #Remove NAs from results vector as these correspond to times where there is
#' #rainfall at certain points in the 6 hour period.
#' no.rain = na.omit(no.rain)
#' #Reset missing values in the rainfall record back to NA
#' S13_Rainfall$Rainfall[z] = NA
#' #Find the start and end times of the 500 events.
#' start = rep(NA,length(S13.Rainfall.Declust$EventsMax))
#' end = rep(NA,length(S13.Rainfall.Declust$EventsMax))
#' for(i in 1:length(S13.Rainfall.Declust$EventsMax)){
#'  start[i] = max(no.rain[which(no.rain<S13.Rainfall.Declust$EventsMax[i])])
#'  end[i] = min(no.rain[which(no.rain>S13.Rainfall.Declust$EventsMax[i])])
#' }
#' start = start + 1
#' end = end - 6
#' d = end - start + 1 #Duration
#' #Simulate some peaks by sampling observed peaks with replacement
#' #I.e., applying the method exactly as in Serinaldi and Kilsby (2013)
#' sim.peak = sample(S13.Rainfall.Declust$EventsMax,size=500,replace=TRUE)
#' #Derive the hyetographs
#' S13.oswl.sample = U_Sample(Data=S13_Rainfall$Rainfall,
#'                            Cluster_Max=S13.Rainfall.Declust$EventsMax,
#'                            D=d,Start=start,End=end,
#'                            Xp=sim.peak)
U_Sample<-function(Data,Cluster_Max,D,Start,End,Xp){

  # Input validation
  if (any(is.na(Data))) {
    warning("Data contains NA values")
  }

  if (any(is.na(c(Cluster_Max, D, Start, End, Xp)))) {
    stop("Input vectors cannot contain NA values")
  }

  if (length(Cluster_Max) != length(D) || length(D) != length(Start) ||
      length(Start) != length(End)) {
    stop("Cluster_Max, D, Start, and End must have the same length")
  }

  if (any(Start > End)) {
    stop("Start indices must be less than or equal to End indices")
  }

  if (any(Cluster_Max < 1) || any(Cluster_Max > length(Data))) {
    stop("Cluster_Max indices are out of bounds for Data")
  }

  if (any(Start < 1) || any(End > length(Data))) {
    stop("Start/End indices are out of bounds for Data")
  }

  if (length(Xp) <= 0) {
    stop("Xp must contain at least one value")
  }

  #Sample size
  Nsim = length(Xp)

  #Vectors for storing the results
  Samp = rep(NA,Nsim)
  v.net.csample = rep(NA,Nsim)
  v.csample = rep(NA,Nsim)
  s = rep(NA,Nsim)
  e = rep(NA,Nsim)

  #Sample a duration for each peak
  d.csample = sample(D,Nsim,replace=T)

  for(i in 1:Nsim){

    #Conditionally sample a cluster maximum event with the sampled duration (d.csample)
    di = which(D==d.csample[i])

    # Check if any matches found
    if (length(di) == 0) {
      warning(paste("No matches found for duration", d.csample[i], "at iteration", i))
      next
    }

    Samp[i] = ifelse(length(di)>1,sample(di,1),di)

    #Indexes of rainfall time series associated with sampled cluster maximum event
    s[i] = Start[Samp[i]]
    e[i] = End[Samp[i]]

    #Calculate net volume
    v.net.csample[i] = sum(Data[(Start[Samp[i]]):(End[Samp[i]])])-Data[Cluster_Max[Samp[i]]]

    #Calculate volume
    v.csample[i] = sum(Data[(Start[Samp[i]]):(End[Samp[i]])])-Data[Cluster_Max[Samp[i]]]+Xp[i]
  }

  #Calculate intensity
  I.csample = v.csample/d.csample

  #Calculate net intensity
  I.net = v.net.csample/d.csample

  #Put results in a data frame
  res = data.frame(Xp,d.csample,Samp,v.csample,v.net.csample,I.csample,I.net,s,e)
  colnames(res) = c("Xp","D","Samp","V","Vn","I","In","Start","End")

  #Output data frame containing the results
  return(res)
}
