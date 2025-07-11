#' Derive water level curves
#'
#' Generates water level curves for simulated extreme water levels based on a simulated "intensity".
#'
#' @param Data A data frame of the time series with the column containing ocean-side water levels labeled \code{"OsWL"}.
#' @param Cluster_Max Numeric vector containing indexes of peaks in the O-sWL column of \code{Data}. If analyzing a sample conditioned on O-sWL derived using \code{Con_Sample_2D()} set equal to the \code{$x.con} output.
#' @param Pre_Low Numeric vector of the indexes of the O-sWL column in \code{Data} containing the preceding low water level.
#' @param Fol_Low Numeric vector of the indexes of the O-sWL column in \code{Data} containing the following low water level.
#' @param Thres Numeric vector of length one, specifying threshold above which to apply the method. Below the threshold an observed curve with an intensity less than \code{limit} is randomly sampled.
#' @param Base_Line Numeric vector of length one, specifying water level about which to calculate the intensity. Default is the mean O-sWL.
#' @param Limit Numeric vector of length one, specifying an upper limit on the observed water level curve intensities to sample for simulated peaks less than \code{Thres}.
#' @param Peak Numeric vector of simulated peak water levels.
#' @param Intensity Numeric vector of the intensity associated with each simulated \code{Peak}.
#' @param Length Numeric vector of length one, specifying the length of time over which the water level curve is simulated before (and after) the time of the simulated peak. Total duration of the water level curve is 2*\code{Length}+1. Minimum is \code{5}. Default is \code{144}.
#' @return A data frame, where each row contains the water level curve generated for corresponding simulated peak in the \code{Peak} input. A vector of the intensity \code{Intensity} of the generated water level curve.
#' @seealso \code{\link{Surge_Criterion}} \code{\link{OsWL_Intensity}}
#' @export
#' @examples
#' #Declustering O-sWL series
#' S13.OsWL.Declust = Decluster(Data=S13.Detrend.df$OsWL,
#'                              SepCrit=24*7, u=0.99667)
#' #Use O-sWL intensity function to obtain index of preceding and following low water levels
#' intensity = OsWL_Intensity(Data=S13.Detrend.df,Cluster_Max=S13.OsWL.Declust$EventsMax)
#'
#' #Four synthetic events
#' sim.peaks = c(3.4,4,4.2,5)
#' sim.intensity = c(38,48,120,140)
#'
#' #Generating the water level curves
#' oswl_ts_oswl = WL_Curve(Data = S13.Detrend.df,
#'                         Cluster_Max = S13.OsWL.Declust$EventsMax,
#'                         Pre_Low = intensity$Pre.Low,
#'                         Fol_Low = intensity$Fol.Low,
#'                         Thres = S13.OsWL.Declust$Threshold, Limit = 45,
#'                         Peak = sim.peaks,
#'                         Intensity = sim.intensity)
#'
#' #Plot the water level curves of the observed peaks
#' plot(1:289,
#'      S13.Detrend.df$OsWL[(S13.OsWL.Declust$EventsMax[1]-144):
#'                          (S13.OsWL.Declust$EventsMax[1]+144)],
#'      type='l',ylim=c(1,5))
#' for(i in 2:length(S13.OsWL.Declust$EventsMax-144)){
#'   lines(1:289,
#'         S13.Detrend.df$OsWL[(S13.OsWL.Declust$EventsMax[i]-144):
#'                             (S13.OsWL.Declust$EventsMax[i]+144)])
#' }
#' #Superimpose the curves generated for the four synthetic events
#' for(i in 1:4){
#'   lines(1:289,oswl_ts_oswl$Series[i,],col=2)
#' }
WL_Curve<-function(Data,Cluster_Max,Pre_Low,Fol_Low,Thres,Base_Line=mean(Data$OsWL,na.rm=T),Limit,Peak,Intensity,Length=144){

  if (is.null(Data) | !is.data.frame(Data)){
    stop("Data must be a data frame.")
  }

  if (!"OsWL" %in% names(Data)){
    stop("Data must contain 'OsWL' column.")
  }

  if (any(Cluster_Max < 1 | Cluster_Max > length(Data$OsWL))){
    stop("Cluster_Max contains invalid indices.")
  }

  if (any(Pre_Low < 1 | Fol_Low > length(Data$OsWL))){
     stop("Pre_Low or Fol_Low contains out-of-bound indices.")
  }

  if (length(Peak) != length(Intensity)){
    stop("Peak and Intensity vectors must be of equal length.")
  }

  # Validate Base_Line parameter
  if (length(Base_Line) != 1) {
    stop("Base_Line must be a single value")
  }

  if ((!is.character(Base_Line) && !is.numeric(Base_Line)) | (is.character(Base_Line) && Base_Line != "Mean")) {
    stop("Base_Line must be either 'Mean' or a numeric value")
  }

  # Remove Cluster_Max values that would cause out-of-bounds issues
  valid_indices <- which((Cluster_Max - max(Length, 5) >= 1) &
                           (Cluster_Max + max(Length, 5) <= length(Data$OsWL)))

  if (length(valid_indices) == 0) {
    stop("No Cluster_Max values are within valid bounds given the specified Length.")
  }

  # Filter all relevant vectors to match only valid indices
  Cluster_Max <- Cluster_Max[valid_indices]
  Pre_Low     <- Pre_Low[valid_indices]
  Fol_Low     <- Fol_Low[valid_indices]

  #Vectors for storing the results
  intensity.scaled = numeric(length(Cluster_Max))
  intensity.event = Intensity
  series = matrix(NA_real_,nrow=length(Intensity),ncol=(2*Length+1))
  e = numeric(length(Peak))

  #Repeat this loop for every simulated peak
  for(k in 1:length(Peak)){

    #For simulated peaks less than the threshold
    if(Peak[k]<Thres){

      #Repeat for each observed peak
      for(j in 1:length(Cluster_Max)){

        #Calculating intensity of the observed peaks
        ts = Data$OsWL[Pre_Low[j]:Fol_Low[j]]
        intensity.scaled[j] = sum(ts[which(ts>Base_Line)]-Base_Line)
      }
      #intensity.scaled[20] = 1000

      #Sample an intensity less than the specified limit
      valid_ce <- which(intensity.scaled < Limit)
      if (length(valid_ce) == 0) stop("No events found with intensity below 'Limit'.")
      ce <- sample(valid_ce, 1)

      #Rescale the sampled event
      new = (Peak[k]-Data$OsWL[Cluster_Max[ce]])+Data$OsWL[(Cluster_Max[ce]-Length):(Cluster_Max[ce]+Length)]

      #Event is defined as the water level curve from preceding low water level to following low water level
      event = (Peak[k]-Data$OsWL[Cluster_Max[ce]])+Data$OsWL[Pre_Low[ce]:Fol_Low[ce]]

      #Compute the intensity
      intensity.event[k] = sum(event[event>Base_Line]-Base_Line)

    } else{
      #For simulated peaks above the threshold
      #Calculating the intensity of the shifted events and choosing the event with the closest intensity less than the simulated intensity
      for(j in 1:length(Cluster_Max)){

        #Re-scale observed water level curve around the time of the peak to coincide with simulated peak
        ts = Data$OsWL[Pre_Low[j]:Fol_Low[j]]
        ts[Cluster_Max[j]-4] = 0.3*(Peak[k]-Data$OsWL[(Cluster_Max[j]-5)])+Data$OsWL[(Cluster_Max[j]-5)]
        ts[Cluster_Max[j]-3] = 0.5*(Peak[k]-Data$OsWL[(Cluster_Max[j]-5)])+Data$OsWL[(Cluster_Max[j]-5)]
        ts[Cluster_Max[j]-2] = 0.7*(Peak[k]-Data$OsWL[(Cluster_Max[j]-5)])+Data$OsWL[(Cluster_Max[j]-5)]
        ts[Cluster_Max[j]-1] = 0.9*(Peak[k]-Data$OsWL[(Cluster_Max[j]-5)])+Data$OsWL[(Cluster_Max[j]-5)]
        ts[Cluster_Max[j]] = Peak[k]
        ts[Cluster_Max[j]+1] = 0.9*(Peak[k]-Data$OsWL[(Cluster_Max[j]+5)])+Data$OsWL[(Cluster_Max[j]+5)]
        ts[Cluster_Max[j]+2] = 0.7*(Peak[k]-Data$OsWL[(Cluster_Max[j]+5)])+Data$OsWL[(Cluster_Max[j]+5)]
        ts[Cluster_Max[j]+3] = 0.5*(Peak[k]-Data$OsWL[(Cluster_Max[j]+5)])+Data$OsWL[(Cluster_Max[j]+5)]
        ts[Cluster_Max[j]+4] = 0.3*(Peak[k]-Data$OsWL[(Cluster_Max[j]+5)])+Data$OsWL[(Cluster_Max[j]+5)]

        #Calculate intensity of re-scaled curve
        intensity.scaled[j] = sum(ts[which(ts>Base_Line)]-Base_Line)

      }

      #Calculate difference between simulated intensity (input) and the intensity of the re-scaled observed events
      d.oswl = intensity.scaled-Intensity[k]

      #
      #d.oswl[20] =1000

      #Select the re-scaled curve with the largest "intensity" that's smaller than the simulated "intensity"
      if (all(d.oswl >= 0)) {
        ce <- which.min(d.oswl)
      } else {
        ce <- which.max(d.oswl[d.oswl < 0])
        ce <- which(d.oswl == d.oswl[ce])[1]
      }

      #Decrease from peak expressed as a proportion of peak O-sWL
      prop = - (Data$OsWL[Pre_Low[ce]:Fol_Low[ce]]-Data$OsWL[Cluster_Max[ce]])/
        Data$OsWL[Cluster_Max[ce]]

      #Set the decreases around the peak to zero as they are re-scaled differently
      prop[(Cluster_Max[ce]-Pre_Low[ce]+1-4):(Cluster_Max[ce]-Pre_Low[ce]+1+4)] = 0

      #Decreases sum to 1
      if (sum(prop) == 0) stop("Cannot normalize zero-sum proportions.")
      prop.stan = prop / sum(prop)

      #Difference in intensity between closest (re-scaled) observed event and simulated event
      OsWL.diff = Intensity[k]-intensity.scaled[ce]

      #Intensity units that need to be added are distributed in proportion with
      #the decreases from the peak OsWL. Larger the distance from the peak the more units added.
      #To create the new OsWL series

      #Add new units to the selected observed event from -144 to 144 time intervals
      new = Data$OsWL[(Cluster_Max[ce]-Length):(Cluster_Max[ce]+Length)]

      #Rescale occurs from time of preceding til time of following low water
      new[(Length-(Cluster_Max[ce]-Pre_Low[ce])+1):(Length+Fol_Low[ce]-Cluster_Max[ce]+1)] = (OsWL.diff*prop.stan)+Data$OsWL[Pre_Low[ce]:Fol_Low[ce]]

      #Water levels about the peak are re-scaled differently to try to maintain smoothness
      new[Length-3] = 0.3*(Peak[k]-Data$OsWL[(Cluster_Max[ce]-5)])+Data$OsWL[(Cluster_Max[ce]-5)]
      new[Length-2] = 0.5*(Peak[k]-Data$OsWL[(Cluster_Max[ce]-5)])+Data$OsWL[(Cluster_Max[ce]-5)]
      new[Length-1] = 0.7*(Peak[k]-Data$OsWL[(Cluster_Max[ce]-5)])+Data$OsWL[(Cluster_Max[ce]-5)]
      new[Length] = 0.9*(Peak[k]-Data$OsWL[(Cluster_Max[ce]-5)])+Data$OsWL[(Cluster_Max[ce]-5)]
      new[Length+1] = Peak[k]
      new[Length+2] = 0.9*(Peak[k]-Data$OsWL[(Cluster_Max[ce]+5)])+Data$OsWL[(Cluster_Max[ce]+5)]
      new[Length+3] = 0.7*(Peak[k]-Data$OsWL[(Cluster_Max[ce]+5)])+Data$OsWL[(Cluster_Max[ce]+5)]
      new[Length+4] = 0.5*(Peak[k]-Data$OsWL[(Cluster_Max[ce]+5)])+Data$OsWL[(Cluster_Max[ce]+5)]
      new[Length+5] = 0.3*(Peak[k]-Data$OsWL[(Cluster_Max[ce]+5)])+Data$OsWL[(Cluster_Max[ce]+5)]

      }

    #Assign water level curve to the appropriate row in the results data frame
    series[k,] = new
    e[k] = ce
  }

  #Put result objects into a list
  res = list("Series" = series, "Intensity" = intensity.event, "Event"=e)

  #Output data frame containing the results
  return(res)
}
