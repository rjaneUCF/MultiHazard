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
#' intensity = Intensity(Data=S13.Detrend.df[,c(1,3)],Cluster_Max=S13.OsWL.Declust$EventsMax)
#' #Plot O-sWL series identifying cluster maximum (in red) and print "intensity" above each maximum
#' plot(as.Date(S13.Detrend.df$Date_Time),
#'      S13.Detrend.df$OsWL)
#' points(as.Date(S13.Detrend.df$Date_Time[S13.OsWL.Declust$EventsMax]),
#'        S13.Detrend.df$OsWL[S13.OsWL.Declust$EventsMax],pch=16,col=2)
#' text(as.Date(S13.Detrend.df$Date_Time[S13.OsWL.Declust$EventsMax]),
#'      S13.Detrend.df$OsWL[S13.OsWL.Declust$EventsMax]+0.2,
#'      round(intensity$Intensity,0),cex=0.5)
Intensity<-function(Data,Cluster_Max,Base_Line="Mean"){

  # Check if Data is provided
  if (is.null(Data) | !is.data.frame(Data)) {
    stop("Data must be a data frame")
  }

  # Check if Data has rows and columns
  if (nrow(Data) == 0 | ncol(Data) == 0) {
    stop("Data must have at least one row and one column")
  }

  # Check if Cluster_Max is numeric
  if (is.null(Cluster_Max) | !is.numeric(Cluster_Max) | any(is.na(Cluster_Max)) | length(Cluster_Max) == 0 | any(Cluster_Max <= 0) | any(Cluster_Max > nrow(Data)) | any(Cluster_Max != round(Cluster_Max))) {
    stop("Cluster_Max must only contain positive integer values that do not exceed length of time series")
  }


  if(any(class(Data[,1]) %in% c("Date", "factor", "POSIXct", "character"))){
    Data <- Data[,-1]
  }

  if (is.data.frame(Data) || is.matrix(Data)) {
   # If multiple columns, use first one and warn
   if (ncol(Data) > 1) {
    warning("Data has multiple columns. Using first column only.")
    Data <- Data[, 1]
   } else {
    Data <- Data[, 1]  # Convert to vector
   }
  } else {
    Data <- as.vector(Data)
  }

  # Validate Base_Line parameter
  if (length(Base_Line) != 1) {
    stop("Base_Line must be a single value")
  }

  if ((!is.character(Base_Line) && !is.numeric(Base_Line)) | (is.character(Base_Line) && Base_Line != "Mean")) {
    stop("Base_Line must be either 'Mean' or a numeric value")
  }

  # Check if local_maximum and local_minimum functions exist
  if (!exists("local_maximum")) {
    stop("Function 'local_maximum' not found. Please ensure it is loaded or defined.")
  }

  if (!exists("local_minimum")) {
    stop("Function 'local_minimum' not found. Please ensure it is loaded or defined.")
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

 #Loop repeated for each of the cluster maximum
 for(i in 1:length(Cluster_Max)){

  #Preceding high water level
  pre_high = x.max[x.max<Cluster_Max[i]]
  if (length(pre_high) == 0 || any(is.infinite(pre_high))) {
   stop("No preceding non-infinite high water level found for event ", i)
  }
  pre.high[i] = max(pre_high, na.rm=T)

  #Following high water level
  fol_high = x.max[x.max>Cluster_Max[i]]
  if (length(fol_high) == 0 || any(is.infinite(fol_high))) {
    stop("No following non-infinite high water level found for event ", i)
  }
  fol.high[i] = min(fol_high, na.rm=T)

  #Preceding low water level
  pre_low = x.min[x.min<pre.high[i]]
  if (length(pre_low) == 0 || any(is.infinite(pre_low))) {
    stop("No preceding non-infinite low water level found for event ", i)
  }
  pre.low[i] = max(pre_low, na.rm=T)

  #Following low water level
  fol_low = x.min[x.min>fol.high[i]]
  if (length(fol_low) == 0 || any(is.infinite(fol_low))) {
    stop("No following non-infinite low water level found for event ", i)
  }
  fol.low[i] = min(fol_low, na.rm=T)


  #Identify which values are above Base_Line
  event = Data[(pre.low[i]):(fol.low[i])]
  above_baseline = which(event>Base_Line)

  #Calculate surge "intensity"
  if (!any(above_baseline)) {
    warning("No data above baseline for event ", i, ". Intensity will be 0.")
    intensity[i] = 0
  } else {
    intensity[i] = sum(event[above_baseline] - Base_Line)
  }

 }

 #Put results in a data frame
 res = data.frame(pre.high,fol.high,pre.low,fol.low,intensity)
 colnames(res) = c("Pre.High","Fol.High","Pre.Low","Fol.Low","Intensity")

 #Output data frame containing the results
 return(res)
}
