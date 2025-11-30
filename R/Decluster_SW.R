#' Declusters a time series using a storm window approach
#'
#' Find peaks with a moving window. The code is based on the \code{IDEVENT} function provided by Sebastian Solari.
#'
#' @param Data Data frame containing two columns. In column: \itemize{
#' \item \code{1} A \code{"Date"} object of equally spaced discrete time steps.
#' \item \code{2} Numeric vector containing corresponding time series values.
#' }
#' @param Window_Width Numeric vector of length one specifying the width, in days, of the window used to ensure events are independent.
#' @return List comprising vectors containing the original time series \code{Detrended}, independent (declustered) events \code{Declustered} and the elements of the original series containing the declustered events \code{EventID}.
#' @export
#' @examples
#' #Declustering the O-sWL at site S22 using a 3-day window.
#' v<-Decluster_SW(Data=S22.Detrend.df[,c(1:2)],Window_Width=7)
#' plot(as.Date(S22.Detrend.df$Date),S22.Detrend.df$Rainfall,pch=16)
#' points(as.Date(S22.Detrend.df$Date)[v$EventID],v$Event,col=2,pch=16)
Decluster_SW<-function(Data,Window_Width){

  # Input validation
  if (missing(Data)) {
    stop("Data parameter is required")
  }

  if (missing(Window_Width)) {
    stop("Window_Width parameter is required")
  }

  # Check if Data is a data frame
  if (!is.data.frame(Data)) {
    stop("Data must be a data frame")
  }

  # Check if Data has exactly 2 columns
  if (ncol(Data) != 2) {
    stop("Data must have exactly 2 columns (Date and numeric values)")
  }

  # Check first column is Date-like
  if (!inherits(Data[,1], c("Date", "POSIXct", "POSIXt"))) {
    # Try to convert to Date if it's character
    if (is.character(Data[,1])) {
      tryCatch({
        Data[,1] <- as.Date(Data[,1])
      }, error = function(e) {
        stop("First column must be a Date object or convertible to Date")
      })
    } else {
      stop("First column must be a Date object")
    }
  }

  # Check for NA values in Date column
  if (any(is.na(Data[,1]))) {
    stop("First column (Date) contains NA values")
  }

  # Check if dates are in ascending order
  if (class(Data[,1])[1]=="Date"){
   if(any(diff(Data[,1]) <= 0)) {
    stop("Dates in first column must be in ascending order with no duplicates")
   }
  }

  # Check if second column contains numeric data
  if (!is.numeric(Data[,2])) {
    stop("Second column of Data must be numeric")
  }

  # Check if Window_Width is numeric
  if (!is.numeric(Window_Width)) {
    stop("Window_Width must be numeric")
  }

  # Check if Window_Width is positive
  if (Window_Width <= 0) {
    stop("Window_Width must be positive")
  }

  if (length(Window_Width) != 1 || Window_Width != round(Window_Width)) {
    stop("Window_Width must be a single integer value")
  }

  N<-1:length(Data[,1])
  EVENT.INDEX<-rep(0,length(Data[,1]))
  EVENT<-rep(0,length(Data[,1]))
  ID1   = 1
  for(i in 1:(length(Data[,2]))){
    ID  = which(N>=max(1,N[i]-(Window_Width-1)/2) & N<=min(max(N),N[i]+(Window_Width-1)/2))
    AUX = Data[ID,2]
    MAX = if(any(!is.na(AUX))) max(AUX, na.rm=T) else NA
    IDMAX = if(!is.na(MAX)) which(AUX==MAX)[1] else 1
    if(is.na(MAX)==F & N[ID[IDMAX]]==N[i]){
      EVENT.INDEX[ID1]<-N[i]
      EVENT[ID1]<-MAX
      ID1 = ID1+1
    }
  }
  EVENT<-EVENT[EVENT.INDEX>0]
  EVENT.INDEX<-EVENT.INDEX[EVENT.INDEX>0]
  Declustered<-rep(NA,nrow(Data))
  Declustered[EVENT.INDEX]<-Data[EVENT.INDEX,2]
  res<-list("Detrend"=Data[,2],"Declustered"=Declustered,"EventID"=EVENT.INDEX)
  return(res)
}
