#' Declusters a Summed time series using a moving (Storm) Window approach
#'
#' Finds the sum of a time series within a moving window then declusters the summed series using another moving window.
#'
#' @param Data Data frame containing two columns. In column: \itemize{
#' \item \code{1} A \code{"Date"} object of equally spaced discrete time steps.
#' \item \code{2} Numeric vector containing corresponding time series values.
#' }
#' @param Window_Width_Sum Numeric vector of length one specifying the window width over which to sum the data.
#' @param Window_Width Numeric vector of length one specifying the width, in days, of the window used to ensure events are independent.
#' @return List comprising vectors containing the original time series \code{Detrend}, the summed series \code{Totals}, independent (declustered) events \code{Declustered}, the elements of the original series containing the start (\code{Event_Start}), center \code{EventID}, and end (\code{Event_End}) of the declustered events. Note for \code{Window_Width_Sum_Type="End"}, \code{Event_End} and \code{EventID} are identical.
#' @export
#' @examples
#' #Declustering 24 hour rainfall totals at site S13 using a 7-day window for declustering the events.
#' plot(S13_Precip$Date,S13_Rainfall$Rainfall)
#' S13_Rainfall_Totals_Declust<-Decluster_S_SW(Data=S13_Rainfall, Window_Width_Sum=24,
#'                                             Window_Width=7*24)
#' plot(S13_Rainfall[,1],
#'      S13_Rainfall_Totals_Declust$Totals,
#'      pch=16,ylim=c(0,10))
#' points(S13_Rainfall[S13_Rainfall_Totals_Declust$EventID,1],
#'        S13_Rainfall_Totals_Declust$Totals[S13_Rainfall_Totals_Declust$EventID],
#'        col=2,pch=16)
Decluster_S_SW<-function(Data, Window_Width_Sum, Window_Width) {

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
  if (any(diff(Data[,1]) <= 0)) {
    stop("Dates in first column must be in ascending order with no duplicates")
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

  # Check if Window_Width is numeric
  if (!is.numeric(Window_Width_Sum)) {
    stop("Window_Width_Sum must be numeric")
  }

  # Check if Window_Width is positive
  if (Window_Width_Sum <= 0) {
    stop("Window_Width_Sum must be positive")
  }

  if (length(Window_Width_Sum) != 1 || Window_Width_Sum != round(Window_Width_Sum)) {
    stop("Window_Width_Sum must be a single integer value")
  }

  #Index of NA elements will be added to z
  z<-0

  #Assigning any NA values to a very small number relative to the data
  if (length(which(is.na(Data[,2]) == T)) > 0) {
    z <- which(is.na(Data[,2]) == T)
    Data[z,2] <- min(Data[,2], na.rm = T) - 1000
  }

  #Summing observations using a centered window
  if(Window_Width_Sum>1){
    Sum<-c(rep(NA,round((Window_Width_Sum - Window_Width_Sum %% 2)/2,0)),c(cumsum(Data[,2])[-c(1:(Window_Width_Sum -1))]) -
             c(0,cumsum(Data[,2])[-c((length(Data[,2])-Window_Width_Sum+1):length(Data[,2]))]),
           rep(NA,round((Window_Width_Sum + Window_Width_Sum %% 2)/2-1,0)))
  }
  if(Window_Width_Sum==1){
    Sum<-Data[,2]
  }

  #Declustering
  Sum_DF<-data.frame(Data[,1],Sum)
  Decl<-Decluster_SW(Data=Sum_DF, Window_Width)

  #Results vector
  Event_Start<-rep(NA,length(Decl$EventID))
  Event_End<-rep(NA,length(Decl$EventID))

  ##Finding the start and end of the declustered events
  if(Window_Width_Sum>1){
    Event_Start=Decl$EventID-round((Window_Width_Sum - Window_Width_Sum %% 2)/2,0)
    Event_End=Decl$EventID+round((Window_Width_Sum + Window_Width_Sum %% 2)/2-1,0)
  }
  if(Window_Width_Sum==1){
    Event_Start=NA
    Event_End=NA
  }

  if (min(z) > 0) {
    Data[z,2] <- NA
    z<-z+rep((-round((Window_Width_Sum-1)/2,0)):(round(Window_Width_Sum/2,0)),each=length(z))
    z<-z[z>1]
    z<-z[z<nrow(Data)]
    Sum[z]<-NA
    Decl$Declustered[z]<-NA
  }

  #Create a list of outputs
  res <- list("Detrend"= Data[, 2], "Totals" = Sum, "Declustered" = Decl$Declustered, "EventID" = Decl$EventID, "Event_Start"=Event_Start, "Event_End"=Event_End)
  return(res)
}
