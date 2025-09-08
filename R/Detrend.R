#' Detrends a time series.
#'
#' Detrends a time series using either a linear fit covering the entire dataset or moving average trend correction with a user-specified window width.
#'
#' @param Data Data frame containing two columns. In column: \itemize{
#' \item \code{1} A \code{"Date"} object of equally spaced discrete time steps.
#' \item \code{2} Numeric vector containing corresponding time series values. No \code{NA}s allowed.
#' }
#' @param Method Character vector of length one specifying approach used to detrend the data. Options are moving average \code{"window"} (default) and \code{"linear"}.
#' @param Window_Width Numeric vector of length one specifying length of the moving average window. Default is \code{89}, window comprises the observation plus 44 days either side, which for daily data corresponds to an approximate 3 month window.
#' @param End_Length Numeric vector of length one specifying number of observations at the end of the time series used to calculate the present day average. Default is \code{1826}, which for daily data corresponds to the final five years of observations.
#' @param PLOT Logical; whether to plot original and detrended series. Default is \code{"FALSE"}.
#' @param x_lab Character vector of length one specifying x-axis label. Default is \code{"Date"}.
#' @param y_lab Character vector of length one specifying y-axis label. Default is \code{"Data"}.
#' @return Numeric vector of the detrended time series.
#' @export
#' @examples
#' #Detrending ocean-side water level at site S22 using a 3 month moving average window and the last
#' #five years of observations to calculate the present day average.
#' Detrend(S22_T_MAX_Daily_Completed,Method = "window",Window_Width= 89,
#'         End_Length = 1826, PLOT=FALSE,x_lab="Data",y_lab="Data")
Detrend<-function(Data, Method = "window",Window_Width= 89, End_Length = 1826, PLOT=FALSE,x_lab="Date",y_lab="Data"){
  data_Detrend<-Data[,2]

  if (!Method %in% c("window", "linear")) {
    stop("Method must be either 'window' or 'linear', got: '", Method, "'")
  }

  if(Method=="window"){
    for(i in 1:(Window_Width/2)){
      data_Detrend[i]<- Data[i,2] - mean(Data[i:(i+(Window_Width/2)),2])
    }

    for(i in ((Window_Width/2)+1):(nrow(Data)-(Window_Width/2))){
      data_Detrend[i]<- Data[i,2] - mean(Data[(i-(Window_Width/2)):(i+(Window_Width/2)),2])
    }
    data_Detrend[1:(length(Data[,2])-(Window_Width/2))] <- data_Detrend[1:(length(Data[,2])-(Window_Width/2))] + mean(Data[(length(Data[,2])-(End_Length+1)):length(Data[,2]),2])
  }

  if(Method=="linear"){
    x<-seq(1,length(Data[,1]),1)
    model <- lm(Data[,2] ~ x)
    residuals <-  Data[,2] - predict(model, data.frame(x))
    data_Detrend <- residuals + mean(Data[(length(Data[,2])-(End_Length+1)):(length(Data[,2])),2], na.rm=T)
  }


  if(PLOT==TRUE){
    plot(as.Date(Data[,1]),Data[,2],type='l',col=1,xlab=x_lab,ylab=y_lab)
    lines(as.Date(Data[,1]),data_Detrend,col=2)
  }
  return(data_Detrend)
}
