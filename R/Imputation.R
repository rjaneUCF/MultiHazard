#' Imputing missing values through linear regression
#'
#' Fits a simple linear regression model, to impute missing values of the dependent variable.
#'
#' @param Data Data frame containing two at least partially concurrent time series. First column may be a \code{"Date"} object. Can be \code{Dataframe_Combine} output.
#' @param Variable Character vector of length one specifying the (column) name of the variable to be imputed i.e. dependent variable in the fitted regression.
#' @param x_lab Character vector of length one specifying the name of the independent variable to appear as the x-axis label on a plot showing the data, imputed values and the linear regression model.
#' @param y_lab Character vector of length one specifying the name of the dependent variable to appear as the y-axis label on plot showing the data, imputed values and the linear regression model.
#' @return List comprising a\itemize{
#' \item \code{Data} data frame containing the original data plus an additional column named \code{Value} where the \code{NA} values of the \code{Variable} of interest have been imputed where possible.
#' \item \code{Model} linear regression model parameters including its coefficient of determination
#' } and a scatter plot of the data (black points), linear regression model (red line) and fitted (imputed) values (blue points).
#' @export
#' @examples
#' ####Objective: Fill in missing values at groundwater well G_3356 using record at G_3355
#' ##Viewing first few rows of G_3356
#' head(G_3356)
#' #Converting date column to a "Date" object
#' G_3356$Date<-seq(as.Date("1985-10-23"), as.Date("2019-05-29"), by="day")
#' #Converting readings to numeric object
#' G_3356$Value<-as.numeric(as.character(G_3356$Value))
#'
#' ##Viewing first few rows of G_3355
#' head(G_3355)
#' #Converting date column to a "Date" object
#' G_3355$Date<-seq(as.Date("1985-08-20"), as.Date("2019-06-02"), by="day")
#' #Converting readings to numeric object
#' G_3355$Value<-as.numeric(as.character(G_3355$Value))
#'
#' ##Merge the two dataframes by date
#' library('dplyr')
#' GW_S20<-merge(G_3356,G_3355,by="Date")
#' colnames(GW_S20)<-c("Date","G3356","G3355")
#' #Carrying out imputation
#' Imputation(Data=GW_S20,Variable="G3356",
#'            x_lab="Groundwater level (ft NGVD 29)",
#'            y_lab="Groundwater level (ft NGVD 29)")
Imputation<-function(Data,Variable,x_lab,y_lab){
  if(class(Data[,1])=="Date" | class(Data[,1])=="factor" | class(Data[,1])=="POSIXct" | class(Data[,1])=="character"){
    data <- Data[,-1]
  } else {
    data <- Data
  }
  variable<-which(names(data)==Variable)
  Other.variable<-c(1:ncol(data))[-which(names(data)==Variable)]
  data.NA<-which(is.na(data[,variable])==TRUE)
  Model<-lm(data[,variable] ~ data[,Other.variable])
  Data[,(ncol(Data)+1)]<-data[,variable]
  Data[data.NA,ncol(Data)]<-coef(Model)[1]+data[data.NA,Other.variable]*coef(Model)[2]
  names(Data)[ncol(Data)] <- 'ValuesFilled'

  #Plot linear model plus predicted points
  plot(data[-data.NA,Other.variable],Data[-data.NA,ncol(Data)],xlab=x_lab,ylab=y_lab,pch=16)
  lines(seq(0,9,0.1),coef(Model)[1]+seq(0,9,0.1)*coef(Model)[2],col=2)
  points(data[data.NA,Other.variable],Data[data.NA,ncol(Data)],col=4,pch=16)
  return(list(Data = Data, Model = summary(Model)))
}


