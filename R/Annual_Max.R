#' Generate annual maximum series
#'
#' Extract annual maximum in years with over a user-defined proportion of non-missing values.
#'
#' @param Data_Detrend Data frame containing two columns. In column: \itemize{
#' \item \code{1} A \code{"Date"} object of equally spaced discrete time steps.
#' \item \code{2} Numeric vector containing corresponding time series values.
#' }
#' @param Complete_Prop Minimum proportion of non-missing values in an annual record for the annual maximum to be extracted. Default is \code{0.8}.
#' @return List comprising the index of the annual maximum \code{Event} and the annual maximum values \code{AM}.
#' @export
#' @examples
#' Annual_Max(Data_Detrend=S20_T_MAX_Daily_Completed_Detrend$Detrend)
Annual_Max<-function(Data_Detrend, Complete_Prop=0.8){

  # Input validation
  if (!is.data.frame(Data_Detrend) && !is.matrix(Data_Detrend)) {
    stop("Data_Detrend must be a data frame or matrix")
  }

  if (ncol(Data_Detrend) != 2) {
    stop("Data_Detrend must comprise two columns, got: ", ncol(Data_Detrend))
  }

  if (!is.numeric(Complete_Prop) || length(Complete_Prop) != 1) {
    stop("Complete_Prop must be a single numeric value, got: ", class(Complete_Prop))
  }

  if (Complete_Prop > 1 || Complete_Prop < 0) {
    stop("Complete_Prop must be between 0 and 1, got: ", Complete_Prop)
  }

  date<-as.numeric(which(sapply(Data_Detrend,class)=="Date"))
  val<-c(1:ncol(Data_Detrend))[-date]

  years<-unique(year(Data_Detrend[,date]))

  prop.complete<-rep(NA,length(years))
  for(i in 1:length(years)){
    prop.complete[i]<-length(which(is.na(Data_Detrend[which(year(Data_Detrend[,date])==years[i]),val])==FALSE))/length(which(year(Data_Detrend[,date])==years[i]))
  }

  year.complete<-which(prop.complete>Complete_Prop)

  x.val<-rep(NA,length(year.complete))
  for(i in 1:length(year.complete)){
    x.val[i] <- (min(which(year(Data_Detrend[,date])==years[year.complete[i]]))-1)+which(Data_Detrend[which(year(Data_Detrend[,date])==years[year.complete[i]]),val] == max(Data_Detrend[which(year(Data_Detrend[,date])==years[year.complete[i]]),val],na.rm=T))
  }
  return(list("Event" = x.val, "AM" = Data_Detrend[x.val,val]))
}
