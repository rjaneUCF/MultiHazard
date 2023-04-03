#' Conditionally sampling a two-dimensional dataset
#'
#' Creates a data frame where the declustered excesses of a (conditioning) variable are paired with co-occurences of another variable.
#'
#' @param Data_Detrend Data frame containing two at least partially concurrent time series, detrended if necessary. Time steps must be equally spaced, with missing values assigned \code{NA}. First column may be a \code{"Date"} object. Can be \code{Dataframe_Combine} output.
#' @param Data_Declust Data frame containing two (independently) declustered at least partially concurrent time series. Time steps must be equally spaced, with missing values assigned \code{NA}. Columns must be in the same order as in \code{Data_Detrend}. First column may be a \code{"Date"} object. Can be \code{Dataframe_Combine} output.
#' @param Con_Variable Column number (1 or 2) or the column name of the conditioning variable. Default is \code{1}.
#' @param u Threshold, as a quantile of the observations of the conditioning variable. Default is \code{0.97}.
#' @param Thres Threshold expressed on the original scale of the observations. Only one of \code{u} and \code{Thres} should be supplied. Default is \code{NA}.
#' @return List comprising the specified \code{Threshold} as the quantile of the conditioning variable above which declustered excesses are paired with co-occurences of the other variable, the resulting two-dimensional sample \code{data} and \code{name} of the conditioning variable. The index of the input dataset that correspond to the events in the conditonal sample \code{x.con} are also provided.
#' @export
#' @examples
#' S20.Rainfall<-Con_Sampling_2D(Data_Detrend=S20.Detrend.df[,-c(1,4)],
#'                               Data_Declust=S20.Detrend.Declustered.df[,-c(1,4)],
#'                               Con_Variable="Rainfall",u=0.97)
Con_Sampling_2D<-function(Data_Detrend,Data_Declust,Con_Variable,u=0.97,Thres=NA){
  if(class(Data_Detrend[,1])[1]=="Date" | class(Data_Detrend[,1])[1]=="factor" | class(Data_Detrend[,1])[1]=="POSIXct" | class(Data_Detrend[,1])[1]=="character"){
    Data_Detrend<-Data_Detrend[,-1]
  }
  if(class(Data_Declust[,1])[1]=="Date" | class(Data_Declust[,1])[1]=="factor" | class(Data_Declust[,1])[1]=="POSIXct" | class(Data_Declust[,1])[1]=="character"){
    Data_Declust<-Data_Declust[,-1]
  }
  if(is.numeric(Con_Variable)==FALSE){
    con<-which(names(Data_Detrend)==Con_Variable)
    noncon <- c(1,2)[-con]
  } else{
    con<-Con_Variable
    noncon<-c(1,2)[-con]
  }
  if(is.na(Thres)==T){
    Thres<-quantile(na.omit(Data_Detrend[,con]), u)
  }
  x<-which(Data_Declust[,con]>=Thres)
  Sample_df<-matrix(c(0),nrow=length(x),ncol=2)
  for(i in 1:length(x)){
    Sample_df[i,con]<-Data_Declust[x[i],con]
    Sample_df[i,noncon]<-Data_Detrend[x[i],noncon]
  }
  if(length(which(is.na(Sample_df[,1])==TRUE | is.na(Sample_df[,2])==TRUE))>0){
    z<-which(is.na(Sample_df[,1])==TRUE | is.na(Sample_df[,2])==TRUE)  #unique(which(is.na(Sample_df[,1] | Sample_df[,2])==TRUE))
    Sample_df<-Sample_df[-z,]
    x<-x[-z]
  }
  Sample_df<-data.frame(Sample_df)
  colnames(Sample_df)<-c(names(Data_Detrend))
  return(list("Threshold"=Thres,"Data"=Sample_df, "Con_Variable"= names(Data_Detrend)[con],"x.con"=x))
}

