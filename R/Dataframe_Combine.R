#' Creates a data frame containing up to five time series
#'
#' Combines up to five time series, detrended where necessary, into a single data frame.
#'
#' @param data.1 Data frames with two columns containing in column \itemize{
#' \item 1 - Continuous sequence of times spanning from the first to the final recorded observations.
#' \item 2 - Corresponding values detrended where necessary.
#' }
#' @param data.2 As for \code{data.1}.
#' @param data.3 As for \code{data.1}.
#' @param data.4 As for \code{data.1}.
#' @param data.5 As for \code{data.1}.
#' @param n Integer \code{1-5} specifying the number of time series. Default is \code{2}.
#' @param names Character vector giving the column names excluding the first column which is labelled as "Date".
#' @return A data frame containing all times from the first to the most up to date reading of any of the variables.
#' @seealso \code{\link{Detrend}}
#' @export
#' @examples
#' #Formatting data
#' S20.Rainfall.df<-Perrine_df
#' S20.Rainfall.df$Date<-as.Date(S20.Rainfall.df$Date)
#' S20.OsWL.df<-S20_T_MAX_Daily_Completed_Detrend_Declustered[,c("Date","ValueFilled")]
#' S20.OsWL.df$Date<-as.Date(S20.OsWL.df$Date)
#' #Detrending O-sWL series at Site S20
#' S20.OsWL.Detrend<-Detrend(Data=S20.OsWL.df,Method = "window",PLOT=FALSE,
#'                          x_lab="Date",y_lab="O-sWL (ft NGVD 29)")
#' #Creating a dataframe with the date alongside the detrended OsWL series
#' S20.OsWL.Detrend.df<-data.frame(as.Date(S20.OsWL.df$Date),S20.OsWL.Detrend)
#' colnames(S20.OsWL.Detrend.df)<-c("Date","OsWL")
#' #Combining the two datasets by Date argument
#' S20.Detrend.df<-Dataframe_Combine(data.1<-S20.Rainfall.df,
#'                                  data.2<-S20.OsWL.Detrend.df,
#'                                  names=c("Rainfall","OsWL"))
Dataframe_Combine<-function(data.1,data.2,data.3=NA,data.4=NA,data.5=NA,n=2,names){

  data_Detrend_1_df<-data.frame(data.1[,1],data.1[,2])
  colnames(data_Detrend_1_df)<-c("Date",colnames(data.1)[2])

  data_Detrend_2_df<-data.frame(data.2[,1],data.2[,2])
  colnames(data_Detrend_2_df)<-c("Date",colnames(data.2)[2])

  if(n==2){
    Detrend_df<-full_join(data_Detrend_1_df, data_Detrend_2_df, by="Date")
  }

  if(n==3){
    data_Detrend_3_df<-data.frame(data.3[,1],data.3[,2])
    colnames(data_Detrend_3_df)<-c("Date",colnames(data.3)[2])

    data_Detrend_df1<-full_join(data_Detrend_1_df, data_Detrend_2_df, by="Date")
    Detrend_df<-full_join(data_Detrend_df1, data_Detrend_3_df, by="Date")
  }

  if(n==4){
    data_Detrend_3_df<-data.frame(data.3[,1],data.3[,2])
    colnames(data_Detrend_3_df)<-c("Date",colnames(data.3)[2])

    data_Detrend_4_df<-data.frame(data.4[,1],data.4[,2])
    colnames(data_Detrend_4_df)<-c("Date",colnames(data.4)[2])

    data_Detrend_df1<-full_join(data_Detrend_1_df, data_Detrend_2_df, by="Date")
    data_Detrend_df2<-full_join(data_Detrend_df1, data_Detrend_3_df, by="Date")
    Detrend_df<-full_join(data_Detrend_df2, data_Detrend_4_df, by="Date")
  }

  if(n==5){
    data_Detrend_3_df<-data.frame(data.3[,1],data.3[,2])
    colnames(data_Detrend_3_df)<-c("Date",colnames(data.3)[2])

    data_Detrend_4_df<-data.frame(data.4[,1],data.4[,2])
    colnames(data_Detrend_4_df)<-c("Date",colnames(data.4)[2])

    data_Detrend_5_df<-data.frame(data.5[,1],data.5[,2])
    colnames(data_Detrend_5_df)<-c("Date",colnames(data.5)[2])

    data_Detrend_df1<-full_join(data_Detrend_1_df, data_Detrend_2_df, by="Date")
    data_Detrend_df2<-full_join(data_Detrend_df1, data_Detrend_3_df, by="Date")
    data_Detrend_df3<-full_join(data_Detrend_df2, data_Detrend_4_df, by="Date")
    Detrend_df<-full_join(data_Detrend_df3, data_Detrend_5_df, by="Date")
  }

  Detrend_df[order(as.Date(Detrend_df$Date, format="%d/%m/%Y")),]
  colnames(Detrend_df)<-c("Date",names)
  return(Detrend_df)
}

