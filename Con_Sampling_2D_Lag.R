#' Conditionally sampling a two dimensional dataset
#'
#' Creates a data frame where the declustered excesses of a (conditioning) variable are paired with the maximum value of a second variable over a specified lag.
#'
#' @param Data_Detrend Data frame containing two at least partially concurrent time series, detrended if necessary. Time steps must be equally spaced, with missing values assigned \code{NA}. First object may be a \code{"Date"} object. Can be \code{Dataframe_Combine} output.
#' @param Data_Declust Data frame containing two (independently) declustered at least partially concurrent time series. Time steps must be equally spaced, with missing values assigned \code{NA}. Columns must be in the same order as in \code{Data_Detrend}. First object may be a \code{"Date"} object. Can be \code{Dataframe_Combine} output.
#' @param Con_Variable Column number (1 or 2) or the column name of the conditioning variable. Default is \code{1}.
#' @param Thres Threshold, as a quantile of the observations of the conditioning variable. Default is \code{0.97}.
#' @param Lag_Backward Positieve lag applied to variable not assigned as the \code{Con_Variable}. Default is \code{0}
#' @param Lag_Forward Negative lag to variable not assigned as the \code{Con_Variable}. Default is \code{0}
#' @return List comprising the specifyied \code{Threshold} as the quantile of the conditioning variable above which declustered excesses are paired with co-occurences of the other variable, the resulting two-dimensional sample \code{data} and \code{name} of the conditioning variable.
#' @export
#' @examples
#' S20.Rainfall<-Con_Sampling_2D(Data_Detrend=S20.Detrend.df[,-c(1,4)],
#'                               Data_Declust=S20.Detrend.Declustered.df[,-c(1,4)],
#'                               Con_Variable="Rainfall",Thres=0.97)
Con_Sampling_2D_Lag<-
  function (Data_Detrend, Data_Declust, Con_Variable, Thres = 0.97, Lag_Backward = 0, Lag_Forward = 0)
  {
    if (class(Data_Detrend[, 1]) == "Date" | class(Data_Detrend[,1]) == "factor" | class(Data_Detrend[,1])=="POSIXct") {
      Data_Detrend <- Data_Detrend[, -1]
    }
    if (class(Data_Declust[, 1]) == "Date" | class(Data_Declust[,1]) == "factor" | class(Data_Declust[,1])=="POSIXct") {
      Data_Declust <- Data_Declust[, -1]
    }
    if (is.numeric(Con_Variable) == FALSE) {
      con <- which(names(Data_Detrend) == Con_Variable)
      noncon <- c(1, 2)[-con]
    }
    else {
      con <- Con_Variable
      noncon <- c(1, 2)[-con]
    }
    x.con <- which(Data_Declust[, con] > quantile(na.omit(Data_Detrend[,con]), Thres))
    x.noncon <- numeric(length(x.con))
    Sample_df <- matrix(c(0), nrow = length(x.con), ncol = 2)
    for (i in 1:length(x.con)) {
      Sample_df[i, con] <- Data_Declust[x.con[i], con]
      Sample_df[i, noncon] <- max(Data_Detrend[c(max((x.con[i]-Lag_Backward),1):(min((x.con[i]+Lag_Forward),nrow(Data_Declust)))), noncon],na.rm=T)[1]
      position<- which(Data_Detrend[c(max((x.con[i]-Lag_Backward),1):(min((x.con[i]+Lag_Forward),nrow(Data_Declust)))), noncon] ==
                         max(Data_Detrend[c(max((x.con[i]-Lag_Backward),1):(min((x.con[i]+Lag_Forward),nrow(Data_Declust)))), noncon]))
      position<-ifelse(length(position)<2,position,ifelse(position==2,2,sample(position,1)))
      x.noncon[i] <- x.con[i] + ((-Lag_Backward):Lag_Forward)[position]
    }
    if(length(which(Sample_df[,1]==-Inf | Sample_df[,2]==-Inf))>0){
      z<-unique(which(Sample_df[,1]==-Inf | Sample_df[,2]==-Inf))
      Sample_df[z,]<-NA
      x.con[z]<-NA
      x.noncon[z]<-NA
    }
    if (length(which(is.na(Sample_df[, 1]) == TRUE | is.na(Sample_df[,2]) == TRUE)) > 0) {
      z <- which(is.na(Sample_df[, 1]) == TRUE | is.na(Sample_df[,2]) == TRUE)
      Sample_df <- Sample_df[-z, ]
      x.con <- x.con[-z]
      x.noncon <- x.noncon[-z]
    }
    Sample_df <- data.frame(Sample_df)
    colnames(Sample_df) <- c(names(Data_Detrend))
    return(list(Threshold = quantile(na.omit(Data_Detrend[, con]),Thres), Data = Sample_df, Con_Variable = names(Data_Detrend)[con], x.con = x.con, x.noncon = x.noncon))
  }


