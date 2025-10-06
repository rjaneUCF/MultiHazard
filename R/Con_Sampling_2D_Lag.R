#' Conditionally sampling a two dimensional dataset
#'
#' Creates a data frame where the declustered excesses of a (conditioning) variable are paired with the maximum value of a second variable over a specified time-lag.
#'
#' @param Data_Detrend Data frame containing two at least partially concurrent time series, detrended if necessary. Time steps must be equally spaced, with missing values assigned \code{NA}. First object may be a \code{"Date"} object. Can be \code{Dataframe_Combine} output.
#' @param Data_Declust Data frame containing two (independently) declustered at least partially concurrent time series. Time steps must be equally spaced, with missing values assigned \code{NA}. Columns must be in the same order as in \code{Data_Detrend}. First object may be a \code{"Date"} object. Can be \code{Dataframe_Combine} output.
#' @param Con_Variable Column number (1 or 2) or the column name of the conditioning variable. Default is \code{1}.
#' @param u Threshold, as a quantile of the observations of the conditioning variable. Default is \code{0.97}.
#' @param Thres Threshold expressed on the original scale of the observations. Only one of \code{u} and \code{Thres} should be supplied. Default is \code{NA}.
#' @param Lag_Backward Positive lag applied to variable not assigned as the \code{Con_Variable}. Default is \code{3}
#' @param Lag_Forward Negative lag to variable not assigned as the \code{Con_Variable}. Default is \code{3}
#' @return List comprising the specified \code{Threshold} as the quantile of the conditioning variable above which declustered excesses are paired with co-occurences of the other variable, the resulting two-dimensional sample \code{data} and \code{Con_Variable} the name of the conditioning variable. The index of the input dataset that correspond to the events of the conditioning variable \code{x.con} and the non-conditioning variable \code{x.noncon} in the conditonal sample are also provided.
#' @export
#' @examples
#' S20.Rainfall<-Con_Sampling_2D_Lag(Data_Detrend=S20.Detrend.df[,-c(1,4)],
#'                                   Data_Declust=S20.Detrend.Declustered.df[,-c(1,4)],
#'                                   Con_Variable="Rainfall",u=0.97)
Con_Sampling_2D_Lag<-
  function (Data_Detrend, Data_Declust, Con_Variable, u = 0.97, Thres=NA, Lag_Backward = 3, Lag_Forward = 3)
  {

    # Input validation
    if (!is.data.frame(Data_Detrend) && !is.matrix(Data_Detrend)) {
      stop("Data_Detrend must be a data frame or matrix")
    }

    if (!is.data.frame(Data_Declust) && !is.matrix(Data_Declust)) {
      stop("Data_Declust must be a data frame or matrix")
    }

    if (ncol(Data_Detrend) != 2) {
      stop("Data_Detrend must comprise two columns, got: ", ncol(Data_Detrend))
    }

    if (ncol(Data_Declust) != 2) {
      stop("Data_Declust must comprise two columns, got: ", ncol(Data_Declust))
    }

    if (!is.numeric(u) || length(u) != 1) {
      stop("u must be a single numeric value, got: ", class(u))
    }

    if (!Con_Variable %in% colnames(Data_Detrend)) {
      stop("Con_Variable must the name of a column in Data_Detrend")
    }

    if (!Con_Variable %in% colnames(Data_Declust)) {
      stop("Con_Variable must be the name of a column in Data_Declust")
    }

    if (u > 1 || u < 0) {
      stop("u must be between 0 and 1, got: ", u)
    }

    if (!is.numeric(Lag_Backward)) {
      stop("Lag_Backward must be a single numeric value, got: ", class(Lag_Backward))
    }

    if (!is.numeric(Lag_Forward)) {
      stop("Lag_Forward must be a single numeric value, got: ", class(Lag_Forward))
    }

    if (class(Data_Detrend[, 1])[1] == "Date" | class(Data_Detrend[,1])[1] == "factor" | class(Data_Detrend[,1])[1]=="POSIXct" | class(Data_Detrend[,1])[1] == "character") {
      Data_Detrend <- Data_Detrend[, -1]
    }
    if (class(Data_Declust[, 1])[1] == "Date" | class(Data_Declust[,1])[1] == "factor" | class(Data_Declust[,1])[1]=="POSIXct" | class(Data_Declust[,1])[1] == "character") {
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
    if(is.na(Thres)){
      Thres<-quantile(na.omit(Data_Detrend[,con]), u)
    }
    x.con <- which(Data_Declust[, con] >= Thres)
    x.noncon <- numeric(length(x.con))
    Sample_df <- matrix(c(0), nrow = length(x.con), ncol = 2)
    for (i in 1:length(x.con)) {
      Sample_df[i, con] <- Data_Declust[x.con[i], con]
      lag_values <- Data_Detrend[c(max((x.con[i]-Lag_Backward),1):(min((x.con[i]+Lag_Forward),nrow(Data_Declust)))), noncon]

      if(all(is.na(lag_values))) {
        Sample_df[i, noncon] <- NA
        x.noncon[i] <- NA
      } else {
        max_val <- max(lag_values, na.rm = TRUE)
        Sample_df[i, noncon] <- max_val
        position <- which(lag_values == max_val)
        position <- ifelse(length(position) < 2, position, ifelse(position == 2, 2, sample(position, 1)))
        x.noncon[i] <- x.con[i] + ((-Lag_Backward):Lag_Forward)[position]
      }
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
    return(list("Threshold" = Thres, "Data" = Sample_df, "Con_Variable" = names(Data_Detrend)[con], "x.con" = x.con, "x.noncon" = x.noncon))
  }

