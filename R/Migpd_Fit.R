#' Fits Multiple independent generalized Pareto models - Fit
#'
#' Fit multiple independent generalized Pareto models to each column of a data frame. Edited version of the \code{migpd} function in \code{texmex}, to allow for \code{NA}s in a time series.
#'
#' @param Data A data frame with \code{n} columns, each comprising a declustered and if necessary detrended time series to be modelled.
#' @param Data_Full A data frame with \code{n} columns, each comprising the original (detrended if necessary) time series to be modelled. Only required if threshold is specified using \code{mqu}.
#' @param mth Marginal thresholds, above which generalized Pareto models are fitted. Numeric vector of length \code{n}.
#' @param mqu Marginal quantiles, above which generalized Pareto models are fitted. \strong{Only one of \code{mth} and \code{mqu} should be supplied.} Numeric vector of length \code{n}.
#' @param penalty See \code{\link{ggplot.migpd}}.
#' @param maxit See \code{\link{ggplot.migpd}}.
#' @param trace See \code{\link{ggplot.migpd}}.
#' @param verbose See \code{\link{ggplot.migpd}}.
#' @param priorParameters See \code{\link{ggplot.migpd}}.
#' @return An object of class \code{"migpd"}. There are \code{coef}, \code{print}, \code{plot}, \code{ggplot} and \code{summary} functions available.
#' @seealso \code{\link{Decluster}} \code{\link{Detrend}} \code{\link{Dataframe_Combine}}
#' @export
#' @examples
#' #With date as first column
#' S22.GPD<-Migpd_Fit(Data=S22.Detrend.Declustered.df, mqu =c(0.99,0.99,0.99))
#' #Without date as first column
#' S22.GPD<-Migpd_Fit(Data=S22.Detrend.Declustered.df[,-1], mqu =c(0.99,0.99,0.99))
#' #Same GPDs fit as above but thresholds given on the original scale
#' S22.Rainfall.Quantile<-quantile(na.omit(S22.Detrend.Declustered.df$Rainfall),0.99)
#' S22.OsWL.Quantile<-quantile(na.omit(S22.Detrend.Declustered.df$OsWL),0.99)
#' S22.GW.Quantile<-quantile(na.omit(S22.Detrend.Declustered.df$Groundwater),0.99)
#' S22.GPD<-Migpd_Fit(Data=S22.Detrend.Declustered.df[,-1],
#'                    mqu =c(S22.Rainfall.Quantile,S22.OsWL.Quantile,S22.GW.Quantile))
Migpd_Fit<-function (Data, Data_Full=NA, mth, mqu, penalty = "gaussian", maxit = 10000,
                      trace = 0, verbose = FALSE, priorParameters = NULL){

  if(class(Data[,1])=="Date" | class(Data[,1])=="factor"){
  data <- Data[,-1]
  } else {
  data <- Data
  }

  if(missing(mth)){
  if(class(Data_Full[,1])=="Date" | class(Data_Full[,1])=="factor"){
    data_full <- Data_Full[,-1]
  } else {
    data_full <- Data_Full
  }
  }

  theCall <- match.call()
  if (is.null(colnames(data))) {
    colnames(data) <- paste(rep("Column", ncol(data)), 1:ncol(data),
                            sep = "")
  }
  d <- dim(data)[2]
  if (missing(mth) & missing(mqu))
    stop("you must provide one of mth or mqu")
  if (!missing(mth) & !missing(mqu))
    stop("you must provide precisely one of mth or mqu")
  if (!missing(mth))
    mth <- rep(mth, length = d)
  if (!missing(mqu))
    mqu <- rep(mqu, length = d)
  if (missing(mqu))
    mqu <- sapply(1:d, function(i, x, mth) 1 - mean(x[, i] > mth[i]), x = data, mth = mth)
  if (missing(mth)){
    mth<-numeric(d)
  for(i in 1:d){
    prob <- mqu[i]
    mth[i] <- quantile(na.omit(data_full[,i]), prob)
  }
  }
  if (penalty %in% c("quadratic", "gaussian") & is.null(priorParameters)) {
    gp = list(c(0, 0), matrix(c(100^2, 0, 0, 0.25), nrow = 2))
    priorParameters <- vector("list", length = length(mth))
    for (i in 1:length(mth)) priorParameters[[i]] <- gp
    names(priorParameters) <- dimnames(data)[[2]]
  }
  else if (penalty %in% c("quadratic", "gaussian")) {
    nm <- names(priorParameters)
    if (is.null(nm))
      stop("priorParameters must be a named list")
    else if (any(!is.element(nm, dimnames(data)[[2]])))
      stop("the names of priorParameters must match the column names of the data")
  }
  wrapgpd <- function(i,x, mth, penalty, maxit, verbose, trace, priorParameters){
    if (verbose)
      cat("Fitting model", i, "\n")
    if (!is.null(priorParameters))
      priorParameters <- priorParameters[[(1:length(priorParameters))[names(priorParameters) ==
                                                                        dimnames(x)[[2]][i]]]]
    x <- c(na.omit(x[, i]))
    mth_ <- mth[i]
    evm(x, th = mth_, penalty = penalty, priorParameters = priorParameters,
        maxit = maxit, trace = trace)
  }
  modlist <- lapply(1:d,wrapgpd, x = data, penalty = penalty,
                    mth = mth, verbose = verbose, priorParameters = priorParameters,
                    maxit = maxit, trace = trace)
  if (length(dimnames(data)[[2]]) == dim(data)[[2]]) {
    names(modlist) <- dimnames(data)[[2]]
  }
  names(mth) <- names(mqu) <- dimnames(data)[[2]]
  res <- list(call = theCall, models = modlist, data = data,
              mth = mth, mqu = mqu, penalty = penalty, priorParameters = priorParameters)
  oldClass(res) <- "migpd"
  invisible(res)
}


