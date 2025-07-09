#' Fit an Archimedean/elliptic copula model - Fit
#'
#' Fit a n-dimensional Archimedean or elliptic copula model. Function is simply a repackaging of the \code{fitCopula} function in the \code{copula} package.
#'
#' @param Data Data frame containing \code{n} at least partially concurrent time series. First column may be a \code{"Date"} object. Can be \code{Dataframe_Combine} output.
#' @param Copula_Type Type of elliptical copula to be fitted, options are \code{"Gaussian"} (Default), \code{"tcopula"}, \code{"Gumbel"}, \code{"Clayton"} and \code{"Frank"}.
#' @return List comprising the \code{Copula_Type} and the fitted copula \code{Model} object.
#' @seealso \code{\link{Dataframe_Combine}} \code{\link{Standard_Copula_Sel}}
#' @export
#' @examples
#' cop<-Standard_Copula_Fit(Data=S20.Detrend.df,Copula_Type="Gaussian")
#' cop<-Standard_Copula_Fit(Data=S20.Detrend.df,Copula_Type="tcopula")
#' cop<-Standard_Copula_Fit(Data=S20.Detrend.df,Copula_Type="Gumbel")
#' cop<-Standard_Copula_Fit(Data=S20.Detrend.df,Copula_Type="Clayton")
#' cop<-Standard_Copula_Fit(Data=S20.Detrend.df,Copula_Type="Frank")
Standard_Copula_Fit<-function(Data,Copula_Type="Gaussian"){

  # Input validation
  if(missing(Data) || is.null(Data)) {
    stop("Error: Data is empty")
  }

  if(!is.data.frame(Data) && !is.matrix(Data)) {
    stop("Error: Data must be a data.frame or matrix")
  }

  if(nrow(Data) < 10) {
    stop("Error: Data must have at least 10 rows")
  }

  if(ncol(Data) < 2) {
    stop("Error: Data must have at least 2 columns for copula fitting")
  }

  valid_copula_types <- c("Gaussian", "tcopula", "Gumbel", "Clayton", "Frank")
  if(!Copula_Type %in% valid_copula_types) {
    stop(paste("Error: Copula_Type must be one of:", paste(valid_copula_types, collapse=", ")))
  }

  # Remove first column if it's Date, factor, POSIXct, or character
  if(class(Data[,1])[1]=="Date" | class(Data[,1])[1]=="factor" | class(Data[,1])[1]=="POSIXct" | class(Data[,1])[1] == "character"){
    Data <- Data[,-1]
  }

  # Check if we still have enough columns after potential removal
  if(ncol(Data) < 2) {
    stop("Error: After removing non-numeric first column, Data must still have at least 2 columns")
  }

  # Check if all remaining columns are numeric
  numeric_cols <- sapply(Data, is.numeric)
  if(!all(numeric_cols)) {
    stop("Error: All data columns must be numeric for copula fitting")
  }

  # Check for sufficient non-missing data
  Data <- na.omit(Data[,1:ncol(Data)])
  if(nrow(Data) < 10) {
    stop("Error: Insufficient non-missing data (need at least 10 complete observations)")
  }

  # Check for constant columns (no variance)
  if(any(apply(Data, 2, var) == 0)) {
    stop("Error: One or more columns have zero variance (constant values)")
  }

  # Convert to pseudo-observations
  pseudo_obs <- pobs(Data)


  if(Copula_Type=="Gaussian"){
    norm.cop <- normalCopula(rep(0.6,ncol(Data)*(ncol(Data)-1)/2), dim=ncol(Data), dispstr="un")
    Gaussian_Copula_Fit <- fitCopula(norm.cop,data=pseudo_obs)
    Model <- normalCopula(c(as.numeric(coef(Gaussian_Copula_Fit))),dim=ncol(Data), dispstr="un")
  }

  if(Copula_Type=="tcopula"){
    t.cop<-tCopula(rep(0.6,ncol(Data)*(ncol(Data)-1)/2), dim=ncol(Data), dispstr="un",df=5)
    t_Copula_Fit <- fitCopula(t.cop,data=pseudo_obs)
    Model <- tCopula(c(as.numeric(coef(t_Copula_Fit)[1:(length(as.numeric(coef(t_Copula_Fit)))-1)])),dim=ncol(Data), dispstr="un",df=as.numeric(coef(t_Copula_Fit)[length(as.numeric(coef(t_Copula_Fit)))]))
  }

  if(Copula_Type=="Gumbel"){
    Gumbel_Copula_Fit <- fitCopula(gumbelCopula(dim=ncol(Data)),pseudo_obs)
    Model <- gumbelCopula(as.numeric(coef(Gumbel_Copula_Fit)[1]),dim=ncol(Data))
  }

  if(Copula_Type=="Clayton"){
    Clayton_Copula_Fit <- fitCopula(claytonCopula(dim=ncol(Data)),pseudo_obs)
    Model <- claytonCopula(as.numeric(coef(Clayton_Copula_Fit)[1]),dim=ncol(Data))
  }

  if(Copula_Type=="Frank"){
    Frank_Copula_Fit <- fitCopula(frankCopula(dim=ncol(Data)),pseudo_obs)
    Model <- frankCopula(as.numeric(coef(Frank_Copula_Fit)[1]),dim=ncol(Data))
  }

  return(Model)
}
