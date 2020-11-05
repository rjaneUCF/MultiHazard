#' Fit an Archimedean/elliptic copula model - Fit
#'
#' Fit a n-dimensional Archimedean or elliptic copula model. Function is simply a repackaging of the \code{fitCopula} function in the \code{copula} package.
#'
#' @param Data Data frame containing \code{n} at least partially concurrent time series. First column may be a \code{"Date"} object. Can be \code{Dataframe_Combine} output.
#' @param Copula_Type Type of elliptical copula to be fitted, options are \code{"Gaussian"} (Default), \code{"tcopula"}, \code{"Gumbel"}, \code{"Clayton"} and \code{"Frank"}.
#' @return List comprising the \code{Copula_Type} and the fitted copula \code{Model} object.
#' @seealso \code{\link{Dataframe_Combine}} \code{\link{Standard_Copula_Sel}} \code{\link{CDVineCopSelect}}  \code{\link{BiCopSelect}}
#' @export
#' @examples
#' cop<-Standard_Copula_Fit(Data=S20.Detrend.df,Copula_Type="Gaussian")
#' cop<-Standard_Copula_Fit(Data=S20.Detrend.df,Copula_Type="tcopula")
#' cop<-Standard_Copula_Fit(Data=S20.Detrend.df,Copula_Type="Gumbel")
#' cop<-Standard_Copula_Fit(Data=S20.Detrend.df,Copula_Type="Clayton")
#' cop<-Standard_Copula_Fit(Data=S20.Detrend.df,Copula_Type="Frank")
Standard_Copula_Fit<-function(Data,Copula_Type="Gaussian"){

  if(class(Data[,1])=="Date" | class(Data[,1])=="factor"){
    Data <- Data[,-1]
  }

  if(Copula_Type=="Gaussian"){
  norm.cop <- normalCopula(rep(0.6,ncol(Data)*(ncol(Data)-1)/2), dim=ncol(Data), dispstr="un")
  Gaussian_Copula_Fit <- fitCopula(norm.cop,data=pobs(na.omit(Data[,1:ncol(Data)])))
  Model <- normalCopula(c(as.numeric(coef(Gaussian_Copula_Fit))),dim=ncol(Data), dispstr="un")
  }

  if(Copula_Type=="tcopula"){
    t.cop<-tCopula(rep(0.6,ncol(Data)*(ncol(Data)-1)/2), dim=ncol(Data), dispstr="un",df=5)
    t_Copula_Fit <- fitCopula(t.cop,data=pobs(na.omit(Data[,1:ncol(Data)])))
    Model <- tCopula(c(as.numeric(coef(t_Copula_Fit)[1:(length(as.numeric(coef(t_Copula_Fit)))-1)])),dim=ncol(Data), dispstr="un",df=as.numeric(coef(t_Copula_Fit)[length(as.numeric(coef(t_Copula_Fit)))]))
  }

  if(Copula_Type=="Gumbel"){
  Gumbel_Copula_Fit <- fitCopula(gumbelCopula(dim=ncol(Data)),pobs(na.omit(Data[,1:ncol(Data)])))
  Model <- gumbelCopula(as.numeric(coef(Gumbel_Copula_Fit)[1]),dim=ncol(Data))
  }

  if(Copula_Type=="Clayton"){
    Clayton_Copula_Fit <- fitCopula(claytonCopula(dim=ncol(Data)),pobs(na.omit(Data[,1:ncol(Data)])))
    Model <- claytonCopula(as.numeric(coef(Clayton_Copula_Fit)[1]),dim=ncol(Data))
  }

  if(Copula_Type=="Frank"){
    Frank_Copula_Fit <- fitCopula(frankCopula(dim=ncol(Data)),pobs(na.omit(Data[,1:ncol(Data)])))
    Model <- frankCopula(as.numeric(coef(Frank_Copula_Fit)[1]),dim=ncol(Data))
  }

  return(Model)
}

