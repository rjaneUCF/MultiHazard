#' C and D-vine Copula - Fitting
#'
#' Fit either a C- or D-vine copula model. Function is a repackaging the \code{RVineStructureSelect} and \code{RVineCopSelect} functions from the \code{RVine} package into a single function.
#'
#' @param Data Data frame containing \code{n} at least partially concurrent time series. First column may be a \code{"Date"} object. Can be \code{Dataframe_Combine} output.
#' @return List comprising the vine copula \code{Structure}, pair-copula families composing the C- or D-vine copula \code{Family}, its parameters \code{Par} and \code{Par2}.
#' @seealso \code{\link{Dataframe_Combine}} \code{\link{Vine_Copula_Sim}}
#' @export
#' @examples
#' S20.Vine<-Vine_Copula_Fit(Data=S20.Detrend.df)
Vine_Copula_Fit<-function(Data){
  if(class(Data[,1])[1]=="Date" | class(Data[,1])[1]=="factor" | class(Data[,1])[1]=="POSIXct" | class(Data[,1])[1] == "character"){
    M<-RVineStructureSelect(pobs(na.omit(Data[,2:ncol(Data)])))
    Model <- RVineCopSelect(pobs(na.omit(Data[,2:ncol(Data)])),Matrix=M$Matrix)
  } else {
    M<-RVineStructureSelect(pobs(na.omit(Data[,1:ncol(Data)])))
    Model <- RVineCopSelect(pobs(na.omit(Data[,1:ncol(Data)])),Matrix=M$Matrix)
  }
  res<- list(Structure = M$Matrix, Family=Model$family, Par = Model$par, Par2 = Model$par2)
  return(res)
}
