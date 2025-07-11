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

  # Check: Input is a data frame or matrix
  if (!is.data.frame(Data) && !is.matrix(Data)) {
    stop("Error: Data must be a data frame or matrix.")
  }

  # Check: Data must have numeric columns (for pobs)
  # If first column is non-numeric, check the rest
  if (inherits(Data[,1], c("Date", "POSIXct", "factor", "character"))) {
    if (!all(sapply(Data[,-1], is.numeric))) {
      stop("Error: All columns except the first must be numeric.")
    }
  } else {
    if (!all(sapply(Data, is.numeric))) {
      stop("Error: All columns must be numeric.")
    }
  }

  # Check: Enough rows
  if (nrow(Data) < 5) {
    stop("Error: Data must contain at least 5 rows.")
  }

  #Proceeds to fit copula model
  if (inherits(Data[,1], c("Date", "factor", "POSIXct", "character"))){
    M<-RVineStructureSelect(pobs(na.omit(Data[,2:ncol(Data)])))
    Model <- RVineCopSelect(pobs(na.omit(Data[,2:ncol(Data)])),Matrix=M$Matrix)
  } else {
    M<-RVineStructureSelect(pobs(na.omit(Data[,1:ncol(Data)])))
    Model <- RVineCopSelect(pobs(na.omit(Data[,1:ncol(Data)])),Matrix=M$Matrix)
  }

  #Output
  res<- list(Structure = M$Matrix, Family=Model$family, Par = Model$par, Par2 = Model$par2)

  #Returns output
  return(res)
}
