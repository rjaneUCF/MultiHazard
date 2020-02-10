#' C and D-vine Copula - Fitting
#'
#' Fit either a C- or D-vine copula model. Function is a repackaging of the \code{CDVineCopSelect} function in the \code{CDVine} package.
#'
#' @param Data Dataframe containing \code{n} at least partially concurrent time series. First column may be a \code{"Date"} object. Can be \code{Dataframe_Combine} output.
#' @param FamilySet Integer vector whuch must include at least one pair-copula family that allows for positive and one that allows for negative dependence. If \code{familyset = NA} (default), selection among all possible families is performed. The coding of pair-copula families is shown below. See help file of the \code{CDVineSim} function to find out copula represented by integers \code{0-40}.
#' @param Type Type of the vine model:\itemize{
#' \item 1 or "CVine" = C-vine
#' \item 2 or "DVine" = D-vine
#' }
#' @param SelCrit Character vector specifying the criterion for choosing among the competing pair-copula. Possible choices: \code{"AIC"} (default) or \code{"BIC"}.
#' @param Indeptest Logical; whether a hypothesis test for the independence of \code{u1} and \code{u2} is performed before bivariate copula selection (default: \code{Indeptest = FALSE}; cp. BiCopIndTest). The independence copula is chosen for a (conditional) pair if the null hypothesis of independence cannot be rejected.
#' @param level Numeric; significance level of the independence test (default: level = 0.05).
#' @return List comprising the pair-copula families composing the C- or D-vine copula \code{Family}, its parameters \code{Par} and \code{Par2} as well as whether it is a C or D-vine \code{Type}.
#' @seealso \code{\link{Detrend_Declustered_Combine}} \code{\link{CDVineCopSelect}}  \code{\link{BiCopSelect}}
#' @export
#' @examples
#' S22.Vine<-Vine_Copula_Fit(Data=S22.Detrend.df, FamilySet=NA, Type="DVine", SelCrit="AIC",Indeptest=FALSE, Level=0.05)
Vine_Copula_Fit<-function(Data, FamilySet=NA, Type="DVine", SelCrit="AIC",Indeptest=FALSE, Level=0.05){
  if(class(Data[,1])=="Date" | class(Data[,1])=="factor"){
  Model <- CDVineCopSelect(pobs(na.omit(Data[,2:ncol(Data)])), familyset=FamilySet, type=Type, selectioncrit=SelCrit,indeptest=Indeptest, level=Level)
  } else {
  Model <- CDVineCopSelect(pobs(na.omit(Data[,2:ncol(Data)])), familyset=FamilySet, type=Type, selectioncrit=SelCrit,indeptest=Indeptest, level=Level)
  }
  res<- list(Family = Model$family, Par = Model$par, Par2 = Model$par2, Type = Type)
  return(res)
}
