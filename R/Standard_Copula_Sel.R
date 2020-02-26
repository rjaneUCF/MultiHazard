#' Selecting best fitting standard (elliptical and Archimeadean) copula
#'
#' Fits five n-dimensional standard copula to a dataset and returns their corresponding AIC values.
#'
#' @param Data Data frame containing n at least partially concurrent time series, detrended if necessary. Time steps must be equally spaced, with missing values assigned \code{NA}. First object may be a \code{"Date"} object. Can be \code{Dataframe_Combine} output.
#' @return Data frame containing copula name in column 1 and associated AIC in column 2.
#' Parameters are estimated using the \code{fitCopula()} function in \code{copula} package using maximum pseudo-likelihood estimator \code{"mpl"}. See \code{\link{fitCopula}} for a more thorough explanation.
#' @export
#' @seealso \code{\link{Dataframe_Combine}} \code{\link{Standard_Copula_Fit}}
#' @examples
#' Standard_Copula_Sel(Data_Detrend=S20.Detrend.df)
Standard_Copula_Sel<-function(Data){

if(class(Data[,1])=="Date" | class(Data[,1])=="factor"){
  Data <- Data[,-1]
}

#Fitting Gaussian
norm.cop <- normalCopula(rep(0.6,(ncol(Data))), dim=(ncol(Data)), dispstr="un")
Gaussian_Copula_Fit <- fitCopula(norm.cop,data=pobs(na.omit(Data[,1:ncol(Data)])))
AIC_Gaussian <- ncol(Data)*2 - logLik(Gaussian_Copula_Fit)[1]

#Fitting student-t
t.cop<-tCopula(rep(0.6,ncol(Data)), dim=ncol(Data), dispstr="un",df=5)
t_Copula_Fit <- fitCopula(t.cop,data=pobs(na.omit(Data[,1:ncol(Data)])))
AIC_tCopula <- (ncol(Data)+1)*2 - logLik(t_Copula_Fit)[1]

#Fitting Gumbel
Gumbel_Copula_Fit <- fitCopula(gumbelCopula(dim=ncol(Data)),pobs(na.omit(Data[,1:ncol(Data)])))
AIC_Gumbel <- 1*2 - logLik(Gumbel_Copula_Fit)[1]

#Fitting Clayton
Clayton_Copula_Fit <- fitCopula(claytonCopula(dim=ncol(Data)),pobs(na.omit(Data[,1:ncol(Data)])))
AIC_Clayton <- 1*2 - logLik(Clayton_Copula_Fit)[1]

#Fitting Frank
Frank_Copula_Fit <- fitCopula(frankCopula(dim=ncol(Data)),pobs(na.omit(Data[,1:ncol(Data)])))
AIC_Frank <- 1*2 - logLik(Frank_Copula_Fit)[1]

res<-data.frame(c("Gaussian","t-cop","Gumbel","Clayton","Frank"),c(AIC_Gaussian,AIC_tCopula,AIC_Gumbel,AIC_Clayton,AIC_Frank))
colnames(res)<-c("Copula","AIC")
return(res)
}



