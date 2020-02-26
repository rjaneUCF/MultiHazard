#' Fits a single generalised Pareto distribution - Fit
#'
#' Fit a Generalized Pareto Distribution (GPD) to a declstered dataset.
#'
#' @param Data Numeric vector containing the declusted data.
#' @param Data_Full Numeric vector containing the non-declustered data.
#' @param u GPD threshold; as a quantile \code{[0,1]} of \code{Data} vector. Default is \code{0.95}.
#' @param Plot Logical; indicating whether to plot diagnostics. Default is \code{FALSE}.
#' @param xlab_hist Character vector of length one. Histogram x-axis label. Default is \code{"Data"}.
#' @param y_lab Character vector of length one. Histogram x-axis label. Default is \code{"Data"}.
#' @section Details:
#' The fitted GPD model, is following parameterised as follows: \eqn{P( X > x| X > u)}
#' @return List comprising the GPD \code{Threshold}, shape parameter \code{xi} and scale parameters \code{sigma} along with their standard errors \code{sigma.SE} and \code{xi.SE}.
#' @export
#' @examples
#' Decluster(Data=S20_T_MAX_Daily_Completed_Detrend$Detrend)
GPD_Fit<-function(Data,Data_Full,u=0.95,PLOT=FALSE,xlab_hist="Data",y_lab="Data"){
  gpd<-evm(na.omit(Data), th=quantile(Data_Full,u),penalty = "gaussian",priorParameters = list(c(0, 0), matrix(c(100^2, 0, 0, 0.25), nrow = 2)))
  if(PLOT==TRUE){
  GPD_diag_HT04(Data=na.omit(Data),
                Data_Full=Data_Full,
                model=gpd,
                param=c(exp(summary(gpd)$coef[[1]]), summary(gpd)$coef[[2]]),
                thres=quantile(Data_Full,u),
                xlab.hist=xlab_hist,
                y.lab=y_lab)
  }
  res<-list("Threshold"=quantile(Data_Full,u),"sigma" = exp(summary(gpd)$coef[[1]]), "xi" = summary(gpd)$coef[[2]],"sigma.SE" = summary(gpd)$coef[[3]],"xi.SE" = summary(gpd)$coef[[4]])
  return(res)
}

