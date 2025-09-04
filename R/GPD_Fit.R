#' Fits a single generalized Pareto distribution - Fit
#'
#' Fit a Generalized Pareto Distribution (GPD) to a declustered dataset.
#'
#' @param Data Numeric vector containing the declusted data.
#' @param Data_Full Numeric vector containing the non-declustered data.
#' @param u GPD threshold expressed as a quantile \code{[0,1]} of \code{Data} vector. Default is \code{0.95}.
#' @param Thres GPD threshold expressed on the original scale of the \code{"Data"}. Only one of \code{u} and \code{Thres} should be supplied. Default is \code{NA}.
#' @param mu Numeric vector of length one specifying (average) occurrence frequency of events in the \code{Data_Full} input. Default is \code{365.25}.
#' @param GPD_Bayes Logical; indicating whether to use a Bayesian approach to estimate GPD parameters. This involves applying a penalty to the likelihood to aid in the stability of the optimization procedure. Default is \code{TRUE}.
#' @param Method Character vector of length one specifying the method of choosing the threshold. \code{"Standard"} (default) chooses the exact threshold specified as either \code{"u"} or \code{"th"}, whereas \code{"Solari"} selects the minimum exceedence of the \code{"Data"} above the user-specified threshold.
#' @param min.RI Numeric vector of length one specifying the minimum return period in the return level plot. Default is \code{1}.
#' @param max.RI Numeric vector of length one specifying the maximum return period in the return level plot. Default is \code{100}.
#' @param PLOT Logical; indicating whether to plot diagnostics. Default is \code{FALSE}.
#' @param xlab_hist Character vector of length one. Histogram x-axis label. Default is \code{"Data"}.
#' @param y_lab Character vector of length one. Histogram x-axis label. Default is \code{"Data"}.
#' @section Details:
#' For excesses of a variable X over a suitably high threshold u the fitted GPD model is parameterized as follows: \deqn{P( X > x| X > u) = \left[1 + \xi \frac{(x-u)}{\sigma}\right]^{-\frac{1}{\xi}}_{+}}
#' where \eqn{\xi} and \eqn{\sigma>0} are the shape and scale parameters  of the GPD and \eqn{[y]_{+}=max(y,0)}.
#' @return List comprising the GPD \code{Threshold}, shape parameter \code{xi} and scale parameters \code{sigma} along with their standard errors \code{sigma.SE} and \code{xi.SE}.
#' @export
#' @examples
#' #Decluster time series
#' S20T_decl = Decluster(Data=S20_T_MAX_Daily_Completed_Detrend_Declustered$Detrend)
#' #Fit GPD
#' GPD_Fit(Data=S20T_decl$Declustered,
#'         Data_Full=S20_T_MAX_Daily_Completed_Detrend_Declustered$Detrend)
GPD_Fit<-function(Data,Data_Full,u=0.95,Thres=NA,mu=365.25,GPD_Bayes=TRUE,Method="Standard",min.RI=1,max.RI=100,PLOT=FALSE,xlab_hist="Data",y_lab="Data"){

  #Checking inputs are valid

  #Data
  if (missing(Data) || !is.numeric(Data)) stop("Data must be a numeric vector.")

  #Data_Full
  if (missing(Data_Full) || !is.numeric(Data_Full)) stop("Data_Full must be a numeric vector.")

  #Quantile threshold
  if (!is.na(u) && (!is.numeric(u) || u <= 0 || u >= 1)) {
    stop("u must be a numeric value strictly between 0 and 1.")
  }

  #Threshold on original scale
  Data_Full<-na.omit(Data_Full)
  if (is.na(Thres)) {
    Thres <- ifelse(is.na(u), Thres, quantile(Data_Full, u))
  }
  if (is.na(Thres) || !is.numeric(Thres)) stop("Thres must be numeric or computable from u.")
  if (all(Data < Thres)) stop("No values in Data exceed the threshold. Cannot fit GPD.")

  #mu
  if (!is.numeric(mu) || mu <= 0) stop("mu must be a positive number.")

  #Method
  if (!Method %in% c("Standard", "Solari")) stop("Method must be Standard or Solari.")
  if (!is.logical(GPD_Bayes) || length(GPD_Bayes) != 1) stop("GPD_Bayes must be TRUE or FALSE.")

  #Plot?
  if (!is.logical(PLOT) || length(PLOT) != 1) stop("PLOT must be TRUE or FALSE.")

  #Plot and estimation limits
  if (!is.numeric(min.RI) || !is.numeric(max.RI) || min.RI <= 0 || max.RI <= 0) stop("min.RI and max.RI must be positive.")
  if (min.RI >= max.RI) stop("min.RI must be less than max.RI.")
  if (!is.character(xlab_hist) || !is.character(y_lab)) stop("xlab_hist and y_lab must be character strings.")

  Data<-na.omit(Data)

  if(Method=="Standard"){
    if(GPD_Bayes==T){
      gpd<-evm(Data, th = Thres,penalty = "gaussian",priorParameters = list(c(0, 0), matrix(c(100^2, 0, 0, 0.25), nrow = 2)))
    } else{
      gpd<-evm(Data, th = Thres)
    }
    gpd$rate<-length(Data[which(Data>=Thres)])/(length(Data_Full)/mu)
  }

  if(Method=="Solari"){
    Exceedence<-Data[which(Data>=Thres)]
    if(GPD_Bayes==T){
      gpd <- evm(Exceedence, th = min(Exceedence), penalty = "gaussian",
                 priorParameters = list(c(0, 0), matrix(c(100^2, 0, 0, 0.25), nrow = 2)))
    } else{
      gpd <- evm(Exceedence, th = min(Exceedence))
    }
    Thres=min(Exceedence)
    gpd$rate<-length(Exceedence)/(length(Data_Full)/mu)
  }

  if(PLOT==TRUE){
    GPD_diag_HT04(Data=na.omit(Data),
                  Data_Full=Data_Full,
                  model=gpd,
                  param=c(exp(gpd$par[1]), gpd$par[2]),
                  thres=Thres,
                  min.RI=min.RI,
                  max.RI=max.RI,
                  xlab.hist=xlab_hist,
                  y.lab=y_lab)
  }
  res<-list("Threshold"=as.numeric(Thres),
            "Rate"=as.numeric(gpd$rate),
            "sigma" = as.numeric(exp(gpd$par[1])),
            "xi" = as.numeric(gpd$par[2]),
            "sigma.SE" = as.numeric(gpd$se[1]),
            "xi.SE" = as.numeric(gpd$se[2]))
  return(res)
}
