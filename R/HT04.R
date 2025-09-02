#' Fits and simulates from the conditional multivariate approach of Heffernan and Tawn (2004)
#'
#' Fits the conditional multivariate approach of Heffernan and Tawn (2004) to a dataset and simulates realizations from the fitted model. Function utilizes the \code{mexDependence} and \code{predict.mex.conditioned} functions from the \code{texmex} package.
#'
#' @param data_Detrend_Dependence_df  A data frame with (n+1) columns, containing in column \itemize{
#' \item 1 - Continuous sequence of dates spanning the first to the final time of any of the variables are recorded.
#' \item 2:(n+1) - Values, detrended where necessary, of the variables to be modelled.
#' }
#' @param data_Detrend_Declustered_df A data frame with (n+1) columns, containing in column \itemize{
#' \item 1 - Continuous sequence of dates spanning the first to the final time of any of the variables are recorded.
#' \item 2:(n+1) - Declustered and if necessary detrended values of the variables to be modelled.
#' }
#' @param Migpd An \code{Migpd} object, containing the generalized Pareto models fitted (independently) to each of the variables.
#' @param mu Numeric vector of length one specifying the (average) occurrence frequency of events in \code{data_Detrend_Dependence_df}. Default is \code{365.25}, daily data.
#' @param N Numeric vector of length one specifying the number of years worth of extremes to simulate. Default is \code{100} years.
#' @param u_Dependence Dependence quantile. Specifies the (sub-sample of) data to which the dependence model is fitted, that for which the conditioning variable exceeds the threshold associated with the prescribed quantile. Default is \code{0.7}, thus the dependence parameters are estimated using the data with the highest \code{30\%} of values of the conditioning variables.
#' @param Margins Character vector specifying the form of margins to which the data are transformed for carrying out dependence estimation. Default is \code{"gumbel"}, alternative is \code{"laplace"}. Under Gumbel margins, the estimated parameters \code{a} and \code{b} describe only positive dependence, while \code{c} and \code{d} describe negative dependence in this case. For Laplace margins, only parameters \code{a} and \code{b} are estimated as these capture both positive and negative dependence.
#' @param V See documentation for mexDependence.
#' @param Maxit See documentation for mexDependence.
#' @return List comprising the fitted HT04 models \code{Models}, proportion of the time each variable is most extreme, given at least one variable is extreme \code{Prop}, residuals \code{z}, as well as the simulated values on the transformed \code{u.sim} and original \code{x.sim} scales.
#' @seealso \code{\link{Dataframe_Combine}} \code{\link{Migpd_Fit}}
#' @export
#' @examples
#' #Fitting and simulating from the Heffernan and Tawn (2004) model
#' S20.HT04<-HT04(data_Detrend_Dependence_df=S20.Detrend.df,
#'                data_Detrend_Declustered_df=S20.Detrend.Declustered.df,
#'                u_Dependence=0.995,Migpd=S20.Migpd,mu=365.25,N=1000)
#' #View model conditioning on rainfall
#' S20.HT04$Model$Rainfall
#' #Assigning simulations (transformed back to the original scale) a name
#' S20.HT04.Sim<-S20.HT04$x.sim
#' #Plotting observed (black) and simulated (red) values
#' S20.Pairs.Plot.Data<-data.frame(rbind(na.omit(S20.Detrend.df[,-1]),S20.HT04.Sim),
#'                                 c(rep("Observation",nrow(na.omit(S20.Detrend.df))),
#'                                   rep("Simulation",nrow(S20.HT04.Sim))))
#' colnames(S20.Pairs.Plot.Data)<-c(names(S20.Detrend.df)[-1],"Type")
#' pairs(S20.Pairs.Plot.Data[,1:3],
#'       col=ifelse(S20.Pairs.Plot.Data$Type=="Observation","Black","Red"),
#'       upper.panel=NULL,pch=16)
HT04<-function(data_Detrend_Dependence_df,data_Detrend_Declustered_df,u_Dependence,Migpd,mu=365.25,N=100,Margins="gumbel",V=10,Maxit=10000){

  # Input Validation

  # Check if data frames exist and are valid
  if(missing(data_Detrend_Dependence_df) || is.null(data_Detrend_Dependence_df) || !is.data.frame(data_Detrend_Dependence_df)){
    stop("data_Detrend_Dependence_df must be a data frame.")
  }

  if(missing(data_Detrend_Declustered_df) || is.null(data_Detrend_Declustered_df) || !is.data.frame(data_Detrend_Declustered_df)){
    stop("data_Detrend_Declustered_df must be a data frame.")
  }

  # Check if data frames have matching dimensions after removing date/factor columns
  temp_dep <- data_Detrend_Dependence_df
  temp_declust <- data_Detrend_Declustered_df

  if(class(temp_dep[,1])[1] %in% c("Date", "factor", "POSIXct", "character")){
    temp_dep <- temp_dep[,-1]
  }
  if(class(temp_declust[,1])[1] %in% c("Date", "factor", "POSIXct", "character")){
    temp_declust <- temp_declust[,-1]
  }

  if(ncol(temp_dep) != ncol(temp_declust)){
    stop("Data frames must have the same number of numeric columns after removing date/factor columns.")
  }

  if(!all(colnames(temp_dep) == colnames(temp_declust))){
    stop("Column names must match between data frames after removing date/factor columns.")
  }

  # Validate u_Dependence
  if(missing(u_Dependence) || is.null(u_Dependence) || !is.numeric(u_Dependence) || length(u_Dependence) != 1 || u_Dependence <= 0 || u_Dependence >= 1){
    stop("u_Dependence must be a single numeric value between 0 and 1 (exclusive).")
  }

  # Validate Migpd
  if(missing(Migpd) || is.null(Migpd) || !is.list(Migpd)){
    stop("Migpd must be a list object (migpd class).")
  }

  if(is.null(Migpd$models)){
    stop("Migpd$models is NULL. Please ensure Migpd object contains fitted models.")
  }

  if(length(Migpd$models) != ncol(temp_dep)){
    stop("Number of models in Migpd$models (", length(Migpd$models),
         ") must match number of data columns (", ncol(temp_dep), ").")
  }

  # Validate numeric parameters
  if(!is.numeric(mu)){
    stop("mu must be numeric.")
  }

  if(length(mu) != 1 || mu <= 0){
     stop("mu must be a positive numeric value.")
  }

  if(!is.numeric(N)){
    stop("N must be a numeric.")
  }

  if(length(N) != 1 || N <= 0 || N != round(N)){
     stop("mu must be a positive numeric value.")
  }

  if(!is.numeric(V)){
    stop("V must be a numeric.")
  }

  if(length(V) != 1 || V <= 0 || V != round(V)){
    stop("V must be a positive integer.")
  }

  if(!is.numeric(Maxit)){
    stop("Maxit must be a numeric value.")
  }

  if(length(Maxit) != 1 || Maxit <= 0 || Maxit != round(Maxit)){
    stop("Maxit must be a positive integer.")
  }

  # Validate Margins parameter
  valid_margins <- c("gumbel", "laplace", "exponential")
  if(!is.character(Margins) || length(Margins) != 1 || !Margins %in% valid_margins){
    stop("Margins must be one of: ", paste(valid_margins, collapse=", "), ".")
  }

  if(inherits(data_Detrend_Dependence_df[,1], c("Date", "factor", "POSIXct", "character"))){
    data_Detrend_Dependence_df <- data_Detrend_Dependence_df[,-1]
  }

  if(inherits(data_Detrend_Declustered_df[,1], c("Date", "factor", "POSIXct", "character"))){
    data_Detrend_Declustered_df <- data_Detrend_Declustered_df[,-1]
  }

  #Output 'vectors'
  HT04_Model<-vector('list',ncol(data_Detrend_Declustered_df))
  u<-array(NA,dim=c(nrow(na.omit(data_Detrend_Dependence_df)),ncol(data_Detrend_Declustered_df)))
  Gumbel_df<-array(NA,dim=c(nrow(na.omit(data_Detrend_Dependence_df)),ncol(data_Detrend_Declustered_df)))
  Gumbel_df.Threshold<-rep(NA,ncol(data_Detrend_Declustered_df))
  u.extremes<-array(0,dim=c(nrow(na.omit(data_Detrend_Dependence_df)),ncol(data_Detrend_Declustered_df)))
  Prop<-rep(NA,ncol(data_Detrend_Declustered_df))
  HT04.Predict<-vector('list',ncol(data_Detrend_Declustered_df))
  HT04_z<-vector('list',ncol(data_Detrend_Declustered_df))

  #Fitting the model
  for(i in 1:ncol(data_Detrend_Declustered_df)){
    HT04.Dataset<-array(NA,dim=c(nrow(data_Detrend_Declustered_df),ncol(data_Detrend_Declustered_df)))
    HT04.Dataset[,i]<-data_Detrend_Declustered_df[,i]
    for(j in 1:length((1:ncol(data_Detrend_Declustered_df))[-i])){
      HT04.Dataset[,(1:ncol(data_Detrend_Declustered_df))[-i][j]]<-data_Detrend_Dependence_df[,(1:ncol(data_Detrend_Declustered_df))[-i][j]]
    }
    colnames(HT04.Dataset)<-colnames(data_Detrend_Dependence_df)
    Migpd$data<-na.omit(HT04.Dataset)
    HT04_Model[[i]]=mexDependence(Migpd,which=colnames(data_Detrend_Dependence_df)[i],dqu=u_Dependence, margins = "laplace",
                                  constrain = FALSE,v = V, maxit = Maxit)
    names(HT04_Model) <- colnames(data_Detrend_Dependence_df)

    Migpd$data<-na.omit(data_Detrend_Dependence_df)
    u[,i]<-as.numeric(transFun.HT04(x=Migpd$data[,i],mod=Migpd$models[[i]]))
    Gumbel_df[,i]<- u[,i] #-log(-log(u[,i]))
    colnames(Gumbel_df)<-colnames(data_Detrend_Dependence_df)
    Gumbel_df.Threshold[i]<-quantile(Gumbel_df[,i],u_Dependence)
    u.extremes[which(Gumbel_df[,i]>Gumbel_df.Threshold[i]),i]<-1
  }
  Gumbel_df_Extremes<-Gumbel_df[which(apply(u.extremes, 1, sum)>0),]
  colnames(Gumbel_df_Extremes)<-colnames(data_Detrend_Dependence_df)

  for(i in 1:ncol(data_Detrend_Declustered_df)){
    Prop[i]<-length(which(as.numeric(unlist(apply(Gumbel_df_Extremes, 1, function(x) which(x==max(x)))))==i))/nrow(Gumbel_df_Extremes)
  }

  #Simulations
  n<-rep(NA,ncol(data_Detrend_Declustered_df))
  for(i in 1:ncol(data_Detrend_Declustered_df)){
    HT04.Predict[[i]]<-predict.mex.conditioned(HT04_Model[[i]], which=colnames(data_Detrend_Dependence_df)[i], pqu = u_Dependence, nsim = round((1-u_Dependence)*mu*N*Prop[i],0), trace=10)
    n[i]<-nrow(HT04.Predict[[i]]$data$transformed)
  }

  HT04_Predict_Transformed<-HT04.Predict[[1]]$data$transformed
  HT04_Predict_Simulated<-HT04.Predict[[1]]$data$simulated
  for(i in 2:ncol(data_Detrend_Declustered_df)){
    HT04_Predict_Transformed<-rbind(HT04_Predict_Transformed,HT04.Predict[[i]]$data$transformed)
    HT04_Predict_Simulated<-rbind(HT04_Predict_Simulated,HT04.Predict[[i]]$data$simulated)
  }

  for(i in 1:ncol(data_Detrend_Dependence_df)){
    HT04_z[[i]]<-HT04.Predict[[i]]$data$z
  }
  names(HT04_z) <- colnames(data_Detrend_Dependence_df)

  res<-list("Model" = HT04_Model, "Prop" = Prop, "z" = HT04_z, "u.sim" = HT04_Predict_Transformed,"x.sim" = HT04_Predict_Simulated)
  return(res)
}
