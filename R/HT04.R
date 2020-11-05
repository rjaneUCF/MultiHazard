#' Fits and simulates from the conditional multivariate approach of Heffernan and Tawn (2004)
#'
#' Fitting and simulating the conditional multivariate approach of Heffernan and Tawn (2004) to a dataset comprising 3 variables. Function utilizes the \code{mexDependence} and \code{predict.mex.conditioned} functions from the \code{texmex} package.
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
#' @param u_Dependence Dependence quantile. Specifies the (sub-sample of) data to which the dependence model is fitted, that for which the conditioning variable exceeds the threshold associated with the prescribed quantile. Default is \code{0.7}, thus the dependence parameters are estimated using the data with the highest \code{30\%} of values of the conditioning variables.
#' @param Margins Character vector specifying the form of margins to which the data are transformed for carrying out dependence estimation. Default is \code{"gumbel"}, alternative is \code{"laplace"}. Under Gumbel margins, the estimated parameters \code{a} and \code{b} describe only positive dependence, while \code{c} and \code{d} describe negative dependence in this case. For Laplace margins, only parameters \code{a} and \code{b} are estimated as these capture both positive and negative dependence.
#' @param V See documentation for mexDependence.
#' @param Maxit See documentation for mexDependence.
#' @return List comprising the fitted HT04 models \code{Models}, proportion of the time each variable is most extreme, given at least one variable is extreme \code{Prop}, as well as the simulated values on the transformed \code{u.sim} and original {x.sim} scales.
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
  
  #Output 'vectors'
  HT04_Model<-vector('list',ncol(data_Detrend_Declustered_df))
  u<-array(NA,dim=c(nrow(na.omit(data_Detrend_Dependence_df)),ncol(data_Detrend_Declustered_df)))
  Gumbel_df<-array(NA,dim=c(nrow(na.omit(data_Detrend_Dependence_df)),ncol(data_Detrend_Declustered_df)))
  Gumbel_df.Threshold<-rep(NA,ncol(data_Detrend_Declustered_df))
  u.extremes<-array(0,dim=c(nrow(na.omit(data_Detrend_Dependence_df)),ncol(data_Detrend_Declustered_df)))
  Prop<-rep(NA,ncol(data_Detrend_Declustered_df))
  HT04.Predict<-vector('list',ncol(data_Detrend_Declustered_df))
  
  #Fitting the model
  for(i in 1:ncol(data_Detrend_Declustered_df)){
    HT04.Dataset<-array(NA,dim=c(nrow(data_Detrend_Declustered_df),ncol(data_Detrend_Declustered_df)))
    HT04.Dataset[,i]<-data_Detrend_Declustered_df[,i]
    for(j in 1:length((1:ncol(data_Detrend_Declustered_df))[-i])){
      HT04.Dataset[,(1:ncol(data_Detrend_Declustered_df))[-i][j]]<-data_Detrend_Dependence_df[,(1:ncol(data_Detrend_Declustered_df))[-i][j]]
    }
    colnames(HT04.Dataset)<-colnames(data_Detrend_Dependence_df)
    Migpd$data<-na.omit(HT04.Dataset)
    HT04_Model[[i]]=mexDependence(Migpd,which=colnames(data_Detrend_Dependence_df)[i],dqu=u_Dependence, margins = "gumbel",
                                  constrain = FALSE,v = V, maxit = Maxit)
    names(HT04_Model) <- colnames(data_Detrend_Dependence_df)
    
    Migpd$data<-na.omit(data_Detrend_Dependence_df)
    u[,i]<-as.numeric(transFun.HT04(x=Migpd$data[,i],mod=Migpd$models[[i]]))
    Gumbel_df[,i]<- -log(-log(u[,i]))
    colnames(Gumbel_df)<-colnames(data_Detrend_Dependence_df)
    Gumbel_df.Threshold[i]<-quantile(Gumbel_df[,i],u_Dependence)
    u.extremes[which(Gumbel_df[,i]>Gumbel_df.Threshold[i]),i]<-1
  }
  
  Gumbel_df_Extremes<-Gumbel_df[which(apply(Gumbel_df_Extremes, 1, sum)>0),]
  colnames(Gumbel_df_Extremes)<-colnames(data_Detrend_Dependence_df)
  
  
  for(i in 1:ncol(data_Detrend_Declustered_df)){
    Prop[i]<-length(which(apply(Gumbel_df_Extremes, 1, function(x) which(x==max(x)))==i))/nrow(Gumbel_df_Extremes)
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
  res<-list("Model" = HT04_Model,"u.sim" = HT04_Predict_Transformed,"x.sim" = HT04_Predict_Simulated)
  return(res)
}