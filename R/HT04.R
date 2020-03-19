#' Fits and simulates from the conditional multivariate approach of Heffernan and Tawn (2004)
#'
#' Fitting and simulating the conditional multivariate approach of Heffernan and Tawn (2004) to a dataset comprising 3 variables. Function utilizes the \code{mexDependence} and \code{} functions from the \code{texmex} package.
#'
#' @param data_Detrend_Dependence_df  A dataframe with (n+1) columns, containing in column \itemize{
#' \item 1 - Continuous sequence of dates spanning the first to the final time of any of the variables are recorded.
#' \item 2:(n+1) - Values, detrended where necessary, of the variables to be modelled.
#' }
#' @param data_Detrend_Declustered_df A dataframe with (n+1) columns, containing in column \itemize{
#' \item 1 - Continuous sequence of dates spanning the first to the final time of any of the variables are recorded.
#' \item 2:(n+1) - Declustered and if necessary detrended values of the variables to be modelled.
#' }
#' @param Migpd An \code{Migpd} object, containing the generalised Pareto models fitted (independently) to each of the variables.
#' @param u_Dependence Dependence quantile. Specifies the (sub-sample of) data to which the dependence model is fitted, that for which the conditioning variable exceeds the threshold associated with the presecribed quantile. Default is \code{0.7}, thus the dependence parameters are estimated using the data with the highest \code{30\%} of values of the conditioning variables.
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
#' #Assigning simulations (transfomed back to the origional scale) a name
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

  if(class(data_Detrend_Declustered_df[,1])=="Date" | class(data_Detrend_Declustered_df[,1])=="factor"){
    data_Detrend_Declustered_df<-data_Detrend_Declustered_df[,-1]
  }
  if(class(data_Detrend_Dependence_df[,1])=="Date" | class(data_Detrend_Dependence_df[,1])=="factor"){
    data_Detrend_Dependence_df<-data_Detrend_Dependence_df[,-1]
  }

  if(ncol(data_Detrend_Dependence_df) == 3){
    HT04.Rainfall.Dataset<-data.frame(data_Detrend_Declustered_df[,1],data_Detrend_Dependence_df[,2],data_Detrend_Dependence_df[,3])
    HT04.OsWL.Dataset<-data.frame(data_Detrend_Dependence_df[,1],data_Detrend_Declustered_df[,2],data_Detrend_Dependence_df[,3])
    HT04.Groundwater.Dataset<-data.frame(data_Detrend_Dependence_df[,1],data_Detrend_Dependence_df[,2],data_Detrend_Declustered_df[,3])
    colnames(HT04.Rainfall.Dataset)<-colnames(data_Detrend_Dependence_df)
    colnames(HT04.OsWL.Dataset)<-colnames(data_Detrend_Dependence_df)
    colnames(HT04.Groundwater.Dataset)<-colnames(data_Detrend_Dependence_df)

    Migpd$data<-na.omit(HT04.Rainfall.Dataset)

    HT04.1=mexDependence(Migpd,which=colnames(data_Detrend_Dependence_df)[1],dqu=u_Dependence, margins = Margins,
                         constrain = FALSE,v = V, maxit = Maxit)


    Migpd$data<-na.omit(HT04.OsWL.Dataset)
    HT04.2=mexDependence(Migpd,which=colnames(data_Detrend_Dependence_df)[2],dqu=u_Dependence, margins = Margins,
                         constrain = FALSE,v = V, maxit = Maxit)

    Migpd$data<-na.omit(HT04.Groundwater.Dataset)
    HT04.3=mexDependence(Migpd,which=colnames(data_Detrend_Dependence_df)[3],dqu=u_Dependence, margins = Margins,
                         constrain = FALSE,v = V, maxit = Maxit)


    ###Proportion of time each variable is most extreme
    Migpd$data<-na.omit(data_Detrend_Dependence_df)

    u.1<-transFun.HT04(x=Migpd$data[,1],mod=Migpd$models[[1]])
    Gumbel.1<- -log(-log(u.1))
    Gumbel.1.Threshold<-quantile(Gumbel.1,u_Dependence)

    u.2<-transFun.HT04(x=Migpd$data[,2],mod=Migpd$models[[2]])
    Gumbel.2<- -log(-log(u.2))
    Gumbel.2.Threshold<-quantile(Gumbel.2,u_Dependence)

    u.3<-transFun.HT04(x=Migpd$data[,3],mod=Migpd$models[[3]])
    Gumbel.3<- -log(-log(u.3))
    Gumbel.3.Threshold<-quantile(Gumbel.3,u_Dependence)

    Gumbel_df<-data.frame(Gumbel.1,Gumbel.2,Gumbel.3)
    colnames(Gumbel_df)<-colnames(data_Detrend_Dependence_df)

    u<-which(Gumbel.1>Gumbel.1.Threshold | Gumbel.2>Gumbel.2.Threshold | Gumbel.3>Gumbel.3.Threshold)
    Gumbel_df_Extremes<-Gumbel_df[u,]
    colnames(Gumbel_df_Extremes)<-colnames(data_Detrend_Dependence_df)

    Prop.1<-length(which((Gumbel_df_Extremes[,1] >= Gumbel_df_Extremes[,2]) & (Gumbel_df_Extremes[,1] >= Gumbel_df_Extremes[,3])))/nrow(Gumbel_df_Extremes)
    Prop.2<-length(which((Gumbel_df_Extremes[,2] >= Gumbel_df_Extremes[,1]) & (Gumbel_df_Extremes[,2] >= Gumbel_df_Extremes[,3])))/nrow(Gumbel_df_Extremes)
    Prop.3<-length(which((Gumbel_df_Extremes[,3] >= Gumbel_df_Extremes[,1]) & (Gumbel_df_Extremes[,3] >= Gumbel_df_Extremes[,2])))/nrow(Gumbel_df_Extremes)

    HT04.Predict.1<-predict.mex.conditioned(HT04.1, which=colnames(data_Detrend_Dependence_df)[1], pqu = u_Dependence, nsim = round((1-u_Dependence)*mu*N*Prop.1,0), trace=10)

    HT04.Predict.2<-predict.mex.conditioned(HT04.2, which=colnames(data_Detrend_Dependence_df)[2], pqu = u_Dependence, nsim =round((1-u_Dependence)*mu*N*Prop.2,0), trace=10)

    HT04.Predict.3<-predict.mex.conditioned(HT04.3, which=colnames(data_Detrend_Dependence_df)[3],pqu = u_Dependence, nsim = round((1-u_Dependence)*mu*N*Prop.3,0), trace=10)


    HT04_Predict_Transformed<-cbind(c(HT04.Predict.1$data$transformed[,1],HT04.Predict.2$data$transformed[,2],HT04.Predict.3$data$transformed[,2]),
                                    c(HT04.Predict.1$data$transformed[,2],HT04.Predict.2$data$transformed[,1],HT04.Predict.3$data$transformed[,3]),
                                    c(HT04.Predict.1$data$transformed[,3],HT04.Predict.2$data$transformed[,3],HT04.Predict.3$data$transformed[,1]))
    colnames(HT04_Predict_Transformed)<-colnames(data_Detrend_Dependence_df)
    HT04_Predict_Simulated<-cbind(c(HT04.Predict.1$data$simulated[,1],HT04.Predict.2$data$simulated[,2],HT04.Predict.3$data$simulated[,2]),
                                  c(HT04.Predict.1$data$simulated[,2],HT04.Predict.2$data$simulated[,1],HT04.Predict.3$data$simulated[,3]),
                                  c(HT04.Predict.1$data$simulated[,3],HT04.Predict.2$data$simulated[,3],HT04.Predict.3$data$simulated[,1]))
    colnames(HT04_Predict_Simulated)<-colnames(data_Detrend_Dependence_df)
    Model<-list(HT04.1,HT04.2,HT04.3)
    names(Model) <- colnames(data_Detrend_Dependence_df)
  }

  if(ncol(data_Detrend_Dependence_df) == 2){
    HT04.Rainfall.Dataset<-data.frame(data_Detrend_Declustered_df[,1],data_Detrend_Dependence_df[,2])
    HT04.OsWL.Dataset<-data.frame(data_Detrend_Dependence_df[,1],data_Detrend_Declustered_df[,2])

    colnames(HT04.Rainfall.Dataset)<-colnames(data_Detrend_Dependence_df)
    colnames(HT04.OsWL.Dataset)<-colnames(data_Detrend_Dependence_df)

    Migpd$data<-na.omit(HT04.Rainfall.Dataset)

    HT04.1=mexDependence(Migpd,which=colnames(data_Detrend_Dependence_df)[1],dqu=u_Dependence, margins = Margins,
                         constrain = FALSE,v = V, maxit = Maxit)


    Migpd$data<-na.omit(HT04.OsWL.Dataset)
    HT04.2=mexDependence(Migpd,which=colnames(data_Detrend_Dependence_df)[2],dqu=u_Dependence, margins = Margins,
                         constrain = FALSE,v = V, maxit = Maxit)

    ###Proportion of time each variable is most extreme
    Migpd$data<-na.omit(data_Detrend_Dependence_df)

    u.1<-transFun.HT04(x=Migpd$data[,1],mod=Migpd$models[[1]])
    Gumbel.1<- -log(-log(u.1))
    Gumbel.1.Threshold<-quantile(Gumbel.1,u_Dependence)

    u.2<-transFun.HT04(x=Migpd$data[,2],mod=Migpd$models[[2]])
    Gumbel.2<- -log(-log(u.2))
    Gumbel.2.Threshold<-quantile(Gumbel.2,u_Dependence)

    Gumbel_df<-data.frame(Gumbel.1,Gumbel.2)
    colnames(Gumbel_df)<-colnames(data_Detrend_Dependence_df)

    u<-which(Gumbel.1>Gumbel.1.Threshold | Gumbel.2>Gumbel.2.Threshold)
    Gumbel_df_Extremes<-Gumbel_df[u,]
    colnames(Gumbel_df_Extremes)<-colnames(data_Detrend_Dependence_df)

    Prop.1<-length(which((Gumbel_df_Extremes[,1] >= Gumbel_df_Extremes[,2])))/nrow(Gumbel_df_Extremes)
    Prop.2<-length(which((Gumbel_df_Extremes[,2] >= Gumbel_df_Extremes[,1])))/nrow(Gumbel_df_Extremes)

    HT04.Predict.1<-predict.mex.conditioned(HT04.1, which=colnames(data_Detrend_Dependence_df)[1], pqu = u_Dependence, nsim = round((1-u_Dependence)*mu*N*Prop.1,0), trace=10)

    HT04.Predict.2<-predict.mex.conditioned(HT04.2, which=colnames(data_Detrend_Dependence_df)[2], pqu = u_Dependence, nsim =round((1-u_Dependence)*mu*N*Prop.2,0), trace=10)


    HT04_Predict_Transformed<-cbind(c(HT04.Predict.1$data$transformed[,1],HT04.Predict.2$data$transformed[,2]),
                                    c(HT04.Predict.1$data$transformed[,2],HT04.Predict.2$data$transformed[,1]))

    colnames(HT04_Predict_Transformed)<-colnames(data_Detrend_Dependence_df)
    HT04_Predict_Simulated<-cbind(c(HT04.Predict.1$data$simulated[,1],HT04.Predict.2$data$simulated[,2]),
                                  c(HT04.Predict.1$data$simulated[,2],HT04.Predict.2$data$simulated[,1]))
    colnames(HT04_Predict_Simulated)<-colnames(data_Detrend_Dependence_df)
    Model<-list(HT04.1,HT04.2)
    names(Model) <- colnames(data_Detrend_Dependence_df)
  }

  res<-list("Model" = Model,"u.sim" = HT04_Predict_Transformed,"x.sim" = HT04_Predict_Simulated)
  return(res)
}
