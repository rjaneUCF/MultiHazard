#' Fits and simulates from the conditional multivariate approach of Heffernan and Tawn (2004)
#'
#' Fitting and simulating the conditional multivariate approach of Heffernan and Tawn (2004) to a dataset. Function utilizes the \code{mexDependence} and \code{predict.mex.conditioned} functions from the \code{texmex} package.
#'
#' @param data_Detrend_Dependence_df  A data frame with (n+1) columns, containing in column \itemize{
#' \item 1 - Continuous sequence of dates spanning the first to the final time of any of the variables are recorded.
#' \item 2:(n+1) - Values, detrended where necessary, of the variables to be modelled.
#' }
#' @param data_Detrend_Declustered_df A data frame with (n+1) columns, containing in column \itemize{
#' \item 1 - Continuous sequence of dates spanning the first to the final time of any of the variables are recorded.
#' \item 2:(n+1) - Declustered and if necessary detrended values of the variables to be modelled.
#' }
#' @param Migpd An \code{Migpd} object, containing the parameterized Pareto models fitted (independently) to each of the variables.
#' @param u_Dependence Dependence quantile. Specifies the (sub-sample of) data to which the dependence model is fitted, that for which the conditioning variable exceeds the threshold associated with the prescribed quantile. Default is \code{0.7}, thus the dependence parameters are estimated using the data with the highest \code{30\%} of values of the conditioning variables.
#' @param Lag Matrix specifying the lags. The no lag i.e. \code{0} lag cases need to be specified. Row n denotes the lags applied to the variable in the nth column of \code{data_Detrend_Dependence_df}. Column n corresponds to the nth largest lag applied to any variable. \code{NA}. Default is \code{matrix(c(0,1,0,NA),nrow=2,byrow = T)}, which corresponds to a lag of 1 being applied to variable in the first column of \code{data_Detrend_Dependence_df} and no lag being applied to the variable in the second column of \code{data_Detrend_Dependence_df}.
#' @param Margins Character vector specifying the form of margins to which the data are transformed for carrying out dependence estimation. Default is \code{"gumbel"}, alternative is \code{"laplace"}. Under Gumbel margins, the estimated parameters \code{a} and \code{b} describe only positive dependence, while \code{c} and \code{d} describe negative dependence in this case. For Laplace margins, only parameters \code{a} and \code{b} are estimated as these capture both positive and negative dependence.
#' @param V See documentation for mexDependence.
#' @param Maxit See documentation for mexDependence.
#' @return List comprising the fitted HT04 models \code{Models}, proportion of the time each variable is most extreme, given at least one variable is extreme \code{Prop}, as well as the simulated values on the transformed \code{u.sim} and original {x.sim} scales.
#' @seealso \code{\link{Dataframe_Combine}} \code{\link{Decluster}} \code{\link{GPD_Fit}} \code{\link{Migpd_Fit}}
#' @export
#' @examples
#' HT04(data_Detrend_Dependence_df = S22.Detrend.df,data_Detrend_Declustered_df = S22.Detrend.Declustered.df ,Migpd = S22_GPD, u_Dependence=0.7,Margins = "gumbel")
HT04_Lag<-function (data_Detrend_Dependence_df, data_Detrend_Declustered_df, Lags, u_Dependence, Migpd, mu = 365.25, N = 100, Margins = "gumbel",V = 10, Maxit = 10000){
  
  HT04_Model<-vector('list',ncol(data_Detrend_Declustered_df))
  u<-array(NA,dim=c(nrow(na.omit(data_Detrend_Dependence_df)),ncol(data_Detrend_Declustered_df)))
  Gumbel_df<-array(NA,dim=c(nrow(na.omit(data_Detrend_Dependence_df)),ncol(data_Detrend_Declustered_df)))
  Gumbel_df.Threshold<-rep(NA,ncol(data_Detrend_Declustered_df))
  u.extremes<-array(0,dim=c(nrow(na.omit(data_Detrend_Dependence_df)),ncol(data_Detrend_Declustered_df)))
  Prop<-rep(NA,ncol(data_Detrend_Declustered_df))
  HT04.Predict<-vector('list',ncol(data_Detrend_Declustered_df))
  HT04_z<-vector('list',ncol(data_Detrend_Declustered_df))
  
  Lag <- function(x, k) {
    if (k > 0) {
      return(c(rep(NA, k), x)[1:length(x)])
    }
    else {
      return(c(x[(-k + 1):length(x)], rep(NA, -k)))
    }
  }
  
  if (class(data_Detrend_Declustered_df[, 1]) == "Date" | 
      class(data_Detrend_Declustered_df[, 1]) == "factor") {
    data_Detrend_Declustered_df <- data_Detrend_Declustered_df[,-1]
  }
  if (class(data_Detrend_Dependence_df[, 1]) == "Date" | 
      class(data_Detrend_Dependence_df[, 1]) == "factor") {
    data_Detrend_Dependence_df <- data_Detrend_Dependence_df[,-1]
  }

  DF<-array(NA, dim=c(nrow(data_Detrend_Dependence_df),length(which(is.na(Lags)==F))))
  colnames_DF<-rep(NA,length(which(is.na(Lags)==F)))
  colnames_DF_1<-rep(NA,length(which(is.na(Lags)==F)))
  
  L<-0
  for(i in 1:ncol(data_Detrend_Dependence_df)){
    for (j in 1:(length(which(is.na(Lags[i,])==F)))){
     L<-L+1
     lags<-na.omit(Lags[i,])
     if(lags[j]==0){
       DF[, L] <- data_Detrend_Dependence_df[, i]
     } else{
       DF[, L] <- Lag(data_Detrend_Dependence_df[, i], lags[j]) 
     }
     colnames_DF[L]<-ifelse(lags[j]==0,names(data_Detrend_Dependence_df)[i],
                                         paste(colnames(data_Detrend_Dependence_df)[i],"_",lags[j],sep=""))
     colnames_DF_1[L]<-names(data_Detrend_Dependence_df)[i]
   }
  }
  colnames(DF)<-colnames_DF

  for(i in 1:ncol(data_Detrend_Dependence_df)){
    HT04.Dataset<-array(NA,dim=c(nrow(data_Detrend_Declustered_df),length((1:ncol(DF))[-which(colnames_DF_1==colnames(data_Detrend_Declustered_df)[i])])+1))
    HT04.Dataset[,1]<-data_Detrend_Declustered_df[,i]
    Migpd_1<-Migpd
    Migpd_1$models[1]<-Migpd_1$models[i]
    for(k in 1:(length((1:ncol(DF))[-which(colnames_DF_1==colnames(data_Detrend_Declustered_df)[i])]))){
      HT04.Dataset[,k+1]<-DF[,(1:ncol(DF))[-which(colnames_DF_1==colnames(data_Detrend_Declustered_df)[i])][k]]
      Migpd_1$models[k+1]<-Migpd$models[which(colnames(data_Detrend_Declustered_df)==colnames_DF_1[-which(colnames_DF_1==colnames(data_Detrend_Declustered_df)[i])][k])]
      Migpd_1$mth[k+1]<-Migpd$mth[which(colnames(data_Detrend_Declustered_df)==colnames_DF_1[-which(colnames_DF_1==colnames(data_Detrend_Declustered_df)[i])][k])]
      Migpd_1$mqu[k+1]<-Migpd$mqu[which(colnames(data_Detrend_Declustered_df)==colnames_DF_1[-which(colnames_DF_1==colnames(data_Detrend_Declustered_df)[i])][k])]
    }
    colnames(HT04.Dataset)<-c(colnames(data_Detrend_Declustered_df)[i],colnames_DF[-which(colnames_DF_1==colnames(data_Detrend_Declustered_df)[i])])
 
       names(Migpd_1$models)<-c(colnames(data_Detrend_Declustered_df)[i],colnames(DF)[-which(colnames_DF_1==colnames(data_Detrend_Declustered_df)[i])])
    names(Migpd_1$mth)<-c(colnames(data_Detrend_Declustered_df)[i],colnames(DF)[-which(colnames_DF_1==colnames(data_Detrend_Declustered_df)[i])])
    names(Migpd_1$mqu)<-c(colnames(data_Detrend_Declustered_df)[i],colnames(DF)[-which(colnames_DF_1==colnames(data_Detrend_Declustered_df)[i])])

    Migpd_1$data<-na.omit(HT04.Dataset)
    HT04_Model[[i]] = mexDependence(Migpd_1, which = colnames(data_Detrend_Dependence_df)[i],dqu = u_Dependence, margins = Margins, constrain = FALSE, v = V, maxit = Maxit)
  
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

  HT04_Predict_Transformed<-array(NA,dim=c(nrow(HT04.Predict[[1]]$data$transformed),ncol(data_Detrend_Declustered_df)))
  HT04_Predict_Simulated<-array(NA,dim=c(nrow(HT04.Predict[[1]]$data$transformed),ncol(data_Detrend_Declustered_df)))
  
  HT04_Predict_Transformed[,1]<-HT04.Predict[[1]]$data$transformed[,1]
  HT04_Predict_Simulated[,1]<-HT04.Predict[[1]]$data$simulated[,1]
    for(j in 2:ncol(data_Detrend_Declustered_df)){
      L<-length(na.omit(Lags[j,])==FALSE)
      if(j==2){
       if(L>1){
         HT04_Predict_Transformed[,j]<-apply(HT04.Predict[[1]]$data$transformed[,2:(2+(length(which(is.na(Lags[j,])==F))-1))], 1, max)
         HT04_Predict_Simulated[,j]<-apply(HT04.Predict[[1]]$data$simulated[,2:(2+(length(which(is.na(Lags[j,])==F))-1))], 1, max)
        } else{
         HT04_Predict_Transformed[,j]<-HT04.Predict[[1]]$data$transformed[,which(colnames(HT04.Predict[[1]]$data$transformed)==paste(colnames(data_Detrend_Dependence_df)[j],".trans",sep=""))]
         HT04_Predict_Simulated[,j]<-HT04.Predict[[1]]$data$simulated[,which(colnames(HT04.Predict[[1]]$data$simulated)==paste(colnames(data_Detrend_Dependence_df)[j]))]
        }
      } else{
       if(L>1){
         HT04_Predict_Transformed[,j]<-apply(HT04.Predict[[1]]$data$transformed[,(length(which(is.na(Lags[2:(j-1),])==F))+2):(length(which(is.na(Lags[2:(j),])==F))+1)], 1, max)
         HT04_Predict_Simulated[,j]<-apply(HT04.Predict[[1]]$data$simulated[,(length(which(is.na(Lags[2:(j-1),])==F))+2):(length(which(is.na(Lags[2:(j),])==F))+1)], 1, max)
        } else{
         HT04_Predict_Transformed[,j]<-HT04.Predict[[1]]$data$transformed[,which(colnames(HT04.Predict[[1]]$data$transformed)==paste(colnames(data_Detrend_Dependence_df)[j],".trans",sep=""))]
         HT04_Predict_Simulated[,j]<-HT04.Predict[[1]]$data$simulated[,which(colnames(HT04.Predict[[1]]$data$simulated)==paste(colnames(data_Detrend_Dependence_df)[j]))]
        }
      }
    }
  colnames(HT04_Predict_Transformed)<-colnames(data_Detrend_Dependence_df)
  colnames(HT04_Predict_Simulated)<-colnames(data_Detrend_Dependence_df)

  for(i in 2:ncol(data_Detrend_Declustered_df)){
   HT04_Predict_Transformed_1<-array(NA,dim=c(nrow(HT04.Predict[[i]]$data$transformed),ncol(data_Detrend_Declustered_df)))
   HT04_Predict_Simulated_1<-array(NA,dim=c(nrow(HT04.Predict[[i]]$data$transformed),ncol(data_Detrend_Declustered_df)))
   HT04_Predict_Transformed_1[,1]<-HT04.Predict[[i]]$data$transformed[,i]
   HT04_Predict_Simulated_1[,1]<-HT04.Predict[[i]]$data$simulated[,i]
  for(j in 2:ncol(data_Detrend_Declustered_df)){
  L<-length(na.omit(Lags[j,])==FALSE)
  if(j==2){
    if(L>1){
      HT04_Predict_Transformed_1[,j]<-apply(HT04.Predict[[i]]$data$transformed[,2:(2+(length(which(is.na(Lags[j,])==F))-1))], 1, max)
      HT04_Predict_Simulated_1[,j]<-apply(HT04.Predict[[i]]$data$simulated[,2:(2+(length(which(is.na(Lags[j,])==F))-1))], 1, max)
    } else{
      HT04_Predict_Transformed_1[,j]<-HT04.Predict[[i]]$data$transformed[,which(colnames(HT04.Predict[[i]]$data$transformed)==paste(colnames(data_Detrend_Dependence_df)[j],".trans",sep=""))]
      HT04_Predict_Simulated_1[,j]<-HT04.Predict[[i]]$data$simulated[,which(colnames(HT04.Predict[[i]]$data$simulated)==paste(colnames(data_Detrend_Dependence_df)[j]))]
    }
  } else{
    if(L>1){
      HT04_Predict_Transformed_1[,j]<-apply(HT04.Predict[[i]]$data$transformed[,(length(which(is.na(Lags[2:(j-1),])==F))+2):(length(which(is.na(Lags[2:(j),])==F))+1)], 1, max)
      HT04_Predict_Simulated_1[,j]<-apply(HT04.Predict[[i]]$data$simulated[,(length(which(is.na(Lags[2:(j-1),])==F))+2):(length(which(is.na(Lags[2:(j),])==F))+1)], 1, max)
    } else{
      HT04_Predict_Transformed_1[,j]<-HT04.Predict[[i]]$data$transformed[,which(colnames(HT04.Predict[[i]]$data$transformed)==paste(colnames(data_Detrend_Dependence_df)[j],".trans",sep=""))]
      HT04_Predict_Simulated_1[,j]<-HT04.Predict[[i]]$data$simulated[,which(colnames(HT04.Predict[[i]]$data$simulated)==paste(colnames(data_Detrend_Dependence_df)[j]))]
    }
  }
  }
  colnames(HT04_Predict_Transformed_1)<-c(colnames(data_Detrend_Dependence_df)[i],colnames(data_Detrend_Dependence_df)[-i])
  colnames(HT04_Predict_Simulated_1)<-c(colnames(data_Detrend_Dependence_df)[i],colnames(data_Detrend_Dependence_df)[-i])
  HT04_Predict_Transformed<-rbind(HT04_Predict_Transformed,HT04_Predict_Transformed_1)
  HT04_Predict_Simulated<-rbind(HT04_Predict_Simulated,HT04_Predict_Simulated_1)
  }
  
  for(i in 1:ncol(data_Detrend_Dependence_df)){
   HT04_z[[i]]<-HT04.Predict[[i]]$data$z
  }
  
  res<-list("Model" = HT04_Model,"z" = HT04_z,"u.sim" = HT04_Predict_Transformed,"x.sim" = HT04_Predict_Simulated)
  return(res)
}

