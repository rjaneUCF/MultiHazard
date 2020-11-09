#' Copula Selection With Threshold 2D - Fit
#'
#' Declustered excesses of a (conditioning) variable are paired with co-occurences of the other variable before the best fitting bivariate copula is selected, using \code{BiCopSelect} function in the \code{VineCopula} package, for a single or range of thresholds. The procedure is automatically repeated with the variables switched.
#'
#' @param Data_Detrend Data frame containing two at least partially concurrent time series, detrended if necessary. Time steps must be equally spaced, with missing values assigned \code{NA}.
#' @param Data_Declust Data frame containing two (independently) declustered at least partially concurrent time series. Time steps must be equally spaced, with missing values assigned \code{NA}.
#' @param Thres A single or sequence of thresholds, given as a quantile of the observations of the conditioning variable. Default, sequence from \code{0.9} to \code{0.99} at intervals of \code{0.01}.
#' @param Lag_Backward_Var1 Numeric vector of length one specifying the negative lag applied to variable in the first column of \code{Data_Detrend}. Default \code{0}.
#' @param Lag_Forward_Var1 Numeric vector of length one specifying positive lag applied to variable in the first column of \code{Data_Detrend}. Default \code{0}.
#' @param Lag_Backward_Var2 Numeric vector of length one specifying negative lag applied to variable in the second column of \code{Data_Detrend}. Default \code{0}.
#' @param Lag_Forward_Var2 Numeric vector of length one specifying  positive lag applied to variable in the second column of \code{Data_Detrend}. Default \code{0}.
#' @param x_lim_min Numeric vector of length one specifying x-axis minimum. Default is the maximum argument in \code{Thres}.
#' @param x_lim_max Numeric vector of length one specifying x-axis maximum. Default is the minimum argument in \code{Thres}.
#' @param y_lim_min Numeric vector of length one specifying y-axis minimum. Default \code{-1.0}.
#' @param y_lim_max Numeric vector of length one specifying y-axis maximum. Default \code{1.0}.
#' @param Upper Numeric vector specifying the element number of the \code{Thres} argument for which the copula family name label to appear above the corresponding point on the Kendall's tau coefficient vs threshold plot, when conditioning on the variable in column 1. Default is \code{0}.
#' @param Lower Numeric vector specifying the element number of the \code{Thres} argument for which the copula family name label to appear below the corresponding point on the Kendall's tau coefficient vs threshold plot, when conditioning on the variable in column 2. Default is \code{0}.
#' @param GAP Numeric vector of length one specifying the distance above or below the copula family name label appears the corresponding point on the Kendall's tau coefficient vs threshold plot. Default is \code{0.05}.
#' @param Legend Logic vector of length one specifying whether a legend should be plotted. Default is \code{TRUE}.
#' @return List comprising: \itemize{
#' \item \code{Kendalls_Tau_Var1}
#' Kendall's tau of a sample
#' \item \code{p_value_Var1}
#' p-value when testing the null hypothesis \code{H_0=0} i.e. that there is no correlation between the variables
#' \item \code{N_Var1}
#' size of the dataset
#' \item \code{Copula_Family_Var1}
#' best fitting copula for the specified thresholds
#' }
#' when the dataset is conditioned on the variable in column 1.
#' Analogous vector \code{Kendalls_Tau_Var2}, \code{p_value_Var2}, \code{N_Var2} and \code{Copula_Family_Var2} for the specified thresholds when the dataset is conditioned on the variable in column 2.
#' @seealso \code{\link{Dataframe_Combine}}
#' @export
#' @examples
#' Copula_Threshold_2D(Data_Detrend=S20.Detrend.df[,-c(1,4)],
#'                     Data_Declust=S20.Detrend.Declustered.df[,-c(1,4)],
#'                     y_lim_min=-0.075, y_lim_max =0.25,
#'                     Upper=c(6,8), Lower=c(6,8),GAP=0.1)
Copula_Threshold_2D_Lag<-function(Data_Detrend,Data_Declust,Thres1=seq(0.9,0.99,0.01),Thres2=seq(0.9,0.99,0.01),Lag_Backward_Var1,Lag_Forward_Var1,Lag_Backward_Var2,Lag_Forward_Var2,x_lim_min=min(c(Thres1,Thres2)),x_lim_max=max(c(Thres1,Thres2)),y_lim_min=-1,y_lim_max=1,Upper=0,Lower=0,GAP=0.05,Legend=TRUE){
  y_lim=y_lim_max-y_lim_min
  copula_table<-data.frame(c(seq(0,40,1)[-(c(11,12,15,21,22,25,31,32,35)+1)],104,114,124,134,204,214,224,234),c("Ind.","Gaussian", "t-copula", "Clayton", "Gumbel","Frank","Joe","BB1","BB6","BB7","BB8","Sur. Clayton","Sur. Gumbel","Sur. Joe",
                                                                                                                "Sur. BB1","Sur. BB6","Sur. BB7","Sur. BB8","Rot. Clayton","Rot. Gumbel","Rot. Joe","Rot. BB1", "Rot. BB6",
                                                                                                                "Rot. BB7","Rot. BB8","Rot. Clayton","Rot. Gumbel","Rot. Joe","Rot. BB1","Rot. BB6","Rot. BB7","Rot. BB8",
                                                                                                                "Tawn type 1","Rot. Tawn type 1","Rot. Tawn type 1","Rot. Tawn type 1","Tawn type 2","Rot. Tawn type 2","Rot. Tawn type 2","Rot. Tawn type 2"))
  colnames(copula_table)<-c("Number","Family")

  correlation_Var1_Value<-numeric(length(Thres1))
  correlation_Var1_Test<-numeric(length(Thres1))
  correlation_Var1_N<-numeric(length(Thres1))
  copula_Var1_Family<-numeric(length(Thres1))
  copula_Var1_Family_Name<-numeric(length(Thres1))

  correlation_Var2_Value<-numeric(length(Thres2))
  correlation_Var2_Test<-numeric(length(Thres2))
  correlation_Var2_N<-numeric(length(Thres2))
  copula_Var2_Family<-numeric(length(Thres2))
  copula_Var2_Family_Name<-numeric(length(Thres2))

  if(class(Data_Detrend[,1])=="Date" | class(Data_Detrend[,1])=="factor"){
    Data_Detrend<-Data_Detrend[,-1]
  }
  if(class(Data_Declust[,1])=="Date" | class(Data_Declust[,1])=="factor"){
    Data_Declust<-Data_Declust[,-1]
  }
  
  #Conditional on Var1
  for(j in 1:length(Thres1)){
    Threshold<-Thres1[j]
    Var1_Var1_x<-which(Data_Declust[,1]>quantile(na.omit(Data_Detrend[,1]),Threshold))
    #Var1_Var1<-numeric(length(Var1_Var1_x))
    Var1_df<-array(0,dim=c(length(Var1_Var1_x),2))
    for(i in 1:length(Var1_Var1_x)){
      Var1_df[i,1]<-Data_Declust[Var1_Var1_x[i],1]
      Var1_df[i,2]<-max(Data_Detrend[c(max((Var1_Var1_x[i]-Lag_Backward_Var1),1):(min((Var1_Var1_x[i]+Lag_Forward_Var1),nrow(Data_Declust)))),2],na.rm=T)[1]
    }
    colnames(Var1_df)<-c(names(Data_Detrend))

    if(length(which(Var1_df[,1]==-Inf | Var1_df[,2]==-Inf))>0){
      z<-unique(which(Var1_df[,1]==-Inf | Var1_df[,2]==-Inf))
      Var1_df[z,]<-NA
      Var1_Var1_x[z]<-NA
    }

    if(length(which(is.na(Var1_df[,1])==TRUE | is.na(Var1_df[,2])==TRUE))>0){
      z<-unique(which(is.na(Var1_df[,1])==TRUE | is.na(Var1_df[,2])==TRUE))
      Var1_df<-Var1_df[-z,]
      Var1_Var1_x<-Var1_Var1_x[-z]
    }
    correlation_Var1_Value[j]<-cor(pobs(Var1_df[,1]), pobs(Var1_df[,2]),method="kendall")
    correlation_Var1_Test[j]<-cor.test(pobs(Var1_df[,1]), pobs(Var1_df[,2]))$p.value

    correlation_Var1_N[j]<-nrow(Var1_df)
    copula_Var1_Family[j]<-BiCopSelect(pobs(Var1_df[,1]), pobs(Var1_df[,2]), familyset = NA, selectioncrit = "AIC",
                                       indeptest = FALSE, level = 0.05, weights = NA, rotations = TRUE,
                                       se = FALSE, presel = TRUE, method = "mle")$family

    copula_Var1_Family_Name[j]<-as.character(copula_table$Family[which(copula_table$Number==copula_Var1_Family[j])])
  }

  #Conditional on Var2
  for(k in 1:length(Thres2)){
    Threshold<-Thres2[k]
    Var2_Var2_x<-which(Data_Declust[,2]>quantile(na.omit(Data_Detrend[,2]),Threshold))
    Var2_df<-array(0,dim=c(length(Var2_Var2_x),2))
    for(i in 1:length(Var2_Var2_x)){
      Var2_df[i,2]<-Data_Declust[Var2_Var2_x[i],2]
      Var2_df[i,1]<-max(Data_Detrend[c(max((Var2_Var2_x[i]-Lag_Backward_Var2),1):(min((Var2_Var2_x[i]+Lag_Forward_Var2),nrow(Data_Declust)))),1],na.rm=T)[1]
    }
    Var2_df<-data.frame(Var2_df)
    colnames(Var2_df)<-c(names(Data_Detrend))

    if(length(which(Var2_df[,1]==-Inf | Var2_df[,2]==-Inf))>0){
      z<-unique(which(Var2_df[,1]==-Inf | Var2_df[,2]==-Inf))
      Var2_df[z,]<-NA
      Var2_Var2_x[z]<-NA
    }

    if(length(which(is.na(Var2_df[,1])==TRUE | is.na(Var2_df[,2])==TRUE))>0){
      z<-unique(which(is.na(Var2_df[,1])==TRUE | is.na(Var2_df[,2])==TRUE))
      Var2_df<-Var2_df[-z,]
      Var2_Var2_x<-Var2_Var2_x[-z]
    }
    correlation_Var2_Value[k]<-cor(pobs(Var2_df[,1]), pobs(Var2_df[,2]),method="kendall")
    correlation_Var2_Test[k]<-cor.test(pobs(Var2_df[,1]), pobs(Var2_df[,2]))$p.value

    correlation_Var2_N[k]<-nrow(Var2_df)
    copula_Var2_Family[k]<-BiCopSelect(pobs(Var2_df[,1]), pobs(Var2_df[,2]), familyset = NA, selectioncrit = "AIC",
                                       indeptest = FALSE, level = 0.05, weights = NA, rotations = TRUE,
                                       se = FALSE, presel = TRUE, method = "mle")$family
    copula_Var2_Family_Name[k]<-as.character(copula_table$Family[which(copula_table$Number==copula_Var2_Family[k])])
   }

  plot(Thres1,correlation_Var1_Value,xlab="Threshold",ylab=expression("Kendall's "*tau*" correlation coefficient"),type='l',lwd=3,xlim=c(x_lim_min,x_lim_max),ylim=c(y_lim_min,y_lim_max),col="Blue")
  points(Thres1,correlation_Var1_Value,pch=ifelse(correlation_Var1_Test<0.05,16,16),col=ifelse(correlation_Var1_Test<0.05,"Blue","White"),cex=5)
  points(Thres1,correlation_Var1_Value,pch=ifelse(correlation_Var1_Test<0.05,16,1),cex=5,col="Blue")
  text(Thres1,correlation_Var1_Value,as.character(correlation_Var1_N),col=ifelse(correlation_Var1_Test<0.05,"White","Black"))

  if(sum(Upper)<1){
    text(Thres1,(correlation_Var1_Value-y_lim*GAP),copula_Var1_Family_Name,col="Blue")
  }
  if(length(Thres1[Upper])==length(Thres1)){
    text(Thres1[Upper],(correlation_Var1_Value+y_lim*GAP)[Upper],copula_Var1_Family_Name[Upper],col="Blue")
  }
  if(sum(Upper)>0 & length(Thres1[Upper])<length(Thres1)){
   text(Thres1[-Upper],(correlation_Var1_Value-y_lim*GAP)[-Upper],copula_Var1_Family_Name[-Upper],col="Blue")
   text(Thres1[Upper],(correlation_Var1_Value+y_lim*GAP)[Upper],copula_Var1_Family_Name[Upper],col="Blue")
  }

  lines(Thres2,correlation_Var2_Value,col="Red")
  points(Thres2,correlation_Var2_Value,pch=ifelse(correlation_Var2_Test<0.05,16,16),col=ifelse(correlation_Var2_Test<0.05,"Red","White"),cex=5)
  points(Thres2,correlation_Var2_Value,pch=ifelse(correlation_Var2_Test<0.05,16,1),cex=5,col="Red")
  text(Thres2,correlation_Var2_Value,as.character(correlation_Var2_N),col=ifelse(correlation_Var2_Test<0.05,"White","Black"))

  if(sum(Lower)<1){
    text(Thres2,(correlation_Var2_Value+y_lim*GAP),copula_Var2_Family_Name,col="Red")
  }
  if(length(Thres2[Lower])==length(Thres2)){
    text(Thres2[Lower],(correlation_Var2_Value-y_lim*GAP)[Lower],copula_Var2_Family_Name[Lower],col="Red")
  }
  if(sum(Lower)>0 & length(Thres2[Lower])<length(Thres2)){
    text(Thres2[-Lower],(correlation_Var2_Value+y_lim*GAP)[-Lower],copula_Var2_Family_Name[-Lower],col="Red")
    text(Thres2[Lower],(correlation_Var2_Value-y_lim*GAP)[Lower],copula_Var2_Family_Name[Lower],col="Red")
  }

  if(Legend==TRUE){
    legend("bottomleft",c(paste("Conditioning on ",colnames(Data_Detrend)[1],sep=""),
                          paste("Conditioning on ",colnames(Data_Detrend)[2],sep="")),
           bty="n",lwd=3,col=c("Blue","Red"))
  }

  res<-list("Kendalls_Tau_Var1" = correlation_Var1_Value,
            "p_value_Var1" = correlation_Var1_Test,
            "N_Var1" = correlation_Var1_N,
            "Copula_Family_Var1" =copula_Var1_Family,
            "Kendalls_Tau_Var2" = correlation_Var2_Value,
            "p_value_Var2" = correlation_Var2_Test,
            "N_Var2" = correlation_Var2_N,
            "Copula_Family_Var2" = copula_Var2_Family)
  return(res)
}
