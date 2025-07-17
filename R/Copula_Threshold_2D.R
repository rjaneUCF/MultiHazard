#' Copula Selection With threshold 2D - Fit
#'
#' Declustered excesses of a (conditioning) variable are paired with co-occurences of the other variable before the best fitting bivariate copula is selected, using \code{BiCopSelect} function in the \code{VineCopula} package, for a single or range of thresholds. The procedure is automatically repeated with the variables switched.
#'
#' @param Data_Detrend Data frame containing two at least partially concurrent time series, detrended if necessary. Time steps must be equally spaced, with missing values assigned \code{NA}.
#' @param Data_Declust Data frame containing two (independently) declustered at least partially concurrent time series. Time steps must be equally spaced, with missing values assigned \code{NA}.
#' @param u1 A single or sequence of thresholds, given as a quantile of the observations of the variable in the first column of \code{Data_Detrend} when it is used as the conditioning variable. Default, sequence from \code{0.9} to \code{0.99} at intervals of \code{0.01}.
#' @param u2 A single or sequence of thresholds, given as a quantile of the observations of the variable in the second column of \code{Data_Detrend} when it is used as the conditioning variable. Default, sequence from \code{0.9} to \code{0.99} at intervals of \code{0.01}.
#' @param PLOT Logical; whether to plot the results. Default is \code{"TRUE"}.
#' @param x_lim_min Numeric vector of length one specifying x-axis minimum. Default is \code{NA}.
#' @param x_lim_max Numeric vector of length one specifying x-axis maximum. Default is \code{NA}.
#' @param y_lim_min Numeric vector of length one specifying y-axis minimum. Default \code{-1.0}.
#' @param y_lim_max Numeric vector of length one specifying y-axis maximum. Default \code{1.0}.
#' @param Upper Numeric vector specifying the element number of the \code{u1} argument for which the copula family name label to appear above the corresponding point on the Kendall's tau coefficient vs threshold plot, when conditioning on the variable in column 1. Default is \code{0}.
#' @param Lower Numeric vector specifying the element number of the \code{u2} argument for which the copula family name label to appear below the corresponding point on the Kendall's tau coefficient vs threshold plot, when conditioning on the variable in column 2. Default is \code{0}.
#' @param GAP Numeric vector of length one specifying the distance above or below the copula family name label appears the corresponding point on the Kendall's tau coefficient vs threshold plot. Default is \code{0.05}.
#' @param Legend Logic vector of length one specifying whether a legend should be plotted. Default is \code{TRUE}.
#' @param Cex_Legend Numeric vector of length one specifying the font size of the legend. Default is \code{1}.
#' @param Cex_Axis Numeric vector of length one specifying the font size of the axes. Default is \code{1}.
#' @param Cex_Axis_Original Numeric vector of length one specifying the font size of the values of the quantiles on the original (data) scale (i.e. second x-axis). Default is \code{1}.
#' @return List comprising: \itemize{
#' \item \code{Kendalls_Tau1}
#' Kendall's tau of a sample
#' \item \code{p_value_Var1}
#' p-value when testing the null hypothesis \code{H_0: tau=0} i.e. that there is no correlation between the variables
#' \item \code{N_Var1}
#' Size of the dataset
#' \item \code{Copula_Family_Var1}
#' Best fitting copula for the specified thresholds
#' }
#' when the dataset is conditioned on the variable in column 1.
#' Analogous vectors \code{Kendalls_Tau2},\code{p_value_Var2}, \code{N_Var2} and \code{Copula_Family_Var2} for the specified thresholds when the data set is conditioned on the variable in column 2.
#' If \code{PLOT=TRUE} then a plot of the Kendall's tau correlation coefficient versus quantile threshold is also returned.
#' Filled circles denote statistically significant correlation at a \code{5\%} significance level. Numbers inside the circles correspond to the sample size while the best fitting copula family is printed above.
#' Numbers below x-axis are the values of the corresponding quantiles on the original (data) scale.
#' @seealso \code{\link{Dataframe_Combine}}
#' @export
#' @examples
#' Copula_Threshold_2D(Data_Detrend=S20.Detrend.df[,-c(1,4)],
#'                     Data_Declust=S20.Detrend.Declustered.df[,-c(1,4)],
#'                     y_lim_min=-0.075, y_lim_max =0.25,
#'                     Upper=c(6,8), Lower=c(6,8),GAP=0.1)
Copula_Threshold_2D<-function(Data_Detrend,Data_Declust,u1=seq(0.9,0.99,0.01),u2=seq(0.9,0.99,0.01),PLOT=TRUE,x_lim_min=NA,x_lim_max=NA,y_lim_min=-1,y_lim_max=1,Upper=NA,Lower=NA,GAP=0.05,Legend=TRUE,Cex_Legend=1,Cex_Axis=1,Cex_Axis_Original=1){

  # Input validation
  if (!is.data.frame(Data_Detrend) && !is.matrix(Data_Detrend)) {
    stop("Data_Detrend must be a data frame or matrix")
  }

  if (!is.data.frame(Data_Declust) && !is.matrix(Data_Declust)) {
    stop("Data_Declust must be a data frame or matrix")
  }

  if (ncol(Data_Detrend) != 2) {
    stop("Data_Detrend must comprise two columns, got: ", ncol(Data_Detrend))
  }

  if (ncol(Data_Declust) != 2) {
    stop("Data_Declust must comprise two columns, got: ", ncol(Data_Declust))
  }


  if (!is.na(u1) & (any(u1 > 1) || any(u1 < 0))) {
    stop("u1 must be between 0 and 1, got values in range: ", min(u1), "to", max(u1))
  }

  if (!is.na(u2) &  (any(u2 > 1) || any(u2 < 0))) {
    stop("u2 must be between 0 and 1, got values in range: ", min(u2), "to", max(u2))
  }

  if (y_lim_min >= y_lim_max) {
    stop("y_lim_min must be less than y_lim_max, given: y_lim_min = ", y_lim_min, ", y_lim_max = ", y_lim_max)
  }

  #Axes limits for plots
  x_lim_min<-ifelse(is.na(x_lim_min)==T,min(u1,u2,na.rm=T),x_lim_min)
  x_lim_max<-ifelse(is.na(x_lim_max)==T,max(u1,u2,na.rm=T),x_lim_max)
  y_lim=y_lim_max-y_lim_min

  #Table of copula codes in Vine copula package
  copula_table<-data.frame(c(seq(0,40,1)[-(c(11,12,15,21,22,25,31,32,35)+1)],104,114,124,134,204,214,224,234),
                           c("Ind.","Gaussian", "t-copula", "Clayton", "Gumbel","Frank","Joe","BB1","BB6","BB7","BB8","Sur. Clayton","Sur. Gumbel","Sur. Joe",
                             "Sur. BB1","Sur. BB6","Sur. BB7","Sur. BB8","Rot. Clayton","Rot. Gumbel","Rot. Joe","Rot. BB1", "Rot. BB6",
                             "Rot. BB7","Rot. BB8","Rot. Clayton","Rot. Gumbel","Rot. Joe","Rot. BB1","Rot. BB6","Rot. BB7","Rot. BB8",
                             "Tawn type 1","Rot. Tawn type 1","Rot. Tawn type 1","Rot. Tawn type 1","Tawn type 2","Rot. Tawn type 2","Rot. Tawn type 2","Rot. Tawn type 2"))
  colnames(copula_table)<-c("Number","Family")

  #Removing date column
  if(class(Data_Detrend[,1])[1]=="Date" | class(Data_Detrend[,1])[1]=="factor" | class(Data_Detrend[,1])[1]=="POSIXct" | class(Data_Detrend[,1])[1]=="character"){
    Data_Detrend<-Data_Detrend[,-1]
  }
  if(class(Data_Declust[,1])[1]=="Date" | class(Data_Declust[,1])[1]=="factor" | class(Data_Declust[,1])[1]=="POSIXct" | class(Data_Detrend[,1])[1]=="character"){
    Data_Declust<-Data_Declust[,-1]
  }

  #Conditional on Var1
  if(is.na(u1[1])==FALSE){
    correlation_Var1_Value<-numeric(length(u1))
    correlation_Var1_Test<-numeric(length(u1))
    correlation_Var1_N<-numeric(length(u1))
    copula_Var1_Family<-numeric(length(u1))
    copula_Var1_Family_Name<-numeric(length(u1))
    for(j in 1:length(u1)){
      u<-u1[j]
      Var1_Var1_x<-which(Data_Declust[,1]>quantile(na.omit(Data_Detrend[,1]),u))
      Var1_Var1<-numeric(length(Var1_Var1_x))
      Var1_df<-array(0,dim=c(length(Var1_Var1_x),2))
      for(i in 1:length(Var1_Var1_x)){
        Var1_df[i,1]<-Data_Declust[Var1_Var1_x[i],1]
        Var1_df[i,2]<-Data_Detrend[Var1_Var1_x[i],2]
      }
      colnames(Var1_df)<-c(names(Data_Detrend))
      if(length(which(is.na(Var1_df[,1])==TRUE | is.na(Var1_df[,2])==TRUE))>0){
        z<-unique(which(is.na(Var1_df[,1])==TRUE | is.na(Var1_df[,2])==TRUE))
        Var1_df<-Var1_df[-z,]
        Var1_Var1_x<-Var1_Var1_x[-z]
      }
      correlation_Var1_Value[j]<-cor(pobs(Var1_df[,1]), pobs(Var1_df[,2]),method="kendall")
      correlation_Var1_Test[j]<-suppressWarnings(cor.test(pobs(Var1_df[,1]), pobs(Var1_df[,2]),method="kendall")$p.value)
      correlation_Var1_N[j]<-nrow(Var1_df)
      copula_Var1_Family[j]<-suppressWarnings(BiCopSelect(pobs(Var1_df[,1]), pobs(Var1_df[,2]), familyset = NA, selectioncrit = "AIC",
                                                          indeptest = FALSE, level = 0.05, weights = NA, rotations = TRUE,
                                                          se = FALSE, presel = TRUE, method = "mle")$family)
      copula_Var1_Family_Name[j]<-as.character(copula_table$Family[which(copula_table$Number==copula_Var1_Family[j])])
    }
  }

  #Conditional on Var2
  if(is.na(u2[1])==FALSE){
    correlation_Var2_Value<-numeric(length(u2))
    correlation_Var2_Test<-numeric(length(u2))
    correlation_Var2_N<-numeric(length(u2))
    copula_Var2_Family<-numeric(length(u2))
    copula_Var2_Family_Name<-numeric(length(u2))
    for(k in 1:length(u2)){
      u<-u2[k]
      Var2_Var2_x<-which(Data_Declust[,2]>quantile(na.omit(Data_Detrend[,2]),u))
      Var2_df<-array(0,dim=c(length(Var2_Var2_x),2))
      for(i in 1:length(Var2_Var2_x)){
        Var2_df[i,2]<-Data_Declust[Var2_Var2_x[i],2]
        Var2_df[i,1]<-Data_Detrend[Var2_Var2_x[i],1]
      }
      Var2_df<-data.frame(Var2_df)
      colnames(Var2_df)<-c(names(Data_Detrend))
      if(length(which(is.na(Var2_df[,1])==TRUE | is.na(Var2_df[,2])==TRUE))>0){
        z<-unique(which(is.na(Var2_df[,1])==TRUE | is.na(Var2_df[,2])==TRUE))
        Var2_df<-Var2_df[-z,]
        Var2_Var2_x<-Var2_Var2_x[-z]
      }
      correlation_Var2_Value[k]<-cor(pobs(Var2_df[,1]), pobs(Var2_df[,2]),method="kendall")
      correlation_Var2_Test[k]<-suppressWarnings(cor.test(pobs(Var2_df[,1]), pobs(Var2_df[,2]),method="kendall")$p.value)
      correlation_Var2_N[k]<-nrow(Var2_df)
      copula_Var2_Family[k]<-suppressWarnings(BiCopSelect(pobs(Var2_df[,1]), pobs(Var2_df[,2]), familyset = NA, selectioncrit = "AIC",
                                                          indeptest = FALSE, level = 0.05, weights = NA, rotations = TRUE,
                                                          se = FALSE, presel = TRUE, method = "mle")$family)
      copula_Var2_Family_Name[k]<-as.character(copula_table$Family[which(copula_table$Number==copula_Var2_Family[k])])
    }
  }

  if(PLOT==TRUE){
    if(is.na(u1[1])==FALSE & is.na(u2[1])==FALSE){
      plot(u1,correlation_Var1_Value,xlab="Threshold",ylab=expression("Kendall's "*tau*" correlation coefficient"),type='l',lwd=3,xlim=c(x_lim_min,x_lim_max),ylim=c(y_lim_min,y_lim_max),col="Blue",cex=Cex_Axis)
      mtext(round(quantile(na.omit(Data_Detrend[,1]),u1),2),at=u1,side=1,line=2,col="Blue",cex=Cex_Axis_Original)
      points(u1,correlation_Var1_Value,pch=ifelse(correlation_Var1_Test<0.05,16,16),col=ifelse(correlation_Var1_Test<0.05,"Blue","White"),cex=5)
      points(u1,correlation_Var1_Value,pch=ifelse(correlation_Var1_Test<0.05,16,1),cex=5,col="Blue")
      text(u1,correlation_Var1_Value,as.character(correlation_Var1_N),col=ifelse(correlation_Var1_Test<0.05,"White","Black"))

      if(is.na(Upper[1])==TRUE){
        text(u1,(correlation_Var1_Value-y_lim*GAP),copula_Var1_Family_Name,col="Blue")
      }
      if(length(u1[Upper])==length(u1)){
        text(u1[Upper],(correlation_Var1_Value+y_lim*GAP)[Upper],copula_Var1_Family_Name[Upper],col="Blue")
      }
      if(sum(Upper)>0 & length(u1[Upper])<length(u1)){
        text(u1[-Upper],(correlation_Var1_Value-y_lim*GAP)[-Upper],copula_Var1_Family_Name[-Upper],col="Blue")
        text(u1[Upper],(correlation_Var1_Value+y_lim*GAP)[Upper],copula_Var1_Family_Name[Upper],col="Blue")
      }

      mtext(round(quantile(na.omit(Data_Detrend[,2]),u2),2),at=u2,side=3,line=1,col="Red")
      lines(u2,correlation_Var2_Value,col="Red")
      points(u2,correlation_Var2_Value,pch=ifelse(correlation_Var2_Test<0.05,16,16),col=ifelse(correlation_Var2_Test<0.05,"Red","White"),cex=5)
      points(u2,correlation_Var2_Value,pch=ifelse(correlation_Var2_Test<0.05,16,1),cex=5,col="Red")
      text(u2,correlation_Var2_Value,as.character(correlation_Var2_N),col=ifelse(correlation_Var2_Test<0.05,"White","Black"))

      if(is.na(Lower[1])==TRUE){
        text(u2,(correlation_Var2_Value+y_lim*GAP),copula_Var2_Family_Name,col="Red")
      }
      if(length(u2[Lower])==length(u2)){
        text(u2[Lower],(correlation_Var2_Value-y_lim*GAP)[Lower],copula_Var2_Family_Name[Lower],col="Red")
      }
      if(sum(Lower)>0 & length(u2[Lower])<length(u2)){
        text(u2[-Lower],(correlation_Var2_Value+y_lim*GAP)[-Lower],copula_Var2_Family_Name[-Lower],col="Red")
        text(u2[Lower],(correlation_Var2_Value-y_lim*GAP)[Lower],copula_Var2_Family_Name[Lower],col="Red")
      }

      if(Legend==TRUE){
        legend("bottomleft",c(paste("Conditioning on ",colnames(Data_Detrend)[1],sep=""),
                              paste("Conditioning on ",colnames(Data_Detrend)[2],sep="")),
               bty="n",lwd=1,col=c("Blue","Red"),cex=Cex_Legend)
      }
    }

    if(is.na(u2[1])==TRUE){
      plot(u1,correlation_Var1_Value,xlab="Threshold",ylab=expression("Kendall's "*tau*" correlation coefficient"),type='l',lwd=3,xlim=c(x_lim_min,x_lim_max),ylim=c(y_lim_min,y_lim_max),col="Blue",cex=Cex_Axis)
      mtext(round(quantile(na.omit(Data_Detrend[,1]),u1),2),at=u1,side=1,line=2,col="Blue",cex=Cex_Axis_Original)
      points(u1,correlation_Var1_Value,pch=ifelse(correlation_Var1_Test<0.05,16,16),col=ifelse(correlation_Var1_Test<0.05,"Blue","White"),cex=5)
      points(u1,correlation_Var1_Value,pch=ifelse(correlation_Var1_Test<0.05,16,1),cex=5,col="Blue")
      text(u1,correlation_Var1_Value,as.character(correlation_Var1_N),col=ifelse(correlation_Var1_Test<0.05,"White","Black"))

      if(is.na(Upper[1])==TRUE){
        text(u1,(correlation_Var1_Value-y_lim*GAP),copula_Var1_Family_Name,col="Blue")
      }
      if(length(u1[Upper])==length(u1)){
        text(u1[Upper],(correlation_Var1_Value+y_lim*GAP)[Upper],copula_Var1_Family_Name[Upper],col="Blue")
      }
      if(sum(Upper)>0 & length(u1[Upper])<length(u1)){
        text(u1[-Upper],(correlation_Var1_Value-y_lim*GAP)[-Upper],copula_Var1_Family_Name[-Upper],col="Blue")
        text(u1[Upper],(correlation_Var1_Value+y_lim*GAP)[Upper],copula_Var1_Family_Name[Upper],col="Blue")
      }

      if(Legend==TRUE){
        legend("bottomleft",paste("Conditioning on ",colnames(Data_Detrend)[1],sep=""),
               bty="n",lwd=3,col="Blue",cex=Cex_Legend)
      }
    }

    if(is.na(u1[1])==TRUE){
      plot(u2,correlation_Var2_Value,xlab="Threshold",ylab=expression("Kendall's "*tau*" correlation coefficient"),type='l',lwd=3,xlim=c(x_lim_min,x_lim_max),ylim=c(y_lim_min,y_lim_max),col="Red",cex=Cex_Axis)
      mtext(round(quantile(na.omit(Data_Detrend[,2]),u2),2),at=u2,side=1,line=2,col="Red",cex=Cex_Axis_Original)
      points(u2,correlation_Var2_Value,pch=ifelse(correlation_Var2_Test<0.05,16,16),col=ifelse(correlation_Var2_Test<0.05,"Red","White"),cex=5)
      points(u2,correlation_Var2_Value,pch=ifelse(correlation_Var2_Test<0.05,16,1),cex=5,col="Red")
      text(u2,correlation_Var2_Value,as.character(correlation_Var2_N),col=ifelse(correlation_Var2_Test<0.05,"White","Black"))

      if(is.na(Lower[1])==TRUE){
        text(u2,(correlation_Var2_Value+y_lim*GAP),copula_Var2_Family_Name,col="Red")
      }
      if(length(u2[Lower])==length(u2)){
        text(u2[Lower],(correlation_Var2_Value-y_lim*GAP)[Lower],copula_Var2_Family_Name[Lower],col="Red")
      }
      if(sum(Lower)>0 & length(u2[Lower])<length(u2)){
        text(u2[-Lower],(correlation_Var2_Value+y_lim*GAP)[-Lower],copula_Var2_Family_Name[-Lower],col="Red")
        text(u2[Lower],(correlation_Var2_Value-y_lim*GAP)[Lower],copula_Var2_Family_Name[Lower],col="Red")
      }

      if(Legend==TRUE){
        legend("bottomleft",paste("Conditioning on ",colnames(Data_Detrend)[2],sep=""),
               bty="n",lwd=3,col="Red",cex=Cex_Legend)
      }
    }
  }

  if(is.na(u1[1])==FALSE & is.na(u2[1])==FALSE){
    res<-list("Kendalls_Tau1" = correlation_Var1_Value,"p_value_Var1" = correlation_Var1_Test,
              "N_Var1" = correlation_Var1_N,"Copula_Family_Var1" =copula_Var1_Family,
              "Kendalls_Tau2" = correlation_Var2_Value,"p_value_Var2" = correlation_Var2_Test,
              "N_Var2" = correlation_Var2_N,"Copula_Family_Var2" = copula_Var2_Family)
  }

  if(is.na(u1[1])==TRUE){
    res<-list("Kendalls_Tau2" = correlation_Var2_Value,"p_value_Var2" = correlation_Var2_Test,
              "N_Var2" = correlation_Var2_N,"Copula_Family_Var2" = copula_Var2_Family)
  }

  if(is.na(u2[1])==TRUE){
    res<-list("Kendalls_Tau1" = correlation_Var1_Value,"p_value_Var1" = correlation_Var1_Test,
              "N_Var1" = correlation_Var1_N,"Copula_Family_Var1" =copula_Var1_Family)
  }
  return(res)
}
