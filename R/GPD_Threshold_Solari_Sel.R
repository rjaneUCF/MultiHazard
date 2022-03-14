#' Goodness-of-fit for the GPD
#'
#' A nonparametric bootstrapping procedure is undertaken to assess the uncertainty in the GPD parameters and associated return levels for a GPD fit to observations above a user specified threshold. The estimates are compared with those obtained at other thresholds by running the \code{GPD_Threshold_Solari} function beforehand, and using its output as an input of this function. The code is based on the \code{AUTOMATICO_MLE_BOOT} function provided by Sebastian Solari.
#'
#' @param Event Numeric vector containing independent events declustered using a moving window approach.
#' @param Data Original time series. Dataframe containing two columns. In column: \itemize{
#' \item \code{1} A \code{"Date"} object of equally spaced discrete time steps.
#' \item \code{2} Numeric vector containing corresponding time series values.
#' }
#' @param Solari_Output Output of the \code{GPD_Threshold_Solari} function.
#' @param Thres Numeric vector of length one specifying the threshold to analyze, chosen by the user based on plots from the \code{GPD_Threshold_Solari} function.
#' @param Alpha Numeric vector of length one specifying the level of confidence associated with the confidence interval i.e., the probability that the interval contains the true value of the parameter is \eqn{1-\frac{Alpha}{2}}. The interval is referred to as the \eqn{100(1-\frac{Alpha}{2})\%} confidence interval. Default is \code{0.1}.
#' @param N_Sim Numeric vector of length one specifying the number of bootstrap samples. Default is \code{10^4}.
#' @param RP_Max Numeric vector of length one specifying the maximum return level to be calculated. Default is \code{1000}.
#' @param RP_Plot Numeric vector of length one specifying the return level in the lower right plot. Default is \code{100}.
#' @param mu (average) occurrence frequency of events in the original time series \code{Data}. Numeric vector of length one. Default is \code{365.25}, daily data.
#' @param y_lab Character vector specifying the y-axis label of the return level plot.
#' @return List containing three objects: \code{Estimate}, \code{CI_Upper} and \code{CI_Lower}. The \code{Estimate} dataframe comprises \itemize{
#' \item \code{xi}
#' GPD shape parameter estimate for the threshold is \code{Thres}.
#' \item \code{sigma}
#' GPD scale parameter estimate for the threshold is \code{Thres}.
#' \item \code{Thres}
#' GPD location parameter estimate for the threshold is \code{Thres}.
#' \item \code{rate}
#' GPD rate parameter i.e., number of independent excesses per year for a threshold of \code{Thres}.
#' \item The remaining columns are \code{RL}
#' Return level estimates from the GPD using a threshold of \code{Thres}.
#' }
#' \code{CI_Upper} and \code{CI_Lower} give the upper and lower bounds of the \eqn{100(1-\frac{Alpha}{2})\%} confidence interval for the corresponding element in \code{Estimate}.
#' Top row: Histograms of the GPD parameter estimates based on a nonparametric bootstrapping simulation. Grey bars correspond to the estimates obtained as the threshold (\code{Thres}) is varied, found by running the function a necessary input of the function. Continuous black lines correspond to results obtained by fixing the threshold at \code{Thres}. Dashed blue lines correspond to the expected values for the fixed threshold.
#' Lower left: Return level plot. Return levels of the observations estimated from the empirical distribution. Grey bars correspond to the maximum of the upper and lower bounds of the \eqn{100(1-\frac{Alpha}{2})\%} confidence intervals as the threshold is varied. Continuous black lines correspond to results obtained by fixing the threshold at \code{Thres}. Dashed blue lines correspond to the expected values for the fixed threshold.
#' Lower right: As in the top row but for the 100 years return period quantile.
#' @export
#' @examples
#  Declustering the O-sWL at site S22 using a 3-day window.
#' Rainfall_Declust_SW<-Decluster_SW(Data=S22.Detrend.df[,c(1:2)],Window_Width=7)
#' Finding an appropriate threshold for the declustered series
#' S22_OsWL_Solari<-GPD_Threshold_Solari(Event=Rainfall_Declust_SW$Declustered,
#'                                   Data=na.omit(S22.Detrend.df[,c(1:2)]))
#'
GPD_Threshold_Solari_Sel<-function(Event,Data,Solari_Output,Thres,Alpha=0.1,N_Sim=10^4,RP_Max=1000,RP_Plot=100,mu=365.25,y_lab="Data"){

  # Auxiliary variables
  Data    = na.omit(Data)
  N_Years = length(Data)/mu
  RP      = sort(unique(c(1:10,seq(20,100,10),seq(200,1000,100),seq(2e3,1e4,1e3),seq(2e4,1e5,1e4),seq(2e5,1e6,1e5),RP_Plot)))
  RP      = RP[RP<=RP_Max]
  Event   = na.omit(Event)

  # POT with MLE and C.I. with bootstrapping ##########################
  # Sort Events
  Event = sort(Event)

  for(J in 1:length(Thres)){

    Exceedence = Event[Event>Thres[J]]
    Rate = length(Exceedence)/N_Years

    Estimate   = unlist(GPD_MLE_Boot(Exceedence,RP,N_Years))
    BOOT   = array(0,dim=c(N_Sim,length(Estimate)))
    for(I in 1:N_Sim){
      BOOT[I,] = unlist(GPD_MLE_Boot(rgpd(length(Exceedence),xi=as.numeric(Estimate[1]),sigma=as.numeric(Estimate[2]),u=as.numeric(Estimate[3])),RP,N_Years))
      while( BOOT[I,1]==0){
        try({
          BOOT[I,] = unlist(GPD_MLE_Boot(sample(Exceedence,length(Exceedence),replace=T),RP,N_Years))
        }, silent = FALSE)
      }

      CI.Upper = apply(BOOT, 2,  quantile, 1-Alpha/2)
      CI.Lower = apply(BOOT, 2,  quantile, Alpha/2)
    }
  }
  #Define the layout of the plots
  layout_matrix <- matrix(c(1, 4, 2, 4, 3, 5), nrow = 2)
  layout(layout_matrix)

  # Histogram of k (GPD shape)
  z<-which(Solari_Output$GPD_MLE[,1]>-0.5 & Solari_Output$GPD_MLE[,1]<0.5 & Solari_Output$GPD_MLE[,6]<(length(na.omit(Data)/mu)))
  h<-hist(Solari_Output$GPD_MLE[z,1],xlab="GP shape",ylab="Frequency",col="Grey",
          xlim=c(min(Solari_Output$GPD_MLE[z,1])-diff(range(Solari_Output$GPD_MLE[z,1]))/4,
                 max(Solari_Output$GPD_MLE[z,1])+diff(range(Solari_Output$GPD_MLE[z,1]))/4),
          main="",boarder = "Grey")
  abline(v=Estimate[1],col="Blue",lwd=2)
  par(new=TRUE)
  K<-BOOT[,1][BOOT[,1]>min(h$breaks) & BOOT[,1]<max(h$breaks)]
  v<-hist(K,breaks=h$breaks,plot=FALSE)
  plot(c(v$breaks[1],v$breaks),c(0,v$counts,0),type="s",
       xlim=c(min(Solari_Output$GPD_MLE[z,1])-diff(range(Solari_Output$GPD_MLE[z,1]))/4,
              max(Solari_Output$GPD_MLE[z,1])+diff(range(Solari_Output$GPD_MLE[z,1]))/4),
       lwd=2,xlab="",ylab="",xaxt='n',yaxt='n')

  # Histogram of sigma (GPD scale)
  h<-hist(Solari_Output$GPD_MLE[z,2],xlab="GP scale",ylab="Frequency",col="Grey",
          xlim=c(min(Solari_Output$GPD_MLE[z,2])-diff(range(Solari_Output$GPD_MLE[z,2]))/4,
                 max(Solari_Output$GPD_MLE[z,2])+diff(range(Solari_Output$GPD_MLE[z,2]))/4),
          main="",boarder = "Grey")
  abline(v=Estimate[2],col="Blue",lwd=2)
  par(new=TRUE)
  K<-BOOT[,2][BOOT[,2]>min(h$breaks) & BOOT[,2]<max(h$breaks)]
  v<-hist(K,breaks=h$breaks,plot=FALSE)
  plot(c(v$breaks[1],v$breaks),c(0,v$counts,0),type="s",
       xlim=c(min(Solari_Output$GPD_MLE[z,2])-diff(range(Solari_Output$GPD_MLE[z,2]))/4,
              max(Solari_Output$GPD_MLE[z,2])+diff(range(Solari_Output$GPD_MLE[z,2]))/4),
       lwd=2,xlab="",ylab="",xaxt='n',yaxt='n')

  # Histogram of u (GPD location)
  h<-hist(Solari_Output$GPD_MLE[z,3],
          xlab="GP position",ylab="Frequency",col="Grey",
          xlim=c(min(Solari_Output$GPD_MLE[z,3])-diff(range(Solari_Output$GPD_MLE[z,3]))/4,
                 max(Solari_Output$GPD_MLE[z,3])+diff(range(Solari_Output$GPD_MLE[z,3]))/4),
          main="",boarder = "Grey")
  abline(v=Estimate[3],col="Blue",lwd=2)
  par(new=TRUE)
  K<-BOOT[,3][BOOT[,3]>min(h$breaks) & BOOT[,3]<max(h$breaks)]
  v<-hist(K,breaks=h$breaks,plot=FALSE)
  plot(c(v$breaks[1],v$breaks),c(0,v$counts,0),type="s",
       xlim=c(min(Solari_Output$GPD_MLE[z,3])-diff(range(Solari_Output$GPD_MLE[z,3]))/4,
              max(Solari_Output$GPD_MLE[z,3])+diff(range(Solari_Output$GPD_MLE[z,3]))/4),
       lwd=2,xlab="",ylab="",xaxt='n',yaxt='n')

  # Histogram for a return level of 100 years
  Solari_Output_RPs<-as.numeric(colnames(Solari_Output$GPD_MLE)[7:ncol(Solari_Output$GPD_MLE)])
  Solari_Output_RPs<-Solari_Output_RPs[which(Solari_Output_RPs<=RP_Max)]
  plot(log10(1/((1-(1:length(Exceedence))/(length(Exceedence)+1))*Rate)),sort(Event[which(Event>Thres)]),
       xlab=expression('log'[10]*'(Return period) [years]'),ylab=y_lab,
       xlim=c(log10(1),log10(RP_Max)),
       ylim=c(min(Event[which(Event>Thres)],min(apply(Solari_Output$CI_Lower[z,c(7:(6+length(Solari_Output_RPs)))],2,min))),max(apply(Solari_Output$CI_Upper[z,rev(c(7:(6+length(Solari_Output_RPs))))],2,max))),
       pch=16,col="Green")
  polygon(c(log10(as.numeric(colnames(Solari_Output$GPD_MLE)[7:(6+length(Solari_Output_RPs))])),log10(as.numeric(colnames(Solari_Output$GPD_MLE)[rev(c(7:(6+length(Solari_Output_RPs))))]))),
          c(apply(Solari_Output$CI_Lower[z,c(7:(6+length(Solari_Output_RPs)))],2,min),apply(Solari_Output$CI_Upper[z,rev(c(7:(6+length(Solari_Output_RPs))))],2,max)),
          col="Grey",boarder = FALSE)
  lines(log10(RP),Estimate[-(1:4)],col="Blue",lwd=2)
  lines(log10(RP),CI.Upper[-(1:4)],lwd=2)
  lines(log10(RP),CI.Lower[-(1:4)],lwd=2)
  points(log10(1/((1-(1:length(Exceedence))/(length(Exceedence)+1))*Rate)),sort(Event[which(Event>Thres)]),col="Green",pch=16)

  h<-hist(Solari_Output$GPD_MLE[z,which(colnames(Solari_Output$GPD_MLE)==RP_Plot)],
          xlab=paste(RP_Plot,'year return level',y_lab),ylab="Frequency",col="Grey",
          xlim=c(min(Solari_Output$GPD_MLE[z,which(colnames(Solari_Output$GPD_MLE)==RP_Plot)])-diff(range(Solari_Output$GPD_MLE[z,which(colnames(Solari_Output$GPD_MLE)==RP_Plot)]))/4,
                 max(Solari_Output$GPD_MLE[z,which(colnames(Solari_Output$GPD_MLE)==RP_Plot)])+diff(range(Solari_Output$GPD_MLE[z,which(colnames(Solari_Output$GPD_MLE)==RP_Plot)]))/4),
          main="",boarder = "Grey")
  abline(v=Estimate[-(1:4)][which(RP==RP_Plot)],col="Blue",lwd=2)
  par(new=TRUE)
  K<-BOOT[,which(RP==RP_Plot)][BOOT[,which(RP==RP_Plot)]>min(h$breaks) & BOOT[,which(RP==RP_Plot)]<max(h$breaks)]
  v<-hist(K,breaks=h$breaks,plot=FALSE)
  plot(c(v$breaks[1],v$breaks),c(0,v$counts,0),type="s",
       xlim=c(min(Solari_Output$GPD_MLE[z,which(colnames(Solari_Output$GPD_MLE)==RP_Plot)])-diff(range(Solari_Output$GPD_MLE[z,which(colnames(Solari_Output$GPD_MLE)==RP_Plot)]))/4,
              max(Solari_Output$GPD_MLE[z,which(colnames(Solari_Output$GPD_MLE)==RP_Plot)])+diff(range(Solari_Output$GPD_MLE[z,which(colnames(Solari_Output$GPD_MLE)==RP_Plot)]))/4),
       lwd=2,xlab="",ylab="",xaxt='n',yaxt='n')

  names(Estimate)<-c("xi","sigma","Thres","rate",as.character(RP))
  res<-list(Estimate=Estimate,CI_Lower=CI.Lower,CI_Upper=CI.Upper,BOOT=BOOT)
  return(res)
}
