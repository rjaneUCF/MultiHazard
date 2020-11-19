#' Demonstrate the goodness of fit of the selected non-extreme marginal distribution
#'
#' Plots demonstrating the goodness of fit of a selected (not truncated) non-extreme marginal distribution to a dataset.
#'
#' @param Data Numeric vector containing realizations of the variable of interest.
#' @param x_lab Numeric vector of length one specifyingLabel on the x-axis of histogram and cummulative distribution plot.
#' @param y_lim_min Numeric vector of length one specifying the lower y-axis limit of the histogram.
#' @param y_lim_max Numeric vector of length one specifying the upper y-axis limit of the histogram.
#' @param Selected Charactor vector of length one specifying the chosen distribution, options are the Gaussian \code{"Gaus"} and logistic \code{"Logis"}.
#' @return Panel consisting of three plots. Upper plot: Plots depicting the AIC of the two fitted distributions. Middle plot: Probabilty Density Functions (PDFs) of the \code{selected} distribtions superimposed on a histgram of the data. Lower plot: Cummulative distribution function (CDFs) of the \code{selected} distribution overlaid on a plot of the empirical CDF.
#' @seealso \code{\link{Diag_Non_Con}}
#' @export
#' @examples
#' S20.Rainfall<-Con_Sampling_2D(Data_Detrend=S20.Detrend.df[,-c(1,4)],
#'                               Data_Declust=S20.Detrend.Declustered.df[,-c(1,4)],
#'                               Con_Variable="Rainfall",Thres=0.97)
#' Diag_Non_Con(Data=S20.Rainfall$Data$OsWL,x_lab="O-sWL (ft NGVD 29)",
#'              y_lim_min=0,y_lim_max=1.5)
#' Diag_Non_Con_Sel(Data=S20.Rainfall$Data$OsWL,x_lab="O-sWL (ft NGVD 29)",
#'                  y_lim_min=0,y_lim_max=1.5,Selected="Twe")
Diag_Non_Con_Sel<-function(Data,x_lab = "Data",y_lim_min = 0,y_lim_max = 1,Selected){
  par(mfrow=c(3,1))
  par(mar=c(4.2,4.2,1,1))

  fit<-fitdistr(Data, "normal")
  AIC.Normal<-2*length(fit$estimate)-2*fit$loglik

  fit<-fitdistr(Data,"logistic")
  AIC.Logistic<-2*length(fit$estimate)-2*fit$loglik

  mypalette<-brewer.pal(9,"Set1")
  plot(0,xlim=c(0,2),ylim=c(min(0,AIC.Normal,AIC.Logistic),max(0,AIC.Normal,AIC.Logistic)),type='n',xlab="Probability Distribution",ylab="AIC",xaxt='n',cex.axis=1,cex.lab=1,las=1)
  axis(1,seq(0.5,1.5,1),c("Gaus","Logis"),cex.axis=0.71)
  colors <- rep("white",2)
  dist<-c("Gaus","Logis")
  j <-  which(dist==Selected)
  colors[j] <- mypalette[c(4,6)][j]
  #rect(0.25,0,0.75,AIC.BS,col=colors[1])
  #rect(1.25,0,1.75,AIC.Exp,col=colors[2])
  #rect(2.25,0,2.75,AIC.Gamma,col=colors[3])
  rect(0.25,0,0.75,AIC.Normal,col=colors[1])
  #rect(4.25,0,4.75,AIC.InverseNormal,col=colors[5])
  rect(1.25,0,1.75,AIC.Logistic,col=colors[2])
  #rect(6.25,0,6.75,AIC.logNormal,col=colors[7])
  #rect(7.25,0,7.75,AIC.Tweedie,col=colors[8])
  #rect(8.25,0,8.75,AIC.Weib,col=colors[9])
  #text(0.175,2225,"(a)",font=2,cex=1.75)

  hist(Data, freq=FALSE,xlab=x_lab,ylim=c(y_lim_min,y_lim_max),col="white",main="",cex.lab=1,cex.axis=1,las=1)
  x<-seq(min(0,min(Data)),max(Data),0.01)

  legend("topright",paste("Fitted",Selected,"dist."),lty=1,col=colors[j],cex=0.9,bty='n')
  #text(5.4,0.4,"(b)",font=2,cex=1.75)

  if(Selected=="Gaus"){
    fit<-fitdistr(Data, "normal")
    lines(x,dnorm(x,fit$estimate[1],fit$estimate[2]),col=mypalette[4],lwd=2)
  }


  if(Selected=="Logis"){
    fit<-fitdistr(Data,"logistic")
    lines(x,dlogis(x,fit$estimate[1],fit$estimate[2]),col=mypalette[6],lwd=2)
  }

  plot(sort(Data),seq(1,length(Data),1)/(length(Data)),xlim=c(min(Data),max(Data)),ylim=c(0,1),xlab=x_lab,ylab="P(X<x)",main="",pch=16,cex.lab=1,cex.axis=1,las=1)
  x<-seq(min(0,min(Data)),max(Data),0.01)
  eta<-sqrt((1/length(Data))*log(2/0.95))
  lines(sort(Data),ifelse(seq(1,length(Data),1)/(length(Data))+eta>1,1,seq(1,length(Data),1)/(length(Data))+eta),col=1,lty=2)
  lines(sort(Data),ifelse(seq(1,length(Data),1)/(length(Data))-eta<0,0,seq(1,length(Data),1)/(length(Data))-eta),col=1,lty=2)
  legend("bottomright",c("95% Conf. Interval",paste("Fitted",Selected, "dist.")),lty=c(2,1),col=c(1,colors[j]),cex=1,bty='n',border = "white")
  #text(5.3,1,"(c)",font=2,cex=1.75)

   if(Selected=="Gaus"){
    fit<-fitdistr(Data,"normal")
    lines(x,pnorm(x,fit$estimate[1],fit$estimate[2]),col=mypalette[4],lwd=2)
  }


  if(Selected=="Logis"){
    fit<-fitdistr(Data,"logistic")
    lines(x,plogis(x,fit$estimate[1],fit$estimate[2]),col=mypalette[6],lwd=2)
  }
}
