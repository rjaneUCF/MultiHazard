#' @noRd
GPD_diag_HT04<-function(Data,Data_Full,model,param,thres,min.RI=min.RI,max.RI=100,xlab.hist="Data",y.lab="Return level"){
  xi<-param[2]
  sigma<-param[1]

  par(mfrow=c(2,2))
  par(mar=c(4.5,4.2,1,1))
  #Probability plot
  Excess<-Data[Data>thres]-thres
  plot( 1:length(Data[Data>thres])/(length(Data[Data>thres])+1), pgpd(sort(Excess),sigma,xi,u=0), xlab = "Empirical", ylab = "Model",main="Probability Plot",pch=16)
  lines(seq(0,1,0.01),seq(0,1,0.01),col=4)

  #Quantile plot
  Emp<-qgpd(1:length(Data[Data>thres])/(length(Data[Data>thres])+1),sigma,xi,u=thres)
  plot(thres+sort(Excess),Emp,xlim=c(min(Emp,thres+sort(Data[Data>thres])),max(Emp,thres+sort(Data[Data>thres]))),ylim=c(min(Emp,thres+sort(Data[Data>thres])),max(Emp,thres+sort(Data[Data>thres]))),xlab = "Empirical", ylab = "Model",main="Quantile Plot",pch=16)
  lines(c(min(Emp,thres+sort(Data[Data>thres])),max(Emp,thres+sort(Data[Data>thres]))),c(min(Emp,thres+sort(Data[Data>thres])),max(Emp,thres+sort(Data[Data>thres]))),col=4)

  #Histogram
  hist(Data[Data>thres],freq=F,xlab=xlab.hist,main="Density Plot")
  density<-dgpd(seq(thres,max(Data),0.1),sigma,xi,u=thres)
  lines(seq(thres,max(Data),0.1),density,col=4)

  #Return Level
  Plot.RI.evm(Data,model,min.RI,max.RI,y.lab=y.lab,main.RI="Return Level Plot")
}

