#' Mean excess plot - GPD threshold selection
#'
#' The empirical mean excess function is linear in the case of a GPD.
#'
#' @param Data A vector comprising a declustered and if necessary detrended time series to be modelled.
#' @return Plot of the empirical mean excess function (black line), average of all observations exceeding a threshold decreased by the threshold, for thresholds spanning the range of the observations. Also provided are \code{95\%} confidence intervals (blue dotted lines) and the observations (black dots).
#' @seealso \code{\link{Decluster}} \code{\link{Detrend}}
#' @export
#' @examples
#' Mean_Excess_Plot(Data=S20.Detrend.Declustered.df$Rainfall)
Mean_Excess_Plot<-function(Data){
  if(length(which(is.na(Data)==TRUE))>0){
    thres<-seq(sort(na.omit(Data),decreasing = T)[length(na.omit(Data))],sort(na.omit(Data),decreasing = T)[2],0.1)
    Data<-na.omit(Data)
  } else{
    thres<-seq(sort(Data,decreasing = T)[length(Data)],sort(Data,decreasing = T)[2],0.1)
  }
  MeanExcess<-numeric(length(thres))
  MeanExcess.Lower<-numeric(length(thres))
  MeanExcess.Upper<-numeric(length(thres))
  for(i in 1:length(MeanExcess)){
    MeanExcess[i]<-mean(Data[Data>=thres[i]]-thres[i])
    MeanExcess.Lower[i]<-mean(Data[Data>=thres[i]]-thres[i])-1.96*sqrt(sd(Data[Data>=thres[i]])/length(which(Data>=thres[i])))
    MeanExcess.Upper[i]<-mean(Data[Data>=thres[i]]-thres[i])+1.96*sqrt(sd(Data[Data>=thres[i]])/length(which(Data>=thres[i])))
  }
  y.range<-range(c(min(MeanExcess.Lower),max(MeanExcess.Upper)))
  y.lim<-y.range+c(-1,1)*diff(y.range)/10
  plot(thres,MeanExcess,type='l',xlim=c(min(thres),max(thres)),ylim=y.lim)
  lines(thres,MeanExcess.Lower,lty=2,col=4)
  lines(thres,MeanExcess.Upper,lty=2,col=4)
  points(Data,rep(min(y.lim),length(Data)),pch=16,cex=0.75)
}
