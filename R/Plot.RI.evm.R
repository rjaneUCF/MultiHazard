#' @noRd
Plot.RI.evm<-function(Data,model,mu,min.RI=1,max.RI=1000,y.lab="Return Level",main.RI="",Cex.axis=1,Cex.lab=1){
  RI.Estimate<-numeric(length(seq(min.RI,max.RI,0.1)))
  RI.Lower<-numeric(length(seq(min.RI,max.RI,0.1)))
  RI.Upper<-numeric(length(seq(min.RI,max.RI,0.1)))
  for(i in 1:length(seq(min.RI,max.RI,0.1))){
    RI.Estimate[i]<-rl(model, M = mu*seq(min.RI,max.RI,0.1)[i],ci.fit=TRUE)[[1]][1]
    RI.Lower[i]<-rl(model, M = mu*seq(min.RI,max.RI,0.1)[i],ci.fit=TRUE)[[1]][2]
    RI.Upper[i]<-rl(model, M = mu*seq(min.RI,max.RI,0.1)[i],ci.fit=TRUE)[[1]][3]
  }
  plot(log10(seq(min.RI,max.RI,0.1)),RI.Estimate,type='l',xlim=c(log10(min.RI),log10(max.RI)),ylim=c(min(RI.Lower),max(RI.Upper)),xlab=expression('log'[10]*'(Return Period)'),ylab=y.lab,main=main.RI,cex.axis=Cex.axis,cex.lab=Cex.lab)
  lines(log10(seq(min.RI,max.RI,0.1)),RI.Upper,col=4)
  lines(log10(seq(min.RI,max.RI,0.1)),RI.Lower,col=4)
  Data<-Data[which(Data>model$threshold)]
  points(log10((length(Data) + 1) / ((1:length(Data) * model$rate * mu)))[1:length(which((length(Data) + 1) / ((1:length(Data) * model$rate* mu))>min.RI))],
         rev(sort(Data))[1:length(which((length(Data) + 1) / ((1:length(Data) * model$rate * mu))>min.RI))],pch=16)
}