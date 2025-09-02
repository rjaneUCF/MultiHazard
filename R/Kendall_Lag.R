#' Kendall's tau correlation coefficient between pairs of variables over a range of lags
#'
#' Kendall's tau correlation coefficient between pairs of up to three variables over a range of lags
#'
#' @param Data A data frame with 3 columns, containing concurrent observations of three time series.
#' @param Lags Integer vector giving the lags over which to calculate coefficient. Default is a vector from \code{-6} to \code{6}.
#' @param PLOT Logical; whether to show plot of Kendall's coefficient vs lag. Default is \code{TRUE}.
#' @param GAP Numeric vector of length one. Length of y-axis above and below max and min Kendall's tau values.
#' @return List comprising Kendall's tau coefficients between the variables pairs composing columns of Data with the specified lags applied to the second named variable \code{Values} and the p-values \code{Test} when testing the null hypothesis H_0: tau=0 i.e. there is no correlation between a pair of variables. Plot of the coefficient with a filled point of hypothesis test (p-value<0.05). Lag applied to variable named second in the legend.
#' @seealso \code{\link{Dataframe_Combine}}
#' @export
#' @examples
#' Kendall_Lag(Data=S20.Detrend.df,GAP=0.1)
Kendall_Lag<-function(Data,Lags=seq(-6,6,1),PLOT=TRUE,GAP=0.1){
  Lag<-function(x,k){
    if(k>0){
      return(c(rep(NA,k),x)[1:length(x)])
    } else{
      return(c(x[(-k+1):length(x)],rep(NA,-k)))
    }
  }

  correlation<-function(Data_Cor,lag){
    for(i in 2:(ncol(Data_Cor))){
      Data_Cor[,i]<-Lag(Data_Cor[,i],-lag[i-1])
    }
    Data_Cor<-na.omit(Data_Cor)
    return(cor(Data_Cor[2:ncol(Data_Cor)],method = "kendall")[which(lower.tri(cor(Data_Cor[2:ncol(Data_Cor)]))==T)])
  }
  n<-ncol(Data)-1
  if(n==3){
    Var1_Var2<-numeric(length(Lags))
    Var2_Var3<-numeric(length(Lags))
    Var1_Var3<-numeric(length(Lags))

    for(i in 1:length(Lags)){
      Var1_Var2[i]<-correlation(Data_Cor=Data,c(Lags[i],0,0))[1]
      Var2_Var3[i]<-correlation(Data_Cor=Data,c(0,0,Lags[i]))[2]
      Var1_Var3[i]<-correlation(Data_Cor=Data,c(0,0,Lags[i]))[3]
    }

    correlation.test<-function(Data_Cor,lag){
      for(i in 2:(ncol(Data_Cor))){
        Data_Cor[,i]<-Lag(Data_Cor[,i],-lag[i-1])
      }
      Data_Cor<-na.omit(Data_Cor)
      return(cor.test(Data_Cor[,2],Data_Cor[,3],method = "kendall")$p.value)
    }

    Var1_Var2_Test<-numeric(length(Lags))
    Var2_Var3_Test<-numeric(length(Lags))
    Var1_Var3_Test<-numeric(length(Lags))

    for(i in 1:length(Lags)){
      Var1_Var2_Test[i]<-correlation.test(Data_Cor=Data[,-4],c(Lags[i],0))
      Var2_Var3_Test[i]<-correlation.test(Data_Cor=Data[,-3],c(0,Lags[i]))
      Var1_Var3_Test[i]<-correlation.test(Data_Cor=Data[,-2],c(0,Lags[i]))
    }

    if(PLOT==TRUE){
      yx<-max(c(Var1_Var2,Var2_Var3,Var1_Var3))-min(c(Var1_Var2,Var2_Var3,Var1_Var3))
      plot(Lags,Var1_Var2,ylim=c(min(c(Var1_Var2,Var2_Var3,Var1_Var3))-GAP*yx,max(c(Var1_Var2,Var2_Var3,Var1_Var3)))+GAP*yx,type='l',xlab="Lag (days)",ylab=expression(paste("Kendall's "*tau*' coefficient')),cex.lab=1.65,cex.axis=1.65,lwd=2.5)
      abline(h=0,lty=2)
      lines(Lags,Var2_Var3,col=2)
      lines(Lags,Var1_Var3,col=3)

      points(Lags,Var1_Var2,pch=ifelse(Var1_Var2_Test<0.05,16,21),bg=ifelse(Var1_Var2_Test<0.05,1,"White"),cex=1.5)
      points(Lags,Var2_Var3,pch=ifelse(Var2_Var3_Test<0.05,16,21),col=2,bg=ifelse(Var2_Var3_Test<0.05,2,"White"),cex=1.5)
      points(Lags,Var1_Var3,pch=ifelse(Var1_Var3_Test<0.05,16,21),col=3,bg=ifelse(Var1_Var3_Test<0.05,3,"White"),cex=1.5)
      legend("topright",c(paste(colnames(Data)[2],"_",colnames(Data)[3],sep=""),
                          paste(colnames(Data)[3],"_",colnames(Data)[4],sep=""),
                          paste(colnames(Data)[2],"_",colnames(Data)[4],sep="")),
             bty="n",lwd=2.5,col=c(1,2,3))
    }
    Value<-list()
    Value[[paste(names(Data)[2],'_',names(Data)[3],sep="")]]= Var1_Var2
    Value[[paste(names(Data)[3],'_',names(Data)[4],sep="")]]= Var2_Var3
    Value[[paste(names(Data)[2],'_',names(Data)[4],sep="")]]= Var1_Var3

    Test<-list()
    Test[[paste(names(Data)[2],'_',names(Data)[3],'_Test',sep="")]]= Var1_Var2_Test
    Test[[paste(names(Data)[3],'_',names(Data)[4],'_Test',sep="")]]= Var2_Var3_Test
    Test[[paste(names(Data)[2],'_',names(Data)[4],'_Test',sep="")]]= Var1_Var3_Test
  }
  if(n==2){
    Var1_Var2 <- numeric(length(Lags))
    for (i in 1:length(Lags)) {
      Var1_Var2[i] <- correlation(Data_Cor = Data, c(Lags[i],0))[1]
    }
    correlation.test <- function(Data_Cor, lag) {
      for (i in 2:(ncol(Data_Cor))) {
        Data_Cor[, i] <- Lag(Data_Cor[, i], -lag[i - 1])
      }
      Data_Cor <- na.omit(Data_Cor)
      return(cor.test(Data_Cor[, 2], Data_Cor[, 3], method = "kendall")$p.value)
    }
    Var1_Var2_Test <- numeric(length(Lags))
    for (i in 1:length(Lags)) {
      Var1_Var2_Test[i] <- correlation.test(Data_Cor = Data[,-4], c(Lags[i], 0))
    }
    if (PLOT == TRUE) {
      yx <- max(c(Var1_Var2)) - min(c(Var1_Var2))
      plot(Lags, Var1_Var2, ylim = c(min(c(Var1_Var2)) - GAP * yx, max(c(Var1_Var2)) + GAP * yx), type = "l", xlab = "Lag (days)",
           ylab = expression(paste("Kendall's " * tau * " coefficient")),
           cex.lab = 1.65, cex.axis = 1.65, lwd = 2.5)
      abline(h = 0, lty = 2)
      points(Lags, Var1_Var2, pch = ifelse(Var1_Var2_Test <  0.05, 16, 21), bg = ifelse(Var1_Var2_Test < 0.05,1, "White"), cex = 1.5)
      legend("topright", c(paste(colnames(Data)[2], "_", colnames(Data)[3],sep = "")), bty = "n", lwd = 2.5, col = c(1, 2, 3))
    }
    Value<-list()
    Value[[paste(names(Data)[2],'_',names(Data)[3],sep="")]]= Var1_Var2

    Test<-list()
    Test[[paste(names(Data)[2],'_',names(Data)[3],'_Test',sep="")]]= Var1_Var2_Test
  }
  res<-list("Value" = Value,"Test" = Test)
  return(res)
}
