#' Goodness of fit of non-extreme marginal distributions
#'
#' Fits ten (truncated) non-extreme marginal distributions to a dataset and returns three plots demonstrating their relative goodness of fit.
#' The distributions are the Birnbaum-Saunders \code{"BS"}, exponential \code{"Exp"}, two-parameter gamma \code{"Gam(2)"}, three-parameter gamma \code{"Gam(3)"}, mixed two-parameter gamma \code{"GamMix(2)"}, mixed three-parameter gamma \code{"GamMix(3)"}, lognormal \code{"LNorm"}, truncated normal \code{"TNorm"}, Tweedie \code{"Twe"} and the Weibull \code{"Weib"}.
#' @param Data Numeric vector containing realizations of the variable of interest.
#' @param Omit Character vector specifying any distributions that are not to be tested. Default \code{"NA"}, all distributions are fit.
#' @param x_lab Character vector of length one specifying the label on the x-axis of histogram and cumulative distribution plot.
#' @param y_lim_min Numeric vector of length one specifying the lower y-axis limit of the histogram. Default is \code{0}.
#' @param y_lim_max Numeric vector of length one specifying the upper y-axis limit of the histogram. Default is \code{1}.
#' @return Dataframe \code{$AIC} giving the AIC associated with each distribution and the name of the best fitting distribution \code{$Best_fit}. Panel consisting of three plots. Upper plot: Plot depicting the AIC of the ten fitted distributions. Middle plot: Probability Density Functions (PDFs) of the fitted distributions superimposed on a histogram of the data. Lower plot: Cumulative Distribution Functions (CDFs) of the fitted distributions overlaid on a plot of the empirical CDF.
#' @seealso \code{\link{Copula_Threshold_2D}}
#' @export
#' @examples
#' S20.OsWL<-Con_Sampling_2D(Data_Detrend=S20.Detrend.df[,-c(1,4)],
#'                          Data_Declust=S20.Detrend.Declustered.df[,-c(1,4)],
#'                          Con_Variable="OsWL",Thres=0.97)
#' Diag_Non_Con_Trunc(Data=S20.OsWL$Data$Rainfall,x_lab="Rainfall (Inches)",
#'                    y_lim_min=0,y_lim_max=2)
Diag_Non_Con_Trunc<-function(Data,Omit=NA,x_lab="Data",y_lim_min=0,y_lim_max=1){
  
  #Colors for plots
  mypalette<-c("Black",brewer.pal(9,"Set1"))
  
  #Distributions to test
  Dist<-c("BS","Exp","Gam(2)","Gam(3)","GamMix(2)","GamMix(3)","LNorm","TNorm","Twe","Weib")
  Test<-1:10
  if(is.na(Omit)==F){
    Test<-Test[-which(Dist %in% Omit)]
  }
  
  #AIC result objects
  AIC.BS<-NA
  AIC.Exp<-NA
  AIC.Gam2<-NA
  AIC.Gam3<-NA
  AIC.GamMix2<-NA
  AIC.GamMix3<-NA
  AIC.logNormal<-NA
  AIC.TNormal<-NA
  AIC.Tweedie<-NA
  AIC.Weib<-NA
  
  #AIC
  if(any(Test==1)){
    bdata2 <- data.frame(shape = exp(-0.5), scale = exp(0.5))
    bdata2 <- transform(bdata2, y = Data)
    fit.BS <- vglm(y ~ 1, bisa, data = bdata2, trace = FALSE)
    AIC.BS<-2*length(coef(fit.BS))-2*logLik(fit.BS)
  }
  if(any(Test==2)){
    fit.Exp<-fitdistr(Data,"exponential")
    AIC.Exp<-2*length(fit.Exp$estimate)-2*fit.Exp$loglik
  }
  if(any(Test==3)){
    fit.Gam2<-fitdistr(Data, "gamma")
    AIC.Gam2<-2*length(fit.Gam2$estimate)-2*fit.Gam2$loglik
  }
  data.gamlss <- data.frame(X=Data)
  if(any(Test==4)){
    ### 3-parameter gamma dist.
    for(i in 1:100){
      fit.Gamma3 <- tryCatch(gamlss(X~1, data=data.gamlss, family=GG),
                             error = function(e) "error")
      if( is.character(fit.Gamma3) ) next
      if( !is.character(fit.Gamma3) ) break
    }
    if( is.character(fit.Gamma3) ){
      #AIC.Gamma3 <- -9999
      Test <- Test[-which(Test==4)]
    }else{
      AIC.Gam3 <- fit.Gamma3$aic
    }
  }
  if(any(Test==5)){
    ### 2 mixture-gamma dist.
    for(i in 1:100){
      fit.GamMIX2_GA <- tryCatch(gamlssMX(X~1, data=data.gamlss, family=GA, K=2),
                                 error = function(e) "error")
      if( is.character(fit.GamMIX2_GA) ) next
      if( !is.character(fit.GamMIX2_GA) ) break
    }
    if( is.character(fit.GamMIX2_GA) ){
      #AIC.GamMIX2_GA <- -9999
      Test <- Test[-which(Test==5)]
    }else{
      AIC.GamMix2 <- fit.GamMIX2_GA$aic
    }
  }
  if(any(Test==6)){
    ### 3 mixture-gamma dist.
    for(i in 1:100){
      fit.GamMIX3_GA <- tryCatch(gamlssMX(X~1, data=data.gamlss, family=GA, K=3),
                                 error = function(e) "error")
      if( is.character(fit.GamMIX3_GA) ) next
      if( !is.character(fit.GamMIX3_GA) ) break
    }
    if( is.character(fit.GamMIX3_GA) ){
      #AIC.GamMIX3_GA <- -9999
      Test <- Test[-which(Test==6)]
    }else{
      AIC.GamMix3 <- fit.GamMIX3_GA$aic
    }
  }
  #fit<-fitdist(Data, "invgauss", start = list(mean = 5, shape = 1))
  #AIC.InverseNormal<-2*length(fit$estimate)-2*fit$loglik
  if(any(Test==7)){
    fit.LNorm<-fitdistr(Data,"lognormal")
    AIC.logNormal<-2*length(fit.LNorm$estimate)-2*fit.LNorm$loglik
  }
  if(any(Test==8)){
    fit.TNorm <- fitdistr(Data, "normal")
    AIC.TNormal <- 2 * length(fit.TNorm$estimate) - 2 * fit.TNorm$loglik
  }
  if(any(Test==9)){
    fit.Twe <- tweedie.profile(Data ~ 1,
                               p.vec=seq(1.5, 2.5, by=0.2), do.plot=FALSE)
    AIC.Tweedie<-2*3-2*fit.Twe$L.max
  }
  if(any(Test==10)){
    fit.Weib<-fitdistr(Data,"weibull")
    AIC.Weib<-2*length(fit.Weib$estimate)-2*fit.Weib$loglik
  }
  
  #Plot layout
  par(mfrow=c(3,1))
  par(mar=c(4.2,4.2,1,1))
  
  #Plotting AIC
  plot(0,xlim=c(0,length(Test)),ylim=c(min(0,AIC.BS,AIC.Exp,AIC.Gam2,AIC.Gam3,AIC.GamMix2,AIC.GamMix3,AIC.logNormal,AIC.TNormal,AIC.Tweedie,AIC.Weib,na.rm=T),max(0,AIC.BS,AIC.Exp,AIC.Gam2,AIC.Gam3,AIC.GamMix2,AIC.GamMix3,AIC.TNormal,AIC.logNormal,AIC.Tweedie,AIC.Weib,na.rm=T)),type='n',xlab="Probability Distribution",ylab="AIC",xaxt='n',cex.axis=1,cex.lab=1,las=1)
  axis(1,seq(0.5,length(Test)-0.5,1),c("Birn-S","Exp","Gam(2)","Gam(3)","GamMix(2)","GamMix(3)","LogN","TNorm","Twe","Weib")[Test],cex.axis=0.71)
  if(any(Test==1)){rect(which(Test==1)-0.5-length(Test)/40,0,which(Test==1)-0.5+length(Test)/40,AIC.BS,col=mypalette[1])}
  if(any(Test==2)){rect(which(Test==2)-0.5-length(Test)/40,0,which(Test==2)-0.5+length(Test)/40,AIC.Exp,col=mypalette[2])}
  if(any(Test==3)){rect(which(Test==3)-0.5-length(Test)/40,0,which(Test==3)-0.5+length(Test)/40,AIC.Gam2,col=mypalette[3])}
  if(any(Test==4)){rect(which(Test==4)-0.5-length(Test)/40,0,which(Test==4)-0.5+length(Test)/40,AIC.Gam3,col=mypalette[4])}
  if(any(Test==5)){rect(which(Test==5)-0.5-length(Test)/40,0,which(Test==5)-0.5+length(Test)/40,AIC.GamMix2,col=mypalette[5])}
  if(any(Test==6)){rect(which(Test==6)-0.5-length(Test)/40,0,which(Test==6)-0.5+length(Test)/40,AIC.GamMix3,col=mypalette[6])}
  if(any(Test==7)){rect(which(Test==7)-0.5-length(Test)/40,0,which(Test==7)-0.5+length(Test)/40,AIC.logNormal,col=mypalette[7])}
  if(any(Test==8)){rect(which(Test==8)-0.5-length(Test)/40,0,which(Test==8)-0.5+length(Test)/40,AIC.TNormal,col=mypalette[8])}
  if(any(Test==9)){rect(which(Test==9)-0.5-length(Test)/40,0,which(Test==9)-0.5+length(Test)/40,AIC.Tweedie,col=mypalette[9])}
  if(any(Test==10)){rect(which(Test==10)-0.5-length(Test)/40,0,which(Test==10)-0.5+length(Test)/40,AIC.Weib,col=mypalette[10])}
  
  #Colors
  
  hist(Data, freq=FALSE,xlab=x_lab,col="white",main="",cex.lab=1,cex.axis=1,ylim=c(y_lim_min,y_lim_max),las=1)
  x<-seq(min(Data),max(Data),0.01)
  #text(5.35,0.1,"(f)",font=2,cex=1.75)
  
  if(any(Test==1)){
    lines(x,dbisa(x,Coef(fit.BS)[1],Coef(fit.BS)[2]),col=mypalette[1],lwd=2)
  }
  if(any(Test==2)){
    lines(x,dexp(x,fit.Exp$estimate[1]),col=mypalette[2],lwd=2)
  }
  if(any(Test==3)){
    lines(x,dgamma(x,fit.Gam2$estimate[1],fit.Gam2$estimate[2]),col=mypalette[3],lwd=2)
  }
  
  if(any(Test==4)){
    lines(x,dGG(x, mu=exp(fit.Gamma3$mu.coefficients), sigma=exp(fit.Gamma3$sigma.coefficients), nu=fit.Gamma3$nu.coefficients),col=mypalette[4],lwd=2)
  }
  
  if(any(Test==5)){
    prob.MX1 <- round(fit.GamMIX2_GA$prob[1],3)
    prob.MX2 <- 1 - prob.MX1
    lines(x,dMX(x, mu=list(mu1=exp(fit.GamMIX2_GA$models[[1]]$mu.coefficients), mu2=exp(fit.GamMIX2_GA$models[[2]]$mu.coefficients)),
                sigma=list(sigma1=exp(fit.GamMIX2_GA$models[[1]]$sigma.coefficients), sigma2=exp(fit.GamMIX2_GA$models[[2]]$sigma.coefficients)),
                pi = list(pi1=prob.MX1, pi2=prob.MX2), family=list(fam1="GA", fam2="GA")),col=mypalette[5],lwd=2)
  }
  
  if(any(Test==6)){
    prob.MX1 <- round(fit.GamMIX3_GA$prob[1],3)
    prob.MX2 <- round(fit.GamMIX3_GA$prob[2],3)
    prob.MX3 <- 1 - prob.MX1 - prob.MX2 
    lines(x,dMX(x, mu=list(mu1=exp(fit.GamMIX3_GA$models[[1]]$mu.coefficients), mu2=exp(fit.GamMIX3_GA$models[[2]]$mu.coefficients), mu3=exp(fit.GamMIX3_GA$models[[3]]$mu.coefficients)),
                sigma=list(sigma1=exp(fit.GamMIX3_GA$models[[1]]$sigma.coefficients), sigma2=exp(fit.GamMIX3_GA$models[[2]]$sigma.coefficients), sigma3=exp(fit.GamMIX3_GA$models[[3]]$sigma.coefficients)),
                pi = list(pi1=prob.MX1, pi2=prob.MX2, pi3=prob.MX3), family=list(fam1="GA", fam2="GA", fam3="GA")),col=mypalette[6],lwd=2)
  }
  
  
  if(any(Test==7)){
    lines(x,dlnorm(x,fit.LNorm$estimate[1],fit.LNorm$estimate[2]),col=mypalette[7],lwd=2)
  }
  if(any(Test==8)){
    lines(x, dtruncnorm(x, a = min(Data), mean = fit.TNorm$estimate[1],
                        sd = fit.TNorm$estimate[2]), col = mypalette[8], lwd = 2)
  }
  if(any(Test==9)){
    lines(x,dtweedie(x,  power=fit.Twe$p.max, mu=mean(Data), phi=fit.Twe$phi.max),col=mypalette[9],lwd=2)
  }
  
  if(any(Test==10)){
    lines(x,dweibull(x,fit.Weib$estimate[1],fit.Weib$estimate[2]),col=mypalette[10],lwd=2)
  }
  
  plot(sort(Data),seq(1,length(Data),1)/(length(Data)),ylim=c(0,1),xlab=x_lab,ylab="P(X<x)",main="",pch=16,cex.lab=1,cex.axis=1,las=1)
  x<-seq(min(Data),max(Data),0.01)
  eta<-sqrt((1/(2*length(Data)))*log(2/0.95))
  lines(sort(Data),ifelse(seq(1,length(Data),1)/(length(Data))+eta>1,1,seq(1,length(Data),1)/(length(Data))+eta),col=1,lty=2)
  lines(sort(Data),ifelse(seq(1,length(Data),1)/(length(Data))-eta<0,0,seq(1,length(Data),1)/(length(Data))-eta),col=1,lty=2)
  legend("bottomright",c("95% Conf. Interval","Fitted distributions"),lty=c(2,1),col=c(1,4),cex=1,bty='n',border = "white")
  #text(2,1,"(g)",font=2,cex=1.75)
  if(any(Test==1)){
    lines(x,pbisa(x,Coef(fit.BS)[1],Coef(fit.BS)[2]),col=mypalette[1],lwd=2)
  }
  if(any(Test==2)){
    lines(x,pexp(x,fit.Exp$estimate[1]),col=mypalette[2],lwd=2)
  }
  if(any(Test==3)){
    lines(x,pgamma(x,fit.Gam2$estimate[1],fit.Gam2$estimate[2]),col=mypalette[3],lwd=2)
  }
  
  if(any(Test==4)){
    lines(x,pGG(x, mu=exp(fit.Gamma3$mu.coefficients), sigma=exp(fit.Gamma3$sigma.coefficients), nu=fit.Gamma3$nu.coefficients),col=mypalette[4],lwd=2)
  }
  
  if(any(Test==5)){
    prob.MX1 <- round(fit.GamMIX2_GA$prob[1],3)
    prob.MX2 <- 1 - prob.MX1
    lines(x,pMX(x, mu=list(mu1=exp(fit.GamMIX2_GA$models[[1]]$mu.coefficients), mu2=exp(fit.GamMIX2_GA$models[[2]]$mu.coefficients)),
                sigma=list(sigma1=exp(fit.GamMIX2_GA$models[[1]]$sigma.coefficients), sigma2=exp(fit.GamMIX2_GA$models[[2]]$sigma.coefficients)),
                pi = list(pi1=prob.MX1, pi2=prob.MX2), family=list(fam1="GA", fam2="GA")),col=mypalette[5],lwd=2)
  }
  
  if(any(Test==6)){
    prob.MX1 <- round(fit.GamMIX3_GA$prob[1],3)
    prob.MX2 <- round(fit.GamMIX3_GA$prob[2],3)
    prob.MX3 <- 1 - prob.MX1 - prob.MX2 
    lines(x,pMX(x, mu=list(mu1=exp(fit.GamMIX3_GA$models[[1]]$mu.coefficients), mu2=exp(fit.GamMIX3_GA$models[[2]]$mu.coefficients), mu3=exp(fit.GamMIX3_GA$models[[3]]$mu.coefficients)),
                sigma=list(sigma1=exp(fit.GamMIX3_GA$models[[1]]$sigma.coefficients), sigma2=exp(fit.GamMIX3_GA$models[[2]]$sigma.coefficients), sigma3=exp(fit.GamMIX3_GA$models[[3]]$sigma.coefficients)),
                pi = list(pi1=prob.MX1, pi2=prob.MX2, pi3=prob.MX3), family=list(fam1="GA", fam2="GA", fam3="GA")),col=mypalette[6],lwd=2)
  }
  
  if(any(Test==7)){
    lines(x,plnorm(x,fit.LNorm$estimate[1],fit.LNorm$estimate[2]),col=mypalette[7],lwd=2)
  }
  if(any(Test==8)){
    lines(x, ptruncnorm(x, a = min(Data), mean = fit.TNorm$estimate[1],
                        sd = fit.TNorm$estimate[2]), col = mypalette[8], lwd = 2)
  }
  if(any(Test==9)){
    lines(x,ptweedie(x,  power=fit.Twe$p.max, mu=mean(Data), phi=fit.Twe$phi.max),col=mypalette[9],lwd=2,pch=16,ylab="P(X<x)")
  }
  if(any(Test==10)){
    lines(x,pweibull(x,fit.Weib$estimate[1],fit.Weib$estimate[2]),col=mypalette[10],lwd=2)
  }
  
  AIC<-data.frame(c("BS","Exp","Gam2","Gam3","GamMix2","GamMix3","LogN","TNorm","Twe","Weib")[Test],c(AIC.BS,AIC.Exp,AIC.Gam2,AIC.Gam3,AIC.GamMix2,AIC.GamMix3,AIC.logNormal,AIC.TNormal,AIC.Tweedie,AIC.Weib)[Test])
  colnames(AIC)<-c("Distribution","AIC")
  Best_fit<-AIC$Distribution[which(AIC$AIC==min(AIC$AIC))]
  res<-list("AIC"=AIC, "Best_fit"=Best_fit)
  return(res)
}