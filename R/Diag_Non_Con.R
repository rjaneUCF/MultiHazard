#' Goodness of fit of non-extreme marginal distributions
#'
#' Fits two (unbounded) non-extreme marginal distributions to a dataset and returns three plots demonstrating their relative goodness of fit.
#' The distributions are the Gaussian \code{"Gaus"}, Gumbel \code{"Gum"}, Laplace \code{"Lapl"}, Logistic \code{"Logis"} and the reverse Gumbel \code{"RGum"}.

#' @param Data Numeric vector containing realizations of the variable of interest.
#' @param Omit Character vector specifying any distributions that are not to be tested. Default \code{"NA"}, all distributions are fit.
#' @param x_lab Character vector of length one specifying the label on the x-axis of histogram and cumulative distribution plot.
#' @param y_lim_min Numeric vector of length one specifying the lower y-axis limit of the histogram. Default is \code{0}.
#' @param y_lim_max Numeric vector of length one specifying the upper y-axis limit of the histogram. Default is \code{1}.
#' @return Dataframe \code{$AIC} giving the AIC associated with each distribution and the name of the best fitting distribution \code{$Best_fit}. Panel consisting of three plots. Upper plot: Plot depicting the AIC of the two fitted distributions. Middle plot: Probability Density Functions (PDFs) of the fitted distributions superimposed on a histogram of the data. Lower plot: Cumulative Distribution Functions (CDFs) of the fitted distributions overlaid on a plot of the empirical CDF.
#' @seealso \code{\link{Copula_Threshold_2D}}
#' @export
#' @examples
#' S20.Rainfall<-Con_Sampling_2D(Data_Detrend=S20.Detrend.df[,-c(1,4)],
#'                               Data_Declust=S20.Detrend.Declustered.df[,-c(1,4)],
#'                               Con_Variable="Rainfall",Thres=0.97)
#' Diag_Non_Con(Data=S20.Rainfall$Data$OsWL,x_lab="O-sWL (ft NGVD 29)",
#'              y_lim_min=0,y_lim_max=1.5)
Diag_Non_Con<-function(Data,Omit=NA,x_lab,y_lim_min=0,y_lim_max=1){

  # Check if Data is provided
  if (missing(Data)) {
    stop("Data parameter is required")
  }

  # Check if Data is numeric
  if (!is.numeric(Data)) {
    stop("Data must be numeric, got: ", class(Data)[1])
  }

  # Check if Data is vector or can be converted to vector
  if (is.matrix(Data) || is.data.frame(Data)) {
    if (ncol(Data) > 1) {
      warning("Data has multiple columns. Using first column only.")
      Data <- Data[, 1]
    } else {
      Data <- as.vector(Data)
    }
  }

  # Check for empty data
  if (length(Data) == 0) {
    stop("Data is empty")
  }

  # Check for all NA values
  if (all(is.na(Data))) {
    stop("Data contains only NA values")
  }

  # Check minimum sample size
  if (length(na.omit(Data)) < 10) {
    stop("Data must have at least 10 non-missing observations, got: ", length(na.omit(Data)))
  }

  # Remove NA values and warn if any exist
  original_length <- length(Data)
  Data <- na.omit(Data)
  if (length(Data) < original_length) {
    warning("Removed ", original_length - length(Data), " NA values from Data")
  }

  # Check for infinite values
  if (any(is.infinite(Data))) {
    inf_count <- sum(is.infinite(Data))
    warning("Data contains ", inf_count, " infinite values. Removing them.")
    Data <- Data[is.finite(Data)]
  }

  # Check if we still have enough data after cleaning
  if (length(Data) < 10) {
    stop("After removing NA/infinite values, data has fewer than 10 observations: ", length(Data))
  }

  # Check for constant data
  if (length(unique(Data)) == 1) {
    stop("Data is constant (all values are the same). Cannot fit distributions.")
  }

  # Check for very low variance
  if (var(Data) < 1e-10) {
    warning("Data has very low variance (", var(Data), "). Distribution fitting may be unreliable.")
  }

  # Validate Omit parameter
  valid_distributions <- c("Gaus", "Gum", "Lapl", "Logis", "RGum")
  if (!is.na(Omit[1])) {
    if (!all(Omit %in% valid_distributions)) {
      invalid_omit <- Omit[!Omit %in% valid_distributions]
      stop("Invalid distribution names in Omit: ", paste(invalid_omit, collapse=", "),
           ". Valid options are: ", paste(valid_distributions, collapse=", "))
    }
    if (length(Omit) >= length(valid_distributions)) {
      stop("Cannot omit all distributions. At least one distribution must be tested.")
    }
  }

  if (y_lim_min >= y_lim_max) {
    stop("y_lim_min must be less than y_lim_max, got: y_lim_min = ", y_lim_min, ", y_lim_max = ", y_lim_max)
  }

  #Plotting colors
  mypalette<-brewer.pal(9,"Set1")

  #Distributions to test
  Dist<-c("Gaus","Gum","Lapl","Logis","RGum")
  Dist.2 <- Dist
  Test<-1:5
  if(is.na(Omit[1])==F){
    Test<-Test[-which(Dist %in% Omit)]
    Dist.2 <- Dist.2[-which(Dist %in% Omit)]
  }

  #Second omit vector for distributions where paramter estimation fails
  Omit.2 = rep(NA,3)

  par(mfrow=c(3,1))
  par(mar=c(4.2,4.2,1,1))

  #AIC result objects
  AIC.Gaus <- NA
  AIC.Gum <- NA
  AIC.Lapl <- NA
  AIC.Logis <- NA
  AIC.RGum <- NA

  #
  if(any(Test==2)){
    fit <- tryCatch(gamlss(Data  ~ 1, family=GU),
                    error = function(e) "error")
   if(fit == "error") {
    Omit.2[1] = "Gum"
   } else {
    Omit.2[1] = ifelse(exp(fit$sigma.coefficients) < 0, "Gum", NA)
   }
  }

  if(any(Test==3)){
   fit <- tryCatch(fitdistr(Data, dlaplace, start=list(location=mean(Data),scale=sd(Data)/sqrt(2))),
                   error = function(e) "error")
   Omit.2[2] = ifelse(fit=="error","Lapl",NA)
  }

  if(any(Test==5)){
   fit <- tryCatch(gamlss(Data ~ 1,family=RG),
                   error = function(e) "error")
   if(fit == "error") {
    Omit.2[3] = "RGum"
   } else {
    Omit.2[3] = ifelse(exp(fit$sigma.coefficients) < 0, "RGum", NA)
   }
  }

  #Distributions to test
  if(any(is.na(Omit.2[1])==F)){
    Test<-Test[-which(Dist.2 %in% Omit.2)]
  }

  #AIC
  if(any(Test==1)){
    fit<-fitdistr(Data, "normal")
    AIC.Gaus<-2*length(fit$estimate)-2*fit$loglik
  }

  if(any(Test==2)){
    fit <- gamlss(Data  ~ 1, family= GU)
    AIC.Gum <-fit$aic
  }

  if(any(Test==3)){
    fit <- fitdistr(Data, dlaplace, start=list(location=mean(Data), scale=sd(Data)/sqrt(2)))
    AIC.Lapl<-2*length(coef(fit))-2*logLik(fit)
  }


  if(any(Test==4)){
    fit<-fitdistr(Data,"logistic")
    AIC.Logis<-2*length(fit$estimate)-2*fit$loglik
  }

  if(any(Test==5)){
    fit <- gamlss(Data ~ 1,family=RG)
    AIC.RGum<-fit$aic
  }


  plot(0,xlim=c(0,length(Test)),ylim=c(min(0,AIC.Gaus,AIC.Gum,AIC.Lapl,AIC.Logis,AIC.RGum,na.rm=T),max(0,AIC.Gaus,AIC.Gum,AIC.Logis,AIC.RGum,na.rm=T)),type='n',xlab="Probability Distribution",ylab="AIC",xaxt='n',cex.axis=1,cex.lab=1,las=1)
  axis(1,seq(0.5,length(Test)-0.5,1),c("Gaussian","Gumbel","Laplace","Logistic","Rev. Gumbel")[Test],cex.axis=0.71)
  if(any(Test==1)){rect(which(Test==1)-0.5-length(Test)/40,0,which(Test==1)-0.5+length(Test)/40,AIC.Gaus,col=mypalette[3])}
  if(any(Test==2)){rect(which(Test==2)-0.5-length(Test)/40,0,which(Test==2)-0.5+length(Test)/40,AIC.Gum,col=mypalette[2])}
  if(any(Test==3)){rect(which(Test==3)-0.5-length(Test)/40,0,which(Test==3)-0.5+length(Test)/40,AIC.Lapl,col=mypalette[4])}
  if(any(Test==4)){rect(which(Test==4)-0.5-length(Test)/40,0,which(Test==4)-0.5+length(Test)/40,AIC.Logis,col=mypalette[6])}
  if(any(Test==5)){rect(which(Test==5)-0.5-length(Test)/40,0,which(Test==5)-0.5+length(Test)/40,AIC.RGum,col=mypalette[1])}

  hist(Data, freq=FALSE,xlab=x_lab,col="white",main="",cex.lab=1,cex.axis=1,ylim=c(y_lim_min,y_lim_max),las=1)
  x<-seq(min(Data),max(Data),0.01)
  #text(5.35,0.1,"(f)",font=2,cex=1.75)

  if(any(Test==1)){
    fit<-fitdistr(Data, "normal")
    lines(x,dnorm(x,fit$estimate[1],fit$estimate[2]),col=mypalette[3],lwd=2)
  }

  if(any(Test==2)){
    fit <- gamlss(Data  ~ 1, family=GU)
    lines(x,dGU(x,fit$mu.coefficients,exp(fit$sigma.coefficients)),col=mypalette[2],lwd=2)
  }

  if(any(Test==3)){
    fit <- fitdistr(Data, dlaplace, start=list(location=mean(Data), scale=sd(Data)/sqrt(2)))
    lines(x,dlaplace(x,fit$estimate[1],fit$estimate[2]),col=mypalette[4],lwd=2)
  }

  if(any(Test==4)){
  fit<-fitdistr(Data,"logistic")
  lines(x,dlogis(x,fit$estimate[1],fit$estimate[2]),col=mypalette[6],lwd=2)
  }

  if(any(Test==5)){
  fit <- gamlss(Data ~ 1,family=RG)
  lines(x,dRG(x,fit$mu.coefficients,exp(fit$sigma.coefficients)),col=mypalette[1],lwd=2)
  }

  plot(sort(Data),seq(1,length(Data),1)/(length(Data)),ylim=c(0,1),xlab=x_lab,ylab="P(X<x)",main="",pch=16,cex.lab=1,cex.axis=1,las=1)
  x<-seq(min(Data),max(Data),0.01)
  eta<-sqrt((1/(2*length(Data)))*log(2/0.95))
  lines(sort(Data),ifelse(seq(1,length(Data),1)/(length(Data))+eta>1,1,seq(1,length(Data),1)/(length(Data))+eta),col=1,lty=2)
  lines(sort(Data),ifelse(seq(1,length(Data),1)/(length(Data))-eta<0,0,seq(1,length(Data),1)/(length(Data))-eta),col=1,lty=2)
  legend("bottomright",c("95% Conf. Interval","Fitted distributions."),lty=c(2,1),col=c(1,1),cex=1,bty='n',border = "white")
  #text(2,1,"(g)",font=2,cex=1.75)

  if(any(Test==1)){
    fit<-fitdistr(Data,"normal")
    lines(x,pnorm(x,fit$estimate[1],fit$estimate[2]),col=mypalette[3],lwd=2)
  }

  if(any(Test==2)){
    fit <- gamlss(Data  ~ 1, family=GU)
  lines(x,pGU(x,fit$mu.coefficients,exp(fit$sigma.coefficients)),col=mypalette[2],lwd=2)
  }

  if(any(Test==3)){
    fit <- fitdistr(Data, dlaplace, start=list(location=mean(Data), scale=sd(Data)/sqrt(2)))
    lines(x,plaplace(x,fit$estimate[1],fit$estimate[2]),col=mypalette[4],lwd=2)
  }

  if(any(Test==4)){
    fit<-fitdistr(Data,"logistic")
    lines(x,plogis(x,fit$estimate[1],fit$estimate[2]),col=mypalette[6],lwd=2)
  }

  if(any(Test==5)){
  fit <- gamlss(Data ~ 1,family=RG)
  lines(x,pRG(x,fit$mu.coefficients,exp(fit$sigma.coefficients)),col=mypalette[1],lwd=2)
  }

  AIC<-data.frame(c("Gaus","Gum","Lapl","Logis","RGum"),c(AIC.Gaus,AIC.Gum,AIC.Lapl,AIC.Logis,AIC.RGum))
  colnames(AIC)<-c("Distribution","AIC")
  Best_fit<-AIC$Distribution[which(AIC$AIC==min(AIC$AIC,na.rm=T))]
  res<-list("AIC"=AIC,"Best_fit"=Best_fit)
  return(res)
}


