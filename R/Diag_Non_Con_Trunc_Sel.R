#' Goodness of fit of the selected non-extreme marginal distribution
#'
#' Plots demonstrating the goodness of fit of a selected (truncated) non-extreme marginal distribution to a dataset.
#'
#' @param Data Numeric vector containing realizations of the variable of interest.
#' @param Selected Character vector of length one specifying the chosen distribution, options are the Birnbaum-Saunders \code{"BS"}, exponential \code{"Exp"}, two-parameter gamma \code{"Gam(2)"}, three-parameter gamma \code{"Gam(3)"}, mixed two-parameter gamma \code{"GamMix(2)"}, mixed three-parameter gamma \code{"GamMix(3)"}, lognormal \code{"LNorm"}, Tweedie \code{"Twe"} and Weibull \code{"Weib"}.
#' @param Omit Character vector specifying any distributions that are not to be tested. Default \code{"NA"}, all distributions are fit.
#' @param x_lab Character vector of length one specifying the label on the x-axis of histogram and cummulative distribution plot.
#' @param y_lim_min Numeric vector of length one specifying the lower y-axis limit of the histogram. Default is \code{0}.
#' @param y_lim_max Numeric vector of length one specifying the upper y-axis limit of the histogram. Default is \code{1}.
#' @return Panel consisting of three plots. Upper plot: Plot depicting the AIC of the eight fitted distributions. Middle plot: Probability Density Functions (PDFs) of the fitted distributions superimposed on a histogram of the data. Lower plot: Cumulative Distribution Functions (CDFs) of the fitted distributions overlaid on a plot of the empirical CDF.
#' @seealso \code{\link{Diag_Non_Con_Trunc}}
#' @export
#' @examples
#' S20.OsWL<-Con_Sampling_2D(Data_Detrend=S20.Detrend.df[,-c(1,4)],
#'                           Data_Declust=S20.Detrend.Declustered.df[,-c(1,4)],
#'                           Con_Variable="OsWL",Thres=0.97)
#' S20.OsWL$Data$Rainfall <- S20.OsWL$Data$Rainfall + runif(length(S20.OsWL$Data$Rainfall),0.001,0.01)
#' Diag_Non_Con_Trunc(Data=S20.OsWL$Data$Rainfall,x_lab="Rainfall (Inches)",
#'                    y_lim_min=0,y_lim_max=2)
#' Diag_Non_Con_Trunc_Sel(Data=S20.OsWL$Data$Rainfall,x_lab="Rainfall (Inches)",
#'                        y_lim_min=0,y_lim_max=2,Selected="Twe")
Diag_Non_Con_Trunc_Sel<-function(Data,Selected,Omit=NA,x_lab="Data",y_lim_min=0,y_lim_max=1){

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

  # Check Selected parameter is specified
  if (missing(Selected)) {
    stop("argument \"Selected\" is missing, with no default")
  }

  # Validate Omit parameter
  valid_distributions <- c("BS","Exp","Gam(2)","Gam(3)","GamMix(2)","GamMix(3)","LNorm","TNorm","Twe","Weib")
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

  # Validate Selected parameter
  valid_distributions <- c("BS","Exp","Gam(2)","Gam(3)","GamMix(2)","GamMix(3)","LNorm","TNorm","Twe","Weib")
  if (!is.na(Selected[1])) {
    if (!all(Selected %in% valid_distributions)) {
      invalid_Selected <- Selected[!Selected %in% valid_distributions]
      stop("Invalid distribution names in Selected: ", paste(invalid_Selected, collapse=", "),
           ". Valid options are: ", paste(valid_distributions, collapse=", "))
    }
    if (length(Selected) >  1){
      stop("Cannot select more than one distribution.")
    }
  }

  if (y_lim_min >= y_lim_max) {
    stop("y_lim_min must be less than y_lim_max, got: y_lim_min = ", y_lim_min, ", y_lim_max = ", y_lim_max)
  }


  #Check gamlss packages are installed
  if (!requireNamespace("gamlss", quietly = TRUE)) {
    stop("The 'gamlss' package is required but not installed.")
  }
  if (!requireNamespace("gamlss.dist", quietly = TRUE)) {
    stop("The 'gamlss.dist' package is required but not installed.")
  }

  # Load density and cumulative functions safely
  GU <- get("GG", envir = asNamespace("gamlss.dist"))
  GA <- get("GA", envir = asNamespace("gamlss.dist"))
  dGG <- get("dGG", envir = asNamespace("gamlss.dist"))
  pGG <- get("pGG", envir = asNamespace("gamlss.dist"))
  dGA <- get("dGA", envir = asNamespace("gamlss.dist"))
  pGA <- get("pGA", envir = asNamespace("gamlss.dist"))

  #Colors for plots
  mypalette<-c("Black",brewer.pal(9,"Set1"))

  #Distributions to test
  Dist<- c("BS","Exp","Gam(2)","Gam(3)","GamMix(2)","GamMix(3)","LNorm","TNorm","Twe","Weib")
  Test<- Dist
  if(!is.na(Omit)){
    Test<-Test[-which(Test %in% Omit)]
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
  if(any(Test=="BS")){
    bdata2 <- data.frame(shape = exp(-0.5), scale = exp(0.5))
    bdata2 <- transform(bdata2, y = Data)
    fit <- vglm(y ~ 1, bisa, data = bdata2, trace = FALSE)
    AIC.BS<-2*length(coef(fit))-2*logLik(fit)
  }

  if(any(Test=="Exp")){
    fit<-fitdistr(Data,"exponential")
    AIC.Exp<-2*length(fit$estimate)-2*fit$loglik
  }

  if(any(Test=="Gam(2)")){
    fit<-fitdistr(Data, "gamma")
    AIC.Gam2<-2*length(fit$estimate)-2*fit$loglik
  }

  if(any(Test=="Gam(3)")){
    data.gamlss <- data.frame(X=Data)
    fit.Gamma3 <- tryCatch(gamlss(X~1, data=data.gamlss, family=GG, trace=FALSE),
                           error = function(e) "error")
    if(is.character(fit.Gamma3)){
      Test <- Test[-which(Test=="Gam(3)")]
    } else {
      AIC.Gam3 <- fit.Gamma3$aic
    }
  }

  data.gamlss <- data.frame(X=Data)

  # Fixed GamMix(2) section - removed duplicate code and fixed logic
  if(any(Test=="GamMix(2)")){
    ### 2 mixture-gamma dist.
    for(i in 1:100){
      fit.GamMIX2_GA <- tryCatch(gamlssMX(X~1, data=data.gamlss, family=GA, K=2, trace=FALSE),
                                 error = function(e) "error")
      if(!is.character(fit.GamMIX2_GA)) break
    }
    if(is.character(fit.GamMIX2_GA)){
      Test <- Test[-which(Test=="GamMix(2)")]
    } else {
      AIC.GamMix2 <- fit.GamMIX2_GA$aic
    }
  }

  if(any(Test=="GamMix(3)")){
    ### 3 mixture-gamma dist.
    for(i in 1:100){
      fit.GamMIX3_GA <- tryCatch(gamlssMX(X~1, data=data.gamlss, family=GA, K=3, trace=FALSE),
                                 error = function(e) "error")
      if( is.character(fit.GamMIX3_GA) ) next
      if( !is.character(fit.GamMIX3_GA) ) break
    }
    if( is.character(fit.GamMIX3_GA) ){
      #AIC.GamMIX3_GA <- -9999
      Test <- Test[-which(Test=="GamMix(3)")]
    }else{
      AIC.GamMix3 <- fit.GamMIX3_GA$aic
    }
  }

  if(any(Test=="LNorm")){
    fit<-fitdistr(Data,"lognormal")
    AIC.logNormal<-2*length(fit$estimate)-2*fit$loglik
  }
  if(any(Test=="TNorm")){
    fit <- fitdistr(Data, "normal")
    AIC.TNormal <- 2 * length(fit$estimate) - 2 * fit$loglik
  }
  if(any(Test=="Twe")){
   capture.output(
    fit <- tweedie.profile(Data ~ 1,p.vec=seq(1.5, 2.5, by=0.2), do.plot=FALSE),
    type = "output"
   )
   AIC.Tweedie<-2*3-2*fit$L.max
  }
  if(any(Test=="Weib")){
    fit<-fitdistr(Data,"weibull")
    AIC.Weib<-2*length(fit$estimate)-2*fit$loglik
  }

  #Plot layout
  par(mfrow=c(3,1))
  par(mar=c(4.2,4.2,1,1))

  #Plotting AIC
  plot(0,xlim=c(0,length(Test)),ylim=c(min(0,AIC.BS,AIC.Exp,AIC.Gam2,AIC.Gam3,AIC.GamMix2,AIC.GamMix3,AIC.logNormal,AIC.TNormal,AIC.Tweedie,AIC.Weib,na.rm=TRUE),max(0,AIC.BS,AIC.Exp,AIC.Gam2,AIC.Gam3,AIC.GamMix2,AIC.GamMix3,AIC.TNormal,AIC.logNormal,AIC.Tweedie,AIC.Weib,na.rm=TRUE)),type='n',xlab="Probability Distribution",ylab="AIC",xaxt='n',cex.axis=1,cex.lab=1,las=1)
  axis(1,seq(0.5,length(Test)-0.5,1),c("Birn-S","Exp","Gam(2)","Gam(3)","GamMix(2)","GamMix(3)","LogN","TNorm","Twe","Weib")[Test],cex.axis=0.71)
  if(any(Test=="BS")){rect(which(Test=="BS")-0.5-length(Test)/40,0,which(Test=="BS")-0.5+length(Test)/40,AIC.BS,col=ifelse(Selected=="BS",mypalette[1],"white"))}
  if(any(Test=="Exp")){rect(which(Test=="Exp")-0.5-length(Test)/40,0,which(Test=="Exp")-0.5+length(Test)/40,AIC.Exp,col=ifelse(Selected=="Exp",mypalette[2],"white"))}
  if(any(Test=="Gam(2)")){rect(which(Test=="Gam(2)")-0.5-length(Test)/40,0,which(Test=="Gam(2)")-0.5+length(Test)/40,AIC.Gam2,col=ifelse(Selected=="Gam(2)",mypalette[3],"white"))}
  if(any(Test=="Gam(3)")){rect(which(Test=="Gam(3)")-0.5-length(Test)/40,0,which(Test=="Gam(3)")-0.5+length(Test)/40,AIC.Gam3,col=ifelse(Selected=="Gam(3)",mypalette[4],"white"))}
  if(any(Test=="GamMix(2)")){rect(which(Test=="GamMix(2)")-0.5-length(Test)/40,0,which(Test=="GamMix(2)")-0.5+length(Test)/40,AIC.GamMix2,col=ifelse(Selected=="GamMix(2)",mypalette[5],"white"))}
  if(any(Test=="GamMix(3)")){rect(which(Test=="GamMix(3)")-0.5-length(Test)/40,0,which(Test=="GamMix(3)")-0.5+length(Test)/40,AIC.GamMix3,col=ifelse(Selected=="GamMix(3)",mypalette[6],"white"))}
  if(any(Test=="LNorm")){rect(which(Test=="LNorm")-0.5-length(Test)/40,0,which(Test=="LNorm")-0.5+length(Test)/40,AIC.logNormal,col=ifelse(Selected=="LNorm",mypalette[7],"white"))}
  if(any(Test=="TNorm")){rect(which(Test=="TNorm")-0.5-length(Test)/40,0,which(Test=="TNorm")-0.5+length(Test)/40,AIC.TNormal,col=ifelse(Selected=="TNorm",mypalette[8],"white"))}
  if(any(Test=="Twe")){rect(which(Test=="Twe")-0.5-length(Test)/40,0,which(Test=="Twe")-0.5+length(Test)/40,AIC.Tweedie,col=ifelse(Selected=="Twe",mypalette[9],"white"))}
  if(any(Test=="Weib")){rect(which(Test=="Weib")-0.5-length(Test)/40,0,which(Test=="Weib")-0.5+length(Test)/40,AIC.Weib,col=ifelse(Selected=="Weib",mypalette[10],"white"))}

  #PDF plot
  hist(Data, freq=FALSE,xlab=x_lab,col="white",main="",cex.lab=1,cex.axis=1,ylim=c(y_lim_min,y_lim_max),las=1)
  x<-seq(min(Data),max(Data),0.01)
  #text(5.35,0.1,"(f)",font=2,cex=1.75)
  legend("topright",paste("Fitted",Selected,"dist."),lty=1,col=mypalette[which(Dist==Selected)],cex=0.9,bty='n')

  if(Selected=="BS"){
    bdata2 <- data.frame(shape = exp(-0.5), scale = exp(0.5))
    bdata2 <- transform(bdata2, y = Data)
    fit <- vglm(y ~ 1, bisa, data = bdata2, trace = FALSE)
    lines(x,dbisa(x,Coef(fit)[1],Coef(fit)[2]),col=mypalette[1],lwd=2)
  }

  if(Selected=="Exp"){
    fit<-fitdistr(Data,"exponential")
    lines(x,dexp(x,fit$estimate[1]),col=mypalette[2],lwd=2)
  }

  if(Selected=="Gam(2)"){
    fit<-fitdistr(Data, "gamma")
    lines(x,dgamma(x,fit$estimate[1],fit$estimate[2]),col=mypalette[3],lwd=2)
  }

  if(Selected=="Gam(3)" && exists("fit.Gamma3")){
    lines(x,dGG(x, mu=exp(fit.Gamma3$mu.coefficients), sigma=exp(fit.Gamma3$sigma.coefficients), nu=fit.Gamma3$nu.coefficients),col=mypalette[4],lwd=2)
  }

  if(Selected=="GamMix(2)" && exists("fit.GamMIX2_GA")){
    prob.MX1 <- round(fit.GamMIX2_GA$prob[1],3)
    prob.MX2 <- 1 - prob.MX1
    #prob.MX2 - round(fit.GamMIX2_GA$prob[2],5) < 0.00001
    lines(x,dMX(x, mu=list(mu1=exp(fit.GamMIX2_GA$models[[1]]$mu.coefficients), mu2=exp(fit.GamMIX2_GA$models[[2]]$mu.coefficients)),
                sigma=list(sigma1=exp(fit.GamMIX2_GA$models[[1]]$sigma.coefficients), sigma2=exp(fit.GamMIX2_GA$models[[2]]$sigma.coefficients)),
                pi = list(pi1=prob.MX1, pi2=prob.MX2), family=list(fam1="GA", fam2="GA")),col=mypalette[5],lwd=2)
  }

  if(Selected=="GamMix(3)" && exists("fit.GamMIX3_GA")){
    prob.MX1 <- round(fit.GamMIX3_GA$prob[1],3)
    prob.MX2 <- round(fit.GamMIX3_GA$prob[2],3)
    prob.MX3 <- 1 - prob.MX1 - prob.MX2
    #prob.MX3 - round(fit.GamMIX3_GA$prob[3],5) < 0.00001
    lines(x,dMX(x, mu=list(mu1=exp(fit.GamMIX3_GA$models[[1]]$mu.coefficients), mu2=exp(fit.GamMIX3_GA$models[[2]]$mu.coefficients), mu3=exp(fit.GamMIX3_GA$models[[3]]$mu.coefficients)),
                sigma=list(sigma1=exp(fit.GamMIX3_GA$models[[1]]$sigma.coefficients), sigma2=exp(fit.GamMIX3_GA$models[[2]]$sigma.coefficients), sigma3=exp(fit.GamMIX3_GA$models[[3]]$sigma.coefficients)),
                pi = list(pi1=prob.MX1, pi2=prob.MX2, pi3=prob.MX3), family=list(fam1="GA", fam2="GA", fam3="GA")),col=mypalette[6],lwd=2)
  }

  if(Selected=="LNorm"){
    fit<-fitdistr(Data,"lognormal")
    lines(x,dlnorm(x,fit$estimate[1],fit$estimate[2]),col=mypalette[7],lwd=2)
  }

  if(Selected=="TNorm"){
    fit<-fitdistr(Data, "normal")
    lines(x,dtruncnorm(x,a=min(Data),mean=fit$estimate[1],sd=fit$estimate[2]),col=mypalette[8],lwd=2)
  }

  if(Selected=="Twe"){
    fit <- tweedie.profile(Data ~ 1,
                           p.vec=seq(1.5, 2.5, by=0.2), do.plot=FALSE)
    lines(x,dtweedie(x,  power=fit$p.max, mu=mean(Data), phi=fit$phi.max),col=mypalette[9],lwd=2)
  }

  if(Selected=="Weib"){
    fit<-fitdistr(Data,"weibull")
    lines(x,dweibull(x,fit$estimate[1],fit$estimate[2]),col=mypalette[10],lwd=2)
  }

  #CDF plot
  plot(sort(Data),seq(1,length(Data),1)/(length(Data)),ylim=c(0,1),xlab=x_lab,ylab="P(X<x)",main="",pch=16,cex.lab=1,cex.axis=1,las=1)
  x<-seq(min(Data),max(Data),0.01)
  eta<-sqrt((1/(2*length(Data)))*log(2/0.95))
  lines(sort(Data),ifelse(seq(1,length(Data),1)/(length(Data))+eta>1,1,seq(1,length(Data),1)/(length(Data))+eta),col=1,lty=2)
  lines(sort(Data),ifelse(seq(1,length(Data),1)/(length(Data))-eta<0,0,seq(1,length(Data),1)/(length(Data))-eta),col=1,lty=2)
  legend("bottomright",c("95% Conf. Interval",paste("Fitted",Selected, "dist.")),lty=c(2,1),col=mypalette[which(Dist==Selected)],cex=1,bty='n',border = "white")

  if(Selected=="BS"){
    bdata2 <- data.frame(shape = exp(-0.5), scale = exp(0.5))
    bdata2 <- transform(bdata2, y = Data)
    fit <- vglm(y ~ 1, bisa, data = bdata2, trace = FALSE)
    lines(x,pbisa(x,Coef(fit)[1],Coef(fit)[2]),col=mypalette[1],lwd=2)
  }

  if(Selected=="Exp"){
    fit<-fitdistr(Data,"exponential")
    lines(x,pexp(x,fit$estimate[1]),col=mypalette[2],lwd=2)
  }

  if(Selected=="Gam(2)"){
    fit<-fitdistr(Data, "gamma")
    lines(x,pgamma(x,fit$estimate[1],fit$estimate[2]),col=mypalette[3],lwd=2)
  }

  if(Selected=="Gam(3)" && exists("fit.GamMIX3_GA")){
    lines(x,pGG(x, mu=exp(fit.Gamma3$mu.coefficients), sigma=exp(fit.Gamma3$sigma.coefficients), nu=fit.Gamma3$nu.coefficients),col=mypalette[4],lwd=2)
  }

  if(Selected=="GamMix(2)" && exists("fit.GamMIX3_GA")){
    prob.MX1 <- round(fit.GamMIX2_GA$prob[1],3)
    prob.MX2 <- 1 - prob.MX1
    lines(x,pMX(x, mu=list(mu1=exp(fit.GamMIX2_GA$models[[1]]$mu.coefficients), mu2=exp(fit.GamMIX2_GA$models[[2]]$mu.coefficients)),
                sigma=list(sigma1=exp(fit.GamMIX2_GA$models[[1]]$sigma.coefficients), sigma2=exp(fit.GamMIX2_GA$models[[2]]$sigma.coefficients)),
                pi = list(pi1=prob.MX1, pi2=prob.MX2), family=list(fam1="GA", fam2="GA")),col=mypalette[5],lwd=2)
  }

  if(Selected=="GamMix(3)" && exists("fit.GamMIX3_GA")){
    prob.MX1 <- round(fit.GamMIX3_GA$prob[1],3)
    prob.MX2 <- round(fit.GamMIX3_GA$prob[2],3)
    prob.MX3 <- 1 - prob.MX1 - prob.MX2
    lines(x,pMX(x, mu=list(mu1=exp(fit.GamMIX3_GA$models[[1]]$mu.coefficients), mu2=exp(fit.GamMIX3_GA$models[[2]]$mu.coefficients), mu3=exp(fit.GamMIX3_GA$models[[3]]$mu.coefficients)),
                sigma=list(sigma1=exp(fit.GamMIX3_GA$models[[1]]$sigma.coefficients), sigma2=exp(fit.GamMIX3_GA$models[[2]]$sigma.coefficients), sigma3=exp(fit.GamMIX3_GA$models[[3]]$sigma.coefficients)),
                pi = list(pi1=prob.MX1, pi2=prob.MX2, pi3=prob.MX3), family=list(fam1="GA", fam2="GA", fam3="GA")),col=mypalette[6],lwd=2)
  }

  if(Selected=="LNorm"){
    fit<-fitdistr(Data,"lognormal")
    lines(x,plnorm(x,fit$estimate[1],fit$estimate[2]),col=mypalette[7],lwd=2)
  }

  if(Selected=="TNorm"){
    fit<-fitdistr(Data,"normal")
    lines(x,ptruncnorm(x,a=min(Data),mean=fit$estimate[1],sd=fit$estimate[2]),col=mypalette[8],lwd=2)
  }

  if(Selected=="Twe"){
    fit <- tweedie.profile(Data ~ 1,p.vec=seq(1.5, 2.5, by=0.2), do.plot=FALSE)
    lines(x,ptweedie(x,  power=fit$p.max, mu=mean(Data), phi=fit$phi.max),col=mypalette[9],lwd=2,pch=16,ylab="P(X<x)")
  }

  if(Selected=="Weib"){
    fit<-fitdistr(Data, "weibull")
    lines(x,pweibull(x,fit$estimate[1],fit$estimate[2]),col=mypalette[10],lwd=2)
  }
}
