#' Calculates joint return periods
#'
#' Calculates joint return periods via simulation. A large number of realizations are simulated from the copulas fit to the conditioned samples, in proportion with the sizes of the conditional samples.
#' The realizations are transformed to the original scale and the relevant probabilities estimated empirically.
#'
#' @param Data Data frame containing two co-occurring time series.
#' @param Data_Con1 Data frame containing the conditional sample (declustered excesses paired with concurrent values of other variable), conditioned on the variable in the first column.
#' @param Data_Con2 Data frame containing the conditional sample (declustered excesses paired with concurrent values of other variable), conditioned on the variable in the second column. Can be obtained using the \code{Con_Sampling_2D} function.
#' @param u1 Numeric vector of length one specifying the (quantile) threshold above which the variable in the first column was sampled in Data_Con1.
#' @param u2 Numeric vector of length one specifying the (quantile) threshold above which the variable in the second column was sampled in Data_Con2.
#' @param Thres1 Numeric vector of length one specifying the threshold above which the variable in the first column was sampled in Data_Con1. Only one of \code{u1} and \code{Thres1} should be supplied. Default is \code{NA}.
#' @param Thres2 Numeric vector of length one specifying the threshold above which the variable in the second column was sampled in Data_Con2. Only one of \code{u2} and \code{Thres2} should be supplied. Default is \code{NA}.
#' @param Copula_Family1 Numeric vector of length one specifying the copula family used to model the \code{Data_Con1} dataset.
#' @param Copula_Family2 Numeric vector of length one specifying the copula family used to model the \code{Data_Con2} dataset. Best fitting of 40 copulas can be found using the \code{Copula_Threshold_2D} function.
#' @param Marginal_Dist1 Character vector of length one specifying (non-extreme) distribution used to model the marginal distribution of the non-conditioned variable in \code{Data_Con1}.
#' @param Marginal_Dist2 Character vector of length one specifying (non-extreme) distribution used to model the marginal distribution of the non-conditioned variable in \code{Data_Con2}.
#' @param Con1 Character vector of length one specifying the name of variable in the first column of \code{Data}.
#' @param Con2 Character vector of length one specifying the name of variable in the second column of \code{Data}.
#' @param mu Numeric vector of length one specifying the (average) occurrence frequency of events in \code{Data}. Default is \code{365.25}, daily data.
#' @param RP_Var1 Numeric vector of length one specifying the univariate return period of variable \code{Con1}.
#' @param RP_Var2 Numeric vector of length one specifying the univariate return period of variable \code{Con2}.
#' @param Var1 Numeric vector specifying the values of the variable in the first column of \code{Data}. Default is \code{NA}.
#' @param Var2 Numeric vector specifying the values of the variable in the second column of \code{Data}. Default is \code{NA}.
#' @param x_lab Character vector specifying the x-axis label.
#' @param y_lab Character vector specifying the y-axis label.
#' @param x_lim_min Numeric vector of length one specifying x-axis minimum. Default is \code{NA}.
#' @param x_lim_max Numeric vector of length one specifying x-axis maximum. Default is \code{NA}.
#' @param y_lim_min Numeric vector of length one specifying y-axis minimum. Default is \code{NA}.
#' @param y_lim_max Numeric vector of length one specifying y-axis maximum. Default is \code{NA}.
#' @param DecP Numeric vector of length one specifying the number of decimal places to round the data in the conditional samples to in order to identify observations in both conditional samples. Default is \code{2}.
#' @param N Numeric vector of length one specifying the size of the sample from the fitted joint distributions used to estimate the density along an isoline. Samples are collected from the two joint distribution with proportions consistent with the total number of extreme events conditioned on each variable. Default is \code{10^6}
#' @return Console output is a vector \code{RP_Copula} of the joint return period of the specified (Var1,Var2) events according to the conditional sampling - copula theory approach.
#' @seealso \code{\link{Design_Event_2D}}
#' @export
#' @examples
#' con.sample.Rainfall<-Con_Sampling_2D(Data_Detrend=S22.Detrend.df[,-c(1,4)],
#'                                      Data_Declust=S22.Detrend.Declustered.df[,-c(1,4)],
#'                                      Con_Variable="Rainfall",u=0.97)
#' con.sample.OsWL<-Con_Sampling_2D(Data_Detrend=S22.Detrend.df[,-c(1,4)],
#'                                  Data_Declust=S22.Detrend.Declustered.df[,-c(1,4)],
#'                                  Con_Variable="OsWL",u=0.97)
#' cop.Rainfall<-Copula_Threshold_2D(Data_Detrend=S22.Detrend.df[,-c(1,4)],
#'                                   Data_Declust=S22.Detrend.Declustered.df[,-c(1,4)],u1 =0.97,
#'                                   y_lim_min=-0.075,y_lim_max=0.25,
#'                                   Upper=c(2,9),Lower=c(2,10),GAP=0.15)$Copula_Family_Var1
#' cop.OsWL<-Copula_Threshold_2D(Data_Detrend=S22.Detrend.df[,-c(1,4)],
#'                               Data_Declust=S22.Detrend.Declustered.df[,-c(1,4)],u2 =0.97,
#'                               y_lim_min=-0.075, y_lim_max =0.25,
#'                               Upper=c(2,9),Lower=c(2,10),GAP=0.15)$Copula_Family_Var2
#' #Joint excedence probability of a 5 inch rainfall and 2 ft O-sWL at S-22.
#' Joint_RP_2D(Data=S22.Detrend.df,
#'             Data_Con1=con.sample.Rainfall$Data, Data_Con2=con.sample.OsWL$Data,
#'             u1=0.98, u2=0.98,
#'             Copula_Family1=cop.Rainfall,Copula_Family2=cop.OsWL,
#'             Marginal_Dist1="Logis", Marginal_Dist2="BS",
#'             Con1 = "Rainfall", Con2 = "OsWL",
#'             mu = 365.25,
#'             RP_Var1 =10, RP_Var2 =10,
#'             x_lab = "Rainfall (Inches)", y_lab = "O-sWL (ft NGVD 29)",
#'             y_lim_max = 10,
#'             N=10^7)
Joint_RP_2D<-function (Data, Data_Con1, Data_Con2, u1, u2,
                       Thres1=NA, Thres2=NA,
                       Copula_Family1, Copula_Family2,
                       Marginal_Dist1, Marginal_Dist2,
                       Con1 = "Rainfall", Con2 = "OsWL", mu = 365.25,
                       RP_Var1, RP_Var2, Var1=NA,Var2=NA,
                       x_lab = "Rainfall (mm)", y_lab = "O-sWL (mNGVD 29)", x_lim_min = NA,
                       x_lim_max = NA, y_lim_min = NA, y_lim_max = NA, DecP = 2, N=10^6)
{
  ###Preliminaries

  #Remove 1st column of Data if it is a Date or factor object.
  if (class(Data[, 1])[1] == "Date" | class(Data[, 1])[1] == "factor" | class(Data[, 1])[1] == "POSIXct" | class(Data[, 1])[1] == "character") {
    Data <- Data[, -1]
  }

  #Find the columns in Data (which should be consistent in terms of column order of the other data input objects)
  #of conditioning variable 1 (Con1) and conditioning variable 2 (Con2).
  var1 <- Var1
  var2 <- Var2
  con1 <- which(names(Data) == Con1)
  con2 <- which(names(Data) == Con2)

  #Axis limits for plots
  x_min <- ifelse(is.na(x_lim_min) == T, min(na.omit(Data[,con1])), x_lim_min)
  x_max <- ifelse(is.na(x_lim_max) == T, max(na.omit(Data[,con1])), x_lim_max)
  y_min <- ifelse(is.na(y_lim_min) == T, min(na.omit(Data[,con2])), y_lim_min)
  y_max <- ifelse(is.na(y_lim_max) == T, max(na.omit(Data[,con2])), y_lim_max)

  #Finding the threshold if specified as a quantile
  if(is.na(Thres1)==T){
    Thres1<-quantile(na.omit(Data[,con1]), u1)
  }

  ##Finding the value of variable con1 associated with a return peroid of RP_Var1
  #Fitting the GPD
  GPD_con1 <- evm(Data_Con1[, con1], th = Thres1, penalty = "gaussian", priorParameters = list(c(0,0), matrix(c(100^2, 0, 0, 0.25), nrow = 2)))
  #Calculate the time period spanned by the original dataset in terms of mu (only including occasions where both variables are observed).
  time.period<-nrow(Data[which(is.na(Data[,con1]) == F & is.na(Data[, con2]) == F),])/mu
  #Calculate the rate of occurrences of extremes (in terms of mu) in Data_Con1.
  rate<-nrow(Data_Con1)/time.period
  #Interarrival time
  EL_Con1<-1/rate
  #Value of con1 with return period RP_Var1
  if(is.na(var1)==T){
   Var1<-as.numeric(u2gpd((1-EL_Con1/RP_Var1), p = 1, th=Thres1 , sigma=exp(GPD_con1$coefficients[1]),xi= GPD_con1$coefficients[2]))
  }
  if(is.na(var2)==F){
   RP_Var1<-1/(1-pgpd(Var1, u=Thres1 , sigma=exp(GPD_con1$coefficients[1]),xi= GPD_con1$coefficients[2]))
  }

  ##Fit the specified marginal distribution (Marginal_Dist1) to the variable con2 in Data_Con1.
  if (Marginal_Dist1 == "BS") {
    bdata2 <- data.frame(shape = exp(-0.5), scale = exp(0.5))
    bdata2 <- transform(bdata2, y = Data_Con1[, con2])
    marginal_non_con1 <- vglm(y ~ 1, bisa, data = bdata2,
                              trace = FALSE)
  }
  if (Marginal_Dist1 == "Exp") {
    marginal_non_con1 <- fitdistr(Data_Con1[, con2], "exponential")
  }
  if (Marginal_Dist1 == "Gam(2)") {
    marginal_non_con1 <- fitdistr(Data_Con1[, con2], "gamma2")
  }
  if(Marginal_Dist1 == "Gam(3)"){
    data.gamlss<-data.frame(X=Data_Con1[,con2])
    marginal_non_con1 <- tryCatch(gamlss(X~1, data=data.gamlss, family=GG),
                                  error = function(e) "error")
  }
  if(Marginal_Dist1 == "GamMix(2)"){
    data.gamlss<-data.frame(X=Data_Con1[,con2])
    marginal_non_con1 <- tryCatch(gamlssMX(X~1, data=data.gamlss, family=GA, K=2),
                                  error = function(e) "error")
  }
  if(Marginal_Dist1 == "GamMix(3)"){
    data.gamlss<-data.frame(X=Data_Con1[,con2])
    marginal_non_con1 <- tryCatch(gamlssMX(X~1, data=data.gamlss, family=GA, K=3),
                                  error = function(e) "error")
  }
  if (Marginal_Dist1 == "Gaus") {
    marginal_non_con1 <- fitdistr(Data_Con1[, con2], "normal")
  }
  if (Marginal_Dist1 == "InvG") {
    marginal_non_con1 <- fitdist(Data_Con1[, con2], "invgauss",
                                 start = list(mean = 5, shape = 1))
  }
  if (Marginal_Dist1 == "Logis") {
    marginal_non_con1 <- fitdistr(Data_Con1[, con2], "logistic")
  }
  if (Marginal_Dist1 == "LogN") {
    marginal_non_con1 <- fitdistr(Data_Con1[, con2], "lognormal")
  }
  if (Marginal_Dist1 == "TNorm") {
    marginal_non_con1 <- fitdistr(Data_Con1[, con2], "normal")
  }
  if (Marginal_Dist1 == "Twe") {
    marginal_non_con1 <- tweedie.profile(Data_Con1[, con2] ~
                                           1, p.vec = seq(1.5, 2.5, by = 0.2), do.plot = FALSE)
  }
  if (Marginal_Dist1 == "Weib") {
    marginal_non_con1 <- fitdistr(Data_Con1[, con2], "weibull")
  }

  #Finding the threshold if specified as a quantile
  if(is.na(Thres2)==T){
    Thres2<-quantile(na.omit(Data[,con2]), u2)
  }

  if (is.na(var2)==F){
    if (Marginal_Dist1 == "BS") {
      RP_Var2_con1 <- 1/(1-pbisa(Var2, as.numeric(Coef(marginal_non_con1)[1]),
                            as.numeric(Coef(marginal_non_con1)[2])))
    }
    if (Marginal_Dist1 == "Exp") {
      RP_Var2_con1 <- 1/(1-pexp(Var2, rate = as.numeric(marginal_non_con1$estimate[1])))
    }
    if (Marginal_Dist1 == "Gam(2)") {
      RP_Var2_con1 <- 1/(1-pgamma(Var2, shape = as.numeric(marginal_non_con1$estimate[1]),
                             rate = as.numeric(marginal_non_con1$estimate[2])))
    }
    if(Marginal_Dist1 == "Gam(3)"){
      RP_Var2_con1<-1/(1-pGG(Var2, mu=exp(marginal_non_con1$mu.coefficients), sigma=exp(marginal_non_con1$sigma.coefficients), nu=marginal_non_con1$nu.coefficients))
    }
    if(Marginal_Dist1=="GamMix(2)"){
      xx <- seq(0, max(Data_Con2[,2])*10, 0.001)
      prob.MX1 <- round(marginal_non_con1$prob[1],3)
      prob.MX2 <- 1 - prob.MX1
      cdf.MX<-pMX(xx, mu=list(mu1=exp(marginal_non_con1$models[[1]]$mu.coefficients), mu2=exp(marginal_non_con1$models[[2]]$mu.coefficients)),
                  sigma=list(sigma1=exp(marginal_non_con1$models[[1]]$sigma.coefficients), sigma2=exp(marginal_non_con1$models[[2]]$sigma.coefficients)),
                  pi = list(pi1=prob.MX1, pi2=prob.MX2), family=list(fam1="GA", fam2="GA"))
      RP_Var2_con1 <- approx(cdf.MX, xx, Var2)$y
    }
    if(Marginal_Dist1=="GamMix(3)"){
      xx <- seq(0, max(Data_Con2[,2])*10, 0.001)
      prob.MX1 <- round(marginal_non_con1$prob[1],3)
      prob.MX2 <- round(marginal_non_con1$prob[2],3)
      prob.MX3 <- 1 - prob.MX1 - prob.MX2
      cdf.MX<-pMX(xx, mu=list(mu1=exp(marginal_non_con1$models[[1]]$mu.coefficients), mu2=exp(marginal_non_con1$models[[2]]$mu.coefficients), mu3=exp(marginal_non_con1$models[[3]]$mu.coefficients)),
                  sigma=list(sigma1=exp(marginal_non_con1$models[[1]]$sigma.coefficients), sigma2=exp(marginal_non_con1$models[[2]]$sigma.coefficients), sigma3=exp(marginal_non_con1$models[[3]]$sigma.coefficients)),
                  pi = list(pi1=prob.MX1, pi2=prob.MX2, pi3=prob.MX3), family=list(fam1="GA", fam2="GA", fam3="GA"))
      RP_Var2_con1 <- approx(cdf.MX, xx, Var2)$y
    }
    if (Marginal_Dist1 == "Gaus") {
      RP_Var2_con1 <- 1/(1-pnorm(Var2, mean = as.numeric(marginal_non_con1$estimate[1]),
                            sd = as.numeric(marginal_non_con1$estimate[2])))
    }
    if (Marginal_Dist1 == "InvG") {
      RP_Var2_con1 <- 1/(1-pinvgauss(Var2, mean = as.numeric(marginal_non_con1$estimate[1]),
                                shape = as.numeric(marginal_non_con1$estimate[2])))
    }
    if (Marginal_Dist1 == "Logis") {
      RP_Var2_con1 <- 1/(1-plogis(Var2, location = as.numeric(marginal_non_con1$estimate[1]),
                             scale = as.numeric(marginal_non_con1$estimate[2])))
    }
    if (Marginal_Dist1 == "LogN") {
      RP_Var2_con1 <- 1/(1-plnorm(Var2, meanlog = as.numeric(marginal_non_con1$estimate[1]),
                             sdlog = as.numeric(marginal_non_con1$estimate[2])))
    }
    if (Marginal_Dist1 == "TNorm") {
      RP_Var2_con1 <- 1/(1-ptruncnorm(Var2, a = min(Data_Con1[,con2]), mean = as.numeric(marginal_non_con1$estimate[1]),
                                 sd = as.numeric(marginal_non_con1$estimate[2])))
    }
    if (Marginal_Dist1 == "Twe") {
      RP_Var2_con1 <- 1/(1-ptweedie(Var2, power = marginal_non_con1$p.max,
                               mu = mean(Data_Con1[, con2]), phi = marginal_non_con1$phi.max))
    }
    if (Marginal_Dist1 == "Weib") {
      RP_Var2_con1 <- 1/(1-pweibull(Var2, shape = as.numeric(marginal_non_con1$estimate[1]),
                               scale = as.numeric(marginal_non_con1$estimate[2])))
    }
  }

  ##Finding the value of variable con2 associated with a return peroid of RP_Var2
  #Fitting the GPD to con2 in Data_Con2
  GPD_con2 <- evm(Data_Con2[, con2], th = Thres2, penalty = "gaussian", priorParameters = list(c(0, 0), matrix(c(100^2, 0, 0, 0.25), nrow = 2)))
  #Calculate the time period spanned by the original dataset in terms of mu (only including occasions where both variables are observed).
  time.period<-nrow(Data[which(is.na(Data[,con1]) == F & is.na(Data[, con2]) == F),])/mu
  #Calculate the rate of occurrences of extremes (in terms of mu) in Data_Con1.
  rate<-nrow(Data_Con2)/time.period
  #Calculate the inter-arrival time of extremes (in terms of mu) in Data_Con1.
  EL_Con2<-1/rate
  #Value of con2 with return period RP_Var2
  if(is.na(var2)==T){
   Var2<-as.numeric(u2gpd((1-EL_Con2/RP_Var2), p = 1, th=Thres2 , sigma=exp(GPD_con2$coefficients[1]),xi= GPD_con2$coefficients[2]))
  }
  if(is.na(var2)==F){
   RP_Var2<-1/(1-pgpd(Var2, u=Thres2 , sigma=exp(GPD_con2$coefficients[1]),xi= GPD_con2$coefficients[2]))
  }

  ##Fit the specified marginal distribution (Marginal_Dist2) to the non-conditioned variable con1 in Data_Con2.
  if (Marginal_Dist2 == "BS") {
    bdata2 <- data.frame(shape = exp(-0.5), scale = exp(0.5))
    bdata2 <- transform(bdata2, y = Data_Con2[, con1])
    marginal_non_con2 <- vglm(y ~ 1, bisa, data = bdata2,
                              trace = FALSE)
  }
  if (Marginal_Dist2 == "Exp") {
    marginal_non_con2 <- fitdistr(Data_Con2[, con1], "exponential")
  }
  if (Marginal_Dist2 == "Gam(2)") {
    marginal_non_con2 <- fitdistr(Data_Con2[, con1], "gamma")
  }
  if(Marginal_Dist2 == "Gam(3)"){
    data.gamlss<-data.frame(X=Data_Con2[,con1])
    marginal_non_con2 <- tryCatch(gamlss(X~1, data=data.gamlss, family=GG),
                                  error = function(e) "error")
  }
  if(Marginal_Dist2 == "GamMix(2)"){
    data.gamlss<-data.frame(X=Data_Con2[,con1])
    marginal_non_con2 <- tryCatch(gamlssMX(X~1, data=data.gamlss, family=GA, K=2),
                                  error = function(e) "error")
  }
  if(Marginal_Dist2 == "GamMix(3)"){
    data.gamlss<-data.frame(X=Data_Con2[,con1])
    marginal_non_con2 <- tryCatch(gamlssMX(X~1, data=data.gamlss, family=GA, K=3),
                                  error = function(e) "error")
  }
  if (Marginal_Dist2 == "Gaus") {
    marginal_non_con2 <- fitdistr(Data_Con2[, con1], "normal")
  }
  if (Marginal_Dist2 == "InvG") {
    marginal_non_con2 <- fitdist(Data_Con2[, con1], "invgauss",
                                 start = list(mean = 5, shape = 1))
  }
  if (Marginal_Dist2 == "Logis") {
    marginal_non_con2 <- fitdistr(Data_Con2[, con1], "logistic")
  }
  if (Marginal_Dist2 == "LogN") {
    marginal_non_con2 <- fitdistr(Data_Con2[, con1], "lognormal")
  }
  if (Marginal_Dist2 == "TNorm") {
    marginal_non_con2 <- fitdistr(Data_Con2[, con1], "normal")
  }
  if (Marginal_Dist2 == "Twe") {
    marginal_non_con2 <- tweedie.profile(Data_Con2[, con1] ~
                                           1, p.vec = seq(1.5, 2.5, by = 0.2), do.plot = FALSE)
  }
  if (Marginal_Dist2 == "Weib") {
    marginal_non_con2 <- fitdistr(Data_Con2[, con1], "weibull")
  }

  if (is.na(var2)==F){
    if (Marginal_Dist2 == "BS") {
      RP_Var1_con2 <- 1/(1-pbisa(Var1, as.numeric(Coef(marginal_non_con2)[1]),
                            as.numeric(Coef(marginal_non_con2)[2])))
    }
    if (Marginal_Dist2 == "Exp") {
      RP_Var1_con2 <- 1/(1-pexp(Var1, rate = as.numeric(marginal_non_con2$estimate[1])))
    }
    if (Marginal_Dist2 == "Gam(2)") {
      RP_Var1_con2 <- 1/(1-pgamma(Var1, shape = as.numeric(marginal_non_con2$estimate[1]),
                             rate = as.numeric(marginal_non_con2$estimate[2])))
    }
    if(Marginal_Dist2=="Gam(3)"){
      RP_Var1_con2<-1/(1-pGG(Var1, mu=exp(marginal_non_con2$mu.coefficients), sigma=exp(marginal_non_con2$sigma.coefficients), nu=marginal_non_con2$nu.coefficients))
    }
    if(Marginal_Dist2=="GamMix(2)"){
      xx <- seq(0, max(Data_Con1[,1])*10, 0.001)
      prob.MX1 <- round(marginal_non_con2$prob[1],3)
      prob.MX2 <- 1 - prob.MX1
      cdf.MX<-pMX(xx, mu=list(mu1=exp(marginal_non_con2$models[[1]]$mu.coefficients), mu2=exp(marginal_non_con2$models[[2]]$mu.coefficients)),
                  sigma=list(sigma1=exp(marginal_non_con2$models[[1]]$sigma.coefficients), sigma2=exp(marginal_non_con2$models[[2]]$sigma.coefficients)),
                  pi = list(pi1=prob.MX1, pi2=prob.MX2), family=list(fam1="GA", fam2="GA"))
      RP_Var1_con2 <- approx(cdf.MX, xx, Var1)$y
    }
    if(Marginal_Dist2=="GamMix(3)"){
      xx <- seq(0, max(Data_Con1[,1])*10, 0.001)
      prob.MX1 <- round(marginal_non_con2$prob[1],3)
      prob.MX2 <- round(marginal_non_con2$prob[2],3)
      prob.MX3 <- 1 - prob.MX1 - prob.MX2
      cdf.MX<-pMX(xx, mu=list(mu1=exp(marginal_non_con2$models[[1]]$mu.coefficients), mu2=exp(marginal_non_con2$models[[2]]$mu.coefficients), mu3=exp(marginal_non_con2$models[[3]]$mu.coefficients)),
                  sigma=list(sigma1=exp(marginal_non_con2$models[[1]]$sigma.coefficients), sigma2=exp(marginal_non_con2$models[[2]]$sigma.coefficients), sigma3=exp(marginal_non_con2$models[[3]]$sigma.coefficients)),
                  pi = list(pi1=prob.MX1, pi2=prob.MX2, pi3=prob.MX3), family=list(fam1="GA", fam2="GA", fam3="GA"))
      RP_Var1_con2 <- approx(cdf.MX, xx, Var1)$y
    }
    if (Marginal_Dist2 == "Gaus") {
      RP_Var1_con2 <- 1/(1-pnorm(Var1, mean = as.numeric(marginal_non_con2$estimate[1]),
                            sd = as.numeric(marginal_non_con2$estimate[2])))
    }
    if (Marginal_Dist2 == "InvG") {
      RP_Var1_con2 <- 1/(1-pinvgauss(Var1, mean = as.numeric(marginal_non_con2$estimate[1]),
                                shape = as.numeric(marginal_non_con2$estimate[2])))
    }
    if (Marginal_Dist2 == "Logis") {
      RP_Var1_con2 <- 1/(1-plogis(Var1, location = as.numeric(marginal_non_con2$estimate[1]),
                             scale = as.numeric(marginal_non_con2$estimate[2])))
    }
    if (Marginal_Dist2 == "LogN") {
      RP_Var1_con2 <- 1/(1-plnorm(Var1, meanlog = as.numeric(marginal_non_con2$estimate[1]),
                             sdlog = as.numeric(marginal_non_con2$estimate[2])))
    }
    if (Marginal_Dist2 == "TNorm") {
      RP_Var1_con2 <- 1/(1-ptruncnorm(Var1, a = min(Data_Con2[,con2]), mean = as.numeric(marginal_non_con2$estimate[1]),
                                 sd = as.numeric(marginal_non_con2$estimate[2])))
    }
    if (Marginal_Dist2 == "Twe") {
      RP_Var1_con2 <- 1/(1-ptweedie(Var1, power = marginal_non_con2$p.max,
                               mu = mean(Data_Con2[, con2]), phi = marginal_non_con2$phi.max))
    }
    if (Marginal_Dist2 == "Weib") {
      RP_Var1_con2 <- 1/(1-pweibull(Var1, shape = as.numeric(marginal_non_con2$estimate[1]),
                               scale = as.numeric(marginal_non_con2$estimate[2])))
    }
  }

  ###Simulating sample from the joint distribution (copula+marginals) fit to the sample conditioned on Con1
  #Fit the specified copula family (Copula_Family1) to the observations in Data_Con1.
  obj1 <- BiCopSelect(pobs(Data_Con1[, 1]), pobs(Data_Con1[,  2]), familyset = Copula_Family1, selectioncrit = "AIC",
                      indeptest = FALSE, level = 0.05, weights = NA, rotations = TRUE,
                      se = FALSE, presel = TRUE, method = "mle")
  #Simulate a sample from the fitted copula. Out of the sample size 'N' the proportion of the sample from the
  #copula assoicated with Data_Con1 is proportional to the size of Data_Con1 relative to Data_Con2.
  sample <- BiCopSim(round(N * nrow(Data_Con1)/(nrow(Data_Con1) + nrow(Data_Con2)), 0), obj1)
  #Transform the realizations of the variable con1 to the original scale using the inverse CDF (quantile function)
  #of the GPD contained in the u2gpd function.
  cop.sample1.con <- u2gpd(sample[, con1], p = 1, th = Thres1, sigma = exp(GPD_con1$coefficients[1]),xi = GPD_con1$coefficients[2])

  #Transform the realizations of variable con2 to the original scale using the inverse CDF (quantile function)
  #of the selected parametric (non-extreme value) distribution (Marginal_Dist1)
  if (Marginal_Dist1 == "BS") {
    cop.sample1.non.con <- qbisa(sample[, con2], as.numeric(Coef(marginal_non_con1)[1]),
                                 as.numeric(Coef(marginal_non_con1)[2]))
  }
  if (Marginal_Dist1 == "Exp") {
    cop.sample1.non.con <- qexp(sample[, con2], rate = as.numeric(marginal_non_con1$estimate[1]))
  }
  if (Marginal_Dist1 == "Gam(2)") {
    cop.sample1.non.con <- qgamma(sample[, con2], shape = as.numeric(marginal_non_con1$estimate[1]),
                                  rate = as.numeric(marginal_non_con1$estimate[2]))
  }
  if(Marginal_Dist1=="Gam(3)"){
    cop.sample1.non.con<-qGG(sample[,con2], mu=exp(marginal_non_con1$mu.coefficients), sigma=exp(marginal_non_con1$sigma.coefficients), nu=marginal_non_con1$nu.coefficients)
  }
  if(Marginal_Dist1=="GamMix(2)"){
    xx <- seq(0, max(Data_Con2[,2])*10, 0.001)
    prob.MX1 <- round(marginal_non_con1$prob[1],3)
    prob.MX2 <- 1 - prob.MX1
    cdf.MX<-pMX(xx, mu=list(mu1=exp(marginal_non_con1$models[[1]]$mu.coefficients), mu2=exp(marginal_non_con1$models[[2]]$mu.coefficients)),
                sigma=list(sigma1=exp(marginal_non_con1$models[[1]]$sigma.coefficients), sigma2=exp(marginal_non_con1$models[[2]]$sigma.coefficients)),
                pi = list(pi1=prob.MX1, pi2=prob.MX2), family=list(fam1="GA", fam2="GA"))
    cop.sample1.non.con <- approx(cdf.MX, xx, sample[,con2])$y
  }
  if(Marginal_Dist1=="GamMix(3)"){
    xx <- seq(0, max(Data_Con2[,2])*10, 0.001)
    prob.MX1 <- round(marginal_non_con1$prob[1],3)
    prob.MX2 <- round(marginal_non_con1$prob[2],3)
    prob.MX3 <- 1 - prob.MX1 - prob.MX2
    cdf.MX<-pMX(xx, mu=list(mu1=exp(marginal_non_con1$models[[1]]$mu.coefficients), mu2=exp(marginal_non_con1$models[[2]]$mu.coefficients), mu3=exp(marginal_non_con1$models[[3]]$mu.coefficients)),
                sigma=list(sigma1=exp(marginal_non_con1$models[[1]]$sigma.coefficients), sigma2=exp(marginal_non_con1$models[[2]]$sigma.coefficients), sigma3=exp(marginal_non_con1$models[[3]]$sigma.coefficients)),
                pi = list(pi1=prob.MX1, pi2=prob.MX2, pi3=prob.MX3), family=list(fam1="GA", fam2="GA", fam3="GA"))
    cop.sample1.non.con <- approx(cdf.MX, xx, sample[,con2])$y
  }
  if (Marginal_Dist1 == "Gaus") {
    cop.sample1.non.con <- qnorm(sample[, con2], mean = as.numeric(marginal_non_con1$estimate[1]),
                                 sd = as.numeric(marginal_non_con1$estimate[2]))
  }
  if (Marginal_Dist1 == "InvG") {
    cop.sample1.non.con <- qinvgauss(sample[, con2], mean = as.numeric(marginal_non_con1$estimate[1]),
                                     shape = as.numeric(marginal_non_con1$estimate[2]))
  }
  if (Marginal_Dist1 == "Logis") {
    cop.sample1.non.con <- qlogis(sample[, con2], location = as.numeric(marginal_non_con1$estimate[1]),
                                  scale = as.numeric(marginal_non_con1$estimate[2]))
  }
  if (Marginal_Dist1 == "LogN") {
    cop.sample1.non.con <- qlnorm(sample[, con2], meanlog = as.numeric(marginal_non_con1$estimate[1]),
                                  sdlog = as.numeric(marginal_non_con1$estimate[2]))
  }
  if (Marginal_Dist1 == "TNorm") {
    cop.sample1.non.con <- qtruncnorm(sample[, con2], a = min(Data_Con1[,
                                                                        con2]), mean = as.numeric(marginal_non_con1$estimate[1]),
                                      sd = as.numeric(marginal_non_con1$estimate[2]))
  }
  if (Marginal_Dist1 == "Twe") {
    cop.sample1.non.con <- qtweedie(sample[, con2], power = marginal_non_con1$p.max,
                                    mu = mean(Data_Con1[, con2]), phi = marginal_non_con1$phi.max)
  }
  if (Marginal_Dist1 == "Weib") {
    cop.sample1.non.con <- qweibull(sample[, con2], shape = as.numeric(marginal_non_con1$estimate[1]),
                                    scale = as.numeric(marginal_non_con1$estimate[2]))
  }
  #Put the realizations that have been transformed to the original scale in a data frame.
  cop.sample1 <- data.frame(cop.sample1.con, cop.sample1.non.con)
  colnames(cop.sample1) <- c("Var1", "Var2")

  ###Simulating sample from the joint distribution (copula+marginals) fit to the sample conditioned on Con2
  #Fit the specified copula family (Copula_Family2) to the observations in Data_Con2.
  obj2 <- BiCopSelect(pobs(Data_Con2[, 1]), pobs(Data_Con2[,2]), familyset = Copula_Family2, selectioncrit = "AIC",
                      indeptest = FALSE, level = 0.05, weights = NA, rotations = TRUE,
                      se = FALSE, presel = TRUE, method = "mle")
  #Simulate a sample from the fitted copula. Out of the sample size 'N' the proportion of the sample from the
  #copula assoicated with Data_Con1 is proportional to the size of Data_Con1 relative to Data_Con2.
  sample <- BiCopSim(round(N * nrow(Data_Con2)/(nrow(Data_Con1) +nrow(Data_Con2)), 0), obj2)

  #Transform the realizations of variable con2 to the original scale using the inverse CDF (quantile function)
  #of the GPD contained in the u2gpd function.
  cop.sample2.con <- u2gpd(sample[, con2], p = 1, th = Thres2, sigma = exp(GPD_con2$coefficients[1]),xi = GPD_con2$coefficients[2])

  #Transform the realizations of variable con1 to the original scale using the inverse CDF (quantile function)
  #of the selected parametric (non-extreme value) distribution (Marginal_Dist2)
  if (Marginal_Dist2 == "BS") {
    cop.sample2.non.con <- qbisa(sample[, con1], as.numeric(Coef(marginal_non_con2)[1]),
                                 as.numeric(Coef(marginal_non_con2)[2]))
  }
  if (Marginal_Dist2 == "Exp") {
    cop.sample2.non.con <- qexp(sample[, con1], rate = as.numeric(marginal_non_con2$estimate[1]))
  }
  if (Marginal_Dist2 == "Gam(2)") {
    cop.sample2.non.con <- qgamma(sample[, con1], shape = as.numeric(marginal_non_con2$estimate[1]),
                                  rate = as.numeric(marginal_non_con2$estimate[2]))
  }
  if(Marginal_Dist2=="Gam(3)"){
    cop.sample2.non.con<-qGG(sample[,con1], mu=exp(marginal_non_con2$mu.coefficients), sigma=exp(marginal_non_con2$sigma.coefficients), nu=marginal_non_con2$nu.coefficients)
  }
  if(Marginal_Dist2=="GamMix(2)"){
    xx <- seq(0, max(Data_Con1[,1])*10, 0.001)
    prob.MX1 <- round(marginal_non_con2$prob[1],3)
    prob.MX2 <- 1 - prob.MX1
    cdf.MX<-pMX(xx, mu=list(mu1=exp(marginal_non_con2$models[[1]]$mu.coefficients), mu2=exp(marginal_non_con2$models[[2]]$mu.coefficients)),
                sigma=list(sigma1=exp(marginal_non_con2$models[[1]]$sigma.coefficients), sigma2=exp(marginal_non_con2$models[[2]]$sigma.coefficients)),
                pi = list(pi1=prob.MX1, pi2=prob.MX2), family=list(fam1="GA", fam2="GA"))
    cop.sample2.non.con <- approx(cdf.MX, xx, sample[,con1])$y
  }
  if(Marginal_Dist2=="GamMix(3)"){
    xx <- seq(0, max(Data_Con1[,1])*10, 0.001)
    prob.MX1 <- round(marginal_non_con2$prob[1],3)
    prob.MX2 <- round(marginal_non_con2$prob[2],3)
    prob.MX3 <- 1 - prob.MX1 - prob.MX2
    cdf.MX<-pMX(xx, mu=list(mu1=exp(marginal_non_con2$models[[1]]$mu.coefficients), mu2=exp(marginal_non_con2$models[[2]]$mu.coefficients), mu3=exp(marginal_non_con2$models[[3]]$mu.coefficients)),
                sigma=list(sigma1=exp(marginal_non_con2$models[[1]]$sigma.coefficients), sigma2=exp(marginal_non_con2$models[[2]]$sigma.coefficients), sigma3=exp(marginal_non_con2$models[[3]]$sigma.coefficients)),
                pi = list(pi1=prob.MX1, pi2=prob.MX2, pi3=prob.MX3), family=list(fam1="GA", fam2="GA", fam3="GA"))
    cop.sample2.non.con <- approx(cdf.MX, xx, sample[,con1])$y
  }
  if (Marginal_Dist2 == "Gaus") {
    cop.sample2.non.con <- qnorm(sample[, con1], mean = as.numeric(marginal_non_con2$estimate[1]),
                                 sd = as.numeric(marginal_non_con2$estimate[2]))
  }
  if (Marginal_Dist2 == "InvG") {
    cop.sample2.non.con <- qinvgauss(sample[, con1], mean = as.numeric(marginal_non_con2$estimate[1]),
                                     shape = as.numeric(marginal_non_con2$estimate[2]))
  }
  if (Marginal_Dist2 == "LogN") {
    cop.sample2.non.con <- qlnorm(sample[, con1], meanlog = as.numeric(marginal_non_con2$estimate[1]),
                                  sdlog = as.numeric(marginal_non_con2$estimate[2]))
  }
  if (Marginal_Dist2 == "Logis") {
    cop.sample2.non.con <- qlogis(sample[, con1], location = as.numeric(marginal_non_con2$estimate[1]),
                                  scale = as.numeric(marginal_non_con2$estimate[2]))
  }
  if (Marginal_Dist2 == "TNorm") {
    cop.sample2.non.con <- qtruncnorm(sample[, con1], a = min(Data_Con2[,
                                                                        con1]), mean = as.numeric(marginal_non_con2$estimate[1]),
                                      sd = as.numeric(marginal_non_con2$estimate[2]))
  }
  if (Marginal_Dist2 == "Twe") {
    cop.sample2.non.con <- qtweedie(sample[, con1], power = marginal_non_con2$p.max,
                                    mu = mean(Data_Con2[, con1]), phi = marginal_non_con2$phi.max)
  }
  if (Marginal_Dist2 == "Weib") {
    cop.sample2.non.con <- qweibull(sample[, con1], shape = as.numeric(marginal_non_con2$estimate[1]),
                                    scale = as.numeric(marginal_non_con2$estimate[2]))
  }
  #Put the realizations that have been transformed to the original scale in a data frame.
  cop.sample2 <- data.frame(cop.sample2.non.con, cop.sample2.con)
  colnames(cop.sample2) <- c("Var1", "Var2")

  #Combine the data frames containg the samples from two joint models (on the original scale)
  cop.sample <- rbind(cop.sample1, cop.sample2)

  Data_Con_Combined <- rbind(round(Data_Con1,DecP),round(Data_Con2,DecP))
  Data_Con_N <- nrow(unique(Data_Con_Combined))

  JRP <- function(x){
  (1 / (Data_Con_N / time.period)) / (length(which(cop.sample[, con1] > x[1] & cop.sample[, con2] > x[2]))/ N)
  }

  y <- data.frame(Var1,Var2)
  RP_Copula <-apply(y,1,JRP)

  return(RP_Copula)
}
