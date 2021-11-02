#' Derives a single or ensemble of bivariate design events
#'
#' Calculates the isoline and relative probability of events on the isoline, given the observational data, for one or more user-specified return periods. Outputs the single "most-likely" design event or an ensemble of possible design events obtained by sampling along the isoline according to these relative probabilities. The design event under the assumption of full dependence is also computed.
#' @param Data Data frame of dimension \code{nx2} containing two co-occurring time series of length \code{n}.
#' @param Data_Con1 Data frame containing the conditional sample (declustered excesses paired with concurrent values of other variable), conditioned on the variable in the first column.
#' @param Data_Con2 Data frame containing the conditional sample (declustered excesses paired with concurrent values of other variable), conditioned on the variable in the second column. Can be obtained using the \code{Con_Sampling_2D} function.
#' @param Thres1 Numeric vector of length one specifying the threshold above which the variable in the first column was sampled in Data_Con1.
#' @param Thres2 Numeric vector of length one specifying the threshold above which the variable in the second column was sampled in Data_Con2.
#' @param Copula_Family1 Numeric vector of length one specifying the copula family used to model the \code{Data_Con1} dataset.
#' @param Copula_Family2 Numeric vector of length one specifying the copula family used to model the \code{Data_Con2} dataset. Best fitting of 40 copulas can be found using the \code{Copula_Threshold_2D} function.
#' @param Marginal_Dist1 Character vector of length one specifying (non-extreme) distribution used to model the marginal distribution of the non-conditioned variable.
#' @param Marginal_Dist2 Character vector of length one specifying (non-extreme) distribution used to model the marginal distribution of the non-conditioned variable.
#' @param Con1 Character vector of length one specifying the name of variable in the first column of \code{Data}.
#' @param Con2 Character vector of length one specifying the name of variable in the second column of \code{Data}.
#' @param mu Numeric vector of length one specifying the (average) occurrence frequency of events in \code{Data}. Default is \code{365.25}, daily data.
#' @param RP Numeric vector specifying the return periods of interest.
#' @param Interval Numeric vector specifying the number of equally spaced points comprising the combined isoline.
#' @param x_lab Character vector specifying the x-axis label.
#' @param y_lab Character vector specifying the y-axis label.
#' @param x_lim_min Numeric vector of length one specifying x-axis minimum. Default is \code{NA}.
#' @param x_lim_max Numeric vector of length one specifying x-axis maximum. Default is \code{NA}.
#' @param y_lim_min Numeric vector of length one specifying y-axis minimum. Default is \code{NA}.
#' @param y_lim_max Numeric vector of length one specifying y-axis maximum. Default is \code{NA}.
#' @param N Numeric vector of length one specifying the size of the sample from the fitted joint distributions used to estimate the density along an isoline. Samples are collected from the two joint distribution with proportions consistent with the total number of extreme events conditioned on each variable. Default is \code{10^6}
#' @param N_Ensemble Numeric vector of length one specifying the number of possible design events sampled along the isoline of interest.
#' @param Sim_Max Numeric vector of length one specifying the maximum value, given as a multiple of the largest observation of each variable, permitted in the sample used to estimate the (relative) probabilities along the isoline.
#' @param Plot_Quantile_Isoline Logical; indicating whether to first plot the quantile isoline. Default is FALSE.
#' @param Isoline_Type Character vector of length one specifying the type of isoline. For isolines obtained using the overlaying method in Bender et al. (2016) use \code{"Combined"} (default). For quantile isoline from the sample conditioned variable \code{Con1}|(\code{Con2}) on use \code{"Con1"}(\code{"Con2"}).
#' @return Plot of all the observations (grey circles) as well as the declustered excesses above Thres1 (blue circles) or Thres2 (blue circles), observations may belong to both conditional samples. Also shown is the isoline associated with \code{RP} contoured according to their relative probability of occurrence on the basis of the sample from the two joint distributions, the "most likely" design event (black diamond), and design event under the assumption of full dependence (black triangle) are also shown in the plot. The function also returns a list comprising the design events assuming full dependence \code{"FullDependence"}, as well as once the dependence between the variables is accounted for the "Most likley" \code{"MostLikelyEvent"} as well as an \code{"Ensemble"} of possible design events.
#' @seealso \code{\link{Copula_Threshold_2D}} \code{\link{Diag_Non_Con}} \code{\link{Diag_Non_Con_Trunc}}
#' @export
#' @examples
#'S22.Rainfall<-Con_Sampling_2D(Data_Detrend=S22.Detrend.df[,-c(1,4)],
#'                              Data_Declust=S22.Detrend.Declustered.df[,-c(1,4)],
#'                              Con_Variable="Rainfall",Thres=0.97)
#'S22.OsWL<-Con_Sampling_2D(Data_Detrend=S22.Detrend.df[,-c(1,4)],
#'                          Data_Declust=S22.Detrend.Declustered.df[,-c(1,4)],
#'                          Con_Variable="OsWL",Thres=0.97)
#'S22.Copula.Rainfall<-Copula_Threshold_2D(Data_Detrend=S22.Detrend.df[,-c(1,4)],
#'                                         Data_Declust=S22.Detrend.Declustered.df[,-c(1,4)],Thres =0.97,
#'                                         y_lim_min=-0.075,y_lim_max=0.25,
#'                                         Upper=c(2,9),Lower=c(2,10),GAP=0.15)$Copula_Family_Var1
#'S22.Copula.OsWL<-Copula_Threshold_2D(Data_Detrend=S22.Detrend.df[,-c(1,4)],
#'                                     Data_Declust=S22.Detrend.Declustered.df[,-c(1,4)],Thres =0.97,
#'                                     y_lim_min=-0.075, y_lim_max =0.25,
#'                                     Upper=c(2,9),Lower=c(2,10),GAP=0.15)$Copula_Family_Var2
#'Design_Event_2D(Data=S22.Detrend.df[,-c(1,4)],
#'                Data_Con1=S22.Rainfall$Data, Data_Con2=S22.OsWL$Data,
#'                Thres1=0.97, Thres2=0.97,
#'                Copula_Family1=S22.Copula.Rainfall, Copula_Family2=S22.Copula.OsWL,
#'                Marginal_Dist1="Logis", Marginal_Dist2="Twe",
#'                RP=100,Interval=10000,N=10,N_Ensemble=10,
#'                Plot_Quantile_Isoline=FALSE)
Design_Event_2D<-function(Data, Data_Con1, Data_Con2, Thres1, Thres2, Copula_Family1, Copula_Family2, Marginal_Dist1, Marginal_Dist2, Con1="Rainfall",Con2="OsWL",mu=365.25, RP,Interval=10000, x_lab="Rainfall (mm)",y_lab="O-sWL (mNGVD 29)",x_lim_min = NA,x_lim_max = NA,y_lim_min = NA,y_lim_max = NA,N=10^6,N_Ensemble=0,Sim_Max=10,Plot_Quantile_Isoline=FALSE,Isoline_Type="Combined"){

  ###Preliminaries

  #Vectors and lists for results
  Quantile_Isoline_1<-vector(mode = "list", length = length(RP))
  names(Quantile_Isoline_1)<-RP
  Quantile_Isoline_2<-vector(mode = "list", length = length(RP))
  names(Quantile_Isoline_2)<-RP
  Isoline<-vector(mode = "list", length = length(RP))
  names(Isoline)<-RP
  Contour<-vector(mode = "list", length = length(RP))
  names(Contour)<-RP
  Ensemble<-vector(mode = "list", length = length(RP))
  names(Ensemble)<-RP
  MostLikelyEvent<-vector(mode = "list", length = length(RP))
  names(MostLikelyEvent)<-RP
  FullDependence<-vector(mode = "list", length = length(RP))
  names(FullDependence)<-RP

  #Remove 1st column of Data if it is a Date or factor object.
  if(class(Data[,1])=="Date" | class(Data[,1])=="factor"){
    Data<-Data[,-1]
  }

  #Find the columns in Data (which should be consistent in terms of column order of the other data input objects) of conditioning variable 1 (Con1) and conditioning variable 2 (Con2).
  con1<-which(names(Data)==Con1)
  con2<-which(names(Data)==Con2)

  ###Fit the 4 marginal distributions (2 GPD and 2 parametric non-extreme value distributions).

  #Fit the GPD to the conditioned variable con1 in Data_Con1.
  GPD_con1<-evm(Data_Con1[,con1], th=quantile(na.omit(Data)[,con1],Thres1) ,penalty = "gaussian",priorParameters = list(c(0, 0), matrix(c(100^2, 0, 0, 0.25), nrow = 2)))

  #Fit the specified marginal distribution (Marginal_Dist1) to the non-conditioned variable con2 in Data_Con1.
  if(Marginal_Dist1 == "BS"){
    bdata2 <- data.frame(shape = exp(-0.5), scale = exp(0.5))
    bdata2 <- transform(bdata2, y = Data_Con1[,con2])
    marginal_non_con1<-vglm(y ~ 1, bisa, data = bdata2, trace = FALSE)
  }
  if(Marginal_Dist1 == "Exp"){
    marginal_non_con1<-fitdistr(Data_Con1[,con2],"exponential")
  }
  if(Marginal_Dist1 == "Gam2"){
    marginal_non_con1<-fitdistr(Data_Con1[,con2], "gamma")
  }
  if(Marginal_Dist1 == "Gam3"){
    data.gamlss<-data.frame(X=Data_Con1[,con2])
    marginal_non_con1 <- tryCatch(gamlss(X~1, data=data.gamlss, family=GG),
                         error = function(e) "error")
  }
  if(Marginal_Dist1 == "GamMix2"){
    data.gamlss<-data.frame(X=Data_Con1[,con2])
    marginal_non_con1 <- tryCatch(gamlssMX(X~1, data=data.gamlss, family=GA, K=2),
                             error = function(e) "error")
  }
  if(Marginal_Dist1 == "GamMix3"){
    data.gamlss<-data.frame(X=Data_Con1[,con2])
    marginal_non_con1 <- tryCatch(gamlssMX(X~1, data=data.gamlss, family=GA, K=3),
                               error = function(e) "error")
  }
  if(Marginal_Dist1 == "Gaus"){
    marginal_non_con1<-fitdistr(Data_Con1[,con2],"normal")
  }
  if(Marginal_Dist1 == "InvG"){
    marginal_non_con1<-fitdist(Data_Con1[,con2], "invgauss", start = list(mean = 5, shape = 1))
  }
  if(Marginal_Dist1 == "Logis"){
    marginal_non_con1<-fitdistr(Data_Con1[,con2], "logistic")
  }
  if(Marginal_Dist1 == "LogN"){
    marginal_non_con1<-fitdistr(Data_Con1[,con2],"lognormal")
  }
  if(Marginal_Dist1 == "TNorm"){
    marginal_non_con1<-fitdistr(Data_Con1[,con2],"normal")
  }
  if(Marginal_Dist1 == "Twe"){
    marginal_non_con1<-tweedie.profile(Data_Con1[,con2] ~ 1,p.vec=seq(1.5, 2.5, by=0.2), do.plot=FALSE)
  }
  if(Marginal_Dist1 == "Weib"){
    marginal_non_con1<-fitdistr(Data_Con1[,con2], "weibull")
  }

  #Fit the GPD to the conditioned variable con2 in Data_Con2.
  GPD_con2<-evm(Data_Con2[,con2], th=quantile(na.omit(Data)[,con2],Thres2) ,penalty = "gaussian",priorParameters = list(c(0, 0), matrix(c(100^2, 0, 0, 0.25), nrow = 2)))

  ##Fit the specified marginal distribution (Marginal_Dist2) to the non-conditioned variable con1 in Data_Con2.
  if(Marginal_Dist2 == "BS"){
    bdata2 <- data.frame(shape = exp(-0.5), scale = exp(0.5))
    bdata2 <- transform(bdata2, y = Data_Con2[,con1])
    marginal_non_con2<-vglm(y ~ 1, bisa, data = bdata2, trace = FALSE)
  }
  if(Marginal_Dist2 == "Exp"){
    marginal_non_con2<-fitdistr(Data_Con2[,con1],"exponential")
  }
  if(Marginal_Dist2 == "Gam2"){
    marginal_non_con2<-fitdistr(Data_Con2[,con1], "gamma")
  }
  if(Marginal_Dist2 == "Gam3"){
    data.gamlss<-data.frame(X=Data_Con2[,con1])
    marginal_non_con2 <- tryCatch(gamlss(X~1, data=data.gamlss, family=GG),
                           error = function(e) "error")
  }
  if(Marginal_Dist2 == "GamMix2"){
    data.gamlss<-data.frame(X=Data_Con2[,con1])
    marginal_non_con2 <- tryCatch(gamlssMX(X~1, data=data.gamlss, family=GA, K=2),
                               error = function(e) "error")
  }
  if(Marginal_Dist2 == "GamMix3"){
    data.gamlss<-data.frame(X=Data_Con2[,con1])
    marginal_non_con2 <- tryCatch(gamlssMX(X~1, data=data.gamlss, family=GA, K=3),
                                  error = function(e) "error")
  }
  if(Marginal_Dist2 == "Gaus"){
    marginal_non_con2<-fitdistr(Data_Con2[,con1],"normal")
  }
  if(Marginal_Dist2 == "InvG"){
    marginal_non_con2<-fitdist(Data_Con2[,con1], "invgauss", start = list(mean = 5, shape = 1))
  }
  if(Marginal_Dist2 == "Logis"){
    marginal_non_con2<-fitdistr(Data_Con2[,con1],"logistic")
  }
  if(Marginal_Dist2 == "LogN"){
    marginal_non_con2<-fitdistr(Data_Con2[,con1],"lognormal")
  }
  if(Marginal_Dist2 == "TNorm"){
    marginal_non_con2<-fitdistr(Data_Con2[,con1],"normal")
  }
  if(Marginal_Dist2 == "Twe"){
    marginal_non_con2<-tweedie.profile(Data_Con2[,con1] ~ 1,p.vec=seq(1.5, 2.5, by=0.2), do.plot=FALSE)
  }
  if(Marginal_Dist2 == "Weib"){
    marginal_non_con2<-fitdistr(Data_Con2[,con1], "weibull")
  }

  ###Generate samples from the copula models to which a kernel density estimate will be applied to estimate relative probabilities along the isoline.

  #Fit the specified copula family (Copula_Family1) to the observations in Data_Con1.
  obj1<-BiCopSelect(pobs(Data_Con1[,1]), pobs(Data_Con1[,2]), familyset=Copula_Family1, selectioncrit = "AIC",
                    indeptest = FALSE, level = 0.05, weights = NA, rotations = TRUE,
                    se = FALSE, presel = TRUE, method = "mle")
  #Simulate a sample from the fitted copula. Out of the sample size 'N' the proportion of the sample from the copula associated with Data_Con1 is proportional to the size of Data_Con1 relative to Data_Con2.
  sample<-BiCopSim(round(N*nrow(Data_Con1)/(nrow(Data_Con1)+nrow(Data_Con2)),0),obj1)
  #Transform the realizations of the conditioned variable con1 to the original scale using inverse cumulative distribution a.k.a. quantile functions (inverse probability integral transform) of the GPD contained in the u2gpd function.
  cop.sample1.con<-u2gpd(sample[,con1], p = 1, th=quantile(na.omit(Data[,con1]),Thres1) , sigma=exp(GPD_con1$coefficients[1]),xi= GPD_con1$coefficients[2])

  #Transform the realizations of the non-conditioned variable con2 to the original scale using the quantile function of the selected parametric (non-extreme value) distribution (Marginal_Dist1).
  if(Marginal_Dist1=="BS"){
    cop.sample1.non.con<-qbisa(sample[,con2], as.numeric(Coef(marginal_non_con1)[1]), as.numeric(Coef(marginal_non_con1)[2]))
  }
  if(Marginal_Dist1=="Exp"){
    cop.sample1.non.con<-qexp(sample[,con2], rate = as.numeric(marginal_non_con1$estimate[1]))
  }
  if(Marginal_Dist1=="Gam2"){
    cop.sample1.non.con<-qgamma(sample[,con2], shape = as.numeric(marginal_non_con1$estimate[1]), rate = as.numeric(marginal_non_con1$estimate[2]))
  }
  if(Marginal_Dist1=="Gam3"){
    cop.sample1.non.con<-qGG(sample[,con2], mu=exp(marginal_non_con1$mu.coefficients), sigma=exp(marginal_non_con1$sigma.coefficients), nu=marginal_non_con1$nu.coefficients)
  }
  if(Marginal_Dist1=="GamMix2"){
    prob.MX1 <- round(marginal_non_con1$prob[1],3)
    prob.MX2 <- 1 - prob.MX1
    cop.sample1.non.con<-qMX(sample[,con2], mu=list(mu1=exp(marginal_non_con1$models[[1]]$mu.coefficients), mu2=exp(marginal_non_con1$models[[2]]$mu.coefficients)),
                             sigma=list(sigma1=exp(marginal_non_con1$models[[1]]$sigma.coefficients), sigma2=exp(marginal_non_con1$models[[2]]$sigma.coefficients)),
                             pi = list(pi1=prob.MX1, pi2=prob.MX2), family=list(fam1="GA", fam2="GA"))
  }
  if(Marginal_Dist1=="GamMix3"){
  prob.MX1 <- round(marginal_non_con1$prob[1],3)
  prob.MX2 <- round(marginal_non_con1$prob[2],3)
  prob.MX3 <- 1 - prob.MX1 - prob.MX2
  cop.sample1.non.con<-qMX(sample[,con2], mu=list(mu1=exp(marginal_non_con1$models[[1]]$mu.coefficients), mu2=exp(marginal_non_con1$models[[2]]$mu.coefficients), mu3=exp(marginal_non_con1$models[[3]]$mu.coefficients)),
                           sigma=list(sigma1=exp(marginal_non_con1$models[[1]]$sigma.coefficients), sigma2=exp(marginal_non_con1$models[[2]]$sigma.coefficients), sigma3=exp(marginal_non_con1$models[[3]]$sigma.coefficients)),
                           pi = list(pi1=prob.MX1, pi2=prob.MX2, pi3=prob.MX3), family=list(fam1="GA", fam2="GA", fam3="GA"))
  }
  if(Marginal_Dist1=="Gaus"){
    cop.sample1.non.con<-qnorm(sample[,con2], mean = as.numeric(marginal_non_con1$estimate[1]), sd = as.numeric(marginal_non_con1$estimate[2]))
  }
  if(Marginal_Dist1=="InvG"){
    cop.sample1.non.con<-qinvgauss(sample[,con2], mean = as.numeric(marginal_non_con1$estimate[1]), shape = as.numeric(marginal_non_con1$estimate[2]))
  }
  if(Marginal_Dist1=="Logis"){
    cop.sample1.non.con<-qlogis(sample[,con2], location = as.numeric(marginal_non_con1$estimate[1]), scale = as.numeric(marginal_non_con1$estimate[2]))
  }
  if(Marginal_Dist1=="LogN"){
    cop.sample1.non.con<-qlnorm(sample[,con2], meanlog = as.numeric(marginal_non_con1$estimate[1]), sdlog = as.numeric(marginal_non_con1$estimate[2]))
  }
  if(Marginal_Dist1=="TNorm"){
    cop.sample1.non.con<-qtruncnorm(sample[,con2], a=min(Data_Con1[,con2]), mean = as.numeric(marginal_non_con1$estimate[1]), sd = as.numeric(marginal_non_con1$estimate[2]))
  }
  if(Marginal_Dist1=="Twe"){
    cop.sample1.non.con<-qtweedie(sample[,con2], power=marginal_non_con1$p.max, mu=mean(Data_Con1[,con2]), phi=marginal_non_con1$phi.max)
  }
  if(Marginal_Dist1=="Weib"){
    cop.sample1.non.con<-qweibull(sample[,con2], shape = as.numeric(marginal_non_con1$estimate[1]), scale=as.numeric(marginal_non_con1$estimate[2]))
  }
  #Put the realizations that have been transformed to the original scale in a data frame
  cop.sample1<-data.frame(cop.sample1.con,cop.sample1.non.con)
  colnames(cop.sample1)<-c("Var1","Var2")

  #Fit the specified copula family (Copula_Family2) to the observations in Data_Con2.
  obj2<-BiCopSelect(pobs(Data_Con2[,1]), pobs(Data_Con2[,2]), familyset=Copula_Family2, selectioncrit = "AIC",
                    indeptest = FALSE, level = 0.05, weights = NA, rotations = TRUE,
                    se = FALSE, presel = TRUE, method = "mle")
  #Simulate a sample from the fitted copula. Out of the sample size 'N' the proportion of the sample from the copula assoicated with Data_Con2 is proportional to the size of Data_Con2 relative to Data_Con1.
  sample<-BiCopSim(round(N*nrow(Data_Con2)/(nrow(Data_Con1)+nrow(Data_Con2)),0),obj2)

  #Transform the realizations of the conditioned variable con2 to the original scale using the inverse CDF (quantile function) of the GPD contained in the u2gpd function.
  cop.sample2.con<-u2gpd(sample[,con2], p = 1, th=quantile(na.omit(Data[,con2]),Thres2) , sigma=exp(GPD_con2$coefficients[1]),xi= GPD_con2$coefficients[2])

  #Transform the realizations of the non-conditioned variable con1 to the original scale using the inverse CDF (quantile function) of the selected parametric (non-extreme value) distribution (Marginal_Dist2).
  if(Marginal_Dist2=="BS"){
    cop.sample2.non.con<-qbisa(sample[,con1], as.numeric(Coef(marginal_non_con2)[1]), as.numeric(Coef(marginal_non_con2)[2]))
  }
  if(Marginal_Dist2=="Exp"){
    cop.sample2.non.con<-qexp(sample[,con1], rate = as.numeric(marginal_non_con2$estimate[1]))
  }
  if(Marginal_Dist2=="Gam2"){
    cop.sample2.non.con<-qgamma(sample[,con1], shape = as.numeric(marginal_non_con2$estimate[1]), rate=as.numeric(marginal_non_con2$estimate[2]))
  }
  if(Marginal_Dist2=="Gam3"){
    cop.sample2.non.con<-qGG(sample[,con1], mu=exp(marginal_non_con2$mu.coefficients), sigma=exp(marginal_non_con2$sigma.coefficients), nu=marginal_non_con2$nu.coefficients)
  }
  if(Marginal_Dist2=="GamMix2"){
    prob.MX1 <- round(marginal_non_con2$prob[1],3)
    prob.MX2 <- 1 - prob.MX1
    cop.sample2.non.con<-qMX(sample[,con1], mu=list(mu1=exp(marginal_non_con2$models[[1]]$mu.coefficients), mu2=exp(marginal_non_con2$models[[2]]$mu.coefficients)),
                             sigma=list(sigma1=exp(marginal_non_con2$models[[1]]$sigma.coefficients), sigma2=exp(marginal_non_con2$models[[2]]$sigma.coefficients)),
                             pi = list(pi1=prob.MX1, pi2=prob.MX2), family=list(fam1="GA", fam2="GA"))
  }
  if(Marginal_Dist2=="GamMix3"){
    prob.MX1 <- round(marginal_non_con2$prob[1],3)
    prob.MX2 <- round(marginal_non_con2$prob[2],3)
    prob.MX3 <- 1 - prob.MX1 - prob.MX2
    cop.sample2.non.con<-qMX(sample[,con1], mu=list(mu1=exp(marginal_non_con2$models[[1]]$mu.coefficients), mu2=exp(marginal_non_con2$models[[2]]$mu.coefficients), mu3=exp(marginal_non_con2$models[[3]]$mu.coefficients)),
                             sigma=list(sigma1=exp(marginal_non_con2$models[[1]]$sigma.coefficients), sigma2=exp(marginal_non_con2$models[[2]]$sigma.coefficients), sigma3=exp(marginal_non_con2$models[[3]]$sigma.coefficients)),
                             pi = list(pi1=prob.MX1, pi2=prob.MX2, pi3=prob.MX3), family=list(fam1="GA", fam2="GA", fam3="GA"))
  }
  if(Marginal_Dist2=="Gaus"){
    cop.sample2.non.con<-qnorm(sample[,con1], mean = as.numeric(marginal_non_con2$estimate[1]), sd=as.numeric(marginal_non_con2$estimate[2]))
  }
  if(Marginal_Dist2=="InvG"){
    cop.sample2.non.con<-qinvgauss(sample[,con1], mean = as.numeric(marginal_non_con2$estimate[1]), shape=as.numeric(marginal_non_con2$estimate[2]))
  }
  if(Marginal_Dist2=="LogN"){
    cop.sample2.non.con<-qlnorm(sample[,con1], meanlog = as.numeric(marginal_non_con2$estimate[1]), sdlog = as.numeric(marginal_non_con2$estimate[2]))
  }
  if(Marginal_Dist2=="Logis"){
    cop.sample2.non.con<-qlogis(sample[,con1], location = as.numeric(marginal_non_con2$estimate[1]), scale=as.numeric(marginal_non_con2$estimate[2]))
  }
  if(Marginal_Dist2=="TNorm"){
    cop.sample1.non.con<-qtruncnorm(sample[,con1], a=min(Data_Con2[,con1]), mean = as.numeric(marginal_non_con2$estimate[1]), sd = as.numeric(marginal_non_con2$estimate[2]))
  }
  if(Marginal_Dist2=="Twe"){
    cop.sample2.non.con<-qtweedie(sample[,con1], power=marginal_non_con2$p.max, mu=mean(Data_Con2[,con1]), phi=marginal_non_con2$phi.max)
  }
  if(Marginal_Dist2=="Weib"){
    cop.sample2.non.con<-qweibull(sample[,con1], shape = as.numeric(marginal_non_con2$estimate[1]), scale=as.numeric(marginal_non_con2$estimate[2]))
  }
  #Put the realizations that have been transformed to the original scale in a data frame.
  cop.sample2<-data.frame(cop.sample2.non.con,cop.sample2.con)
  colnames(cop.sample2)<-c("Var1","Var2")

  #Combine the data frames containg the samples from two copulas (on the original scale)
  cop.sample<-rbind(cop.sample1,cop.sample2)

  ###Deriving the quantile isoline from the sample conditioned on variable 'Con2' i.e. Data_Con1

  for(k in 1:length(RP)){

    #Generate a regular grid on the unit square.
    x<- c(10^(-4),seq(999.9*10^(-4),1-(1*10^(-5)),10^(-3)))
    y<- c(10^(-4),seq(999.9*10^(-4),1-(1*10^(-5)),10^(-3)))
    u<-expand.grid(x,y,10^(-3))
    #Evaluate the copula at each point on the grid.
    u1<-BiCopCDF(u[,1], u[,2], obj1)

    #Calculate the time period spanned by the original dataset in terms of mu (only including occasions where both variables are observed).
    time.period<-nrow(Data[which(is.na(Data[,1])==FALSE & is.na(Data[,2])==FALSE),])/mu
    #Calculate the rate of occurrences of extremes (in terms of mu) in Data_Con1.
    rate<-nrow(Data_Con1)/time.period
    #Calculate the inter-arrival time of extremes (in terms of mu) in Data_Con1.
    EL<-1/rate

    #Define a function which evaluates the return period at a given point (x,y).
    f<-function(x,y){EL/(1-x-y+u1[which(u[,1]==x & u[,2]==y)]) }
    #Evaluate the return period at each point on the grid 'u' (the 'outer' function creates the grid internally using the points on the boundary i.e. the x and y we defined earlier).
    z<- outer(x,y,f)
    #The contourLines function in the grDevices package extracts the isoline with the specified return period - 'RP' in our case.
    xy160<-contourLines(x,y,z,levels= RP[k])
    #Transform the points on the contour to the original scale using the inverse cumulative distribution a.k.a. quantile functions (inverse probability integral transform)
    #Transform the conditioned variable in Data_Con1, Con1 to the original scale using the inverse CDF of the GPD contained in the u2gpd function
    con1.x<-u2gpd(as.numeric(unlist(xy160[[1]][2])), p = 1, th=quantile(na.omit(Data[,con1]),Thres1) , sigma=exp(GPD_con1$coefficients[1]),xi= GPD_con1$coefficients[2] )
    #Transform the non-conditioned variable in Data_Con1, Con2' to the original scale using the quantile function of the selected parameteric (non-extreme value) distributions
    if(Marginal_Dist1=="BS"){
      con1.y<-qbisa(as.numeric(unlist(xy160[[1]][3])),as.numeric(Coef(marginal_non_con1)[1]),as.numeric(Coef(marginal_non_con1)[2]))
    }
    if(Marginal_Dist1=="Exp"){
      con1.y<-qexp(as.numeric(unlist(xy160[[1]][3])),as.numeric(marginal_non_con1$estimate[1]))
    }
    if(Marginal_Dist1=="Gam2"){
      con1.y<-qgamma(as.numeric(unlist(xy160[[1]][3])),as.numeric(marginal_non_con1$estimate[1]),as.numeric(marginal_non_con1$estimate[2]))
    }
    if(Marginal_Dist1=="Gam3"){
      con1.y<-qGG(as.numeric(unlist(xy160[[1]][3])), mu=exp(marginal_non_con1$mu.coefficients), sigma=exp(marginal_non_con1$sigma.coefficients), nu=marginal_non_con1$nu.coefficients)
    }
    if(Marginal_Dist1=="GamMix2"){
      prob.MX1 <- round(marginal_non_con1$prob[1],3)
      prob.MX2 <- 1 - prob.MX1
      con1.y<-qMX(as.numeric(unlist(xy160[[1]][3])), mu=list(mu1=exp(marginal_non_con1$models[[1]]$mu.coefficients), mu2=exp(marginal_non_con1$models[[2]]$mu.coefficients)),
                  sigma=list(sigma1=exp(marginal_non_con1$models[[1]]$sigma.coefficients), sigma2=exp(marginal_non_con1$models[[2]]$sigma.coefficients)),
                  pi = list(pi1=prob.MX1, pi2=prob.MX2), family=list(fam1="GA", fam2="GA"))
    }
    if(Marginal_Dist1=="GamMix3"){
      prob.MX1 <- round(marginal_non_con1$prob[1],3)
      prob.MX2 <- round(marginal_non_con1$prob[2],3)
      prob.MX3 <- 1 - prob.MX1 - prob.MX2
      #prob.MX3 - round(fit.GamMIX3_GA$prob[3],5) < 0.00001
      con1.y<-qMX(as.numeric(unlist(xy160[[1]][3])), mu=list(mu1=exp(marginal_non_con1$models[[1]]$mu.coefficients), mu2=exp(marginal_non_con1$models[[2]]$mu.coefficients), mu3=exp(marginal_non_con1$models[[3]]$mu.coefficients)),
                  sigma=list(sigma1=exp(marginal_non_con1$models[[1]]$sigma.coefficients), sigma2=exp(marginal_non_con1$models[[2]]$sigma.coefficients), sigma3=exp(marginal_non_con1$models[[3]]$sigma.coefficients)),
                  pi = list(pi1=prob.MX1, pi2=prob.MX2, pi3=prob.MX3), family=list(fam1="GA", fam2="GA", fam3="GA"))
    }
    if(Marginal_Dist1=="Gaus"){
      con1.y<-qnorm(as.numeric(unlist(xy160[[1]][3])),as.numeric(marginal_non_con1$estimate[1]),as.numeric(marginal_non_con1$estimate[2]))
    }
    if(Marginal_Dist1=="InvG"){
      con1.y<-qinvgauss(as.numeric(unlist(xy160[[1]][3])),as.numeric(marginal_non_con1$estimate[1]),as.numeric(marginal_non_con1$estimate[2]))
    }
    if(Marginal_Dist1=="Logis"){
      con1.y<-qlogis(as.numeric(unlist(xy160[[1]][3])),as.numeric(marginal_non_con1$estimate[1]),as.numeric(marginal_non_con1$estimate[2]))
    }
    if(Marginal_Dist1=="LogN"){
      con1.y<-qlnorm(as.numeric(unlist(xy160[[1]][3])),as.numeric(marginal_non_con1$estimate[1]),as.numeric(marginal_non_con1$estimate[2]))
    }
    if(Marginal_Dist1=="TNorm"){
      con1.y<-qtruncnorm(as.numeric(unlist(xy160[[1]][3])),a=min(Data_Con1[,con2]),as.numeric(marginal_non_con1$estimate[1]),as.numeric(marginal_non_con1$estimate[2]))
    }
    if(Marginal_Dist1=="Twe"){
      con1.y<-qtweedie(as.numeric(unlist(xy160[[1]][3])), power=marginal_non_con1$p.max, mu=mean(Data_Con1[,con2]), phi=marginal_non_con1$phi.max)
    }
    if(Marginal_Dist1=="Weib"){
      con1.y<-qweibull(as.numeric(unlist(xy160[[1]][3])),as.numeric(marginal_non_con1$estimate[1]),as.numeric(marginal_non_con1$estimate[2]))
    }

    #Linearly interpolate the points at a 0.01 increment  on the x-axis between the smallest and largest x-value on the contour.
    prediction.points<-approx(c(con1.x),c(con1.y),xout=seq(min(con1.x),max(con1.x),0.01))$y
    #Put the results of the interpolation in a data frame.
    prediction.points<-data.frame(seq(min(con1.x),max(con1.x),0.01),prediction.points)

    #Linearly interpolate the points at a 0.01 increment  on the y-axis between the smallest and largest y-value on the contour (as above but with x and y reversed).
    prediction.points.reverse<-approx(c(con1.y),c(con1.x),xout=seq(min(con1.y),max(con1.y),0.01))$y
    #Put the results of the interpolation in a data frame.
    prediction.points.reverse<-data.frame(seq(min(con1.y),max(con1.y),0.01),prediction.points.reverse)

    #Combine the two data frames derived above - ordering the rows in terms of the magnitudes of the x-values.
    con1.prediction.points.ALL<-data.frame(c(prediction.points[,1],prediction.points.reverse[,2])[order((c(prediction.points[,1],prediction.points.reverse[,2])))],c(prediction.points[,2],prediction.points.reverse[,1])[order((c(prediction.points[,1],prediction.points.reverse[,2])))])
    colnames(con1.prediction.points.ALL)<-c(names(Data)[1],names(Data)[2])

    ###Deriving the quantile isoline from the sample conditioned on variable 'Con2' i.e. Data_Con2.

    #Generate a regular grid on the unit square.
    x<- c(10^(-4),seq(999.9*10^(-4),1-(1*10^(-5)),10^(-3)))
    y<- c(10^(-4),seq(999.9*10^(-4),1-(1*10^(-5)),10^(-3)))
    u<-expand.grid(x,y,10^(-3))
    #Evaluate the copula at each point on the grid.
    u1<-BiCopCDF(u[,1], u[,2], obj2)

    #Calculate the rate of occurrences of extremes (in terms of mu) in Data_Con2.
    rate<-nrow(Data_Con2)/time.period
    #Calculate the inter-arrival time of extremes (in terms of mu) in Data_Con2.
    EL<-1/rate

    #Define a function which evaluates the return period at a given point (x,y).
    f<-function(x,y){EL/(1-x-y+u1[which(u[,1]==x & u[,2]==y)]) }
    #Evaluate the return period at each point on the grid 'u' (the 'outer' function creates the grid internally using the points on the boundary i.e. the x and y we defined earlier)
    z<- outer(x,y,f)
    #The contourLines function in the grDevices package extracts the isoline with the specified return period - 'RP' in our case.
    xy160<-contourLines(x,y,z,levels= RP[k])

    #Transform the points on the contour to the original scale using the inverse cumulative distributions a.k.a. quantile functions (i.e. using the inverse probability integral transform).
    #Transforming the conditioned variable in Data_Con2, Con2 to the original scale using the inverse CDF of the GPD contained in the u2gpd function.
    con2.y<-u2gpd(as.numeric(unlist(xy160[[1]][3])), p = 1, th=quantile(na.omit(Data[,con2]),Thres2) , sigma=exp(GPD_con2$coefficients[1]),xi= GPD_con2$coefficients[2] )
    #Transform the non-conditioned variable in Data_Con2, Con1' to the original scale using the quantile function of the selected parameteric (non-extreme value) distributions.
    if(Marginal_Dist2=="BS"){
      con2.x<-qbisa(as.numeric(unlist(xy160[[1]][2])), as.numeric(Coef(marginal_non_con2)[1]),as.numeric(Coef(marginal_non_con2)[2]))
    }
    if(Marginal_Dist2=="Exp"){
      con2.x<-qexp(as.numeric(unlist(xy160[[1]][2])), as.numeric(marginal_non_con2$estimate[1]))
    }
    if(Marginal_Dist2=="Gam2"){
      con2.x<-qgamma(as.numeric(unlist(xy160[[1]][2])), shape = as.numeric(marginal_non_con2$estimate[1]), rate = as.numeric(marginal_non_con2$estimate[2]))
    }
    if(Marginal_Dist2=="Gam3"){
      con2.x<-qGG(as.numeric(unlist(xy160[[1]][2])), mu=exp(marginal_non_con2$mu.coefficients), sigma=exp(marginal_non_con2$sigma.coefficients), nu=marginal_non_con2$nu.coefficients)
    }
    if(Marginal_Dist2=="GamMix2"){
      prob.MX1 <- round(marginal_non_con2$prob[1],3)
      prob.MX2 <- 1 - prob.MX1
      con2.x<-qMX(as.numeric(unlist(xy160[[1]][2])), mu=list(mu1=exp(marginal_non_con2$models[[1]]$mu.coefficients), mu2=exp(marginal_non_con2$models[[2]]$mu.coefficients)),
                  sigma=list(sigma1=exp(marginal_non_con2$models[[1]]$sigma.coefficients), sigma2=exp(marginal_non_con2$models[[2]]$sigma.coefficients)),
                  pi = list(pi1=prob.MX1, pi2=prob.MX2), family=list(fam1="GA", fam2="GA"))
    }
    if(Marginal_Dist2=="GamMix3"){
      prob.MX1 <- round(marginal_non_con1$prob[1],3)
      prob.MX2 <- round(marginal_non_con1$prob[2],3)
      prob.MX3 <- 1 - prob.MX1 - prob.MX2
      #prob.MX3 - round(fit.GamMIX3_GA$prob[3],5) < 0.00001
      con2.x<-qMX(as.numeric(unlist(xy160[[1]][3])), mu=list(mu1=exp(marginal_non_con1$models[[1]]$mu.coefficients), mu2=exp(marginal_non_con1$models[[2]]$mu.coefficients), mu3=exp(marginal_non_con1$models[[3]]$mu.coefficients)),
                  sigma=list(sigma1=exp(marginal_non_con1$models[[1]]$sigma.coefficients), sigma2=exp(marginal_non_con1$models[[2]]$sigma.coefficients), sigma3=exp(marginal_non_con1$models[[3]]$sigma.coefficients)),
                  pi = list(pi1=prob.MX1, pi2=prob.MX2, pi3=prob.MX3), family=list(fam1="GA", fam2="GA", fam3="GA"))
    }
    if(Marginal_Dist2=="Gaus"){
      con2.y<-qnorm(as.numeric(unlist(xy160[[1]][2])), as.numeric(marginal_non_con2$estimate[1]), as.numeric(marginal_non_con2$estimate[2]))
    }
    if(Marginal_Dist2=="InvG"){
      con2.x<-qinvgauss(as.numeric(unlist(xy160[[1]][2])), as.numeric(marginal_non_con2$estimate[1]), as.numeric(marginal_non_con2$estimate[2]))
    }
    if(Marginal_Dist2=="Logis"){
      con2.x<-qlogis(as.numeric(unlist(xy160[[1]][2])),as.numeric(marginal_non_con2$estimate[1]),as.numeric(marginal_non_con2$estimate[2]))
    }
    if(Marginal_Dist2=="LogN"){
      con2.x<-qlnorm(as.numeric(unlist(xy160[[1]][2])), as.numeric(marginal_non_con2$estimate[1]), as.numeric(marginal_non_con2$estimate[2]))
    }
    if(Marginal_Dist1=="TNorm"){
      con2.x<-qtruncnorm(as.numeric(unlist(xy160[[1]][2])),a=min(Data_Con2[,con1]),as.numeric(marginal_non_con2$estimate[1]),as.numeric(marginal_non_con2$estimate[2]))
    }
    if(Marginal_Dist2=="Twe"){
      con2.x<-qtweedie(as.numeric(unlist(xy160[[1]][2])), power=marginal_non_con2$p.max, mu=mean(Data_Con2[,con1]), phi=marginal_non_con2$phi.max)
    }
    if(Marginal_Dist2=="Weib"){
      con2.x<-qweibull(as.numeric(unlist(xy160[[1]][2])), as.numeric(marginal_non_con2$estimate[1]), as.numeric(marginal_non_con2$estimate[2]))
    }

    #Linearly interpolate the points at a 0.01 increment on the x-axis between the smallest and largest x-value on the contour.
    prediction.points<-approx(c(con2.x),c(con2.y),xout=seq(min(con2.x),max(con2.x),0.01))$y
    #Put the results of the interpolation in a data frame.
    prediction.points<-data.frame(seq(min(con2.x),max(con2.x),0.01),prediction.points)

    #Linearly interpolate the points at a 0.01 increment on the y-axis between the smallest and largest y-value on the contour (as above but with x and y reversed).
    prediction.points.reverse<-approx(c(con2.y),c(con2.x),xout=seq(min(con2.y),max(con2.y),0.01))$y
    #Put the results of the interpolation in a data frame.
    prediction.points.reverse<-data.frame(seq(min(con2.y),max(con2.y),0.01),prediction.points.reverse)

    #Combine the two data frames derived above - ordering the rows in terms of the magnitudes of the x-values.
    con2.prediction.points.ALL<-data.frame(c(prediction.points[,1],prediction.points.reverse[,2])[order((c(prediction.points[,1],prediction.points.reverse[,2])))],c(prediction.points[,2],prediction.points.reverse[,1])[order((c(prediction.points[,1],prediction.points.reverse[,2])))])
    colnames(con2.prediction.points.ALL)<-c(names(Data)[1],names(Data)[2])
    ##plot(con1.prediction.points.ALL,col=1,xlim=c(9,11),ylim=c(-2,9))

    x_min<-ifelse(is.na(x_lim_min)==T,min(na.omit(Data[,con1])),x_lim_min)
    x_max<-ifelse(is.na(x_lim_max)==T,max(na.omit(Data[,con1])),x_lim_max)
    y_min<-ifelse(is.na(y_lim_min)==T,min(na.omit(Data[,con2])),y_lim_min)
    y_max<-ifelse(is.na(y_lim_max)==T,max(na.omit(Data[,con2])),y_lim_max)

    if(Plot_Quantile_Isoline==T){
      #Plot
      par(mar=c(4.5,4.2,0.5,0.5))
      plot(Data[, con1], Data[, con2], xlim =c(x_min, x_max), ylim = c(y_min, y_max), col = "Light Grey",xlab = x_lab, ylab = y_lab, cex.lab = 1.5, cex.axis = 1.5)
      points(Data_Con1[,con1],Data_Con1[,con2],col=4,cex=1.5)
      points(Data_Con2[,con1],Data_Con2[,con2],col="Red",pch=4,cex=1.5)
      points(con1.prediction.points.ALL,col=4,pch=16,cex=1)
      points(con2.prediction.points.ALL,col="Red",pch=16,cex=1)
    }

    Quantile_Isoline_1[[k]]<-con1.prediction.points.ALL
    Quantile_Isoline_2[[k]]<-con2.prediction.points.ALL

    if(Isoline_Type==Con1){
      ###Combining the two quantile isolines

      #In the following lines of code the maximum y-values at each x-value from the two quantile isolines are extracted

      #Generate a sequence of x-values at a 0.01 increment starting at the minimum and ending at the maximum points from the two conditional isolines. Round to 2 decimal places.
      x.1<-round(seq(min(Data[,1],con1.prediction.points.ALL[,1],con2.prediction.points.ALL[,1],na.rm = T),max(Data[,1],con1.prediction.points.ALL[,1],con2.prediction.points.ALL[,1],na.rm = T),0.01),2)
      #Round the x-values  from both quantile isolines to 2 decimal places
      con1.prediction.points.ALL.Round<-round(con1.prediction.points.ALL[,1],2)
      #con2.prediction.points.ALL.Round<-round(con2.prediction.points.ALL[,1],2)
      #Find the maximum y-value from the two quantile isolines at each x-value in x.1.
      y.1<-numeric(length(x.1))
      for(i in 1:length(x.1)){
        y.1[i]<-max(con1.prediction.points.ALL[,2][which(con1.prediction.points.ALL.Round==x.1[i])])
      }
      #If any y.1 elements are '-Inf' then remove.
      if(any(y.1==-Inf)==T){
        y.1[which(y.1==-Inf)]<-ifelse(y.1[which(y.1==-Inf)]==max(y.1,na.rm=T),
                                      max(c(con1.prediction.points.ALL[,2])),
                                      NA)
      }

     #If any x.1 or y.1 elements are 'NA' then remove.
      if(length(which(is.na(y.1)==T))>0){
        x.1<-x.1[-which(is.na(y.1)==TRUE)]
        y.1<-y.1[-which(is.na(y.1)==TRUE)]
      }
      ##points(x.1,y.1,col=5)
      #In the following lines of code, the maximum x-values at each y-value from the two quantile isolines are extracted

      #Generate a sequence of y-values at a 0.01 increment starting at the minimum and ending at the maximum points from the two conditional isolines. Round to 2 decimal places.
      y.2<-round(seq(min(Data[,2],con1.prediction.points.ALL[,2],na.rm=T),max(Data[,2],con1.prediction.points.ALL[,2],na.rm=T),0.01),2)
      #Round the y-values from both quantile isolines to 2 decimal places
      con1.prediction.points.ALL.Round<-round(con1.prediction.points.ALL[,2],2)
      #con2.prediction.points.ALL.Round<-round(con2.prediction.points.ALL[,2],2)

      #Find the maximum x-value from the two quantile isolines at each x-value in y.2.
      x.2<-numeric(length(y.2))
      for(i in 1:length(y.2)){
        x.2[i]<-max(round(con1.prediction.points.ALL[,1],2)[which(con1.prediction.points.ALL.Round==y.2[i])])
        #con2.prediction.points.ALL[,1][which(con2.prediction.points.ALL.Round==y.2[i])])
      }

      if(any(x.2==-Inf)==T){
        x.2[which(x.2==-Inf)]<-NA
      }
      ##points(x.2,y.2,col=2)
      #If any x.2 or y.2 elements are 'NA' then remove.
      if(length(which(is.na(x.2)==T))>0){
        y.2<-y.2[-which(is.na(x.2)==TRUE)]
        x.2<-x.2[-which(is.na(x.2)==TRUE)]
      }

      if(Isoline_Type==Con2){
        ###Combining the two quantile isolines

        #In the following lines of code the maximum y-values at each x-value from the two quantile isolines are extracted

        #Generate a sequence of x-values at a 0.01 increment starting at the minimum and ending at the maximum points from the two conditional isolines. Round to 2 decimal places.
        x.1<-round(seq(min(Data[,1],con1.prediction.points.ALL[,1],con2.prediction.points.ALL[,1],na.rm = T),max(Data[,1],con2.prediction.points.ALL[,1],na.rm = T),0.01),2)
        #Round the x-values  from both quantile isolines to 2 decimal places
        con2.prediction.points.ALL.Round<-round(con2.prediction.points.ALL[,1],2)
        #Find the maximum y-value from the two quantile isolines at each x-value in x.1.
        y.1<-numeric(length(x.1))
        for(i in 1:length(x.1)){
          y.1[i]<-con2.prediction.points.ALL[,2][which(con2.prediction.points.ALL.Round==x.1[i])]
        }
        #If any y.1 elements are '-Inf' then remove.
        if(any(y.1==-Inf)==T){
          y.1[which(y.1==-Inf)]<-ifelse(y.1[which(y.1==-Inf)]==max(y.1,na.rm=T),
                                        con2.prediction.points.ALL[,2],
                                        NA)
        }

        #If any x.1 or y.1 elements are 'NA' then remove.
        if(length(which(is.na(y.1)==T))>0){
          x.1<-x.1[-which(is.na(y.1)==TRUE)]
          y.1<-y.1[-which(is.na(y.1)==TRUE)]
        }

        #In the following lines of code, the maximum x-values at each y-value from the two quantile isolines are extracted

        #Generate a sequence of y-values at a 0.01 increment starting at the minimum and ending at the maximum points from the two conditional isolines. Round to 2 decimal places.
        y.2<-round(seq(min(Data[,2],con1.prediction.points.ALL[,2],con2.prediction.points.ALL[,2],na.rm=T),max(Data[,2],con1.prediction.points.ALL[,2],con2.prediction.points.ALL[,2],na.rm=T),0.01),2)
        #Round the y-values from both quantile isolines to 2 decimal places
        con2.prediction.points.ALL.Round<-round(con2.prediction.points.ALL[,2],2)

        #abline(h=con2.prediction.points.ALL.Round,col=4)

        #Find the maximum x-value from the two quantile isolines at each x-value in y.2.
        x.2<-numeric(length(y.2))
        for(i in 1:length(y.2)){
          x.2[i]<-max(round(con1.prediction.points.ALL[,1],2)[which(con1.prediction.points.ALL.Round==y.2[i])])
        }

        if(any(x.2==-Inf)==T){
          x.2[which(x.2==-Inf)]<-NA
        }
        ##points(x.2,y.2,col=2)
        #If any x.2 or y.2 elements are 'NA' then remove.
        if(length(which(is.na(x.2)==T))>0){
          y.2<-y.2[-which(is.na(x.2)==TRUE)]
          x.2<-x.2[-which(is.na(x.2)==TRUE)]
        }

        prediction.points.ALL<-data.frame(c(x.1,x.2),c(y.1,y.2))[-1,]
        colnames(prediction.points.ALL)<-c(names(Data)[1],names(Data)[2])

        #Round the points defining the isoline to 2 decimal places.
        prediction.points.ALL[,1]<-round(prediction.points.ALL[,1],2)
        #To ensure each x is only paired withe y value and vice versa we remove any
        prediction.points.ALL<-prediction.points.ALL[!duplicated(prediction.points.ALL[,1]), ]
        #Order the rows in terms of magnitude of the x-values.
        z<-order(prediction.points.ALL[,1])
        prediction.points.ALL.1<-(prediction.points.ALL[z,1]-min(prediction.points.ALL[z,1]))/(max(prediction.points.ALL[z,1])-min(prediction.points.ALL[z,1]))
        prediction.points.ALL.2<-(prediction.points.ALL[z,2]-min(prediction.points.ALL[z,2]))/(max(prediction.points.ALL[z,2])-min(prediction.points.ALL[z,2]))
        #Calculate the distance between adjacent points on the isoline.
        x.diff<-c(0,diff(prediction.points.ALL.1))
        y.diff<-c(0,diff(prediction.points.ALL.2))
        d<-sqrt(x.diff^2 + y.diff^2)
        #Linearly interpolate the x values with respect to their cumulative distance along the isoline.
        v.x<-approx(x=cumsum(d),y=prediction.points.ALL.1,xout=seq(0,sum(d),length.out=Interval))$y*(max(prediction.points.ALL[z,1])-min(prediction.points.ALL[z,1]))+min(prediction.points.ALL[z,1])
        #Linearly interpolate the y values with respect to their cumulative distance along the isoline.
        v.y<-approx(x=cumsum(d),y=prediction.points.ALL.2,xout=seq(0,sum(d),length.out=Interval))$y*(max(prediction.points.ALL[z,2])-min(prediction.points.ALL[z,2]))+min(prediction.points.ALL[z,2])
        ##plot(prediction.points.ALL.1,cumsum(d))
        Iso<-data.frame(v.x,v.y)
        colnames(Iso)<-c(names(Data)[1],names(Data)[2])
        #Put the points composing the isoline into a data frame to form part of the function's output.
        Isoline[[k]] <- data.frame(x=Iso[,1],y=Iso[,2])
        #colnames(Isoline) <- c(names(Data)[1],names(Data)[2])

        ###Estimate the (relative) probabilty of events along the isoline
        #Estimate the (relative) probability of events along the isoline by applying a KDE to 'cop.sample'
        #i.e. the large sample of events generated from the two fitted copulas (with sample sizes proportional
        #to the size of the two conditional samples) and transformed back to the original scale. These probabilities are
        #used as estimates of the relative probability of the points on the isoline according to the original data.
        remove<-which(cop.sample[,1] > Sim_Max*max(Data[,1],na.rm=T) | cop.sample[,2] > Sim_Max*max(Data[,2],na.rm=T))
        if(length(remove)>1){
          cop.sample<-cop.sample[-remove,]
        }
        prediction<-kde(x=cop.sample, eval.points=Iso)$estimate

        #(relative) Probabilities implied by the data for the points composing the isoline. Probabilities are scaled to [0,1].
        Contour[[k]] <- (prediction-min(prediction))/(max(prediction)-min(prediction))

        ###Extract design event(s)

        #Find the 'most likely' design event and add it to the plot (denoted by a diamond).
        #x.MostLikelyEvent.AND[k]<-as.numeric(Iso[which(prediction==max(prediction,na.rm=T)),1])
        MostLikelyEvent.AND<-data.frame(as.numeric(Iso[which(prediction==max(prediction,na.rm=T)),1]),as.numeric(Iso[which(prediction==max(prediction,na.rm=T)),2]))
        colnames(MostLikelyEvent.AND) <- c(names(Data)[1],names(Data)[2])
        MostLikelyEvent[[k]]<-MostLikelyEvent.AND

        #Find the design event under the assumption of full dependence and add it to the plot (denoted by a triangle).
        FullDependence.AND<-data.frame(max(Iso[,1],na.rm=T),max(Iso[-1,2],na.rm=T))
        colnames(FullDependence.AND)<- c(names(Data)[1],names(Data)[2])
        FullDependence[[k]]<-FullDependence.AND
        #Generate a sample of events along the contour. Sample is weighted according to the probabilities
        #given by the KDE estimate for each point on the isoline. Sample size is N_Ensemble.
        sample.AND <- Iso[sample(1:length(prediction[prediction>0]),size = N_Ensemble, replace = TRUE, prob=prediction[prediction>0]),]
        colnames(sample.AND) <- c(names(Data)[1],names(Data)[2])
        #Put the ensemble of design event into a data frame to form part of the function's output.
        Ensemble[[k]] <- data.frame(sample.AND)
        #colnames(Ensemble) <- c(names(Data)[1],names(Data)[2])
      }
    }

    if(Isoline_Type=="Combined"){
      ###Combining the two quantile isolines

      #In the following lines of code the maximum y-values at each x-value from the two quantile isolines are extracted

      #Generate a sequence of x-values at a 0.01 increment starting at the minimum and ending at the maximum points from the two conditional isolines. Round to 2 decimal places.
      x.1<-round(seq(min(Data[,1],con1.prediction.points.ALL[,1],con2.prediction.points.ALL[,1],na.rm = T),max(Data[,1],con1.prediction.points.ALL[,1],con2.prediction.points.ALL[,1],na.rm = T),0.01),2)
      #Round the x-values  from both quantile isolines to 2 decimal places
      con1.prediction.points.ALL.Round<-round(con1.prediction.points.ALL[,1],2)
      con2.prediction.points.ALL.Round<-round(con2.prediction.points.ALL[,1],2)
      #Find the maximum y-value from the two quantile isolines at each x-value in x.1.
      y.1<-numeric(length(x.1))
      for(i in 1:length(x.1)){
        y.1[i]<-max(con1.prediction.points.ALL[,2][which(con1.prediction.points.ALL.Round==x.1[i])],
                    con2.prediction.points.ALL[,2][which(con2.prediction.points.ALL.Round==x.1[i])])
      }
      #If any y.1 elements are '-Inf' then remove.
      if(any(y.1==-Inf)==T){
        y.1[which(y.1==-Inf)]<-ifelse(y.1[which(y.1==-Inf)]==max(y.1,na.rm=T),
                                      max(c(con1.prediction.points.ALL[,2],con2.prediction.points.ALL[,2])),
                                      NA)
      }

      #If any x.1 or y.1 elements are 'NA' then remove.
      if(length(which(is.na(y.1)==T))>0){
        x.1<-x.1[-which(is.na(y.1)==TRUE)]
        y.1<-y.1[-which(is.na(y.1)==TRUE)]
      }

      #In the following lines of code, the maximum x-values at each y-value from the two quantile isolines are extracted

      #Generate a sequence of y-values at a 0.01 increment starting at the minimum and ending at the maximum points from the two conditional isolines. Round to 2 decimal places.
      y.2<-round(seq(min(Data[,2],con1.prediction.points.ALL[,2],con2.prediction.points.ALL[,2],na.rm=T),max(Data[,2],con1.prediction.points.ALL[,2],con2.prediction.points.ALL[,2],na.rm=T),0.01),2)
      #Round the y-values from both quantile isolines to 2 decimal places
      con1.prediction.points.ALL.Round<-round(con1.prediction.points.ALL[,2],2)
      con2.prediction.points.ALL.Round<-round(con2.prediction.points.ALL[,2],2)

     #Find the maximum x-value from the two quantile isolines at each x-value in y.2.
     x.2<-numeric(length(y.2))
     for(i in 1:length(y.2)){
         x.2[i]<-max(round(con1.prediction.points.ALL[,1],2)[which(con1.prediction.points.ALL.Round==y.2[i])])
     }

     if(any(x.2==-Inf)==T){
         x.2[which(x.2==-Inf)]<-NA
      }

     #If any x.2 or y.2 elements are 'NA' then remove.
        if(length(which(is.na(x.2)==T))>0){
          y.2<-y.2[-which(is.na(x.2)==TRUE)]
          x.2<-x.2[-which(is.na(x.2)==TRUE)]
        }

        prediction.points.ALL<-data.frame(c(x.1,x.2),c(y.1,y.2))[-1,]
        colnames(prediction.points.ALL)<-c(names(Data)[1],names(Data)[2])

        #Round the points defining the isoline to 2 decimal places.
        prediction.points.ALL[,1]<-round(prediction.points.ALL[,1],2)
        #To ensure each x is only paired withe y value and vice versa we remove any
        prediction.points.ALL<-prediction.points.ALL[!duplicated(prediction.points.ALL[,1]), ]
        #Order the rows in terms of magnitude of the x-values.
        z<-order(prediction.points.ALL[,1])
        prediction.points.ALL.1<-(prediction.points.ALL[z,1]-min(prediction.points.ALL[z,1]))/(max(prediction.points.ALL[z,1])-min(prediction.points.ALL[z,1]))
        prediction.points.ALL.2<-(prediction.points.ALL[z,2]-min(prediction.points.ALL[z,2]))/(max(prediction.points.ALL[z,2])-min(prediction.points.ALL[z,2]))
        #Calculate the distance between adjacent points on the isoline.
        x.diff<-c(0,diff(prediction.points.ALL.1))
        y.diff<-c(0,diff(prediction.points.ALL.2))
        d<-sqrt(x.diff^2 + y.diff^2)
        #Linearly interpolate the x values with respect to their cumulative distance along the isoline.
        v.x<-approx(x=cumsum(d),y=prediction.points.ALL.1,xout=seq(0,sum(d),length.out=Interval))$y*(max(prediction.points.ALL[z,1])-min(prediction.points.ALL[z,1]))+min(prediction.points.ALL[z,1])
        #Linearly interpolate the y values with respect to their cumulative distance along the isoline.
        v.y<-approx(x=cumsum(d),y=prediction.points.ALL.2,xout=seq(0,sum(d),length.out=Interval))$y*(max(prediction.points.ALL[z,2])-min(prediction.points.ALL[z,2]))+min(prediction.points.ALL[z,2])
        ##plot(prediction.points.ALL.1,cumsum(d))
        Iso<-data.frame(v.x,v.y)
        colnames(Iso)<-c(names(Data)[1],names(Data)[2])
        #Put the points composing the isoline into a data frame to form part of the function's output.
        Isoline[[k]] <- data.frame(x=Iso[,1],y=Iso[,2])
        #colnames(Isoline) <- c(names(Data)[1],names(Data)[2])

        ###Estimate the (relative) probabilty of events along the isoline
        #Estimate the (relative) probability of events along the isoline by applying a KDE to 'cop.sample'
        #i.e. the large sample of events generated from the two fitted copulas (with sample sizes proportional
        #to the size of the two conditional samples) and transformed back to the original scale. These probabilities are
        #used as estimates of the relative probability of the points on the isoline according to the original data.
        remove<-which(cop.sample[,1] > Sim_Max*max(Data[,1],na.rm=T) | cop.sample[,2] > Sim_Max*max(Data[,2],na.rm=T))
        if(length(remove)>1){
          cop.sample<-cop.sample[-remove,]
        }
        prediction<-kde(x=cop.sample, eval.points=Iso)$estimate

        #(relative) Probabilities implied by the data for the points composing the isoline. Probabilities are scaled to [0,1].
        Contour[[k]] <- (prediction-min(prediction))/(max(prediction)-min(prediction))

        ###Extract design event(s)

        #Find the 'most likely' design event and add it to the plot (denoted by a diamond).
        #x.MostLikelyEvent.AND[k]<-as.numeric(Iso[which(prediction==max(prediction,na.rm=T)),1])
        MostLikelyEvent.AND<-data.frame(as.numeric(Iso[which(prediction==max(prediction,na.rm=T)),1]),as.numeric(Iso[which(prediction==max(prediction,na.rm=T)),2]))
        colnames(MostLikelyEvent.AND) <- c(names(Data)[1],names(Data)[2])
        MostLikelyEvent[[k]]<-MostLikelyEvent.AND

        #Find the design event under the assumption of full dependence and add it to the plot (denoted by a triangle).
        FullDependence.AND<-data.frame(max(Iso[,1],na.rm=T),max(Iso[-1,2],na.rm=T))
        colnames(FullDependence.AND)<- c(names(Data)[1],names(Data)[2])
        FullDependence[[k]]<-FullDependence.AND
        #Generate a sample of events along the contour. Sample is weighted according to the probabilities
        #given by the KDE estimate for each point on the isoline. Sample size is N_Ensemble.
        sample.AND <- Iso[sample(1:length(prediction[prediction>0]),size = N_Ensemble, replace = TRUE, prob=prediction[prediction>0]),]
        colnames(sample.AND) <- c(names(Data)[1],names(Data)[2])
        #Put the ensemble of design event into a data frame to form part of the function's output.
        Ensemble[[k]] <- data.frame(sample.AND)
        #colnames(Ensemble) <- c(names(Data)[1],names(Data)[2])
      }
    }

  ###Plot the isoline

  #Find the minimum and maximum x- and y-axis limits for the plot. If the limits are not specified in the input use the minimum and maximum values of the Data.
  x_min<-ifelse(is.na(x_lim_min)==T,min(na.omit(Data[,con1])),x_lim_min)
  x_max<-ifelse(is.na(x_lim_max)==T,max(na.omit(Data[,con1])),x_lim_max)
  y_min<-ifelse(is.na(y_lim_min)==T,min(na.omit(Data[,con2])),y_lim_min)
  y_max<-ifelse(is.na(y_lim_max)==T,max(na.omit(Data[,con2])),y_lim_max)

  #Plot
  par(mar=c(4.5,4.2,0.5,0.5))
  plot(Data[, con1], Data[, con2], xlim = c(x_min, x_max), ylim = c(y_min, y_max), col = "Light Grey",xlab = x_lab, ylab = y_lab, cex.lab = 1.5, cex.axis = 1.5)
  points(Data_Con1[,con1],Data_Con1[,con2],col=4,cex=1.5)
  points(Data_Con2[,con1],Data_Con2[,con2],col="Red",pch=4,cex=1.5)
  for(k in 1:length(RP)){
    points(Isoline[[k]][,1],Isoline[[k]][,2],col=rev(heat.colors(150))[20:120][1+100*Contour[[k]]],lwd=3,pch=16,cex=1.75)
    if(N_Ensemble>0){
      points(Ensemble[[k]][,1],Ensemble[[k]][,2],col=1,lwd=2,pch=16,cex=1)
    }
    points(MostLikelyEvent[[k]][,1],MostLikelyEvent[[k]][,2],pch=18,cex=1.75)
    text(MostLikelyEvent[[k]][,1],MostLikelyEvent[[k]][,2],paste(RP[k]),col="White",cex=0.5)
    points(FullDependence[[k]][,1],FullDependence[[k]][,2],pch=17,cex=1.75)
    text(FullDependence[[k]][,1],FullDependence[[k]][,2],paste(RP[k]),col="White",cex=0.5)
  }
  #Create a list of outputs.
  res<-list("FullDependence" = FullDependence, "MostLikelyEvent" = MostLikelyEvent, "Ensemble"=Ensemble, "Isoline" = Isoline, "Contour"= Contour, "Quantile_Isoline_1" = Quantile_Isoline_1, "Quantile_Isoline_2" = Quantile_Isoline_2)
  return(res)
}
