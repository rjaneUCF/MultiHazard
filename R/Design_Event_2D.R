#' Derives a single or ensemble of bivariate design events
#'
#' Calculates the single design event under the assumption of full dependence, or once accounting for dependence between variables the single "most-likely" or an ensemble of possible design events.
#'
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
#' @param RP Numeric vector of length one specifying the return period of interest.
#' @param x_lab Character vector specifying the x-axis label.
#' @param y_lab Character vector specifying the y-axis label.
#' @param x_lim_min Numeric vector of length one specifying x-axis minimum. Default is \code{NA}.
#' @param x_lim_max Numeric vector of length one specifying x-axis maximum. Default is \code{NA}.
#' @param y_lim_min Numeric vector of length one specifying y-axis minimum. Default is \code{NA}.
#' @param y_lim_max Numeric vector of length one specifying y-axis maximum. Default is \code{NA}.
#' @param delta Numeric vector of length one specifying of the resolution at which the copula CDF is  evaluated on the \code{[0,1]2} square. Default is \code{0.0001}.
#' @param N Numeric vector of length one specifying the size of the sample from the fitted joint distributions used to estimate the density along an isoline. Samples are collected from the two joint distribution with proportions consistent with the total number of extreme events conditioned on each variable.
#' @param N_Ensemble Numeric vector of length one specifying the number of possible design events sampled along the isoline of interest.
#' @return Plot of all the observations (grey circles) as well as the declustered excesses above Thres1 (blue circles) or Thres2 (blue circles), observations may belong to both conditional samples. Also shown is the isoline associated with \code{RP} contoured according to their relative probability of occurrence on the basis of the sample from the two joint distributions, the "most likely" design event (black diamond), and design event under the assumption of full dependence (black triangle) are also shown in the plot. The function also returns a list comprising the design events assuming full dependence \code{"FullDependence"}, as well as once the dependence between the variables is accounted for the "Most likley" {"MostLikelyEvent"} as well as an {"Ensemble"} of possible design events.
#' @seealso \code{\link{Dataframe_Combine}} \code{\link{Copula_Threshold_2D}} \code{\link{Diag_Non_Con}}  \code{\link{Diag_Non_Con_Trunc}}
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
#'                                    Upper=c(2,9),Lower=c(2,10),GAP=0.15)$Copula_Family_Var2
#'Design_Event_2D(Data=S22.Detrend.df[,-c(1,4)], Data_Con1=S22.Rainfall$Data,
#'                Data_Con2=S22.OsWL$Data, Thres1=0.97, Thres2=0.97,
#'                Copula_Family1=S22.Copula.Rainfall, Copula_Family2=S22.Copula.OsWL,
#'                Marginal_Dist1="Logis", Marginal_Dist2="Twe",RP=100,N=10,N_Ensemble=10)
Design_Event_2D<-function(Data, Data_Con1, Data_Con2, Thres1, Thres2, Copula_Family1, Copula_Family2, Marginal_Dist1, Marginal_Dist2, Con1="Rainfall",Con2="OsWL",mu=365.25, RP,x_lab="Rainfall (mm)",y_lab="O-sWL (mNGVD 29)",x_lim_min = NA,x_lim_max = NA,y_lim_min = NA,y_lim_max = NA,N,N_Ensemble){

  if(class(Data[,1])=="Date" | class(Data[,1])=="factor"){
    Data<-Data[,-1]
  }


  con1<-which(names(Data)==Con1)
  con2<-which(names(Data)==Con2)

  x_min<-ifelse(is.na(x_lim_min)==T,min(na.omit(Data[,con1])),x_lim_min)
  x_max<-ifelse(is.na(x_lim_max)==T,max(na.omit(Data[,con1])),x_lim_max)
  y_min<-ifelse(is.na(y_lim_min)==T,min(na.omit(Data[,con2])),y_lim_min)
  y_max<-ifelse(is.na(y_lim_max)==T,max(na.omit(Data[,con2])),y_lim_max)

  ##Fitting marginal distributions

  if(Marginal_Dist1 == "BS"){
    bdata2 <- data.frame(shape = exp(-0.5), scale = exp(0.5))
    bdata2 <- transform(bdata2, y = Data_Con1[,con2])
    marginal_non_con1<-vglm(y ~ 1, bisa, data = bdata2, trace = FALSE)
  }
  if(Marginal_Dist1 == "Exp"){
    marginal_non_con1<-fitdistr(Data_Con1[,con2],"exponential")
  }
  if(Marginal_Dist1 == "Gam"){
    marginal_non_con1<-fitdistr(Data_Con1[,con2], "gamma")
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

  GPD_con1<-evm(Data_Con1[,con1], th=quantile(na.omit(Data[,con1]),Thres1) ,penalty = "gaussian",priorParameters = list(c(0, 0), matrix(c(100^2, 0, 0, 0.25), nrow = 2)))

  if(Marginal_Dist2 == "BS"){
    bdata2 <- data.frame(shape = exp(-0.5), scale = exp(0.5))
    bdata2 <- transform(bdata2, y = Data_Con2[,con1])
    marginal_non_con2<-vglm(y ~ 1, bisa, data = bdata2, trace = FALSE)
  }
  if(Marginal_Dist2 == "Exp"){
    marginal_non_con2<-fitdistr(Data_Con2[,con1],"exponential")
  }
  if(Marginal_Dist2 == "Gam"){
    marginal_non_con2<-fitdistr(Data_Con2[,con1], "gamma")
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

  GPD_con2<-evm(Data_Con2[,con2], th=quantile(na.omit(Data[,con2]),Thres2) ,penalty = "gaussian",priorParameters = list(c(0, 0), matrix(c(100^2, 0, 0, 0.25), nrow = 2)))

  #Simulating from copulas
  obj1<-BiCopSelect(pobs(Data_Con1[,1]), pobs(Data_Con1[,2]), familyset=Copula_Family1, selectioncrit = "AIC",
                    indeptest = FALSE, level = 0.05, weights = NA, rotations = TRUE,
                    se = FALSE, presel = TRUE, method = "mle")
  sample<-BiCopSim(round(N*nrow(Data_Con1)/(nrow(Data_Con1)+nrow(Data_Con2)),0),obj1)

  if(Marginal_Dist1=="BS"){
    cop.sample1.non.con<-qbisa(sample[,con2], as.numeric(Coef(marginal_non_con1)[1]), as.numeric(Coef(marginal_non_con1)[2]))
  }
  if(Marginal_Dist1=="Exp"){
    cop.sample1.non.con<-qexp(sample[,con2], rate = as.numeric(marginal_non_con1$estimate[1]))
  }
  if(Marginal_Dist1=="Gam"){
    cop.sample1.non.con<-qgamma(sample[,con2], shape = as.numeric(marginal_non_con1$estimate[1]), rate = as.numeric(marginal_non_con1$estimate[2]))
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

  cop.sample1.con<-u2gpd(sample[,con1], p = 1, th=quantile(na.omit(Data[,con1]),Thres1) , sigma=exp(GPD_con1$coefficients[1]),xi= GPD_con1$coefficients[2])
  cop.sample1<-data.frame(cop.sample1.con,cop.sample1.non.con)
  colnames(cop.sample1)<-c("Var1","Var2")

  obj2<-BiCopSelect(pobs(Data_Con2[,1]), pobs(Data_Con2[,2]), familyset=Copula_Family2, selectioncrit = "AIC",
                    indeptest = FALSE, level = 0.05, weights = NA, rotations = TRUE,
                    se = FALSE, presel = TRUE, method = "mle")
  sample<-BiCopSim(round(N*nrow(Data_Con2)/(nrow(Data_Con1)+nrow(Data_Con2)),0),obj2)

  if(Marginal_Dist2=="BS"){
    cop.sample2.non.con<-qbisa(sample[,con1], as.numeric(Coef(marginal_non_con2)[1]), as.numeric(Coef(marginal_non_con2)[2]))
  }
  if(Marginal_Dist2=="Exp"){
    cop.sample2.non.con<-qexp(sample[,con1], rate = as.numeric(marginal_non_con2$estimate[1]))
  }
  if(Marginal_Dist2=="Gam"){
    cop.sample2.non.con<-qgamma(sample[,con1], shape = as.numeric(marginal_non_con2$estimate[1]), rate=as.numeric(marginal_non_con2$estimate[2]))
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

  cop.sample2.con<-u2gpd(sample[,con2], p = 1, th=quantile(na.omit(Data[,con2]),Thres2) , sigma=exp(GPD_con2$coefficients[1]),xi= GPD_con2$coefficients[2])
  cop.sample2<-data.frame(cop.sample2.non.con,cop.sample2.con)
  colnames(cop.sample2)<-c("Var1","Var2")
  cop.sample<-rbind(cop.sample1,cop.sample2)

  #Result vectors
  x.MostLikelyEvent.AND<-numeric(1)
  y.MostLikelyEvent.AND<-numeric(1)
  x.full.dependence<-numeric(1)
  y.full.dependence<-numeric(1)

  #Copula object 1
  u<-expand.grid(c(10^(-4),seq(999.9*10^(-4),1-(1*10^(-5)),10^(-3))),c(10^(-4),seq(999.9*10^(-4),1-(1*10^(-5)),10^(-3))))
  u1<-BiCopCDF(u[,1], u[,2], obj1)

  par(mar=c(4.5,4.2,0.5,0.5))
  plot(Data[, con1], Data[, con2], xlim = c(x_min, x_max), ylim = c(y_min, y_max), col = "Light Grey",xlab = x_lab, ylab = y_lab, cex.lab = 1.5, cex.axis = 1.5)
  points(Data_Con1[,con1],Data_Con1[,con2],col=4,cex=1.5)
  points(Data_Con2[,con1],Data_Con2[,con2],col="Red",pch=4,cex=1.5)

  x<- c(10^(-4),seq(999.9*10^(-4),1-(1*10^(-5)),10^(-3)))
  y<- c(10^(-4),seq(999.9*10^(-4),1-(1*10^(-5)),10^(-3)))
  years<-length(which(is.na(Data[,1])==FALSE & is.na(Data[,2])==FALSE))/mu
  EL<-1/(nrow(Data_Con1)/years)

  f<-function(x,y){EL/(1-x-y+u1[which(u[,1]==x & u[,2]==y)]) }
  z<- outer(x,y,f)

  xy160<-contourLines(x,y,z,levels= RP)

  con1.x<-u2gpd(as.numeric(unlist(xy160[[1]][2])), p = 1, th=quantile(na.omit(Data[,con1]),Thres1) , sigma=exp(GPD_con1$coefficients[1]),xi= GPD_con1$coefficients[2] )

  if(Marginal_Dist1=="BS"){
    con1.y<-qbisa(as.numeric(unlist(xy160[[1]][3])),as.numeric(Coef(marginal_non_con1)[1]),as.numeric(Coef(marginal_non_con2)[2]))
  }
  if(Marginal_Dist1=="Exp"){
    con1.y<-qexp(as.numeric(unlist(xy160[[1]][3])),as.numeric(marginal_non_con1$estimate[1]))
  }
  if(Marginal_Dist1=="Gam"){
    con1.y<-qgamma(as.numeric(unlist(xy160[[1]][3])),as.numeric(marginal_non_con1$estimate[1]),as.numeric(marginal_non_con1$estimate[2]))
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

  prediction.points<-approx(c(con1.x),c(con1.y),xout=seq(min(con1.x),max(con1.x),0.01))$y
  prediction.points<-data.frame(seq(min(con1.x),max(con1.x),0.01),prediction.points)

  prediction.points.reverse<-approx(c(con1.y),c(con1.x),xout=seq(min(con1.y),max(con1.y),0.01))$y
  prediction.points.reverse<-data.frame(seq(min(con1.y),max(con1.y),0.01),prediction.points.reverse)

  con1.prediction.points.ALL<-data.frame(c(prediction.points[,1],prediction.points.reverse[,2])[order((c(prediction.points[,1],prediction.points.reverse[,2])))],c(prediction.points[,2],prediction.points.reverse[,1])[order((c(prediction.points[,1],prediction.points.reverse[,2])))])
  colnames(con1.prediction.points.ALL)<-c(names(Data)[1],names(Data)[2])
  u<-expand.grid(c(10^(-5),seq(999.9*10^(-5),1-(1*10^(-6)),100*10^(-5))),c(10^(-5),seq(999.9*10^(-5),1-(1*10^(-6)),100*10^(-5))))

  #OsWL
  u1<-BiCopCDF(u[,1], u[,2], obj2)
  x<- c(10^(-5),seq(999.9*10^(-5),1-(1*10^(-6)),100*10^(-5)))
  y<- c(10^(-5),seq(999.9*10^(-5),1-(1*10^(-6)),100*10^(-5)))
  EL<-1/(nrow(Data_Con2)/years)

  f<-function(x,y){EL/(1-x-y+u1[which(u[,1]==x & u[,2]==y)]) }
  z<- outer(x,y,f)
  xy160<-contourLines(x,y,z,levels= RP)

  con2.y<-u2gpd(as.numeric(unlist(xy160[[1]][3])), p = 1, th=quantile(na.omit(Data[,con2]),Thres2) , sigma=exp(GPD_con2$coefficients[1]),xi= GPD_con2$coefficients[2] )

  if(Marginal_Dist2=="BS"){
    con2.x<-qbisa(as.numeric(unlist(xy160[[1]][2])), as.numeric(Coef(marginal_non_con2)[1]),as.numeric(Coef(marginal_non_con2)[2]))
  }
  if(Marginal_Dist2=="Exp"){
    con2.x<-qexp(as.numeric(unlist(xy160[[1]][2])), as.numeric(marginal_non_con2$estimate[1]))
  }
  if(Marginal_Dist2=="Gam"){
    con2.x<-qgamma(as.numeric(unlist(xy160[[1]][2])), shape = as.numeric(marginal_non_con2$estimate[1]), rate = as.numeric(marginal_non_con2$estimate[2]))
  }
  if(Marginal_Dist2=="Gaus"){
    con2.x<-qnorm(as.numeric(unlist(xy160[[1]][2])), as.numeric(marginal_non_con2$estimate[1]), as.numeric(marginal_non_con2$estimate[2]))
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
  if(Marginal_Dist2=="Wei"){
    con2.x<-qweibull(as.numeric(unlist(xy160[[1]][2])), as.numeric(marginal_non_con2$estimate[1]), as.numeric(marginal_non_con2$estimate[2]))
  }

  prediction.points<-approx(c(con2.x),c(con2.y),xout=seq(min(con2.x),max(con2.x),0.01))$y
  prediction.points<-data.frame(seq(min(con2.x),max(con2.x),0.01),prediction.points)

  prediction.points.reverse<-approx(c(con2.y),c(con2.x),xout=seq(min(con2.y),max(con2.y),0.01))$y
  prediction.points.reverse<-data.frame(seq(min(con2.y),max(con2.y),0.01),prediction.points.reverse)

  con2.prediction.points.ALL<-data.frame(c(prediction.points[,1],prediction.points.reverse[,2])[order((c(prediction.points[,1],prediction.points.reverse[,2])))],c(prediction.points[,2],prediction.points.reverse[,1])[order((c(prediction.points[,1],prediction.points.reverse[,2])))])
  colnames(con2.prediction.points.ALL)<-c(names(Data)[1],names(Data)[2])


  x<-seq(0,max(con1.prediction.points.ALL[,1],con2.prediction.points.ALL[,1],20),0.01)

  con1.prediction.points.ALL.Round<-round(con1.prediction.points.ALL[,1],2)
  con2.prediction.points.ALL.Round<-round(con2.prediction.points.ALL[,1],2)

  y<-numeric(length(x))
  for(i in 1:length(x)){
    y[i]<-max(con1.prediction.points.ALL[,2][which(con1.prediction.points.ALL.Round==x[i])],
              con2.prediction.points.ALL[,2][which(con2.prediction.points.ALL.Round==x[i])])
  }

  if(any(y==-Inf)==T){
    y[which(y==-Inf)]<-NA
  }

  if(length(which(is.na(y)==T))>0){
    x<-x[-which(is.na(y)==TRUE)]
    y<-y[-which(is.na(y)==TRUE)]
  }

  y.2<-round(seq(min(con1.prediction.points.ALL[,2],con2.prediction.points.ALL[,2]),max(con1.prediction.points.ALL[,2],con2.prediction.points.ALL[,2]),0.01),2)

  con1.prediction.points.ALL.Round<-round(con1.prediction.points.ALL[,2],2)
  con2.prediction.points.ALL.Round<-round(con2.prediction.points.ALL[,2],2)

  x.2<-numeric(length(y.2))
  for(i in 1:length(y.2)){
    x.2[i]<-max(con1.prediction.points.ALL[,1][which(con1.prediction.points.ALL.Round==y.2[i])],
                con2.prediction.points.ALL[,1][which(con2.prediction.points.ALL.Round==y.2[i])])
  }

  if(any(x.2==-Inf)==T){
    x.2[which(x.2==-Inf)]<-NA
  }

  if(length(which(is.na(x.2)==T))>0){
    y.2<-y.2[-which(is.na(x.2)==TRUE)]
    x.2<-x.2[-which(is.na(x.2)==TRUE)]
  }
  x1<-x
  y1<-y

  points(x.2,y.2,col="Green")
  points(x,y,col="Green")

  prediction.points.ALL<-data.frame(c(x,x.2),c(y,y.2))[-1,]
  #prediction.points.ALL<-data.frame(c(x,Data[, con1],na.rm=T)),c(y,min(y)))[-1,]
  colnames(prediction.points.ALL)<-c(names(Data)[1],names(Data)[2])

  prediction.points.ALL<-prediction.points.ALL[!duplicated(prediction.points.ALL[,1:2]), ]

  prediction<-kde(x=cop.sample, eval.points=prediction.points.ALL)$estimate
  k=1
  #lines(prediction.points.ALL[,1],prediction.points.ALL[,2],col=max(rev(heat.colors(150))[20:120][1+100*((prediction-min(prediction))/(max(prediction)-min(prediction)))]),lwd=10)
  points(prediction.points.ALL[,1],prediction.points.ALL[,2],col=rev(heat.colors(150))[20:120][1+100*((prediction-min(prediction))/(max(prediction)-min(prediction)))],lwd=3,pch=16,cex=1.75)

  prediction.points.ALL<-data.frame(x1,y1)[-1,]
  #prediction.points.ALL<-data.frame(c(x,Data[, con1],na.rm=T)),c(y,min(y)))[-1,]
  colnames(prediction.points.ALL)<-c(names(Data)[1],names(Data)[2])

  prediction.points.ALL<-prediction.points.ALL[!duplicated(prediction.points.ALL[,1:2]), ]
  prediction<-kde(x=cop.sample, eval.points=prediction.points.ALL)$estimate

  x.MostLikelyEvent.AND[k]<-as.numeric(prediction.points.ALL[which(prediction==max(prediction)),][1])
  y.MostLikelyEvent.AND[k]<-as.numeric(prediction.points.ALL[which(prediction==max(prediction)),][2])
  points(x.MostLikelyEvent.AND[k],y.MostLikelyEvent.AND[k],pch=18,cex=1.75)

  x.full.dependence[k]<-max(x)
  y.full.dependence[k]<-max(y[-1])
  points(x.full.dependence[k],y.full.dependence[k],pch=17,cex=1.75)

  FullDependence<-data.frame(x.full.dependence,y.full.dependence)
  colnames(FullDependence) <- c(names(Data)[1],names(Data)[2])

  MostLikelyEvent<-data.frame(x.MostLikelyEvent.AND,y.MostLikelyEvent.AND)
  colnames(MostLikelyEvent) <- c(names(Data)[1],names(Data)[2])

  sample.AND <- sample(1:length(prediction[prediction>0]),size = N_Ensemble, replace = TRUE, prob=prediction[prediction>0])
  Ensemble <- data.frame(prediction.points.ALL[sample.AND,])
  colnames(Ensemble) <- c(names(Data)[1],names(Data)[2])

  Isoline <- data.frame(x=prediction.points.ALL[,1],y=prediction.points.ALL[,2])
  colnames(Isoline) <- c(names(Data)[1],names(Data)[2])

  Contour <- (prediction-min(prediction))/(max(prediction)-min(prediction))

  res<-list("FullDependence" = FullDependence, "MostLikelyEvent" = MostLikelyEvent, "Ensemble"=Ensemble, "Isoline" = Isoline, "Contour"= Contour)
  return(res)
}
