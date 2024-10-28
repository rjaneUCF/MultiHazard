#' Derives a single or ensemble of bivariate design events
#'
#' Calculates the isoline and relative probability of events on the isoline, where the data contains events from two populations. Outputs the single "most-likely" design event or an ensemble of possible design events obtained by sampling along the isoline according to these relative probabilities. The design event under the assumption of full dependence is also computed. Isoline is derived by calculating annual exceedance probabilities from both copula models on a user-specified grid rather by overlaying the partial isolines from the two copula models as in \code{Design_Event_2D}.
#' @param Data Data frame of dimension \code{nx2} containing two co-occurring time series of length \code{n}.
#' @param Data_Con1 Data frame containing the conditional sample (declustered excesses paired with concurrent values of other variable), conditioned on the variable in the first column for population 1.
#' @param Data_Con2 Data frame containing the conditional sample (declustered excesses paired with concurrent values of other variable), conditioned on the variable in the second column for population 1.
#' @param Data_Con3 Data frame containing the conditional sample (declustered excesses paired with concurrent values of other variable), conditioned on the variable in the first column for population 2.
#' @param Data_Con4 Data frame containing the conditional sample (declustered excesses paired with concurrent values of other variable), conditioned on the variable in the second column for population 2.
#' @param u1 Numeric vector of length one specifying the threshold, expressed as a quantile, above which the variable in the first column was sampled in \code{Data_Con1}.
#' @param u2 Numeric vector of length one specifying the threshold, expressed as a quantile, above which the variable in the second column was sampled in \code{Data_Con2}.
#' @param u3 Numeric vector of length one specifying the threshold, expressed as a quantile, above which the variable in the first column was sampled in \code{Data_Con3}.
#' @param u4 Numeric vector of length one specifying the threshold, expressed as a quantile, above which the variable in the second column was sampled in \code{Data_Con4}.
#' @param Thres1 Numeric vector of length one specifying the threshold above which the variable in the first column was sampled in \code{Data_Con1}. Only one of \code{u1} and \code{Thres1} should be supplied. Default is \code{NA}.
#' @param Thres2 Numeric vector of length one specifying the threshold above which the variable in the second column was sampled in \code{Data_Con2}. Only one of \code{u2} and \code{Thres2} should be supplied. Default is \code{NA}.
#' @param Thres3 Numeric vector of length one specifying the threshold above which the variable in the first column was sampled in \code{Data_Con3}. Only one of \code{u3} and \code{Thres3} should be supplied. Default is \code{NA}.
#' @param Thres4 Numeric vector of length one specifying the threshold above which the variable in the second column was sampled in \code{Data_Con4}. Only one of \code{u4} and \code{Thres4} should be supplied. Default is \code{NA}.
#' @param N_Both_1 Numeric vector of length one specifying the number of data points in population 1 that feature in both conditional samples.
#' @param N_Both_2 Numeric vector of length one specifying the number of data points in population 1 that feature in both conditional samples.
#' @param Copula_Family1 Numeric vector of length one specifying the copula family used to model the \code{Data_Con1} dataset.
#' @param Copula_Family2 Numeric vector of length one specifying the copula family used to model the \code{Data_Con2} dataset.
#' @param Copula_Family3 Numeric vector of length one specifying the copula family used to model the \code{Data_Con3} dataset.
#' @param Copula_Family4 Numeric vector of length one specifying the copula family used to model the \code{Data_Con4} dataset.
#' @param Marginal_Dist1 Character vector of length one specifying (non-extreme) distribution used to model the marginal distribution of the non-conditioned variable in \code{Data_Con1}.
#' @param Marginal_Dist2 Character vector of length one specifying (non-extreme) distribution used to model the marginal distribution of the non-conditioned variable in \code{Data_Con2}.
#' @param Marginal_Dist3 Character vector of length one specifying (non-extreme) distribution used to model the marginal distribution of the non-conditioned variable in \code{Data_Con3}.
#' @param Marginal_Dist4 Character vector of length one specifying (non-extreme) distribution used to model the marginal distribution of the non-conditioned variable in \code{Data_Con4}.
#' @param Marginal_Dist1_Par Character vector of length one specifying (non-extreme) parameters of  \code{Marginal_Dist1}. Default is \code{NA}, specified distribution is fit within the procedure.
#' @param Marginal_Dist2_Par Character vector of length one specifying (non-extreme) parameters of  \code{Marginal_Dist2}. Default is \code{NA}, specified distribution is fit within the procedure.
#' @param Marginal_Dist3_Par Character vector of length one specifying (non-extreme) parameters of  \code{Marginal_Dist3}. Default is \code{NA}, specified distribution is fit within the procedure.
#' @param Marginal_Dist4_Par Character vector of length one specifying (non-extreme) parameters of  \code{Marginal_Dist4}. Default is \code{NA}, specified distribution is fit within the procedure.
#' @param Con1 Character vector of length one specifying the name of variable in the first column of \code{Data_Con1}.
#' @param Con2 Character vector of length one specifying the name of variable in the second column of \code{Data_Con2}.
#' @param Con1 Character vector of length one specifying the name of variable in the first column of \code{Data_Con3}.
#' @param Con2 Character vector of length one specifying the name of variable in the second column of \code{Data_Con4}.
#' @param GPD1 Output of \code{GPD_Fit} applied to variable \code{con1} i.e., GPD fit \code{con1}. Default \code{NA}. Only one of \code{u1}, \code{Thres1}, \code{GPD1} and \code{Tab1} is required.
#' @param GPD2 Output of \code{GPD_Fit} applied to variable \code{con2} i.e., GPD fit \code{con2}. Default \code{NA}. Only one of \code{u2}, \code{Thres2}, \code{GPD2} and \code{Tab2} is required.
#' @param GPD3 Output of \code{GPD_Fit} applied to variable \code{con3} i.e., GPD fit \code{con3}. Default \code{NA}. Only one of \code{u3}, \code{Thres3}, \code{GPD3} and \code{Tab3} is required.
#' @param GPD4 Output of \code{GPD_Fit} applied to variable \code{con4} i.e., GPD fit \code{con4}. Default \code{NA}. Only one of \code{u4}, \code{Thres4}, \code{GPD4} and \code{Tab4} is required.
#' @param Rate_Con1 Numeric vector of length one specifying the occurrence rate of observations in \code{Data_Con1}. Default is \code{NA}.
#' @param Rate_Con2 Numeric vector of length one specifying the occurrence rate of observations in \code{Data_Con2}. Default is \code{NA}.
#' @param Rate_Con3 Numeric vector of length one specifying the occurrence rate of observations in \code{Data_Con3}. Default is \code{NA}.
#' @param Rate_Con4 Numeric vector of length one specifying the occurrence rate of observations in \code{Data_Con4}. Default is \code{NA}.
#' @param Tab1 Data frame specifying the return periods of variable \code{con1}, when conditioning on \code{con1}. First column specifies the return period and the second column gives the corresponding levels. First row must contain the return level of \code{con1} for the inter-arrival time (1/rate) of the sample.
#' @param Tab2 Data frame specifying the return periods of variable \code{con2}, when conditioning on \code{con2}. First column specifies the return period and the second column gives the corresponding levels. First row must contain the return level of \code{con2} for the inter-arrival time (1/rate) of the sample.
#' @param Tab3 Data frame specifying the return periods of variable \code{con3}, when conditioning on \code{con3}. First column specifies the return period and the second column gives the corresponding levels. First row must contain the return level of \code{con3} for the inter-arrival time (1/rate) of the sample.
#' @param Tab4 Data frame specifying the return periods of variable \code{con4}, when conditioning on \code{con4}. First column specifies the return period and the second column gives the corresponding levels. First row must contain the return level of \code{con4} for the inter-arrival time (1/rate) of the sample.
#' @param mu Numeric vector of length one specifying the (average) occurrence frequency of events in \code{Data}. Default is \code{365.25}, daily data.
#' @param GPD_Bayes Logical; indicating whether to use a Bayesian approach to estimate GPD parameters. This involves applying a penalty to the likelihood to aid in the stability of the optimization procedure. Default is \code{FALSE}.
#' @param RP Numeric vector specifying the return periods of interest.
#' @param Decimal_Palace Numeric vector specifying the number of decimal places to which to specify the isoline. Default is \code{2}
#' @param Grid_x_min Numeric vector of length one specifying the minimum value of the variable in first column of \code{Data} contained in the grid.
#' @param Grid_x_max Numeric vector of length one specifying the maximum value of the variable in first column of \code{Data} contained in the grid.
#' @param Grid_x_min Numeric vector of length one specifying the minimum value of the variable in second column of \code{Data} contained in the grid.
#' @param Grid_x_max Numeric vector of length one specifying the maximum value of the variable in second column of \code{Data} contained in the grid.
#' @param Grid_x_interval Numeric vector of length one specifying the resolution of the grid in terms of the variable in first column of \code{Date}. Default is an interval \code{2} of between consecutive values.
#' @param Grid_y_interval Numeric vector of length one specifying the resolution of the grid in terms of the variable in second column of \code{Date}. Default is an interval \code{0.1} of between consecutive values.
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
#' @param Plot_Quantile_Isoline Logical; indicating whether to first plot the quantile isoline. Default is \code{FALSE}.
#' @param Isoline_Type Character vector of length one specifying the type of isoline. For isolines obtained using the overlaying method in Bender et al. (2016) use \code{"Combined"} (default). For quantile isoline from the sample conditioned on variable \code{Con1}|(\code{Con2}) use \code{"Con1"}(\code{"Con2"}).
#' @return Plot of all the observations (grey circles) as well as the declustered excesses above Thres1 (blue circles) or Thres2 (blue circles), observations may belong to both conditional samples. Also shown is the isoline associated with \code{RP} contoured according to their relative probability of occurrence on the basis of the sample from the two joint distributions, the "most likely" design event (black diamond), and design event under the assumption of full dependence (black triangle) are also shown in the plot. The function also returns a list comprising the design events assuming full dependence \code{"FullDependence"}, as well as once the dependence between the variables is accounted for the "Most likley" \code{"MostLikelyEvent"} as well as an \code{"Ensemble"} of possible design events and relative probabilities of events on the isoline \code{Contour}. The quantile isolines with \code{Quantile_Isoline_1} and \code{Quantile_Isoline_2}, and GPD thresholds with \code{Threshold_1} and \code{Threshold_2}.
#' @seealso \code{\link{Copula_Threshold_2D}} \code{\link{Diag_Non_Con}} \code{\link{Diag_Non_Con_Trunc}}
#' @export
Design_Event_2D_Multi_Pop<-function(Data, Data_Con1, Data_Con2, Data_Con3, Data_Con4, u1, u2, u3, u4, Thres1=NA, Thres2=NA, Thres3=NA, Thres4=NA, N_Both_1, N_Both_2, Copula_Family1, Copula_Family2, Copula_Family3, Copula_Family4, Marginal_Dist1, Marginal_Dist2, Marginal_Dist3, Marginal_Dist4, Marginal_Dist1_Par=NA, Marginal_Dist2_Par=NA, Marginal_Dist3_Par=NA, Marginal_Dist4_Par=NA, Con1="Rainfall",Con2="OsWL", Con3="Rainfall", Con4="OsWL", GPD1=NA, GPD2=NA, GPD3=NA, GPD4=NA, Rate_Con1=NA, Rate_Con2=NA, Rate_Con3=NA, Rate_Con4=NA, Tab1= NA, Tab2 = NA, Tab3 = NA, Tab4 = NA, mu=365.25, GPD_Bayes=FALSE, Decimal_Place=2, Grid_x_min = NA ,Grid_x_max = NA, Grid_y_min = NA, Grid_y_max = NA, Grid_x_interval=NA, Grid_y_interval=NA, RP, Interval=10000, End=F, Resolution="Low", x_lab="Rainfall (mm)",y_lab="O-sWL (mNGVD 29)",x_lim_min = NA,x_lim_max = NA,y_lim_min = NA,y_lim_max = NA,Isoline_Probs="Sample", N=10^6,N_Ensemble=0,Sim_Max=10,Plot_Quantile_Isoline=FALSE){

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
  if(class(Data[,1])=="Date" | class(Data[,1])=="factor" | class(Data[,1])[1]=="POSIXct"){
    Data<-Data[,-1]
  }

  #Define the grid over which to calculate Annual Excedence Probabilities
  Grid_x_min = ifelse(is.na(Grid_x_min),min(Data[,1],na.rm=T),Grid_x_min)
  Grid_x_max = 2*ifelse(is.na(Grid_x_max),max(Data[,1],na.rm=T),Grid_x_max)
  Grid_y_min = ifelse(is.na(Grid_y_min),min(Data[,2],na.rm=T),Grid_y_min)
  Grid_y_max = 2*ifelse(is.na(Grid_y_max),max(Data[,2],na.rm=T),Grid_y_max)
  Grid_x_interval = ifelse(is.na(Grid_x_interval),2,Grid_x_interval)
  Grid_y_interval = ifelse(is.na(Grid_y_interval),0.1,Grid_y_interval)

  # Creating a grid over the event space
  var1<- seq(Grid_x_min,Grid_x_max,Grid_x_interval)
  var2<- seq(Grid_y_min,Grid_y_max,Grid_y_interval)
  Pgrid<-expand.grid(var1,var2)

  #Find the columns in Data (which should be consistent in terms of column order of the other data input objects) of conditioning variable 1 (Con1) and conditioning variable 2 (Con2).
  con1<-which(names(Data)==Con1)
  con2<-which(names(Data)==Con2)
  con3<-which(names(Data)==Con3)
  con4<-which(names(Data)==Con4)

  #Find the occurrence rates of the two conditional samples

  #Calculate the time period spanned by the original dataset
  time.period<-nrow(Data[which(is.na(Data[,1])==FALSE & is.na(Data[,2])==FALSE),])/mu

  #Calculate the rate of occurrences of extremes (in terms of mu) in Data_Con1.
  if(is.na(Rate_Con1)==T){
    Rate_Con1<-(nrow(Data_Con1)+N_Both_1/2)/time.period
  }

  #Calculate the rate of occurrences of extremes (in terms of mu) in Data_Con2.
  if(is.na(Rate_Con2)==T){
    Rate_Con2<-(nrow(Data_Con2)+N_Both_1/2)/time.period
  }

  #Calculate the rate of occurrences of extremes (in terms of mu) in Data_Con3.
  if(is.na(Rate_Con3)==T){
    Rate_Con3<-(nrow(Data_Con3)+N_Both_2/2)/time.period
  }

  #Calculate the rate of occurrences of extremes (in terms of mu) in Data_Con4.
  if(is.na(Rate_Con4)==T){
    Rate_Con4<-(nrow(Data_Con4)+N_Both_2/2)/time.period
  }


  ###Fit the 4 marginal distributions (2 GPD and 2 parametric non-extreme value distributions) to the samples from each population
  ###i.e, 8 distributions in total.

  #Fit the GPD to the conditioned variable con1 in Data_Con1.
  if(is.na(GPD1[[1]][1]) & is.na(Thres1) & is.na(Tab1[[1]][1])){
    Thres1<-quantile(na.omit(Data[,con1]),u1)
  }

  if(is.na(GPD1[[1]][1]) & GPD_Bayes & is.na(Tab1[[1]][1])){
    GPD_con1<-evm(Data_Con1[,con1], th = Thres1,penalty = "gaussian",priorParameters = list(c(0, 0), matrix(c(100^2, 0, 0, 0.25), nrow = 2)))
  }
  if(is.na(GPD1[[1]][1]) & !GPD_Bayes & is.na(Tab1[[1]][1])){
    GPD_con1<-evm(Data_Con1[,con1], th = Thres1)
  }


  #Fit the specified marginal distribution (Marginal_Dist1) to the non-conditioned variable con2 in Data_Con1.
  if(Marginal_Dist1 == "BS"){
    if(is.na(Marginal_Dist1_Par)==T){
      bdata2 <- data.frame(shape = exp(-0.5), scale = exp(0.5))
      bdata2 <- transform(bdata2, y = Data_Con1[,con2])
      marginal_non_con1<-vglm(y ~ 1, bisa, data = bdata2, trace = FALSE)
    }else{
      marginal_non_con1<-Marginal_Dist1_Par
    }
  }
  if(Marginal_Dist1 == "Exp"){
    if(is.na(Marginal_Dist1_Par)==T){
      marginal_non_con1<-fitdistr(Data_Con1[,con2],"exponential")
    }else{
      marginal_non_con1<-Marginal_Dist1_Par
    }
  }
  if(Marginal_Dist1 == "Gam(2)"){
    if(is.na(Marginal_Dist1_Par)==T){
      marginal_non_con1<-fitdistr(Data_Con1[,con2], "gamma")
    }else{
      marginal_non_con1<-Marginal_Dist1_Par
    }
  }
  if(Marginal_Dist1 == "Gam(3)"){
    if(is.na(Marginal_Dist1_Par)==T){
      data.gamlss<-data.frame(X=Data_Con1[,con2])
      marginal_non_con1 <- tryCatch(gamlss(X~1, data=data.gamlss, family=GG),
                                    error = function(e) "error")
    }else{
      marginal_non_con1<-Marginal_Dist1_Par
    }
  }
  if(Marginal_Dist1 == "GamMix(2)"){
    if(is.na(Marginal_Dist1_Par)==T){
      data.gamlss<-data.frame(X=Data_Con1[,con2])
      marginal_non_con1 <- tryCatch(gamlssMX(X~1, data=data.gamlss, family=GA, K=2),
                                    error = function(e) "error")
    }else{
      marginal_non_con1<-Marginal_Dist1_Par
    }
  }
  if(Marginal_Dist1 == "GamMix(3)"){
    if(is.na(Marginal_Dist1_Par)==T){
      data.gamlss<-data.frame(X=Data_Con1[,con2])
      marginal_non_con1 <- tryCatch(gamlssMX(X~1, data=data.gamlss, family=GA, K=3),
                                    error = function(e) "error")
    }else{
      marginal_non_con1<-Marginal_Dist1_Par
    }
  }
  if(Marginal_Dist1 == "Gaus"){
    if(is.na(Marginal_Dist1_Par)==T){
      marginal_non_con1<-fitdistr(Data_Con1[,con2],"normal")
    }else{
      marginal_non_con1<-Marginal_Dist1_Par
    }
  }
  if(Marginal_Dist1 == "Gum"){
    if(is.na(Marginal_Dist1_Par)==T){
      marginal_non_con1 <- gamlss(Data_Con1[,con2]  ~ 1, family= GU)
    }else{
      marginal_non_con1<-Marginal_Dist1_Par
    }
  }
  if(Marginal_Dist1 == "InvG"){
    if(is.na(Marginal_Dist1_Par)==T){
      marginal_non_con1<-fitdist(Data_Con1[,con2], "invgauss", start = list(mean = 5, shape = 1))
    }else{
      marginal_non_con1<-Marginal_Dist1_Par
    }
  }
  if(Marginal_Dist1 == "Lapl"){
    if(is.na(Marginal_Dist1_Par)==T){
      marginal_non_con1<-fitdistr(Data_Con1[,con2], dlaplace, start=list(location=mean(Data_Con1[,con2]), scale=sd(Data_Con1[,con2])/sqrt(2)))
    }else{
    marginal_non_con1<-Marginal_Dist1_Par
    }
  }
  if(Marginal_Dist1 == "Logis"){
    if(is.na(Marginal_Dist1_Par)==T){
      marginal_non_con1<-fitdistr(Data_Con1[,con2], "logistic")
    }else{
      marginal_non_con1<-Marginal_Dist1_Par
    }
  }
  if(Marginal_Dist1 == "LogN"){
    if(is.na(Marginal_Dist1_Par)==T){
      marginal_non_con1<-fitdistr(Data_Con1[,con2],"lognormal")
    }else{
      marginal_non_con1<-Marginal_Dist1_Par
    }
  }
  if(Marginal_Dist1 == "RGum"){
    if(is.na(Marginal_Dist1_Par)==T){
      marginal_non_con1 <- gamlss(Data_Con1[,con2] ~ 1,family=RG)
    }else{
      marginal_non_con1<-Marginal_Dist1_Par
    }
  }
  if(Marginal_Dist1 == "TNorm"){
    if(is.na(Marginal_Dist1_Par)==T){
      marginal_non_con1<-fitdistr(Data_Con1[,con2],"normal")
    }else{
      marginal_non_con1<-Marginal_Dist1_Par
    }
  }
  if(Marginal_Dist1 == "Twe"){
    if(is.na(Marginal_Dist1_Par)==T){
      marginal_non_con1<-tweedie.profile(Data_Con1[,con2] ~ 1,p.vec=seq(1.5, 2.5, by=0.2), do.plot=FALSE)
    }else{
      marginal_non_con1<-Marginal_Dist1_Par
    }
  }
  if(Marginal_Dist1 == "Weib"){
    if(is.na(Marginal_Dist1_Par)==T){
      marginal_non_con1<-fitdistr(Data_Con1[,con2], "weibull")
    }else{
      marginal_non_con1<-Marginal_Dist1_Par
    }
  }

  #Fit the GPD to the conditioned variable con2 in Data_Con2.
  if(is.na(GPD2[[1]][1]) & is.na(Thres2) & is.na(Tab2[[1]][1])){
    Thres2<-quantile(na.omit(Data[,con2]),u2)
  }

  if(is.na(GPD2[[1]][1]) & GPD_Bayes & is.na(Tab2[[1]][1])){
    GPD_con2<-evm(Data_Con2[,con2], th=Thres2 ,penalty = "gaussian",priorParameters = list(c(0, 0), matrix(c(100^2, 0, 0, 0.25), nrow = 2)))
  }
  if(is.na(GPD2[[1]][1]) & !GPD_Bayes & is.na(Tab2[[1]][1])){
    GPD_con2<-evm(Data_Con2[,con2], th= Thres2)
  }

  ##Fit the specified marginal distribution (Marginal_Dist2) to the non-conditioned variable con1 in Data_Con2.
  if(Marginal_Dist2 == "BS"){
    if(is.na(Marginal_Dist2_Par)==T){
      bdata2 <- data.frame(shape = exp(-0.5), scale = exp(0.5))
      bdata2 <- transform(bdata2, y = Data_Con2[,con1])
      marginal_non_con2<-vglm(y ~ 1, bisa, data = bdata2, trace = FALSE)
    }else{
      marginal_non_con2<-Marginal_Dist2_Par
    }
  }
  if(Marginal_Dist2 == "Exp"){
    if(is.na(Marginal_Dist2_Par)==T){
      marginal_non_con2<-fitdistr(Data_Con2[,con1],"exponential")
    }else{
      marginal_non_con2<-Marginal_Dist2_Par
    }
  }
  if(Marginal_Dist2 == "Gam(2)"){
    if(is.na(Marginal_Dist2_Par)==T){
      marginal_non_con2<-fitdistr(Data_Con2[,con1], "gamma")
    }else{
      marginal_non_con2<-Marginal_Dist2_Par
    }
  }
  if(Marginal_Dist2 == "Gam(3)"){
    if(is.na(Marginal_Dist2_Par)==T){
      data.gamlss<-data.frame(X=Data_Con2[,con1])
      marginal_non_con2 <- tryCatch(gamlss(X~1, data=data.gamlss, family=GG),
                                    error = function(e) "error")
    }else{
      marginal_non_con2<-Marginal_Dist2_Par
    }
  }
  if(Marginal_Dist2 == "GamMix(2)"){
    if(is.na(Marginal_Dist2_Par)==T){
      data.gamlss<-data.frame(X=Data_Con2[,con1])
      marginal_non_con2 <- tryCatch(gamlssMX(X~1, data=data.gamlss, family=GA, K=2),
                                    error = function(e) "error")
    }else{
      marginal_non_con2<-Marginal_Dist2_Par
    }
  }
  if(Marginal_Dist2 == "GamMix(3)"){
    if(is.na(Marginal_Dist2_Par)==T){
      data.gamlss<-data.frame(X=Data_Con2[,con1])
      marginal_non_con2 <- tryCatch(gamlssMX(X~1, data=data.gamlss, family=GA, K=3),
                                    error = function(e) "error")
    }else{
      marginal_non_con2<-Marginal_Dist2_Par
    }
  }
  if(Marginal_Dist2 == "Gaus"){
    if(is.na(Marginal_Dist2_Par)==T){
      marginal_non_con2<-fitdistr(Data_Con2[,con1],"normal")
    }else{
      marginal_non_con2<-Marginal_Dist2_Par
    }
  }
  if(Marginal_Dist2 == "Gum"){
    if(is.na(Marginal_Dist2_Par)==T){
      marginal_non_con2 <- gamlss(Data_Con2[,con1]  ~ 1, family= GU)
    }else{
      marginal_non_con2<-Marginal_Dist2_Par
    }
  }
  if(Marginal_Dist2 == "InvG"){
    if(is.na(Marginal_Dist2_Par)==T){
      marginal_non_con2<-fitdist(Data_Con2[,con1], "invgauss", start = list(mean = 5, shape = 1))
    }else{
      marginal_non_con2<-Marginal_Dist2_Par
    }
  }
  if(Marginal_Dist2 == "Lapl"){
    if(is.na(Marginal_Dist2_Par)==T){
      marginal_non_con2<-fitdistr(Data_Con2[,con1], dlaplace, start=list(location=mean(Data_Con2[,con1]), scale=sd(Data_Con2[,con1])/sqrt(2)))
    }else{
      marginal_non_con2<-Marginal_Dist2_Par
    }
  }
  if(Marginal_Dist2 == "Logis"){
    if(is.na(Marginal_Dist2_Par)==T){
      marginal_non_con2<-fitdistr(Data_Con2[,con1],"logistic")
    }else{
      marginal_non_con2<-Marginal_Dist2_Par
    }
  }
  if(Marginal_Dist2 == "LogN"){
    if(is.na(Marginal_Dist2_Par)==T){
      marginal_non_con2<-fitdistr(Data_Con2[,con1],"lognormal")
    }else{
      marginal_non_con2<-Marginal_Dist2_Par
    }
  }
  if(Marginal_Dist2 == "RGum"){
    if(is.na(Marginal_Dist2_Par)==T){
      marginal_non_con2 <- gamlss(Data_Con2[,con1] ~ 1,family=RG)
    }else{
      marginal_non_con2<-Marginal_Dist2_Par
    }
  }
  if(Marginal_Dist2 == "TNorm"){
    if(is.na(Marginal_Dist2_Par)==T){
      marginal_non_con2<-fitdistr(Data_Con2[,con1],"normal")
    }else{
      marginal_non_con2<-Marginal_Dist2_Par
    }
  }
  if(Marginal_Dist2 == "Twe"){
    if(is.na(Marginal_Dist2_Par)==T){
      marginal_non_con2<-tweedie.profile(Data_Con2[,con1] ~ 1,p.vec=seq(1.5, 2.5, by=0.2), do.plot=FALSE)
    }else{
      marginal_non_con2<-Marginal_Dist2_Par
    }
  }
  if(Marginal_Dist2 == "Weib"){
    if(is.na(Marginal_Dist2_Par)==T){
      marginal_non_con2<-fitdistr(Data_Con2[,con1], "weibull")
    }else{
      marginal_non_con2<-Marginal_Dist2_Par
    }
  }

  ###Generate samples from the copula models to which a kernel density estimate will be applied to estimate relative probabilities along the isoline.

  #Fit the specified copula family (Copula_Family1) to the observations in Data_Con1.
  obj1<-BiCopSelect(pobs(Data_Con1[,1]), pobs(Data_Con1[,2]), familyset=Copula_Family1, selectioncrit = "AIC",
                    indeptest = FALSE, level = 0.05, weights = NA, rotations = TRUE,
                    se = FALSE, presel = TRUE, method = "mle")
  #Simulate a sample from the fitted copula. Out of the sample size 'N' the proportion of the sample from the copula associated with Data_Con1 is proportional to the size of Data_Con1 relative to Data_Con2.
  sample<-BiCopSim(round(N*Rate_Con1/(Rate_Con1+Rate_Con2+Rate_Con3+Rate_Con4),0),obj1)

  #Transform the realizations of the conditioned variable con1 to the original scale using inverse cumulative distribution a.k.a. quantile functions (inverse probability integral transform) of the GPD contained in the u2gpd function.
  if(is.na(GPD1[[1]][1])==T & is.na(Tab1[[1]][1])==T){
    cop.sample1.con<-u2gpd(sample[,con1], p = 1, th=Thres1 , sigma=exp(GPD_con1$coefficients[1]),xi= GPD_con1$coefficients[2])
  }
  if(is.na(GPD1[[1]][1])==F){
    cop.sample1.con<-u2gpd(sample[,con1], p = 1, th = GPD1$Threshold, sigma = GPD1$sigma, xi= GPD1$xi)
  }
  if(is.na(Tab1[[1]][1])==F){
    cop.sample1.con = approx(1-1/Tab1[,1],Tab1[,2],xout=u1+sample[,con1]*(1-u1))$y
  }

  #Transform the realizations of the non-conditioned variable con2 to the original scale using the quantile function of the selected parametric (non-extreme value) distribution (Marginal_Dist1).
  if(Marginal_Dist1=="BS"){
    cop.sample1.non.con<-qbisa(sample[,con2], as.numeric(Coef(marginal_non_con1)[1]), as.numeric(Coef(marginal_non_con1)[2]))
  }
  if(Marginal_Dist1=="Exp"){
    cop.sample1.non.con<-qexp(sample[,con2], rate = as.numeric(marginal_non_con1$estimate[1]))
  }
  if(Marginal_Dist1=="Gam(2)"){
    cop.sample1.non.con<-qgamma(sample[,con2], shape = as.numeric(marginal_non_con1$estimate[1]), rate = as.numeric(marginal_non_con1$estimate[2]))
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
  if(Marginal_Dist1=="Gaus"){
    cop.sample1.non.con<-qnorm(sample[,con2], mean = as.numeric(marginal_non_con1$estimate[1]), sd = as.numeric(marginal_non_con1$estimate[2]))
  }
  if(Marginal_Dist1=="Gum"){
    cop.sample1.non.con<-qGU(sample[,con2],as.numeric(marginal_non_con1$mu.coefficients),exp(as.numeric(marginal_non_con1$sigma.coefficients)))
  }
  if(Marginal_Dist1=="InvG"){
    cop.sample1.non.con<-qinvgauss(sample[,con2], mean = as.numeric(marginal_non_con1$estimate[1]), shape = as.numeric(marginal_non_con1$estimate[2]))
  }
  if(Marginal_Dist1=="Lapl"){
    cop.sample1.non.con <- qlaplace(sample[,con2],as.numeric(marginal_non_con1$estimate[1]), as.numeric(marginal_non_con1$estimate[2]))
  }
  if(Marginal_Dist1=="Logis"){
    cop.sample1.non.con<-qlogis(sample[,con2], location = as.numeric(marginal_non_con1$estimate[1]), scale = as.numeric(marginal_non_con1$estimate[2]))
  }
  if(Marginal_Dist1=="LogN"){
    cop.sample1.non.con<-qlnorm(sample[,con2], meanlog = as.numeric(marginal_non_con1$estimate[1]), sdlog = as.numeric(marginal_non_con1$estimate[2]))
  }
  if(Marginal_Dist1=="RGum"){
    cop.sample1.non.con<-qRG(sample[,con2],marginal_non_con1$mu.coefficients,exp(marginal_non_con1$sigma.coefficients))
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
  sample<-BiCopSim(round(N*Rate_Con2/(Rate_Con1+Rate_Con2+Rate_Con3+Rate_Con4),0),obj2)

  #Transform the realizations of the conditioned variable con2 to the original scale using the inverse CDF (quantile function) of the GPD contained in the u2gpd function.
  if(is.na(GPD2[[1]][1]) & is.na(Tab2[[1]][1])){
    cop.sample2.con<-u2gpd(sample[,con2], p = 1, th=Thres2, sigma=exp(GPD_con2$coefficients[1]),xi= GPD_con2$coefficients[2])
  }
  if(!is.na(GPD2[[1]][1])){
    cop.sample2.con<-u2gpd(sample[,con2], p = 1, th = GPD2$Threshold, sigma = GPD2$sigma, xi= GPD2$xi)
  }
  if(!is.na(Tab2[[1]][1])){
    cop.sample2.con = approx(1-1/Tab2[,1],Tab2[,2],xout=u2+sample[,con2]*(1-u2))$y
  }

  #Transform the realizations of the non-conditioned variable con1 to the original scale using the inverse CDF (quantile function) of the selected parametric (non-extreme value) distribution (Marginal_Dist2).
  if(Marginal_Dist2=="BS"){
    cop.sample2.non.con<-qbisa(sample[,con1], as.numeric(Coef(marginal_non_con2)[1]), as.numeric(Coef(marginal_non_con2)[2]))
  }
  if(Marginal_Dist2=="Exp"){
    cop.sample2.non.con<-qexp(sample[,con1], rate = as.numeric(marginal_non_con2$estimate[1]))
  }
  if(Marginal_Dist2=="Gam(2)"){
    cop.sample2.non.con<-qgamma(sample[,con1], shape = as.numeric(marginal_non_con2$estimate[1]), rate=as.numeric(marginal_non_con2$estimate[2]))
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
  if(Marginal_Dist2=="Gaus"){
    cop.sample2.non.con<-qnorm(sample[,con1], mean = as.numeric(marginal_non_con2$estimate[1]), sd=as.numeric(marginal_non_con2$estimate[2]))
  }
  if(Marginal_Dist2=="Gum"){
    cop.sample2.non.con<-qGU(sample[,con1],as.numeric(marginal_non_con2$mu.coefficients),exp(as.numeric(marginal_non_con2$sigma.coefficients)))
  }
  if(Marginal_Dist2=="InvG"){
    cop.sample2.non.con<-qinvgauss(sample[,con1], mean = as.numeric(marginal_non_con2$estimate[1]), shape=as.numeric(marginal_non_con2$estimate[2]))
  }
  if(Marginal_Dist2=="Lapl"){
    cop.sample2.non.con <- qlaplace(sample[,con1],as.numeric(marginal_non_con2$estimate[1]), as.numeric(marginal_non_con2$estimate[2]))
  }
  if(Marginal_Dist2=="LogN"){
    cop.sample2.non.con<-qlnorm(sample[,con1], meanlog = as.numeric(marginal_non_con2$estimate[1]), sdlog = as.numeric(marginal_non_con2$estimate[2]))
  }
  if(Marginal_Dist2=="Logis"){
    cop.sample2.non.con<-qlogis(sample[,con1], location = as.numeric(marginal_non_con2$estimate[1]), scale=as.numeric(marginal_non_con2$estimate[2]))
  }
  if(Marginal_Dist2=="RGum"){
    cop.sample2.non.con<-qRG(sample[,con1],marginal_non_con2$mu.coefficients,exp(marginal_non_con2$sigma.coefficients))
  }
  if(Marginal_Dist2=="TNorm"){
    cop.sample2.non.con<-qtruncnorm(sample[,con1], a=min(Data_Con2[,con1]), mean = as.numeric(marginal_non_con2$estimate[1]), sd = as.numeric(marginal_non_con2$estimate[2]))
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
  cop.sample.pop.1<-rbind(cop.sample1,cop.sample2)

  ###Samples from second population

  #Fit the GPD to the conditioned variable con1 in Data_Con1.
  if(is.na(GPD3[[1]][1]) & is.na(Thres3) & is.na(Tab3[[1]][1])){
    Thres3<-quantile(na.omit(Data[,con3]),u3)
  }

  if(is.na(GPD3[[1]][1]) & GPD_Bayes & is.na(Tab3[[1]][1])){
    GPD_con3<-evm(Data_Con3[,con3], th = Thres3,penalty = "gaussian",priorParameters = list(c(0, 0), matrix(c(100^2, 0, 0, 0.25), nrow = 2)))
  }
  if(is.na(GPD3[[1]][1]) & !GPD_Bayes & is.na(Tab3[[1]][1])){
    GPD_con3<-evm(Data_Con3[,con3], th = Thres3)
  }

  #Fit the specified marginal distribution (Marginal_Dist3) to the non-conditioned variable con4 in Data_Con3.
  if(Marginal_Dist3 == "BS"){
    if(is.na(Marginal_Dist3_Par)==T){
      bdata2 <- data.frame(shape = exp(-0.5), scale = exp(0.5))
      bdata2 <- transform(bdata2, y = Data_Con3[,con4])
      marginal_non_con3<-vglm(y ~ 1, bisa, data = bdata2, trace = FALSE)
    }else{
      marginal_non_con3<-Marginal_Dist3_Par
    }
  }
  if(Marginal_Dist3 == "Exp"){
    if(is.na(Marginal_Dist3_Par)==T){
      marginal_non_con3<-fitdistr(Data_Con3[,con4],"exponential")
    }else{
      marginal_non_con3<-Marginal_Dist3_Par
    }
  }
  if(Marginal_Dist3 == "Gam(2)"){
    if(is.na(Marginal_Dist3_Par)==T){
      marginal_non_con3<-fitdistr(Data_Con3[,con4], "gamma")
    }else{
      marginal_non_con3<-Marginal_Dist3_Par
    }
  }
  if(Marginal_Dist3 == "Gam(3)"){
    if(is.na(Marginal_Dist3_Par)==T){
      data.gamlss<-data.frame(X=Data_Con3[,con4])
      marginal_non_con3 <- tryCatch(gamlss(X~1, data=data.gamlss, family=GG),
                                    error = function(e) "error")
    }else{
      marginal_non_con3<-Marginal_Dist3_Par
    }
  }
  if(Marginal_Dist3 == "GamMix(2)"){
    if(is.na(Marginal_Dist3_Par)==T){
      data.gamlss<-data.frame(X=Data_Con3[,con4])
      marginal_non_con3 <- tryCatch(gamlssMX(X~1, data=data.gamlss, family=GA, K=2),
                                    error = function(e) "error")
    }else{
      marginal_non_con3<-Marginal_Dist3_Par
    }
  }
  if(Marginal_Dist3 == "GamMix(3)"){
    if(is.na(Marginal_Dist3_Par)==T){
      data.gamlss<-data.frame(X=Data_Con3[,con4])
      marginal_non_con3 <- tryCatch(gamlssMX(X~1, data=data.gamlss, family=GA, K=3),
                                    error = function(e) "error")
    }else{
      marginal_non_con3<-Marginal_Dist3_Par
    }
  }
  if(Marginal_Dist3 == "Gaus"){
    if(is.na(Marginal_Dist3_Par)==T){
      marginal_non_con3<-fitdistr(Data_Con3[,con4],"normal")
    }else{
      marginal_non_con3<-Marginal_Dist3_Par
    }
  }
  if(Marginal_Dist3 == "Gum"){
    if(is.na(Marginal_Dist3_Par)==T){
      marginal_non_con3 <- gamlss(Data_Con3[,con4]  ~ 1, family= GU)
    }else{
      marginal_non_con3<-Marginal_Dist3_Par
    }
  }
  if(Marginal_Dist3 == "InvG"){
    if(is.na(Marginal_Dist3_Par)==T){
      marginal_non_con3<-fitdist(Data_Con3[,con4], "invgauss", start = list(mean = 5, shape = 1))
    }else{
      marginal_non_con3<-Marginal_Dist3_Par
    }
  }
  if(Marginal_Dist3 == "Lapl"){
    if(is.na(Marginal_Dist3_Par)==T){
      marginal_non_con3<-fitdistr(Data_Con3[,con4], dlaplace, start=list(location=mean(Data_Con3[,con4]), scale=sd(Data_Con3[,con4])/sqrt(2)))
    }else{
      marginal_non_con3<-Marginal_Dist3_Par
    }
  }
  if(Marginal_Dist3 == "Logis"){
    if(is.na(Marginal_Dist3_Par)==T){
      marginal_non_con3<-fitdistr(Data_Con3[,con4], "logistic")
    }else{
      marginal_non_con3<-Marginal_Dist3_Par
    }
  }
  if(Marginal_Dist3 == "LogN"){
    if(is.na(Marginal_Dist3_Par)==T){
      marginal_non_con3<-fitdistr(Data_Con3[,con4],"lognormal")
    }else{
      marginal_non_con3<-Marginal_Dist3_Par
    }
  }
  if(Marginal_Dist3 == "RGum"){
    if(is.na(Marginal_Dist3_Par)==T){
      marginal_non_con3 <- gamlss(Data_Con3[,con4] ~ 1,family=RG)
    }else{
      marginal_non_con3<-Marginal_Dist3_Par
    }
  }
  if(Marginal_Dist3 == "TNorm"){
    if(is.na(Marginal_Dist3_Par)==T){
      marginal_non_con3<-fitdistr(Data_Con3[,con4],"normal")
    }else{
      marginal_non_con3<-Marginal_Dist3_Par
    }
  }
  if(Marginal_Dist3 == "Twe"){
    if(is.na(Marginal_Dist3_Par)==T){
      marginal_non_con3<-tweedie.profile(Data_Con3[,con4] ~ 1,p.vec=seq(1.5, 2.5, by=0.2), do.plot=FALSE)
    }else{
      marginal_non_con3<-Marginal_Dist3_Par
    }
  }
  if(Marginal_Dist3 == "Weib"){
    if(is.na(Marginal_Dist3_Par)==T){
      marginal_non_con3<fitdistr(Data_Con3[,con4], "weibull")
    }else{
      marginal_non_con3<-Marginal_Dist3_Par
    }
  }

  #Fit the GPD to the conditioned variable con4 in Data_Con4.
  if(is.na(GPD4[[1]][1]) & is.na(Thres4) & is.na(Tab4[[1]][1])){
    Thres4<-quantile(na.omit(Data_4[,con4]),u4)
  }

  if(is.na(GPD4[[1]][1]) & GPD_Bayes & is.na(Tab4[[1]][1])){
    GPD_con4<-evm(Data_Con4[,con4], th=Thres4 ,penalty = "gaussian",priorParameters = list(c(0, 0), matrix(c(100^2, 0, 0, 0.25), nrow = 2)))
  }
  if(is.na(GPD4[[1]][1]) & !GPD_Bayes & is.na(Tab4[[1]][1])){
    GPD_con4<-evm(Data_Con4[,con4], th= Thres4)
  }

  ##Fit the specified marginal distribution (Marginal_Dist4) to the non-conditioned variable con3 in Data_Con4.
  if(Marginal_Dist4 == "BS"){
    if(is.na(Marginal_Dist4_Par)==T){
      bdata2 <- data.frame(shape = exp(-0.5), scale = exp(0.5))
      bdata2 <- transform(bdata2, y = Data_Con4[,con3])
      marginal_non_con4<-vglm(y ~ 1, bisa, data = bdata2, trace = FALSE)
    }else{
      marginal_non_con4<-Marginal_Dist4_Par
    }
  }
  if(Marginal_Dist4 == "Exp"){
    if(is.na(Marginal_Dist4_Par)==T){
      marginal_non_con4<-fitdistr(Data_Con4[,con3],"exponential")
    }else{
      marginal_non_con4<-Marginal_Dist4_Par
    }
  }
  if(Marginal_Dist4 == "Gam(2)"){
    if(is.na(Marginal_Dist4_Par)==T){
      marginal_non_con4<-fitdistr(Data_Con4[,con3], "gamma")
    }else{
      marginal_non_con4<-Marginal_Dist4_Par
    }
  }
  if(Marginal_Dist4 == "Gam(3)"){
    if(is.na(Marginal_Dist4_Par)==T){
      data.gamlss<-data.frame(X=Data_Con4[,con3])
      marginal_non_con4 <- tryCatch(gamlss(X~1, data=data.gamlss, family=GG),
                                    error = function(e) "error")
    }else{
      marginal_non_con4<-Marginal_Dist4_Par
    }
  }
  if(Marginal_Dist4 == "GamMix(2)"){
    if(is.na(Marginal_Dist4_Par)==T){
      data.gamlss<-data.frame(X=Data_Con4[,con3])
      marginal_non_con4 <- tryCatch(gamlssMX(X~1, data=data.gamlss, family=GA, K=2),
                                    error = function(e) "error")
    }else{
      marginal_non_con4<-Marginal_Dist4_Par
    }
  }
  if(Marginal_Dist4 == "GamMix(3)"){
    if(is.na(Marginal_Dist4_Par)==T){
      data.gamlss<-data.frame(X=Data_Con4[,con3])
      marginal_non_con4 <- tryCatch(gamlssMX(X~1, data=data.gamlss, family=GA, K=3),
                                    error = function(e) "error")
    }else{
      marginal_non_con4<-Marginal_Dist4_Par
    }
  }
  if(Marginal_Dist4 == "Gaus"){
    if(is.na(Marginal_Dist4_Par)==T){
      marginal_non_con4<-fitdistr(Data_Con4[,con3],"normal")
    }else{
      marginal_non_con4<-Marginal_Dist4_Par
    }
  }
  if(Marginal_Dist4 == "Gum"){
    if(is.na(Marginal_Dist4_Par)==T){
      marginal_non_con4 <- gamlss(Data_Con4[,con3]  ~ 1, family= GU)
    }else{
      marginal_non_con4<-Marginal_Dist4_Par
    }
  }
  if(Marginal_Dist4 == "InvG"){
    if(is.na(Marginal_Dist4_Par)==T){
      marginal_non_con4<-fitdist(Data_Con4[,con3], "invgauss", start = list(mean = 5, shape = 1))
    }else{
      marginal_non_con4<-Marginal_Dist4_Par
    }
  }
  if(Marginal_Dist4 == "Lapl"){
    if(is.na(Marginal_Dist4_Par)==T){
      marginal_non_con4<-fitdistr(Data_Con4[,con3], dlaplace, start=list(location=mean(Data_Con4[,con3]), scale=sd(Data_Con4[,con3])/sqrt(2)))
    }else{
      marginal_non_con4<-Marginal_Dist4_Par
    }
  }
  if(Marginal_Dist4 == "Logis"){
    if(is.na(Marginal_Dist4_Par)==T){
      marginal_non_con4<-fitdistr(Data_Con4[,con3],"logistic")
    }else{
      marginal_non_con4<-Marginal_Dist4_Par
    }
  }
  if(Marginal_Dist4 == "LogN"){
    if(is.na(Marginal_Dist4_Par)==T){
      marginal_non_con4<-fitdistr(Data_Con4[,con3],"lognormal")
    }else{
      marginal_non_con4<-Marginal_Dist4_Par
    }
  }
  if(Marginal_Dist4 == "RGum"){
    if(is.na(Marginal_Dist4_Par)==T){
      marginal_non_con4 <- gamlss(Data_Con4[,con3] ~ 1,family=RG)
    }else{
      marginal_non_con4<-Marginal_Dist4_Par
    }
  }
  if(Marginal_Dist4 == "TNorm"){
    if(is.na(Marginal_Dist4_Par)==T){
      marginal_non_con4<-fitdistr(Data_Con4[,con3],"normal")
    }else{
      marginal_non_con4<-Marginal_Dist4_Par
    }
  }
  if(Marginal_Dist4 == "Twe"){
    if(is.na(Marginal_Dist4_Par)==T){
      marginal_non_con4<-tweedie.profile(Data_Con4[,con3] ~ 1,p.vec=seq(1.5, 2.5, by=0.2), do.plot=FALSE)
    }else{
      marginal_non_con4<-Marginal_Dist4_Par
    }
  }
  if(Marginal_Dist4 == "Weib"){
    if(is.na(Marginal_Dist4_Par)==T){
      marginal_non_con4<-fitdistr(Data_Con4[,con3], "weibull")
    }else{
      marginal_non_con4<-Marginal_Dist4_Par
    }
  }

  ###Generate samples from the copula models to which a kernel density estimate will be applied to estimate relative probabilities along the isoline.

  #Fit the specified copula family (Copula_Family1) to the observations in Data_Con1.
  obj3<-BiCopSelect(pobs(Data_Con3[,1]), pobs(Data_Con3[,2]), familyset=Copula_Family1, selectioncrit = "AIC",
                    indeptest = FALSE, level = 0.05, weights = NA, rotations = TRUE,
                    se = FALSE, presel = TRUE, method = "mle")
  #Simulate a sample from the fitted copula. Out of the sample size 'N' the proportion of the sample from the copula associated with Data_Con1 is proportional to the size of Data_Con1 relative to Data_Con2.
  sample<-BiCopSim(round(N*Rate_Con3/(Rate_Con1+Rate_Con2+Rate_Con3+Rate_Con4),0),obj3)

  #Transform the realizations of the conditioned variable con1 to the original scale using inverse cumulative distribution a.k.a. quantile functions (inverse probability integral transform) of the GPD contained in the u2gpd function.
  if(is.na(GPD3[[1]][1]) & is.na(Tab3[[1]][1])){
    cop.sample3.con<-u2gpd(sample[,con3], p = 1, th=Thres3 , sigma=exp(GPD_con3$coefficients[1]),xi= GPD_con3$coefficients[2])
  }
  if(!is.na(GPD3[[1]][1])){
    cop.sample3.con<-u2gpd(sample[,con3], p = 1, th = GPD3$Threshold, sigma = GPD3$sigma, xi= GPD3$xi)
  }
  if(!is.na(Tab3[[1]][1])){
    cop.sample3.con = approx(1-1/Tab3[,1],Tab3[,2],xout=u3+sample[,con3]*(1-u3))$y
  }

  #Transform the realizations of the non-conditioned variable con2 to the original scale using the quantile function of the selected parametric (non-extreme value) distribution (Marginal_Dist1).
  if(Marginal_Dist3=="BS"){
    cop.sample3.non.con<-qbisa(sample[,con4], as.numeric(Coef(marginal_non_con3)[1]), as.numeric(Coef(marginal_non_con3)[2]))
  }
  if(Marginal_Dist3=="Exp"){
    cop.sample3.non.con<-qexp(sample[,con4], rate = as.numeric(marginal_non_con3$estimate[1]))
  }
  if(Marginal_Dist3=="Gam(2)"){
    cop.sample3.non.con<-qgamma(sample[,con4], shape = as.numeric(marginal_non_con3$estimate[1]), rate = as.numeric(marginal_non_con3$estimate[2]))
  }
  if(Marginal_Dist3=="Gam(3)"){
    cop.sample3.non.con<-qGG(sample[,con4], mu=exp(marginal_non_con3$mu.coefficients), sigma=exp(marginal_non_con3$sigma.coefficients), nu=marginal_non_con3$nu.coefficients)
  }
  if(Marginal_Dist3=="GamMix(2)"){
    xx <- seq(0, max(Data_Con3[,2])*10, 0.001)
    prob.MX1 <- round(marginal_non_con3$prob[1],3)
    prob.MX2 <- 1 - prob.MX1
    cdf.MX<-pMX(xx, mu=list(mu1=exp(marginal_non_con3$models[[1]]$mu.coefficients), mu2=exp(marginal_non_con3$models[[2]]$mu.coefficients)),
                sigma=list(sigma1=exp(marginal_non_con3$models[[1]]$sigma.coefficients), sigma2=exp(marginal_non_con3$models[[2]]$sigma.coefficients)),
                pi = list(pi1=prob.MX1, pi2=prob.MX2), family=list(fam1="GA", fam2="GA"))
    cop.sample3.non.con <- approx(cdf.MX, xx, sample[,con4])$y
  }
  if(Marginal_Dist1=="GamMix(3)"){
    xx <- seq(0, max(Data_Con2[,2])*10, 0.001)
    prob.MX1 <- round(marginal_non_con3$prob[1],3)
    prob.MX2 <- round(marginal_non_con3$prob[2],3)
    prob.MX3 <- 1 - prob.MX1 - prob.MX2
    cdf.MX<-pMX(xx, mu=list(mu1=exp(marginal_non_con3$models[[1]]$mu.coefficients), mu2=exp(marginal_non_con3$models[[2]]$mu.coefficients), mu3=exp(marginal_non_con3$models[[3]]$mu.coefficients)),
                sigma=list(sigma1=exp(marginal_non_con3$models[[1]]$sigma.coefficients), sigma2=exp(marginal_non_con3$models[[2]]$sigma.coefficients), sigma3=exp(marginal_non_con3$models[[3]]$sigma.coefficients)),
                pi = list(pi1=prob.MX1, pi2=prob.MX2, pi3=prob.MX3), family=list(fam1="GA", fam2="GA", fam3="GA"))
    cop.sample3.non.con <- approx(cdf.MX, xx, sample[,con4])$y
  }
  if(Marginal_Dist3=="Gaus"){
    cop.sample3.non.con<-qnorm(sample[,con4], mean = as.numeric(marginal_non_con3$estimate[1]), sd = as.numeric(marginal_non_con3$estimate[2]))
  }
  if(Marginal_Dist3=="Gum"){
    cop.sample3.non.con<-qGU(sample[,con4],as.numeric(marginal_non_con3$mu.coefficients),exp(as.numeric(marginal_non_con3$sigma.coefficients)))
  }
  if(Marginal_Dist3=="InvG"){
    cop.sample3.non.con<-qinvgauss(sample[,con4], mean = as.numeric(marginal_non_con3$estimate[1]), shape = as.numeric(marginal_non_con3$estimate[2]))
  }
  if(Marginal_Dist3=="Lapl"){
    cop.sample3.non.con <- qlaplace(sample[,con4],as.numeric(marginal_non_con3$estimate[1]), as.numeric(marginal_non_con3$estimate[2]))
  }
  if(Marginal_Dist3=="Logis"){
    cop.sample3.non.con<-qlogis(sample[,con4], location = as.numeric(marginal_non_con3$estimate[1]), scale = as.numeric(marginal_non_con3$estimate[2]))
  }
  if(Marginal_Dist3=="LogN"){
    cop.sample3.non.con<-qlnorm(sample[,con4], meanlog = as.numeric(marginal_non_con3$estimate[1]), sdlog = as.numeric(marginal_non_con3$estimate[2]))
  }
  if(Marginal_Dist3=="RGum"){
    cop.sample3.non.con<-qRG(sample[,con4],marginal_non_con3$mu.coefficients,exp(marginal_non_con3$sigma.coefficients))
  }
  if(Marginal_Dist3=="TNorm"){
    cop.sample3.non.con<-qtruncnorm(sample[,con4], a=min(Data_Con3[,con4]), mean = as.numeric(marginal_non_con3$estimate[1]), sd = as.numeric(marginal_non_con3$estimate[2]))
  }
  if(Marginal_Dist3=="Twe"){
    cop.sample3.non.con<-qtweedie(sample[,con4], power=marginal_non_con3$p.max, mu=mean(Data_Con3[,con4]), phi=marginal_non_con3$phi.max)
  }
  if(Marginal_Dist3=="Weib"){
    cop.sample3.non.con<-qweibull(sample[,con4], shape = as.numeric(marginal_non_con3$estimate[1]), scale=as.numeric(marginal_non_con3$estimate[2]))
  }
  #Put the realizations that have been transformed to the original scale in a data frame
  cop.sample3<-data.frame(cop.sample3.con,cop.sample3.non.con)
  colnames(cop.sample3)<-c("Var1","Var2")

  #Fit the specified copula family (Copula_Family2) to the observations in Data_Con2.
  obj4<-BiCopSelect(pobs(Data_Con4[,1]), pobs(Data_Con4[,2]), familyset=Copula_Family4, selectioncrit = "AIC",
                    indeptest = FALSE, level = 0.05, weights = NA, rotations = TRUE,
                    se = FALSE, presel = TRUE, method = "mle")
  #Simulate a sample from the fitted copula. Out of the sample size 'N' the proportion of the sample from the copula assoicated with Data_Con2 is proportional to the size of Data_Con2 relative to Data_Con1.
  sample<-BiCopSim(round(N*Rate_Con4/(Rate_Con1+Rate_Con2+Rate_Con3+Rate_Con4),0),obj4)

  #Transform the realizations of the conditioned variable con2 to the original scale using the inverse CDF (quantile function) of the GPD contained in the u2gpd function.
  if(is.na(GPD4[[1]][1]) & is.na(Tab4[[1]][1])){
    cop.sample4.con<-u2gpd(sample[,con4], p = 1, th=Thres4, sigma=exp(GPD_con4$coefficients[1]),xi= GPD_con4$coefficients[2])
  }
  if(!is.na(GPD4[[1]][1])){
    cop.sample4.con<-u2gpd(sample[,con4], p = 1, th = GPD4$Threshold, sigma = GPD4$sigma, xi= GPD4$xi)
  }
  if(!is.na(Tab4[[1]][1])){
    cop.sample4.con = approx(1-1/Tab4[,1],Tab4[,2],xout=u4+sample[,con4]*(1-u4))$y
  }

  #Transform the realizations of the non-conditioned variable con1 to the original scale using the inverse CDF (quantile function) of the selected parametric (non-extreme value) distribution (Marginal_Dist2).
  if(Marginal_Dist4=="BS"){
    cop.sample4.non.con<-qbisa(sample[,con3], as.numeric(Coef(marginal_non_con4)[1]), as.numeric(Coef(marginal_non_con4)[2]))
  }
  if(Marginal_Dist4=="Exp"){
    cop.sample4.non.con<-qexp(sample[,con3], rate = as.numeric(marginal_non_con4$estimate[1]))
  }
  if(Marginal_Dist4=="Gam(2)"){
    cop.sample4.non.con<-qgamma(sample[,con3], shape = as.numeric(marginal_non_con4$estimate[1]), rate=as.numeric(marginal_non_con4$estimate[2]))
  }
  if(Marginal_Dist4=="Gam(3)"){
    cop.sample4.non.con<-qGG(sample[,con3], mu=exp(marginal_non_con4$mu.coefficients), sigma=exp(marginal_non_con4$sigma.coefficients), nu=marginal_non_con4$nu.coefficients)
  }
  if(Marginal_Dist4=="GamMix(2)"){
    xx <- seq(0, max(Data_Con4[,con3])*10, 0.001)
    prob.MX1 <- round(marginal_non_con4$prob[1],3)
    prob.MX2 <- 1 - prob.MX1
    cdf.MX<-pMX(xx, mu=list(mu1=exp(marginal_non_con4$models[[1]]$mu.coefficients), mu2=exp(marginal_non_con4$models[[2]]$mu.coefficients)),
                sigma=list(sigma1=exp(marginal_non_con4$models[[1]]$sigma.coefficients), sigma2=exp(marginal_non_con4$models[[2]]$sigma.coefficients)),
                pi = list(pi1=prob.MX1, pi2=prob.MX2), family=list(fam1="GA", fam2="GA"))
    cop.sample4.non.con <- approx(cdf.MX, xx, sample[,con3])$y
  }
  if(Marginal_Dist4=="GamMix(3)"){
    xx <- seq(0, max(Data_Con4[,con3])*10, 0.001)
    prob.MX1 <- round(marginal_non_con4$prob[1],3)
    prob.MX2 <- round(marginal_non_con4$prob[2],3)
    prob.MX3 <- 1 - prob.MX1 - prob.MX2
    cdf.MX<-pMX(xx, mu=list(mu1=exp(marginal_non_con4$models[[1]]$mu.coefficients), mu2=exp(marginal_non_con4$models[[2]]$mu.coefficients), mu3=exp(marginal_non_con4$models[[3]]$mu.coefficients)),
                sigma=list(sigma1=exp(marginal_non_con4$models[[1]]$sigma.coefficients), sigma2=exp(marginal_non_con4$models[[2]]$sigma.coefficients), sigma3=exp(marginal_non_con4$models[[3]]$sigma.coefficients)),
                pi = list(pi1=prob.MX1, pi2=prob.MX2, pi3=prob.MX3), family=list(fam1="GA", fam2="GA", fam3="GA"))
    cop.sample4.non.con <- approx(cdf.MX, xx, sample[,con3])$y
  }
  if(Marginal_Dist4=="Gaus"){
    cop.sample4.non.con<-qnorm(sample[,con3], mean = as.numeric(marginal_non_con4$estimate[1]), sd=as.numeric(marginal_non_con4$estimate[2]))
  }
  if(Marginal_Dist4=="Gum"){
    cop.sample4.non.con<-qGU(sample[,con3],as.numeric(marginal_non_con4$mu.coefficients),exp(as.numeric(marginal_non_con4$sigma.coefficients)))
  }
  if(Marginal_Dist4=="InvG"){
    cop.sample4.non.con<-qinvgauss(sample[,con3], mean = as.numeric(marginal_non_con4$estimate[1]), shape=as.numeric(marginal_non_con4$estimate[2]))
  }
  if(Marginal_Dist4=="Lapl"){
    cop.sample4.non.con <- qlaplace(sample[,con3],as.numeric(marginal_non_con4$estimate[1]), as.numeric(marginal_non_con4$estimate[2]))
  }
  if(Marginal_Dist4=="LogN"){
    cop.sample4.non.con<-qlnorm(sample[,con3], meanlog = as.numeric(marginal_non_con4$estimate[1]), sdlog = as.numeric(marginal_non_con4$estimate[2]))
  }
  if(Marginal_Dist4=="Logis"){
    cop.sample4.non.con<-qlogis(sample[,con3], location = as.numeric(marginal_non_con4$estimate[1]), scale=as.numeric(marginal_non_con4$estimate[2]))
  }
  if(Marginal_Dist4=="RGum"){
    cop.sample4.non.con<-qRG(sample[,con3],marginal_non_con4$mu.coefficients,exp(marginal_non_con4$sigma.coefficients))
  }
  if(Marginal_Dist4=="TNorm"){
    cop.sample4.non.con<-qtruncnorm(sample[,con3], a=min(Data_Con4[,con3]), mean = as.numeric(marginal_non_con4$estimate[1]), sd = as.numeric(marginal_non_con4$estimate[2]))
  }
  if(Marginal_Dist4=="Twe"){
    cop.sample4.non.con<-qtweedie(sample[,con3], power=marginal_non_con4$p.max, mu=mean(Data_Con4[,con3]), phi=marginal_non_con4$phi.max)
  }
  if(Marginal_Dist4=="Weib"){
    cop.sample4.non.con<-qweibull(sample[,con3], shape = as.numeric(marginal_non_con4$estimate[1]), scale=as.numeric(marginal_non_con4$estimate[2]))
  }
  #Put the realizations that have been transformed to the original scale in a data frame.
  cop.sample4<-data.frame(cop.sample4.non.con,cop.sample4.con)
  colnames(cop.sample4)<-c("Var1","Var2")

  #Combine the data frames containg the samples from two copulas (on the original scale)
  cop.sample.pop.2<-rbind(cop.sample3,cop.sample4)

  #Combined sample
  cop.sample <- rbind(cop.sample.pop.1,cop.sample.pop.2)

  ###Deriving the quantile isoline from the sample conditioned on variable 'Con2' i.e. Data_Con1

  for(k in 1:length(RP)){

    #Transform the points on the grid to the (0,1) using the fitted cumulative distribution
    #Transform the values of the conditioned variable of Data_Con1 (Con1) in the grid to the (0,1) scale using the CDF of the GPD contained in the u2gpd function
    #Population 1
    if(is.na(GPD1[[1]][1]) & is.na(Tab1[[1]][1])){
      con1.x.u<-pgpd(Pgrid[,1], u=Thres1 , sigma=exp(GPD_con1$coefficients[1]),xi= GPD_con1$coefficients[2] )
    }

    if(!is.na(GPD1[[1]][1])){
      con1.x.u<-pgpd(Pgrid[,1], u = GPD1$Threshold, sigma = GPD1$sigma, xi = GPD1$xi)
    }

    if(!is.na(Tab1[[1]][1])){
      con1.x = approx(Tab1[,2],1-1/Tab1[,1],xout=Pgrid[,1])$y
    }

    #Transform the values of the non-conditioned variable in Data_Con1 (Con2) to the (0,1) scale using the CDF function of the selected parameteric (non-extreme value) distributions
    if(Marginal_Dist1=="BS"){
      con1.y.u<-pbisa(Pgrid[,2],as.numeric(Coef(marginal_non_con1)[1]),as.numeric(Coef(marginal_non_con1)[2]))
    }
    if(Marginal_Dist1=="Exp"){
      con1.y.u<-pexp(Pgrid[,2],as.numeric(marginal_non_con1$estimate[1]))
    }
    if(Marginal_Dist1=="Gam(2)"){
      con1.y.u<-pgamma(Pgrid[,2],as.numeric(marginal_non_con1$estimate[1]),as.numeric(marginal_non_con1$estimate[2]))
    }
    if(Marginal_Dist1=="Gam(3)"){
      con1.y.u<-pGG(Pgrid[,2], mu=exp(marginal_non_con1$mu.coefficients), sigma=exp(marginal_non_con1$sigma.coefficients), nu=marginal_non_con1$nu.coefficients)
    }
    if(Marginal_Dist1=="GamMix(2)"){
      xx <- seq(0, max(Data_Con2[,2])*10, 0.001)
      prob.MX1 <- round(marginal_non_con1$prob[1],3)
      prob.MX2 <- 1 - prob.MX1
      cdf.MX<-pMX(xx, mu=list(mu1=exp(marginal_non_con1$models[[1]]$mu.coefficients), mu2=exp(marginal_non_con1$models[[2]]$mu.coefficients)),
                  sigma=list(sigma1=exp(marginal_non_con1$models[[1]]$sigma.coefficients), sigma2=exp(marginal_non_con1$models[[2]]$sigma.coefficients)),
                  pi = list(pi1=prob.MX1, pi2=prob.MX2), family=list(fam1="GA", fam2="GA"))
      con1.y.u <- approx(xx, cdf.MX, Pgrid[,2])$y
    }
    if(Marginal_Dist1=="GamMix(3)"){
      xx <- seq(0, max(Data_Con2[,2])*10, 0.001)
      prob.MX1 <- round(marginal_non_con1$prob[1],3)
      prob.MX2 <- round(marginal_non_con1$prob[2],3)
      prob.MX3 <- 1 - prob.MX1 - prob.MX2
      #prob.MX3 - round(fit.GamMIX3_GA$prob[3],5) < 0.00001
      cdf.MX<-pMX(xx, mu=list(mu1=exp(marginal_non_con1$models[[1]]$mu.coefficients), mu2=exp(marginal_non_con1$models[[2]]$mu.coefficients), mu3=exp(marginal_non_con1$models[[3]]$mu.coefficients)),
                  sigma=list(sigma1=exp(marginal_non_con1$models[[1]]$sigma.coefficients), sigma2=exp(marginal_non_con1$models[[2]]$sigma.coefficients), sigma3=exp(marginal_non_con1$models[[3]]$sigma.coefficients)),
                  pi = list(pi1=prob.MX1, pi2=prob.MX2, pi3=prob.MX3), family=list(fam1="GA", fam2="GA", fam3="GA"))
      con1.y.u <- approx(xx, cdf.MX, Pgrid[,2])$y
    }
    if(Marginal_Dist1=="Gaus"){
      con1.y.u<-pnorm(Pgrid[,2],as.numeric(marginal_non_con1$estimate[1]),as.numeric(marginal_non_con1$estimate[2]))
    }
    if(Marginal_Dist1=="Gum"){
      con1.y.u<-pGU(Pgrid[,2],as.numeric(marginal_non_con1$mu.coefficients),exp(as.numeric(marginal_non_con1$sigma.coefficients)))
    }
    if(Marginal_Dist1=="InvG"){
      con1.y.u<-pinvgauss(Pgrid[,2],as.numeric(marginal_non_con1$estimate[1]),as.numeric(marginal_non_con1$estimate[2]))
    }
    if(Marginal_Dist1=="Lapl"){
      con1.y.u <- plaplace(Pgrid[,2],as.numeric(marginal_non_con1$estimate[1]), as.numeric(marginal_non_con1$estimate[2]))
    }
    if(Marginal_Dist1=="Logis"){
      con1.y.u<-plogis(Pgrid[,2],as.numeric(marginal_non_con1$estimate[1]),as.numeric(marginal_non_con1$estimate[2]))
    }
    if(Marginal_Dist1=="LogN"){
      con1.y.u<-plnorm(Pgrid[,2],meanlog = as.numeric(marginal_non_con1$estimate[1]), sdlog = as.numeric(marginal_non_con1$estimate[2]))
    }
    if(Marginal_Dist1=="RGum"){
      con1.y.u<-pRG(Pgrid[,2],marginal_non_con1$mu.coefficients,exp(marginal_non_con1$sigma.coefficients))
    }
    if(Marginal_Dist1=="TNorm"){
      con1.y.u<-ptruncnorm(Pgrid[,2],a=min(Data_Con1[,con2]),as.numeric(marginal_non_con1$estimate[1]),as.numeric(marginal_non_con1$estimate[2]))
    }
    if(Marginal_Dist1=="Twe"){
      con1.y.u<-ptweedie(Pgrid[,2], power=marginal_non_con1$p.max, mu=mean(Data_Con1[,con2]), phi=marginal_non_con1$phi.max)
    }
    if(Marginal_Dist1=="Weib"){
      con1.y.u<-pweibull(Pgrid[,2],as.numeric(marginal_non_con1$estimate[1]),as.numeric(marginal_non_con1$estimate[2]))
    }

    #Evaluate the copula associated with the sample conditioned on Con1 at each grid point in the (0,1) space
    UU1<-BiCopCDF(con1.x.u, con1.y.u, obj1)

    #Calculate the inter-arrival time of extremes (in terms of mu) in Data_Con1.
    EL_con1<-1/Rate_Con1

    #Evaluate the return period at each point on the grid 'u'
    con1_func<-function(x,y){(1-x-y+UU1[which(con1.x.u==x & con1.y.u==y)]) / EL_con1 }
    AEP_con1 <- con1_func(con1.x.u,con1.y.u)

    ###Converting the grid to the unit square using the fitted distributions for the sample conditioned on variable 2.

    #Transforming the values of the conditioning variable in Data_Con2 (Con2) in the grid to the (0,1) scale using the CDF of the GPD.
    if(is.na(GPD2[[1]][1]) & is.na(Tab2[[1]][1])){
      con2.y.u<-pgpd(Pgrid[,2], u=Thres2 , sigma=exp(GPD_con2$coefficients[1]),xi= GPD_con2$coefficients[2] )
    }

    if(!is.na(GPD2[[1]][1])){
      con2.y.u<-pgpd(Pgrid[,2], u = GPD2$Threshold, sigma = GPD2$sigma, xi = GPD2$xi)
    }
    #if(is.na(Tab2[[1]][1])==F){
    #  con2.y = approx(1-1/Tab2[,1],Tab2[,2],xout=u2+as.numeric(unlist(xy160[[1]][3]))*(((1-1/(mu*RP[k]))-u2)/max(as.numeric(unlist(xy160[[1]][3])))))$y
    #}

    #Transform the values of the conditioned variable in Data_Con2 (Con1) to the (0,1) scale using the CDF function of the selected parameteric (non-extreme value) distributions.
    if(Marginal_Dist2=="BS"){
      con2.x.u<-pbisa(Pgrid[,1], as.numeric(Coef(marginal_non_con2)[1]),as.numeric(Coef(marginal_non_con2)[2]))
    }
    if(Marginal_Dist2=="Exp"){
      con2.x.u<-pexp(Pgrid[,1], as.numeric(marginal_non_con2$estimate[1]))
    }
    if(Marginal_Dist2=="Gam(2)"){
      con2.x.u<-pgamma(Pgrid[,1], shape = as.numeric(marginal_non_con2$estimate[1]), rate = as.numeric(marginal_non_con2$estimate[2]))
    }
    if(Marginal_Dist2=="Gam(3)"){
      con2.x.u<-pGG(Pgrid[,1], mu=exp(marginal_non_con2$mu.coefficients), sigma=exp(marginal_non_con2$sigma.coefficients), nu=marginal_non_con2$nu.coefficients)
    }
    if(Marginal_Dist2=="GamMix(2)"){
      xx <- seq(0, max(Data_Con1[,1])*10, 0.001)
      prob.MX1 <- round(marginal_non_con2$prob[1],3)
      prob.MX2 <- 1 - prob.MX1
      cdf.MX<-pMX(xx, mu=list(mu1=exp(marginal_non_con2$models[[1]]$mu.coefficients), mu2=exp(marginal_non_con2$models[[2]]$mu.coefficients)),
                  sigma=list(sigma1=exp(marginal_non_con2$models[[1]]$sigma.coefficients), sigma2=exp(marginal_non_con2$models[[2]]$sigma.coefficients)),
                  pi = list(pi1=prob.MX1, pi2=prob.MX2), family=list(fam1="GA", fam2="GA"))
      con2.x <- approx(xx, cdf.MX, Pgrid[,1])$y
    }
    if(Marginal_Dist2=="GamMix(3)"){
      xx <- seq(0, max(Data_Con1[,1])*10, 0.001)
      prob.MX1 <- round(marginal_non_con2$prob[1],3)
      prob.MX2 <- round(marginal_non_con2$prob[2],3)
      prob.MX3 <- 1 - prob.MX1 - prob.MX2
      #prob.MX3 - round(fit.GamMIX3_GA$prob[3],5) < 0.00001
      cdf.MX<-pMX(xx, mu=list(mu1=exp(marginal_non_con2$models[[1]]$mu.coefficients), mu2=exp(marginal_non_con2$models[[2]]$mu.coefficients), mu3=exp(marginal_non_con2$models[[3]]$mu.coefficients)),
                  sigma=list(sigma1=exp(marginal_non_con2$models[[1]]$sigma.coefficients), sigma2=exp(marginal_non_con2$models[[2]]$sigma.coefficients), sigma3=exp(marginal_non_con2$models[[3]]$sigma.coefficients)),
                  pi = list(pi1=prob.MX1, pi2=prob.MX2, pi3=prob.MX3), family=list(fam1="GA", fam2="GA", fam3="GA"))
      con2.x <- approx(xx, cdf.MX, Pgrid[,1])$y
    }
    if(Marginal_Dist2=="Gaus"){
      con2.x.u<-pnorm(Pgrid[,1], as.numeric(marginal_non_con2$estimate[1]), as.numeric(marginal_non_con2$estimate[2]))
    }
    if(Marginal_Dist2=="Gum"){
      con2.x.u<-qGU(Pgrid[,1],as.numeric(marginal_non_con2$mu.coefficients),exp(as.numeric(marginal_non_con2$sigma.coefficients)))
    }
    if(Marginal_Dist2=="InvG"){
      con2.x.u<-pinvgauss(Pgrid[,1], as.numeric(marginal_non_con2$estimate[1]), as.numeric(marginal_non_con2$estimate[2]))
    }
    if(Marginal_Dist2=="Lapl"){
      con2.x.u <- plaplace(Pgrid[,1],as.numeric(marginal_non_con2$estimate[1]), as.numeric(marginal_non_con2$estimate[2]))
    }
    if(Marginal_Dist2=="Logis"){
      con2.x.u<-plogis(Pgrid[,1],as.numeric(marginal_non_con2$estimate[1]),as.numeric(marginal_non_con2$estimate[2]))
    }
    if(Marginal_Dist2=="LogN"){
      con2.x.u<-plnorm(Pgrid[,1], meanlog = as.numeric(marginal_non_con2$estimate[1]), sdlog = as.numeric(marginal_non_con2$estimate[2]))
    }
    if(Marginal_Dist2=="RGum"){
      con2.x.u<-pRG(Pgrid[,1],marginal_non_con2$mu.coefficients,exp(marginal_non_con2$sigma.coefficients))
    }
    if(Marginal_Dist2=="TNorm"){
      con2.x.u<-ptruncnorm(Pgrid[,1],a=min(Data_Con2[,con1]),as.numeric(marginal_non_con2$estimate[1]),as.numeric(marginal_non_con2$estimate[2]))
    }
    if(Marginal_Dist2=="Twe"){
      con2.x.u<-ptweedie(Pgrid[,1], power=marginal_non_con2$p.max, mu=mean(Data_Con2[,con1]), phi=marginal_non_con2$phi.max)
    }
    if(Marginal_Dist2=="Weib"){
      con2.x.u<-pweibull(Pgrid[,1], as.numeric(marginal_non_con2$estimate[1]), as.numeric(marginal_non_con2$estimate[2]))
    }

    #Evaluate the copula CDF at each grid point
    UU2<-BiCopCDF(con2.x.u, con2.y.u, obj2)

    #Interarrival time of events in sample conditioning on var 2
    EL_con2<-1/Rate_Con2

    #Function for evaluating the AEP at each point in the event space
    con2_func<-function(x,y){ (1-x-y+UU2[which(con2.x.u==x & con2.y.u==y)]) / EL_con2 }
    AEP_con2 <- con2_func(con2.x.u, con2.y.u)

    #Put the two AEPs into a dataframe
    AEP = data.frame(AEP_con1,AEP_con2)

    #Select the maximum of the AEPs from the two samples at each grid point
    AEP = apply(AEP,1,function(x) max(x, na.rm = TRUE))

    #Convert data to matrix to match dimensions of PGrid
    AEP.pop.1 = matrix(AEP,nrow=length(var1))

    ##Population 2

    #Transform the points on the grid to the (0,1) using the fitted cumulative distribution
    #Transform the values of the conditioned variable of Data_Con1 (Con1) in the grid to the (0,1) scale using the CDF of the GPD contained in the u2gpd function

    if(is.na(GPD3[[1]][1]) & is.na(Tab3[[1]][1])){
      con3.x.u<-pgpd(Pgrid[,1], u=Thres3 , sigma=exp(GPD_con3$coefficients[1]),xi= GPD_con3$coefficients[2] )
    }

    if(!is.na(GPD3[[1]][1])){
      con3.x.u<-pgpd(Pgrid[,1], u = GPD3$Threshold, sigma = GPD3$sigma, xi = GPD3$xi)
    }

    if(!is.na(Tab1[[1]][1])){
      con3.x.u = approx(Tab1[,2],1-1/Tab1[,1],xout=Pgrid[,1])$y
    }

    #Transform the values of the non-conditioned variable in Data_Con1 (Con2) to the (0,1) scale using the CDF function of the selected parameteric (non-extreme value) distributions
    if(Marginal_Dist3=="BS"){
      con3.y.u<-pbisa(Pgrid[,2],as.numeric(Coef(marginal_non_con3)[1]),as.numeric(Coef(marginal_non_con3)[2]))
    }
    if(Marginal_Dist3=="Exp"){
      con3.y.u<-pexp(Pgrid[,2],as.numeric(marginal_non_con3$estimate[1]))
    }
    if(Marginal_Dist3=="Gam(2)"){
      con3.y.u<-pgamma(Pgrid[,2],as.numeric(marginal_non_con3$estimate[1]),as.numeric(marginal_non_con3$estimate[2]))
    }
    if(Marginal_Dist3=="Gam(3)"){
      con3.y.u<-pGG(Pgrid[,2], mu=exp(marginal_non_con3$mu.coefficients), sigma=exp(marginal_non_con3$sigma.coefficients), nu=marginal_non_con3$nu.coefficients)
    }
    if(Marginal_Dist3=="GamMix(2)"){
      xx <- seq(0, max(Data_Con3[,con4])*10, 0.001)
      prob.MX1 <- round(marginal_non_con3$prob[1],3)
      prob.MX2 <- 1 - prob.MX1
      cdf.MX<-pMX(xx, mu=list(mu1=exp(marginal_non_con3$models[[1]]$mu.coefficients), mu2=exp(marginal_non_con3$models[[2]]$mu.coefficients)),
                  sigma=list(sigma1=exp(marginal_non_con3$models[[1]]$sigma.coefficients), sigma2=exp(marginal_non_con3$models[[2]]$sigma.coefficients)),
                  pi = list(pi1=prob.MX1, pi2=prob.MX2), family=list(fam1="GA", fam2="GA"))
      con3.y.u <- approx(xx, cdf.MX, Pgrid[,2])$y
    }
    if(Marginal_Dist3=="GamMix(3)"){
      xx <- seq(0, max(Data_Con3[,con4])*10, 0.001)
      prob.MX1 <- round(marginal_non_con3$prob[1],3)
      prob.MX2 <- round(marginal_non_con3$prob[2],3)
      prob.MX3 <- 1 - prob.MX1 - prob.MX2
      #prob.MX3 - round(fit.GamMIX3_GA$prob[3],5) < 0.00001
      cdf.MX<-pMX(xx, mu=list(mu1=exp(marginal_non_con3$models[[1]]$mu.coefficients), mu2=exp(marginal_non_con3$models[[2]]$mu.coefficients), mu3=exp(marginal_non_con3$models[[3]]$mu.coefficients)),
                  sigma=list(sigma1=exp(marginal_non_con3$models[[1]]$sigma.coefficients), sigma2=exp(marginal_non_con3$models[[2]]$sigma.coefficients), sigma3=exp(marginal_non_con3$models[[3]]$sigma.coefficients)),
                  pi = list(pi1=prob.MX1, pi2=prob.MX2, pi3=prob.MX3), family=list(fam1="GA", fam2="GA", fam3="GA"))
      con3.y.u <- approx(xx, cdf.MX, Pgrid[,2])$y
    }
    if(Marginal_Dist3=="Gaus"){
      con3.y.u<-pnorm(Pgrid[,2],as.numeric(marginal_non_con3$estimate[1]),as.numeric(marginal_non_con3$estimate[2]))
    }
    if(Marginal_Dist3=="Gum"){
      con3.x.u<-qGU(Pgrid[,2],as.numeric(marginal_non_con3$mu.coefficients),exp(as.numeric(marginal_non_con3$sigma.coefficients)))
    }
    if(Marginal_Dist3=="InvG"){
      con3.y.u<-pinvgauss(Pgrid[,2],as.numeric(marginal_non_con3$estimate[1]),as.numeric(marginal_non_con3$estimate[2]))
    }
    if(Marginal_Dist3=="Lapl"){
      con3.y.u <- plaplace(Pgrid[,2],as.numeric(marginal_non_con3$estimate[1]), as.numeric(marginal_non_con3$estimate[2]))
    }
    if(Marginal_Dist3=="Logis"){
      con3.y.u<-plogis(Pgrid[,2],as.numeric(marginal_non_con3$estimate[1]),as.numeric(marginal_non_con3$estimate[2]))
    }
    if(Marginal_Dist3=="LogN"){
      con3.y.u<-plnorm(Pgrid[,2],meanlog = as.numeric(marginal_non_con3$estimate[1]), sdlog = as.numeric(marginal_non_con3$estimate[2]))
    }
    if(Marginal_Dist3=="RGum"){
      con3.y.u<-pRG(Pgrid[,2],marginal_non_con3$mu.coefficients,exp(marginal_non_con3$sigma.coefficients))
    }
    if(Marginal_Dist3=="TNorm"){
      con3.y.u<-ptruncnorm(Pgrid[,2],a=min(Data_Con3[,con4]),as.numeric(marginal_non_con3$estimate[1]),as.numeric(marginal_non_con3$estimate[2]))
    }
    if(Marginal_Dist3=="Twe"){
      con3.y.u<-ptweedie(Pgrid[,2], power=marginal_non_con3$p.max, mu=mean(Data_Con3[,con4]), phi=marginal_non_con3$phi.max)
    }
    if(Marginal_Dist3=="Weib"){
      con3.y.u<-pweibull(Pgrid[,2],as.numeric(marginal_non_con3$estimate[1]),as.numeric(marginal_non_con3$estimate[2]))
    }

    #Evaluate the copula associated with the sample conditioned on Con1 at each grid point in the (0,1) space
    UU3<-BiCopCDF(con3.x.u, con3.y.u, obj3)

    #Calculate the inter-arrival time of extremes (in terms of mu) in Data_Con1.
    EL_con3<-1/Rate_Con3

    #Evaluate the return period at each point on the grid 'u'
    con3_func<-function(x,y){(1-x-y+UU3[which(con3.x.u==x & con3.y.u==y)]) / EL_con3 }
    AEP_con3 <- con3_func(con3.x.u,con3.y.u)

    ###Converting the grid to the unit square using the fitted distributions for the sample conditioned on variable 2.

    #Transforming the values of the conditioning variable in Data_Con2 (Con2) in the grid to the (0,1) scale using the CDF of the GPD.
    if(is.na(GPD4[[1]][1]) & is.na(Tab4[[1]][1])){
      con4.y.u<-pgpd(Pgrid[,2], u=Thres4 , sigma=exp(GPD_con4$coefficients[1]),xi= GPD_con4$coefficients[2] )
    }

    if(!is.na(GPD4[[1]][1])){
      con4.y.u<-pgpd(Pgrid[,2], u = GPD4$Threshold, sigma = GPD4$sigma, xi = GPD4$xi)
    }
    #if(is.na(Tab2[[1]][1])==F){
    #  con2.y = approx(1-1/Tab2[,1],Tab2[,2],xout=u2+as.numeric(unlist(xy160[[1]][3]))*(((1-1/(mu*RP[k]))-u2)/max(as.numeric(unlist(xy160[[1]][3])))))$y
    #}

    #Transform the values of the conditioned variable in Data_Con2 (Con1) to the (0,1) scale using the CDF function of the selected parameteric (non-extreme value) distributions.
    if(Marginal_Dist4=="BS"){
      con4.x.u<-pbisa(Pgrid[,1], as.numeric(Coef(marginal_non_con4)[1]),as.numeric(Coef(marginal_non_con4)[2]))
    }
    if(Marginal_Dist4=="Exp"){
      con4.x.u<-pexp(Pgrid[,1], as.numeric(marginal_non_con4$estimate[1]))
    }
    if(Marginal_Dist4=="Gam(2)"){
      con4.x.u<-pgamma(Pgrid[,1], shape = as.numeric(marginal_non_con4$estimate[1]), rate = as.numeric(marginal_non_con4$estimate[2]))
    }
    if(Marginal_Dist4=="Gam(3)"){
      con4.x.u<-pGG(Pgrid[,1], mu=exp(marginal_non_con4$mu.coefficients), sigma=exp(marginal_non_con4$sigma.coefficients), nu=marginal_non_con4$nu.coefficients)
    }
    if(Marginal_Dist4=="GamMix(2)"){
      xx <- seq(0, max(Data_Con4[,con3])*10, 0.001)
      prob.MX1 <- round(marginal_non_con4$prob[1],3)
      prob.MX2 <- 1 - prob.MX1
      cdf.MX<-pMX(xx, mu=list(mu1=exp(marginal_non_con4$models[[1]]$mu.coefficients), mu2=exp(marginal_non_con4$models[[2]]$mu.coefficients)),
                  sigma=list(sigma1=exp(marginal_non_con4$models[[1]]$sigma.coefficients), sigma2=exp(marginal_non_con4$models[[2]]$sigma.coefficients)),
                  pi = list(pi1=prob.MX1, pi2=prob.MX2), family=list(fam1="GA", fam2="GA"))
      con4.x.u <- approx(xx, cdf.MX, Pgrid[,1])$y
    }
    if(Marginal_Dist4=="GamMix(3)"){
      xx <- seq(0, max(Data_Con4[,con3])*10, 0.001)
      prob.MX1 <- round(marginal_non_con4$prob[1],3)
      prob.MX2 <- round(marginal_non_con4$prob[2],3)
      prob.MX3 <- 1 - prob.MX1 - prob.MX2
      #prob.MX3 - round(fit.GamMIX3_GA$prob[3],5) < 0.00001
      cdf.MX<-pMX(xx, mu=list(mu1=exp(marginal_non_con4$models[[1]]$mu.coefficients), mu2=exp(marginal_non_con4$models[[2]]$mu.coefficients), mu3=exp(marginal_non_con4$models[[3]]$mu.coefficients)),
                  sigma=list(sigma1=exp(marginal_non_con4$models[[1]]$sigma.coefficients), sigma2=exp(marginal_non_con4$models[[2]]$sigma.coefficients), sigma3=exp(marginal_non_con4$models[[3]]$sigma.coefficients)),
                  pi = list(pi1=prob.MX1, pi2=prob.MX2, pi3=prob.MX3), family=list(fam1="GA", fam2="GA", fam3="GA"))
      con4.x.u <- approx(xx, cdf.MX, Pgrid[,1])$y
    }
    if(Marginal_Dist4=="Gaus"){
      con4.x.u<-pnorm(Pgrid[,1], as.numeric(marginal_non_con4$estimate[1]), as.numeric(marginal_non_con4$estimate[2]))
    }
    if(Marginal_Dist4=="Gum"){
      con4.x.u<-qGU(Pgrid[,1],as.numeric(marginal_non_con4$mu.coefficients),exp(as.numeric(marginal_non_con4$sigma.coefficients)))
    }
    if(Marginal_Dist4=="InvG"){
      con4.x.u<-pinvgauss(Pgrid[,1], as.numeric(marginal_non_con4$estimate[1]), as.numeric(marginal_non_con4$estimate[2]))
    }
    if(Marginal_Dist4=="Lapl"){
      con4.x.u <- plaplace(Pgrid[,1],as.numeric(marginal_non_con4$estimate[1]), as.numeric(marginal_non_con4$estimate[2]))
    }
    if(Marginal_Dist4=="Logis"){
      con4.x.u<-plogis(Pgrid[,1],as.numeric(marginal_non_con4$estimate[1]),as.numeric(marginal_non_con4$estimate[2]))
    }
    if(Marginal_Dist4=="LogN"){
      con4.x.u<-plnorm(Pgrid[,1], meanlog = as.numeric(marginal_non_con4$estimate[1]), sdlog = as.numeric(marginal_non_con4$estimate[2]))
    }
    if(Marginal_Dist4=="RGum"){
      con4.x.u<-pRG(Pgrid[,1],marginal_non_con4$mu.coefficients,exp(marginal_non_con4$sigma.coefficients))
    }
    if(Marginal_Dist4=="TNorm"){
      con4.x.u<-ptruncnorm(Pgrid[,1],a=min(Data_Con4[,con3]),as.numeric(marginal_non_con4$estimate[1]),as.numeric(marginal_non_con4$estimate[2]))
    }
    if(Marginal_Dist4=="Twe"){
      con4.x.u<-ptweedie(Pgrid[,1], power=marginal_non_con4$p.max, mu=mean(Data_Con4[,con3]), phi=marginal_non_con4$phi.max)
    }
    if(Marginal_Dist4=="Weib"){
      con4.x.u<-pweibull(Pgrid[,1], as.numeric(marginal_non_con4$estimate[1]), as.numeric(marginal_non_con4$estimate[2]))
    }

    #Evaluate the copula CDF at each grid point
    UU4<-BiCopCDF(con4.x.u, con4.y.u, obj4)

    #Interarrival time of events in sample conditioning on var 2
    EL_con4<-1/Rate_Con4

    #Function for evaluating the AEP at each point in the event space
    con4_func<-function(x,y){ (1-x-y+UU4[which(con4.x.u==x & con4.y.u==y)]) / EL_con4 }
    AEP_con4 <- con4_func(con4.x.u, con4.y.u)

    #Put the two AEPs into a dataframe
    AEP = data.frame(AEP_con3,AEP_con4)

    #Select the maximum of the AEPs from the two samples at each grid point
    AEP = apply(AEP,1,function(x) max(x, na.rm = TRUE))

    #Convert data to matrix to match dimensions of PGrid
    AEP.pop.2 = matrix(AEP,nrow=length(var1))

    #Populations are independent so overall non-exceedance probability is product of individual non-exceedence probabilities
    ANEP = (1 - AEP.pop.1) * (1 - AEP.pop.2)

    #Convert non-excedence probabilities to return periods
    RP_Grid = 1 / ( 1 - ANEP )

    #Compute isoline
    iso = contourLines(var1,var2,RP_Grid,levels= RP[k])
    Isoline[[k]] = data.frame(as.numeric(unlist(iso[[1]][2])),as.numeric(unlist(iso[[1]][3])))
    Iso = Isoline[[k]]

    #Remove very extreme values that can exert significant leverage on contours
    remove<-which(cop.sample[,1] > Sim_Max*max(Data[,1],na.rm=T) | cop.sample[,2] > Sim_Max*max(Data[,2],na.rm=T))
    if(Isoline_Probs=="Sample"){
      if(length(remove)>1){
        cop.sample<-cop.sample[-remove,]
      }
      prediction<-kde(x=cop.sample, eval.points=Iso)$estimate
    }
    if(Isoline_Probs=="Observations"){
      prediction<-kde(x=na.omit(Data), eval.points=Iso)$estimate
    }      #(relative) Probabilities implied by the data for the points composing the isoline. Probabilities are scaled to [0,1].
    Contour[[k]] <- (prediction-min(prediction))/(max(prediction)-min(prediction))

    ###Extract design event(s)

    #Find the 'most likely' design event and add it to the plot (denoted by a diamond).
    MostLikelyEvent.AND<-data.frame(as.numeric(Iso[which(prediction==max(prediction,na.rm=T)),1]),as.numeric(Iso[which(prediction==max(prediction,na.rm=T)),2]))
    colnames(MostLikelyEvent.AND) <- c(names(Data)[1],names(Data)[2])
    MostLikelyEvent[[k]]<-MostLikelyEvent.AND

    FullDependence.AND<-data.frame(max(Iso[,1]),max(Iso[,2]))
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

  ###Plot the isoline

 #Find the minimum and maximum x- and y-axis limits for the plot. If the limits are not specified in the input use the minimum and maximum values of the Data.
 x_min<-ifelse(is.na(x_lim_min)==T,min(na.omit(Data[,con1])),x_lim_min)
 x_max<-ifelse(is.na(x_lim_max)==T,max(na.omit(Data[,con1])),x_lim_max)
 y_min<-ifelse(is.na(y_lim_min)==T,min(na.omit(Data[,con2])),y_lim_min)
 y_max<-ifelse(is.na(y_lim_max)==T,max(na.omit(Data[,con2])),y_lim_max)

 #Plot
 par(mar=c(4.5,4.2,0.5,0.5))
 plot(Data[, con1], Data[, con2], xlim = c(x_min, x_max), ylim = c(y_min, y_max), col = "Light Grey",xlab = x_lab, ylab = y_lab, cex.lab = 1.5, cex.axis = 1.5)
 points(Data_Con1[,con1],Data_Con1[,con2],col="Red",pch=4,cex=1.5)
 points(Data_Con2[,con1],Data_Con2[,con2],col="Red",pch=4,cex=1.5)
 points(Data_Con3[,con3],Data_Con3[,con4],col="Blue",pch=1,cex=1.5)
 points(Data_Con4[,con3],Data_Con4[,con4],col="Blue",pch=1,cex=1.5)
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

 RPs = as.numeric(RP_Grid)

 #Create a list of outputs.
 res<-list("FullDependence" = FullDependence, "MostLikelyEvent" = MostLikelyEvent, "Ensemble"=Ensemble, "Isoline" = Isoline, "Contour"= Contour, "Quantile_Isoline_1" = Quantile_Isoline_1, "Quantile_Isoline_2" = Quantile_Isoline_2, "Threshold_1" = Thres1, "Threshold_2"=Thres2)
 return(res)
}
