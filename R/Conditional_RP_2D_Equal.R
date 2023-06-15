#' Calculates joint and conditional return periods
#'
#' A large number of realizations are simulated from the copulas fit to the conditioned samples, in proportion with the sizes of the conditional samples. The realization are transformed to the original scale and the relevant probabilities estimated empirically. The conditional probabilities return period of the conditioning variable equals
#'
#' @param Data Data frame of dimension \code{nx2} containing two co-occurring time series of length \code{n}.
#' @param Data_Con1 Data frame containing the conditional sample (declustered excesses paired with concurrent values of other variable), conditioned on the variable in the first column.
#' @param Data_Con2 Data frame containing the conditional sample (declustered excesses paired with concurrent values of other variable), conditioned on the variable in the second column. Can be obtained using the \code{Con_Sampling_2D} function.
#' @param u1 Numeric vector of length one specifying the (quantile) threshold above which the variable in the first column was sampled in \code{Data_Con1}.
#' @param u2 Numeric vector of length one specifying the (quantile) threshold above which the variable in the second column was sampled in \code{Data_Con2}.
#' @param Thres1 Numeric vector of length one specifying the threshold above which the variable in the first column was sampled in \code{Data_Con1}. Only one of \code{u1} and \code{Thres1} should be supplied. Default is \code{NA}.
#' @param Thres2 Numeric vector of length one specifying the threshold above which the variable in the second column was sampled in \code{Data_Con2}. Only one of \code{u2} and \code{Thres2} should be supplied. Default is \code{NA}.
#' @param Copula_Family1 Numeric vector of length one specifying the copula family used to model the \code{Data_Con1} dataset.
#' @param Copula_Family2 Numeric vector of length one specifying the copula family used to model the \code{Data_Con2} dataset. Best fitting of 40 copulas can be found using the \code{Copula_Threshold_2D} function.
#' @param Marginal_Dist1 Character vector of length one specifying (non-extreme) distribution used to model the marginal distribution of the non-conditioned variable in \code{Data_Con1}.
#' @param Marginal_Dist2 Character vector of length one specifying (non-extreme) distribution used to model the marginal distribution of the non-conditioned variable in \code{Data_Con2}.
#' @param Con1 Character vector of length one specifying the name of variable in the first column of \code{Data}.
#' @param Con2 Character vector of length one specifying the name of variable in the second column of \code{Data}.
#' @param mu Numeric vector of length one specifying the (average) occurrence frequency of events in \code{Data}. Default is \code{365.25}, daily data.
#' @param Con_Var Character vector of length one specifying the (column) name of the conditioning variable.
#' @param RP_Con Numeric vector of length one specifying the return period of the conditioning variable \code{Con_Var}.
#' @param RP_Non_Con Numeric vector of length one specifying the return period of the non-conditioning variable.
#' @param Width Numeric vector of length one specifying the distance above and below the \code{RP_Con} event of \code{Con_Var} the simulated events are used to estimate the conditional probability.
#' @param x_lab Character vector specifying the x-axis label.
#' @param y_lab Character vector specifying the y-axis label.
#' @param x_lim_min Numeric vector of length one specifying x-axis minimum. Default is \code{NA}.
#' @param x_lim_max Numeric vector of length one specifying x-axis maximum. Default is \code{NA}.
#' @param y_lim_min Numeric vector of length one specifying y-axis minimum. Default is \code{NA}.
#' @param y_lim_max Numeric vector of length one specifying y-axis maximum. Default is \code{NA}.
#' @param N Numeric vector of length one specifying the size of the sample from the fitted joint distributions used to estimate the density along an isoline. Samples are collected from the two joint distribution with proportions consistent with the total number of extreme events conditioned on each variable. Default is \code{10^6}
#' @return  Console output: \itemize{
#' \item Con_Var
#' Name of the conditioning variable
#' \item RP_Var1
#' Return period of variable Con1 i.e., variable in second column of \code{Data}
#' \item RP_Var2
#' Return period of variable Con2 i.e., variable in third column of \code{Data}
#' \item Var1
#' Value of Con1 at the return period of interest
#' \item Var2
#' Value of Con2 at the return period of interest
#' \item RP_Full_Dependence
#' Joint return period of the (Var1,Var2) event under full dependence
#' \item RP_Independence
#' Joint return period of the (Var1,Var2) event under independence
#' \item RP_Copula
#' Joint return period of the (Var1,Var2) event according to the two sided conditional sampling - copula theory approach
#' \item Prob
#' Probability associated with \code{RP_Copula}
#' \item N_Sub_Sample
#' Number of realizations of the \code{Con_Var} within +/- width of the value of \code{Con_Var} with retunr period \code{}.
#' \item Non_Con_Var_X
#' Values of the non-conditioned variable of the (conditional) Cummulative Distribution Function (CDF) i.e. x-axis of bottom left plot
#' \item Con_Prob
#' \code{Con_Prob} CDF of the non-conditioned variable given the return period of \code{Con_Var} equals \code{RP_Con}
#' \item Con_Prob_Est
#' Probability the non-conditioned variable is less than or equal to \code{RP_Non_Con} given the return period of \code{Con_Var} equals \code{RP_Con}
#' }
#' Graphical output: \itemize{
#' \item Top Left: Sample conditioned on rainfall (red crosses) and O-sWL (blue circles). Black dot is the event with a marginal return period of the conditioned variable \code{Var_Con} and non-conditioned variable equal to \code{RP_Con} and \code{RP_Non_Con}, respectively. The joint return period of the event using the conditional sampling - copula theory approach and under the assumptions of full dependence and independence between the variables are printed.
#' \item Top Right: Sample used to estimate the joint return period of the event of interest. Black dots denote the \code{N_Excess} sized subset of the sample where the marginal return period of the conditioned variable \code{Var_Con} exceeds \code{RP_Con} (years). The subset is used to estimate the conditional probabilities in part two of the question.
#' \item Bottom Left: Conditional Cumulative Distribution Function (CDF) of the non-conditioned variable given the marginal return period of the conditioned variable \code{Var_Con} exceeds \code{RP_Con} years i.e. the black dots in the top right plot.
#' \item Bottom Right: Conditional return period of the non-conditioned variable given the conditioned variable \code{Var_Con} has a return period longer than \code{RP_Con}.
#' }
#' @seealso \code{\link{Design_Event_2D}} \code{\link{Conditional_RP_2D}}
#' @export
#' @examples
#' #Under a 10yr rainfall event condition, what is the joint probability that a 10yr surge (O-sWL)
#' #event occurs simultaneously?  What is the cumulative probability of events with the frequency
#' #equal to or less than a 10yr surge event?
#' Conditional_RP_2D_Equal(Data=S22.Detrend.df,
#'                         Data_Con1=con.sample.Rainfall$Data, Data_Con2=con.sample.OsWL$Data,
#'                         u1=0.98, u2=0.98,
#'                         Copula_Family1=cop.Rainfall,Copula_Family2=cop.OsWL,
#'                         Marginal_Dist1="Logis", Marginal_Dist2="Twe",
#'                         Con1 = "Rainfall", Con2 = "OsWL",
#'                         mu = 365.25,
#'                         Con_Var="Rainfall",
#'                         RP_Con=10, RP_Non_Con=10,
#'                         x_lab = "Rainfall (Inches)", y_lab = "O-sWL (ft NGVD 29)",
#'                         y_lim_max = 10,
#'                         N=10^8)
Conditional_RP_2D_Equal<-function(Data, Data_Con1, Data_Con2,
                                  Thres1, Thres2, u1, u2,
                                  Copula_Family1,Copula_Family2,
                                  Marginal_Dist1, Marginal_Dist2, Con1 = "Rainfall",
                                  Con2 = "OsWL", mu = 365.25, Con_Var, RP_Con, RP_Non_Con, Width=0.1, x_lab = "Rainfall (mm)",
                                  y_lab = "O-sWL (mNGVD 29)", x_lim_min = NA, x_lim_max = NA,
                                  y_lim_min = NA, y_lim_max = NA, N){
  ###Preliminaries

  #Remove 1st column of Data if it is a Date or factor object.
  if (class(Data[, 1])[1] == "Date" | class(Data[, 1])[1] == "factor" | class(Data[,1])[1]=="POSIXct" | class(Data[,1])[1] == "character") {
    Data <- Data[, -1]
  }

  #Find the columns in Data (which should be consistent in terms of column order of the other data input objects)
  #of conditioning variable 1 (Con1) and conditioning variable 2 (Con2).
  con1 <- which(names(Data) == Con1)
  con2 <- which(names(Data) == Con2)
  con_var <-which(names(Data) == Con_Var)
  RP_Var1<-ifelse(con_var==1,RP_Con,RP_Non_Con)
  RP_Var2<-ifelse(con_var==2,RP_Con,RP_Non_Con)

  #Axis limits for plots
  x_min <- ifelse(is.na(x_lim_min) == T, min(na.omit(Data[, con1])), x_lim_min)
  x_max <- ifelse(is.na(x_lim_max) == T, max(na.omit(Data[, con1])), x_lim_max)
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
  Var1<-as.numeric(u2gpd((1-EL_Con1/RP_Var1), p = 1, th=Thres1 , sigma=exp(GPD_con1$coefficients[1]),xi= GPD_con1$coefficients[2]))

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
    marginal_non_con1 <- fitdistr(Data_Con1[, con2], "gamma")
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
  Var2<-as.numeric(u2gpd((1-EL_Con2/RP_Var2), p = 1, th= Thres2 , sigma=exp(GPD_con2$coefficients[1]),xi= GPD_con2$coefficients[2]))

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

  ###Simulating sample from the joint distribution (copula+marginals) fit to the sample conditioned on Con1
  #Fit the specified copula family (Copula_Family1) to the observations in Data_Con1.
  obj1 <- BiCopSelect(pobs(Data_Con1[, 1]), pobs(Data_Con1[,2]), familyset = Copula_Family1, selectioncrit = "AIC",
                      indeptest = FALSE, level = 0.05, weights = NA, rotations = TRUE,
                      se = FALSE, presel = TRUE, method = "mle")
  #Simulate a sample from the fitted copula. Out of the sample size 'N' the proportion of the sample from the
  #copula assoicated with Data_Con1 is proportional to the size of Data_Con1 relative to Data_Con2.
  sample <- BiCopSim(round(N * nrow(Data_Con1)/(nrow(Data_Con1) +
                                                  nrow(Data_Con2)), 0), obj1)
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
    prob.MX1 <- round(marginal_non_con1$prob[1],3)
    prob.MX2 <- 1 - prob.MX1
    cop.sample1.non.con<-qMX(sample[,con2], mu=list(mu1=exp(marginal_non_con1$models[[1]]$mu.coefficients), mu2=exp(marginal_non_con1$models[[2]]$mu.coefficients)),
                             sigma=list(sigma1=exp(marginal_non_con1$models[[1]]$sigma.coefficients), sigma2=exp(marginal_non_con1$models[[2]]$sigma.coefficients)),
                             pi = list(pi1=prob.MX1, pi2=prob.MX2), family=list(fam1="GA", fam2="GA"))
  }
  if(Marginal_Dist1=="GamMix(3)"){
    prob.MX1 <- round(marginal_non_con1$prob[1],3)
    prob.MX2 <- round(marginal_non_con1$prob[2],3)
    prob.MX3 <- 1 - prob.MX1 - prob.MX2
    cop.sample1.non.con<-qMX(sample[,con2], mu=list(mu1=exp(marginal_non_con1$models[[1]]$mu.coefficients), mu2=exp(marginal_non_con1$models[[2]]$mu.coefficients), mu3=exp(marginal_non_con1$models[[3]]$mu.coefficients)),
                             sigma=list(sigma1=exp(marginal_non_con1$models[[1]]$sigma.coefficients), sigma2=exp(marginal_non_con1$models[[2]]$sigma.coefficients), sigma3=exp(marginal_non_con1$models[[3]]$sigma.coefficients)),
                             pi = list(pi1=prob.MX1, pi2=prob.MX2, pi3=prob.MX3), family=list(fam1="GA", fam2="GA", fam3="GA"))
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
  obj2 <- BiCopSelect(pobs(Data_Con2[, 1]), pobs(Data_Con2[,
                                                           2]), familyset = Copula_Family2, selectioncrit = "AIC",
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
    prob.MX1 <- round(marginal_non_con2$prob[1],3)
    prob.MX2 <- 1 - prob.MX1
    cop.sample2.non.con<-qMX(sample[,con1], mu=list(mu1=exp(marginal_non_con2$models[[1]]$mu.coefficients), mu2=exp(marginal_non_con2$models[[2]]$mu.coefficients)),
                             sigma=list(sigma1=exp(marginal_non_con2$models[[1]]$sigma.coefficients), sigma2=exp(marginal_non_con2$models[[2]]$sigma.coefficients)),
                             pi = list(pi1=prob.MX1, pi2=prob.MX2), family=list(fam1="GA", fam2="GA"))
  }
  if(Marginal_Dist2=="GamMix(3)"){
    prob.MX1 <- round(marginal_non_con2$prob[1],3)
    prob.MX2 <- round(marginal_non_con2$prob[2],3)
    prob.MX3 <- 1 - prob.MX1 - prob.MX2
    cop.sample2.non.con<-qMX(sample[,con1], mu=list(mu1=exp(marginal_non_con2$models[[1]]$mu.coefficients), mu2=exp(marginal_non_con2$models[[2]]$mu.coefficients), mu3=exp(marginal_non_con2$models[[3]]$mu.coefficients)),
                             sigma=list(sigma1=exp(marginal_non_con2$models[[1]]$sigma.coefficients), sigma2=exp(marginal_non_con2$models[[2]]$sigma.coefficients), sigma3=exp(marginal_non_con2$models[[3]]$sigma.coefficients)),
                             pi = list(pi1=prob.MX1, pi2=prob.MX2, pi3=prob.MX3), family=list(fam1="GA", fam2="GA", fam3="GA"))
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
    cop.sample1.non.con <- qtruncnorm(sample[, con1], a = min(Data_Con2[,
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

  ##Calculating the joint probability
  #Joint probability from the model conditioned on Con1
  RP_Con1<-(EL_Con1/(1-(1-EL_Con1/RP_Var1)-(1-EL_Con1/RP_Var2)+BiCopCDF(1-EL_Con1/RP_Var1, 1-EL_Con1/RP_Var2, obj1)))
  #Joint probability from the model conditioned on Con2
  RP_Con2<-(EL_Con2/(1-(1-EL_Con2/RP_Var1)-(1-EL_Con2/RP_Var2)+BiCopCDF(1-EL_Con2/RP_Var1, 1-EL_Con2/RP_Var2, obj2)))
  #Joint return period will be maximum of the retunr periods from the two models
  RP_Copula<-max(RP_Con1,RP_Con2)

  #Plotting the results so far
  par(mfrow = c(2, 2))
  par(mar = c(4.5, 4.2, 0.5, 0.5))
  plot(Data[, con1], Data[, con2], xlim = c(x_min, x_max),
       ylim = c(y_min, y_max), col = "Light Grey", xlab = x_lab,
       ylab = y_lab, cex.lab = 1.5, cex.axis = 1.5)
  points(Data_Con1[, con1], Data_Con1[, con2], col = 4, cex = 1.5)
  points(Data_Con2[, con1], Data_Con2[, con2], col = "Red",
         pch = 4, cex = 1.5)
  points(Var1, Var2, pch = 16, cex = 1.5)
  legend("topright", c(paste("Full dependence RP = ", min(RP_Var1,RP_Var2), " years", sep = ""), paste("Joint RP = ", round(RP_Copula,0), " years", sep = ""), paste("Independence RP = ", RP_Var1 * RP_Var2, " years", sep = "")), bty = "n", cex = 1.25)
  segments(Var1,0,Var1,Var2,lty=2)
  axis(1,Var1,labels=paste(round(Var1,2)),line=1.2)
  segments(0,Var2,Var1,Var2,lty=2)
  axis(2,Var2,labels=paste(round(Var2,2)),line=1.2)

  ##Calculating the conditional probabilities using a simulation approach
  if(con_var==con1){
    #plot(cop.sample[, con1],
    #     cop.sample[, con2],xlim = c(min(cop.sample[, con1]), max(cop.sample[, con1])),
    #     ylim = c(min(cop.sample[, con2]), max(cop.sample[, con2])), col = "Light Grey", xlab = x_lab,
    #     ylab = y_lab, cex.lab = 1.5, cex.axis = 1.5)
    #rect(Var1-Width,min(cop.sample[, con2],na.omit(Data[,con2])),Var1+Width,max(cop.sample[, con2],na.omit(Data[,con2])),col="Light grey")
    #points(cop.sample[which(cop.sample[,con1]>(Var1-Width) & cop.sample[,con1]<(Var1+Width)),con1],
    #       cop.sample[which(cop.sample[,con1]>(Var1-Width) & cop.sample[,con1]<(Var1+Width)),con2],col=1,pch=16)
    #Rate
    rate<-length(which(cop.sample[, con1] > Var1))/N
    #Histogram of the values of the non-conditioned variable when the conditioning variable is in the interval around Var1 i.e. [Var1-width,Var1+width]
    hist(cop.sample[which(cop.sample[,con1]>(Var1-Width) & cop.sample[,con1]<(Var1+Width)),con2],
         cex.lab = 1.5, cex.axis = 1.5,
         freq=F,main="",xlab = y_lab,col="Grey")
    #Cummulative distribution function (CDF) of the non-conditioned variable when the conditioning variable is approximately Var1 i.e. [Var1-width,Var1+width]
    CDF_Var<-approx(seq(1,length(which(cop.sample[,con1]>(Var1-Width) & cop.sample[,con1]<(Var1+Width))),1)/length(which(cop.sample[,con1]>(Var1-Width) & cop.sample[,con1]<(Var1+Width))),
                    cop.sample[which(cop.sample[,con1]>(Var1-Width) & cop.sample[,con1]<(Var1+Width)),con2][order(cop.sample[which(cop.sample[,con1]>(Var1-Width) & cop.sample[,con1]<(Var1+Width)),con2])],
                    seq(round(min(1/length(which(cop.sample[,con1]>(Var1-Width) & cop.sample[,con1]<(Var1+Width)))+0.0005),3),1,0.001))
    #Plotting the conditional CDF
    plot(CDF_Var$y,CDF_Var$x,xlab=y_lab,ylab="Conditional CDF",xlim = c(y_min, y_max),cex.lab = 1.5, cex.axis = 1.5,type='l',lwd=2.5)
    #Probability that the non-conditioned variable is less than Var2 given the conditioned variable has a return period of RP_Con
    Con_Prob_Est<-approx(CDF_Var$y,CDF_Var$x,Var2)$y
    CDF_Var<-approx(CDF_Var$y,CDF_Var$x,seq(round(min(CDF_Var$y),2),round(max(CDF_Var$y),2),0.01))
    #Number of realizations where the conditioned variable is in [Var1-width,Var1+width]
    N_Sub_Sample<-length(which(cop.sample[,con1]>(Var1-Width) & cop.sample[,con1]<(Var1+Width)))
  }

  if(con_var==con2){
    #plot(cop.sample[, con1], cop.sample[, con2], xlim = c(min(cop.sample[, con1]), max(cop.sample[, con1])),
    #     ylim = c(min(cop.sample[, con2]), max(cop.sample[, con2])), col = "Light Grey", xlab = x_lab,
    #     ylab = y_lab, cex.lab = 1.5, cex.axis = 1.5)
    #rect(Var2-Width,min(cop.sample[, con1],na.omit(Data[,con1])),Var2+Width,max(cop.sample[, con1],na.omit(Data[,con1])),col="Light grey")
    #points(cop.sample[which(cop.sample[,con2]>(Var2-Width) & cop.sample[,con2]<(Var2+Width)),con1],
    #       cop.sample[which(cop.sample[,con2]>(Var2-Width) & cop.sample[,con2]<(Var2+Width)),con2],col=1,pch=16)
    #box()
    #Rate
    rate<-length(which(cop.sample[, con2] > Var2))/N
    #Histogram of the values of the non-conditioned variable when the conditioning variable is in the interval around Var2 i.e. [Var2-width,Var2+width]
    hist(cop.sample[which(cop.sample[,con2]>(Var2-Width) & cop.sample[,con2]<(Var2+Width)),con1],
         cex.lab = 1.5, cex.axis = 1.5,freq=F,main="", xlab = y_lab,col="Grey")
    #Cummulative distribution function (CDF) of the non-conditioned variable when the conditioning variable is approximately Var2 i.e. [Var2-width,Var2+width]
    CDF_Var<-approx(seq(1,length(which(cop.sample[,con2]>(Var2-Width) & cop.sample[,con2]<(Var2+Width))),1)/length(which(cop.sample[,con2]>(Var2-Width) & cop.sample[,con2]<(Var2+Width))),
                    cop.sample[which(cop.sample[,con2]>(Var2-Width) & cop.sample[,con2]<(Var2+Width)),con1][order(cop.sample[which(cop.sample[,con2]>(Var2-Width) & cop.sample[,con2]<(Var2+Width)),con1])],
                    seq(round(min(1/length(which(cop.sample[,con2]>(Var2-Width) & cop.sample[,con2]<(Var2+Width)))+0.0005),3),1,0.001))
    #Plotting the conditional CDF
    plot(CDF_Var$y,CDF_Var$x,xlab=x_lab,ylab="Conditional CDF",xlim = c(x_min, x_max),cex.lab = 1.5, cex.axis = 1.5,type='l',lwd=2.5)
    #Probability that the non-conditioned variable is less than Var1 given the conditioned variable has a return period of RP_Con
    Con_Prob_Est<-approx(CDF_Var$y,CDF_Var$x,Var1)$y
    CDF_Var<-approx(CDF_Var$y,CDF_Var$x,seq(round(min(CDF_Var$y),2),round(max(CDF_Var$y),2),0.01))
    #Number of realizations where the conditioned variable is in [Var2-width,Var2+width]
    N_Sub_Sample<-length(which(cop.sample[,con2]>(Var2-Width) & cop.sample[,con2]<(Var2+Width)))
  }

  #Create a list of outputs.
  res<-list(Con_Var=Con_Var,
            RP_Var1=RP_Var1,RP_Var2=RP_Var2,
            Var1=Var1,Var2=Var2,
            RP_Full_Dependence=min(RP_Var1,RP_Var2),
            RP_Independence=RP_Var1*RP_Var2,
            RP_Copula = RP_Copula,
            Prob = 1/RP_Copula,
            N_Sub_Sample=N_Sub_Sample,
            Non_Con_Var_X=CDF_Var$x,Con_Prob=CDF_Var$y,
            Con_RP=1/(rate*(1-CDF_Var$y)),
            Con_Prob_Est=Con_Prob_Est)
  return(res)
}
