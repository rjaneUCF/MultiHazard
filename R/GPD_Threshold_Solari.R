#' Solari et al (2017) automatic GPD threshold selection
#'
#' Automatic threshold selection method in Solari et al. (2017) is implemented to find the threshold above which excesses are follow a GPD. The code is based on the \code{ANALISIS_POT_LNORM} function provided by Sebastian Solari.
#'
#' @param Event Numeric vector containing the declustered events.
#' @param Data Original time series. Dataframe containing two columns. In column: \itemize{
#' \item \code{1} A \code{"Date"} object of equally spaced discrete time steps.
#' \item \code{2} Numeric vector containing corresponding time series values.
#' }
#' @param RPs Numeric vector specifying the return levels calculated from the GPD fits over the thresholds. Default is \code{c(50,100,500,100)} plus the return period associated with the minimum candidate threshold.
#' @param RPs_PLOT Numeric vector of length three specifying which elements of \code{RPs} are plotted in the middle row of the graphical output. Default is \code{c(1,2,3)}.
#' @param Min_Quantile Numeric vector of length one specifying the minimum threshold, expressed as a quantile of the original time series (2nd column of \code{Data}) to be tested. Default \code{0.95}.
#' @param Alpha Numeric vector of length one specifying the level of confidence associated with the confidence interval i.e., the probability that the interval contains the true value of the parameter is \eqn{1-\frac{Alpha}{2}}. The interval is referred to as the \eqn{100(1-\frac{Alpha}{2})\%} confidence interval. Default is \code{0.1}.
#' @param mu (average) occurrence frequency of events in the original time series \code{Data}. Numeric vector of length one. Default is \code{365.25}, daily data.
#' @param N_Sim Numeric vector of length one specifying the number of bootstrap samples. Default is \code{10}.
#' @section Details:
#' EDF-statistics are goodness-of-fit statistics based on a comparison of the Empirical Distribution Function (EDF) \eqn{F_n} and a candidate parametric probability distribution \eqn{F} Stephens et al. (1974). Quadratic EDF test measure the distance between \eqn{F} and \eqn{F_n} by:\deqn{n\int^{\infty}_{-\infty}=(F(x)-F_n(x))^2w(x) dx}
#' where \eqn{n} is the number of elements in the original sample and \eqn{w(x)} is a weighting function.
#' In the Cramer Von Misses statistic \eqn{w(x)=1}, whereas the Anderson-Darling statistic \eqn{A^2}, assigns more weight to the tails of the data by setting \eqn{w(x)=\frac{1}{F(x)(1-F(x))}}. Under the null hypothesis that the sample \eqn{x_1,\dots,x_n} is from a GPD, the transformation \eqn{z=F_1(x)} a sample z uniformly distribution between \eqn{0} and \eqn{1}.  \deqn{A^2=-\frac{1}{n}\sum_{i=1}^{n} \{(2i-1)[log(z_i)+log(1-z_{n+z-i})]\}-n}
#' Sinclair et al. (1990) proposed the right-tail weighted Anderson Darling statistic \eqn{A_R^2} which allocates more weight to the upper tail and less to the lower tail of the distribution than \eqn{A^2} and is given by:
#' \deqn{A_R^2=\frac{n}{2}\sum_{i=1}^{n} \left [2-\frac{(2i-1)}{n}log(1-z_i)+2z_{i} \right]}
#'
#' Solari et al. (2017) formalized EDF statistic - GOF test threshold selection procedures used to test the null hypothesis that a sample is from a GPD distribution. creating an automated approach adopting the \eqn{A_R^2} as the EDF statistic. The authors also proposed combining the approach with a bootstrapping technique to assess the influence of threshold on the uncertainly of higher return period quantiles. The approach in Solari et al. (2017) comprises the following steps: \enumerate{
#'
#' \item Decluster the time series to produce a series of \eqn{n_p} independent cluster maxima \eqn{\{x_i:i=1,\dots,n_p\}} and sort such that \eqn{\{x_1\leq\dots\leq x_p\}}.
#'
#' \item The sorted series defines a series of \eqn{n_u} thresholds after excluding repeated values i.e., \eqn{n_u \leq n_p}.  For each threshold \eqn{\{u_j, j=1,\dots,n_u\}} fit the GPD via L-Moments using only the excesses satisfying \eqn{x>u_j}. Then, calculate the R-AD statistic and its associated p-value for each threshold.
#'
#' \item Select the threshold that minimizes one minus the p-value i.e.,
#' \deqn{u_0 = argmin_{u_j} (1-p(u_j)).}
#' }
#' @return List comprising \itemize{
#' \item \code{Thres_Candidate}
#' Thresholds tested which are the cluster maxima in \code{Events} exceeding the \code{Min_Quantile} quantile of the original time series (given in column 2 of \code{Data}).
#' \item \code{GPD_MLE}
#' GPD parameter estimates, Mean Residual Life Plot (MRLP) values and return level estimates associated with each \code{Thres_Candidate}.
#' \item \code{CI_Upper}
#' Upper limits of the confidence interval for the point estimates of the corresponding element of \code{GPD_MLE}.
#' \item \code{CI_Lower}
#' Lower limits of the confidence interval for the point estimates of the corresponding element of \code{GPD_MLE}.
#' \item \code{AR2}
#' Value of the right-tail weighted Anderson Darling statistic \eqn{A_R^2}, the test statistic used in the Solari et al. (2017) method for each \code{Thres_Candidate}.
#' \item \code{AR2_pValue}
#' p-value asssociated with \eqn{A_R^2}.
#' }
#' To interpret the graphical output. Top row: The GPD exhibits certain threshold stability properties. The guiding principle for threshold choice is to find the lowest value of the threshold such that the parameter estimates stabilize to a constant value which is sustained at all higher thresholds, once the sample uncertainty has been accounted for (typically assessed by pointwise uncertainty intervals). Mean residual life plot (left). If the GPD is a valid model for excesses above a threshold then the mean of these excesses will be a linear function of the threshold. We therefore select the lowest threshold where there is a linear trend in the mean residual life plot. Parameter stability plots for the shape (center) and scale (right) parameters. If the GPD is a suitable model for a threshold then for all higher thresholds it will also be suitable, with the shape and scale parameters being constant. The lowest threshold - to reduce the associated uncertainty - at which the parameter estimates are stable for all higher thresholds should be selected.
#' Middle row: Return levels estimated from the GPD fitted at various thresholds.
#' Lower row: Right-tail weighted Anderson Darling statistic \eqn{A_R^2} associated with the GPD fitted using various thresholds. Lower \eqn{A_R^2} statistic values signify less (quadratic) distance between the empirical distribution and the GPD i.e., GPD is a better fit for these thresholds (left).  \eqn{1-p_{value}} associated with the \eqn{A_R^2} for each threshold. The \eqn{A_R^2} goodness of fit tests, tests the null hypothesis that the observations are from a GPD. At smaller \eqn{1-p_{value}} figure there is less chance of rejecting the null hypothesis i.e., the GPD is more suitable at these thresholds (center). Events per year at each threshold (right).
#' @export
#' @examples
#' #Declustering the rainfall at site S22 using a 7-day window.
#' Rainfall_Declust_SW<-Decluster_SW(Data=S22.Detrend.df[,c(1:2)],Window_Width=7)
#' #Finding an appropriate threshold for the declustered series
#' GPD_Threshold_Solari(Event=Rainfall_Declust_SW$Declustered,
#'                      Data=22.Detrend.df[,2])
GPD_Threshold_Solari<-function(Event,Data,RPs=c(10,50,100,500,1000),RPs_PLOT=c(2,3,4),Min_Quantile=0.95,Alpha=0.1,mu=365.25,N_Sim=10){
  
  # Loads the p-values matrix
  p_val <- array(c(PVAL_AU2_LMOM_1$PVAL,PVAL_AU2_LMOM_2$PVAL,PVAL_AU2_LMOM_3$PVAL,PVAL_AU2_LMOM_4$PVAL),dim=c(24,6,100001))
  
  # Removing NAs in Data
  Data = na.omit(Data)
  
  # Auxiliary variables
  N_Years = length(Data)/mu
  
  # POT with L-moments and with bootstrapping confidence intervals
  Event = na.omit(Event)
  Event = sort(Event[Event>quantile(Data,Min_Quantile)])
  
  # initialize variables
  u_Candidate         = sort(Event)
  u_Candidate         = u_Candidate[1:(length(u_Candidate)-30)]
  u_Candidate         = u_Candidate[diff(u_Candidate)>0]
  N_u_Candidate      = length(u_Candidate)
  GPD.MLE       = array(0,dim=c(N_u_Candidate,6+length(RPs)))
  INCONS_MLE    = array(0,dim=c(N_u_Candidate,1))
  CI.Upper      = array(0,dim=c(N_u_Candidate,6+length(RPs)))
  CI.Lower      = array(0,dim=c(N_u_Candidate,6+length(RPs)))
  AR2           = rep(NA,N_u_Candidate)
  AR2.pValue    = rep(NA,N_u_Candidate)
  
  # Analysis
  for(i in 1:N_u_Candidate){
    GPD.MLE[i,] = unlist(GPD_Eval_MLE(Event[Event>u_Candidate[i]],RPs,N_Years))
    
    # Bootstraping
    BOOT = array(0,dim=c(N_Sim,length(GPD.MLE[i,])))
    for(J in 1:N_Sim){
      while( BOOT[J,1]==0){
        try({
          BOOT[J,] = unlist(GPD_Eval_MLE(sample(Event[Event>u_Candidate[i]],length(Event[Event>u_Candidate[i]]),replace=T),RPs,N_Years))
        }, silent = TRUE)
      }
    }
    
    CI.Upper[i,] =  apply(BOOT, 2,  quantile, 1-Alpha/2)
    CI.Lower[i,] =  apply(BOOT, 2,  quantile, Alpha/2)
    
    # Anderson-Darling Upper without confidence intervals for L-Moments
    AR2[i] = AR2(GPD.MLE[i,1:3],Event[Event>u_Candidate[i]])
    x_AR2 <- c(10,12,14,16,18,20,22,24,26,28,30,35,40,45,50,60,70,80,90,100,200,300,400,500)
    y_AR2 <- c(-0.5,-0.3,-0.1,0.1,0.3,0.5)
    x <- max(10.1,min(499.1,length(Event[Event>u_Candidate[i]])))
    y <- sign(GPD.MLE[i,1])*min(abs(GPD.MLE[i,1]),0.4999)
    m <- max(1,min((1+(round(AR2[i],4)*10^4)),100001))
    #AR2.pValue[i] = interp(x,y,p_val[,,m],
    #                       sign(GPD.MLE[i,1])*min(abs(GPD.MLE[i,1]),0.5),
    #                       max(10,min(500,length(Event[Event>u_Candidate[i]]))))$z
    x.min <- max(x_AR2[x_AR2<=x])
    x.max <- min(x_AR2[x_AR2>x])
    y.min <- max(y_AR2[y_AR2<=y])
    y.max <- min(y_AR2[y_AR2>y])
    f.x.ymin <- ((x-x.min)/(x.max-x.min)) *  p_val[which(x_AR2==x.max),which(y_AR2==y.min),m] + ((x.max-x)/(x.max-x.min)) * p_val[which(x_AR2==x.min),which(y_AR2==y.min),m]
    f.x.ymax <- ((x-x.min)/(x.max-x.min)) *  p_val[which(x_AR2==x.max),which(y_AR2==y.max),m] + ((x.max-x)/(x.max-x.min)) * p_val[which(x_AR2==x.min),which(y_AR2==y.max),m]
    AR2.pValue[i] <- ((y-y.min)/(y.max-y.min)) *  f.x.ymax + ((y.max-y)/(y.max-y.min)) * f.x.ymin
  }
  colnames(GPD.MLE)<-c("xi","sigma","u","MRLP","mod_sigma","rate",paste(RPs))
  
  RPs<-as.character(RPs)
  z<-which(GPD.MLE[,1]>-0.5 & GPD.MLE[,1]<0.5 & GPD.MLE[,6]<(length(na.omit(Data)/365.25)))
  par(mfrow=c(3,3))
  par(mar=c(4.2,4.5,0.5,0.5))
  #MARLP
  plot(u_Candidate[z],GPD.MLE[z,4],type="l",ylim=c(0,max(CI.Upper[z,4])),xlab="Threshold",ylab="Mean residual life plot")
  lines(u_Candidate[z],CI.Upper[z,4],lty=5)
  lines(u_Candidate[z],CI.Lower[z,4],lty=5)
  #Shape
  plot(u_Candidate[z],GPD.MLE[z,1],type="l",ylim=c(-0.5,0.5),xlab="Threshold",ylab="GP shape parameter")
  lines(u_Candidate[z],CI.Upper[z,1],lty=5)
  lines(u_Candidate[z],CI.Lower[z,1],lty=5)
  #Modified Scale
  plot(u_Candidate[z],GPD.MLE[z,5],type="l",ylim=c(min(CI.Lower[z,5]),max(CI.Upper[z,5])),xlab="Threshold",ylab="Modified GP scale parameter")
  lines(u_Candidate[z],CI.Upper[z,5],lty=5)
  lines(u_Candidate[z],CI.Lower[z,5],lty=5)
  #1st return level plot
  plot(u_Candidate[z],GPD.MLE[z,which(colnames(GPD.MLE)==RPs[RPs_PLOT[1]])],type="l",ylim=c(0,max(CI.Upper[z,which(colnames(GPD.MLE)==RPs[RPs_PLOT[1]])])),xlab="Threshold",ylab=paste(RPs[RPs_PLOT[1]],'year return period'))
  lines(u_Candidate[z],CI.Upper[z,which(colnames(GPD.MLE)==RPs[RPs_PLOT[1]])],lty=5)
  lines(u_Candidate[z],CI.Lower[z,which(colnames(GPD.MLE)==RPs[RPs_PLOT[1]])],lty=5)
  #2nd return level plot
  plot(u_Candidate[z],GPD.MLE[z,which(colnames(GPD.MLE)==RPs[RPs_PLOT[2]])],type="l",ylim=c(0,max(CI.Upper[z,which(colnames(GPD.MLE)==RPs[RPs_PLOT[2]])])),xlab="Threshold",ylab=paste(RPs[RPs_PLOT[2]],'year return period'))
  lines(u_Candidate[z],CI.Upper[z,which(colnames(GPD.MLE)==RPs[RPs_PLOT[2]])],lty=5)
  lines(u_Candidate[z],CI.Lower[z,which(colnames(GPD.MLE)==RPs[RPs_PLOT[2]])],lty=5)
  #3rd return level plot
  plot(u_Candidate[z],GPD.MLE[z,which(colnames(GPD.MLE)==RPs[RPs_PLOT[3]])],type="l",ylim=c(0,max(CI.Upper[z,which(colnames(GPD.MLE)==RPs[RPs_PLOT[3]])])),xlab="Threshold",ylab=paste(RPs[RPs_PLOT[3]],'year return period'))
  lines(u_Candidate[z],CI.Upper[z,which(colnames(GPD.MLE)==RPs[RPs_PLOT[3]])],lty=5)
  lines(u_Candidate[z],CI.Lower[z,which(colnames(GPD.MLE)==RPs[RPs_PLOT[3]])],lty=5)
  #AR2
  plot(u_Candidate[z],AR2[z],type="l",ylim=c(min(AR2[z]),max(AR2[z])),xlab="Threshold",ylab=expression('AR'^2*' statistic'))
  #1-p_Value
  plot(u_Candidate[z],AR2.pValue[z],type="l",ylim=c(min(AR2.pValue[z]),max(AR2.pValue[z])),xlab="Threshold",ylab=expression('1-p'[value]))
  #Events per year
  plot(u_Candidate[z],GPD.MLE[z,6],type="l",ylim=c(0,max(GPD.MLE[z,6])),xlab="Threshold",ylab="Events per year")
  par(mfrow=c(1,1))
  ecdf_fun <- function(x,perc) ecdf(x)(perc)
  u_Candidate_Quantile<-ecdf_fun(Data,Event)
  
  Candidate_Thres<-u_Candidate[z][which(AR2.pValue[z]==min(AR2.pValue[z]))]
  res<-list("Thres_Candidate"=u_Candidate,"Thres_Candidate_Quantile"=u_Candidate_Quantile,"GPD_MLE"=GPD.MLE,"CI_Upper"=CI.Upper,"CI_Lower"=CI.Lower,"AR2"=AR2,"AR2_pValue"=AR2.pValue,"Candidate_Thres"=Candidate_Thres)
  return(res)
}