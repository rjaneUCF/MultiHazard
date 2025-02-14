#' Evaluates the goodness of fit of the return curve estimates
#' 
#' The procedure calculates the empirical probability of observing data within the survival regions defined by a subset of points on the return curve. If the curve is a good fit, the empirical probabilities should closely match the probabilities associated with the return level curve.  The procedure which is introduced in Murphy-Barltrop et al. (2023) uses bootstrap resampling of the original data set to obtain confidence intervals for the empirical estimates.
#'
#' @param data Data frame of raw data detrended if necessary. First column should be of class \code{Date}.
#' @param q Numeric vector of length one specifying quantile level for fitting GPDs and the HT04 and WT13 models.
#' @param rp Numeric vector of length one specifying return period of interest.
#' @param mu Numeric vector of length one specifying the (average) occurrence frequency of events in Data. Default is 365.25, daily data.
#' @param n_sim Numeric vector of length one specifying the number of simulations for HT model. Default is \code{50}.
#' @param n_grad Numeric vector of length one specifying number of number of rays along which to compute points on the curve. Default is \code{50}. 
#' @param boot_method Character vector of length one specifying the bootstrap method. Options are \code{"basic"} (default), \code{"block"} or \code{"monthly"}. 
#' @param boot_replace Character vector of length one specifying whether simple bootstrapping is carried out with \code{"T"} or without \code{"F"} replacement. Only required if \code{boot_method = "basic"}. Default is \code{NA}. 
#' @param block_length Numeric vector of length one specifying block length. Only required if \code{boot_method = "block"}. Default is \code{NA}. 
#' @param boot_prop Numeric vector of length one specifying the minimum proportion of non-missing values of at least of the variables for a month to be included in the bootstrap. Only required if \code{boot_method = "monthly"}. Default is \code{0.8}.
#' @param n_boot Numeric vector of length one specifying number of bootstrap samples. Default is \code{100}. 
#' @param decl_method_x Character vector of length one specifying the declustering method to apply to the first variable. Options are the storm window approach \code{"window"} (default) and the runs method \code{"runs"}. 
#' @param decl_method_y Character vector of length one specifying the declustering method to apply to the second variable. Options are the storm window approach \code{"window"} (default) and the runs method \code{"runs"}. 
#' @param window_length_x Numeric vector of length one specifying the storm window length to apply during the declustering of the  first variable if \code{decl_method_x = "window"}. 
#' @param window_length_y Numeric vector of length one specifying the storm window length to apply during the declustering of the second variable if \code{decl_method_y = "window"}. 
#' @param u_x Numeric vector of length one specifying the threshold to adopt in the declustering of the first variable if \code{decl_method_x = "runs"}. Default is \code{NA}.
#' @param u_y Numeric vector of length one specifying the threshold to adopt in the declustering of the second variable if \code{decl_method_y = "runs"}. Default is \code{NA}.
#' @param sep_crit_x Numeric vector of length one specifying the separation criterion to apply during the declustering of the first variable if \code{decl_method_x = "runs"}. Default is \code{NA}.
#' @param sep_crit_y Numeric vector of length one specifying the separation criterion to apply during the declustering of the second variable if \code{decl_method_y = "runs"}. Default is \code{NA}.
#' @param boot_method_all Character vector of length one specifying the bootstrapping procedure to use when estimating the distribution of empirical (survival) probabilities from the original dataset (without any declustering). Options are \code{"basic"} (default) and \code{"block"}.
#' @param boot_replace_all Character vector of length one specifying whether bootstrapping of original dataset (without any declustering) when estimating the distribution of empirical (survival) probabilities is carried out with \code{"T"} or without \code{"F"} replacement. Only required if \code{boot_method_all = "basic"}. Default is \code{NA}. 
#' @param block_length_all Numeric vector of length one specifying block length. Only required if \code{boot_method_all = "block"}. Default is \code{14}.
#' @param alpha Numeric vector of length one specifying the \code{100(1-alpha)\%} confidence interval. Default is \code{0.1}.
#' @param x_lab Character vector specifying the x-axis label.
#' @param y_lab Character vector specifying the y-axis label.
#' @param x_lim_min Numeric vector of length one specifying x-axis minimum. Default is \code{NA}.
#' @param x_lim_max Numeric vector of length one specifying x-axis maximum. Default is \code{NA}.
#' @param y_lim_min Numeric vector of length one specifying y-axis minimum. Default is \code{NA}.
#' @param y_lim_max Numeric vector of length one specifying y-axis maximum. Default is \code{NA}.
#' @section Details:
#' The HT04 model is fit to two conditional samples. One sample comprises the declustered time series of the first variable paired with concurrent values of the other variable. The second sample is obtained in the same way but with the variables reversed. The empirical probabilities are calculated using these two conditional samples and the original dataset (without any declustering).
#' The return period should be chosen to ensure there is sufficient data for estimating empirical probabilities, yet the curve is sufficiently 'extreme'. An example could be to consider the fit using the 1 year return period curve rather than the 100 year return period curve. 
#' @return List comprising the angles \code{"ang_ind"} associated with the points on the curve for which the empirical probability estimates were calculated. 
#' For the HT04 model: Median {"med_x_ht04"}, lower \code{"lb_x_ht04"} and upper \code{"ub_x_ht04"} bounds associated with the probabilities calculated using the sample conditioned on the first variable.  
#' Median {"med_y_ht04"}, lower \code{"lb_y_ht04"} and upper \code{"ub_y_ht04"} bounds associated with the probabilities calculated using the sample conditioned on the second variable.  
#' Median {"med_ht04"}, lower \code{"lb_ht04"} and upper \code{"ub_ht04"} bounds associated with the original dataset (without any declustering).
#' 
#' For the WT13 model: Median {"med_x_wt13"}, lower \code{"lb_x_wt13"} and upper \code{"ub_x_wt13"} bounds associated with the probabilities calculated using the sample conditioned on the first variable. 
#' Median {"med_y_wt13"}, lower \code{"lb_y_wt13"} and upper \code{"ub_y_wt13"} bounds associated with the probabilities calculated using the sample conditioned on the second variable.  
#' Median {"med_wt13"}, lower \code{"lb_wt13"} and upper \code{"ub_wt13"} bounds associated with the original dataset (without any declustering).
#' @export
#' @examples
#' #' #Data starts on first day of 1948
#' head(S22.Detrend.df)
#' 
#' #Dataframe ends on 1948-02-03 
#' tail(S22.Detrend.df)
#'
#' #Adding dates to complete final month of combined records
#' final.month = data.frame(seq(as.Date("2019-02-04"),as.Date("2019-02-28"),by="day"),NA,NA,NA)
#' colnames(final.month) = c("Date","Rainfall","OsWL","Groundwater")
#' S22.Detrend.df.extended = rbind(S22.Detrend.df,final.month)
#' #Derive return curves
#' rc_diag(data=S22.Detrend.df.extended[,1:3],q=0.985,rp=1,mu=365.25,n_sim=100,
#'         n_grad=50,n_boot=100,boot_method="monthly", 
#'         boot_replace=NA, block_length=NA, boot_prop=0.8, 
#'         decl_method_x="runs", decl_method_y="runs", 
#'         window_length_x=NA,window_length_y=NA, 
#'         u_x=0.95, u_y=0.95, 
#'         sep_crit_x=36, sep_crit_y=36, 
#'         alpha=0.1, x_lab=NA, y_lab=NA,
#'         boot_method_all="block", boot_replace_all=NA, block_length_all=14)
rc_diag = function(data,q,rp,mu,n_sim,n_grad,n_boot,boot_method, boot_replace, block_length, boot_prop, decl_method_x, decl_method_y, window_length_x, window_length_y, u_x=NA, u_y=NA, sep_crit_x=NA, sep_crit_y=NA, boot_method_all="block", boot_replace_all=NA, block_length_all=14, boot_prop_all=0.8,alpha=0.1, x_lab=NA, y_lab=NA,x_lim_min=min(data_df[,2],na.rm=T),x_lim_max=max(data_df[,2],na.rm=T)+0.3*diff(range(data[,2],na.rm=T)),y_lim_min=min(data[,3],na.rm=T),y_lim_max=max(data[,3],na.rm=T)+0.3*diff(range(data[,2],na.rm=T))){

  #Convert return period (rp) to probability
  p = (1/mu) / rp
  
  curve = rc_est(data=data,q=q,rp=rp,mu=mu,n_sim=n_sim,
                 n_grad=n_grad,n_boot=n_boot,boot_method=boot_method, 
                 boot_replace=boot_replace, block_length=boot_length, boot_prop=boot_prop, 
                 decl_method_x=decl_method_x, decl_method_y=decl_method_y, 
                 window_length_x=window_length_x,window_length_y=window_length_y, 
                 u_x=u_x, u_y=u_y, 
                 sep_crit_x=sep_crit_x, sep_crit_y=sep_crit_y, 
                 alpha=alpha, x_lab=NA, y_lab=NA,
                 plot=F)
 
  median_diag = curve$median_ht04
  median_diag2 = curve$median_wt13
  
#Computing diagnostic for median curves 

emp_prob_x_1 <- lapply(1:n_grad, function(i) vector())
emp_prob_y_1 <- lapply(1:n_grad, function(i) vector())

emp_prob_x_2 <- lapply(1:n_grad, function(i) vector())
emp_prob_y_2 <- lapply(1:n_grad, function(i) vector())

for(j in 1:n_boot){
  
  #Bootstrap data using user specified method
  if(boot_method=="basic"){
    index = sample(1:nrow(data),replace=boot_replace)
    boot_df = data[index,]
  }
  
  if(boot_method=="block"){
    boot_df = bootstrap_block(data,k=block_length)
  }
  
  if(boot_method=="monthly"){
    boot_df = bootstrap_month(data,boot_prop)
  }
  
  #Decluster bootstrapped series
  #Variable 1
  if(decl_method_x == "window"){
    decl = Decluster_SW(Data=boot_df[,1:2], Window_Width = window_length_x)
    boot_x_dec_df = data.frame(as.Date(data[,1]),decl$Declustered)
  }
  
  if(decl_method_x == "runs"){
    decl = Decluster(Data=boot_df[,2], u=u_x, SepCrit = sep_crit_x, mu = mu)
    boot_x_dec_df = data.frame(as.Date(data[,1]),decl$Declustered)
  }
  
  colnames(boot_x_dec_df) = colnames(data[,1:2])
  boot_x_dec_df$Date = as.Date(boot_x_dec_df$Date)
  
  #Variable 2
  if(decl_method_y == "window"){
    decl = Decluster_SW(Data=boot_df[,c(1,3)], Window_Width = window_length_y)
    boot_y_dec_df = data.frame(as.Date(data[,1]),decl$Declustered)
  }
  
  if(decl_method_y == "runs"){
    decl = Decluster(Data=boot_df[,3], u = u_y, SepCrit = sep_crit_y, mu = mu)
    boot_y_dec_df = data.frame(as.Date(data[,1]),decl$Declustered)
  }
  
  colnames(boot_y_dec_df) = colnames(data[,c(1,3)])
  boot_y_dec_df$Date = as.Date(boot_y_dec_df$Date)
  
  #Dataframes
  df = data.frame("Date"=seq.Date(as.Date(min(data$Date,data$Date)),as.Date(max(data$Date,data$Date)),by="day"))
  df_1 = left_join(df,boot_x_dec_df,by="Date")
  boot_decl_df = left_join(df_1,boot_y_dec_df,by="Date")
  colnames(boot_decl_df) = colnames(data)
  
  #Dataframes to which HT04 model is applied 
  data_x = data.frame(boot_decl_df[,2],boot_df[,3])
  data_y = data.frame(boot_df[,2],boot_decl_df[,3])
  
  colnames(data_x) = colnames(data[,2:3])
  colnames(data_y) = colnames(data[,2:3])
  
  data_x = na.omit(data_x)
  data_y = na.omit(data_y)

  for(i in 1:n_grad){
    emp_prob_x_1[[i]][j] <- mean(data_x[, 1] > median_diag[i, 1] & data_x[, 2] > median_diag[i, 2])
    
    emp_prob_y_1[[i]][j] <- mean(data_y[, 1] > median_diag[i, 1] & data_y[, 2] > median_diag[i, 2])
    
    emp_prob_x_2[[i]][j] <- mean(data_x[, 1] > median_diag2[i, 1] & data_x[, 2] > median_diag2[i, 2])
    
    emp_prob_y_2[[i]][j] <- mean(data_y[, 1] > median_diag2[i, 1] & data_y[, 2] > median_diag2[i, 2])
  }

}

lb_x_1 <- sapply(1:n_grad, function(i) quantile(emp_prob_x_1[[i]], alpha/2))
ub_x_1 <- sapply(1:n_grad, function(i) quantile(emp_prob_x_1[[i]], 1 - alpha/2))
med_x_1 <- sapply(1:n_grad, function(i) quantile(emp_prob_x_1[[i]], 0.5))

lb_y_1 <- sapply(1:n_grad, function(i) quantile(emp_prob_y_1[[i]], alpha/2))
ub_y_1 <- sapply(1:n_grad, function(i) quantile(emp_prob_y_1[[i]], 1 - alpha/2))
med_y_1 <- sapply(1:n_grad, function(i) quantile(emp_prob_y_1[[i]], 0.5))

lb_x_2 <- sapply(1:n_grad, function(i) quantile(emp_prob_x_2[[i]], alpha/2))
ub_x_2 <- sapply(1:n_grad, function(i) quantile(emp_prob_x_2[[i]], 1 - alpha/2))
med_x_2 <- sapply(1:n_grad, function(i) quantile(emp_prob_x_2[[i]], 0.5))

lb_y_2 <- sapply(1:n_grad, function(i) quantile(emp_prob_y_2[[i]], alpha/2))
ub_y_2 <- sapply(1:n_grad, function(i) quantile(emp_prob_y_2[[i]], 1 - alpha/2))
med_y_2 <- sapply(1:n_grad, function(i) quantile(emp_prob_y_2[[i]], 0.5))

par(mfrow=c(1,2),mgp=c(2.2,1,0),mar=c(5,4,4,2)+0.1)

ang_ind = 1:n_grad

plot(ang_ind,med_x_1,xlab="Angle Index",ylab = "Probability",ylim=range(med_x_1,lb_x_1,ub_x_1,p),type='n',sub=paste('Decl.',colnames(data[2])),main="HT isoline",cex.lab=1.3, cex.axis=1.4,cex.main=1.7)
polygon(c(rev(ang_ind), ang_ind), c(rev( lb_x_1), ub_x_1), col = 'grey80', border = NA)
lines(ang_ind,med_x_1,type="l",col=1,lwd=2)
lines(ang_ind, ub_x_1, lty = 'dashed', col = 'blue',lwd=2)
lines(ang_ind, lb_x_1, lty = 'dashed', col = 'blue',lwd=2)
lines(ang_ind,rep(p,length(ang_ind)),type="l",lwd=3,col=2)
legend("topleft",legend = c("True Probability","Median Estimate","95% Confidence Intervals"),lty=c('solid','solid','dashed'),col=c(2,1,"blue"),cex=1.2,lwd=2,bg="white")

plot(ang_ind,med_y_1,xlab="Angle Index",ylab = "Probability",ylim=range(med_y_1,lb_y_1,ub_y_1,p),type='n',sub=paste('Decl.',colnames(data[2])),main="HT isoline",cex.lab=1.3, cex.axis=1.4,cex.main=1.7)
polygon(c(rev(ang_ind), ang_ind), c(rev( lb_y_1), ub_y_1), col = 'grey80', border = NA)
lines(ang_ind,med_y_1,type="l",col=1,lwd=2)
lines(ang_ind, ub_y_1, lty = 'dashed', col = 'blue',lwd=2)
lines(ang_ind, lb_y_1, lty = 'dashed', col = 'blue',lwd=2)
lines(ang_ind,rep(p,length(ang_ind)),type="l",lwd=3,col=2)
legend("topleft",legend = c("True Probability","Median Estimate","95% Confidence Intervals"),lty=c('solid','solid','dashed'),col=c(2,1,"blue"),cex=1.2,lwd=2,bg="white")

plot(ang_ind,med_x_2,xlab="Angle Index",ylab = "Probability",ylim=range(med_x_2,lb_x_2,ub_x_2,p),type='n',sub=paste('Decl.',colnames(data[3])),main="WT isoline",cex.lab=1.3, cex.axis=1.4,cex.main=1.7)
polygon(c(rev(ang_ind), ang_ind), c(rev( lb_x_2), ub_x_2), col = 'grey80', border = NA)
lines(ang_ind,med_x_2,type="l",col=1,lwd=2)
lines(ang_ind, ub_x_2, lty = 'dashed', col = 'blue',lwd=2)
lines(ang_ind, lb_x_2, lty = 'dashed', col = 'blue',lwd=2)
lines(ang_ind,rep(p,length(ang_ind)),type="l",lwd=3,col=2)
legend("topleft",legend = c("True Probability","Median Estimate","95% Confidence Intervals"),lty=c('solid','solid','dashed'),col=c(2,1,"blue"),cex=1.2,lwd=2,bg="white")

plot(ang_ind,med_y_2,xlab="Angle Index",ylab = "Probability",ylim=range(med_y_2,lb_y_2,ub_y_2,p),type='n',sub=paste('Decl.',colnames(data[3])),main="WT isoline",cex.lab=1.3, cex.axis=1.4,cex.main=1.7)
polygon(c(rev(ang_ind), ang_ind), c(rev( lb_y_2), ub_y_2), col = 'grey80', border = NA)
lines(ang_ind,med_y_2,type="l",col=1,lwd=2)
lines(ang_ind, ub_y_2, lty = 'dashed', col = 'blue',lwd=2)
lines(ang_ind, lb_y_2, lty = 'dashed', col = 'blue',lwd=2)
lines(ang_ind,rep(p,length(ang_ind)),type="l",lwd=3,col=2)
legend("topleft",legend = c("True Probability","Median Estimate","95% Confidence Intervals"),lty=c('solid','solid','dashed'),col=c(2,1,"blue"),cex=1.2,lwd=2,bg="white")


#Computing diagnostic for median curves 

emp_prob_1 <- lapply(1:n_grad, function(i) vector())

emp_prob_2 <- lapply(1:n_grad, function(i) vector())

for(j in 1:n_boot){
  
  boot_data = data[!(is.na(data[,2]) | is.na(data[,3])),c(2,3)]
  
  if(j > 1){
   if(boot_method_all=="basic"){
     index = sample(1:nrow(data),replace=boot_replace_all)
     boot_data = data[index,]
   }
   if(boot_method_all=="block"){
    boot_data = bootstrap_block(boot_data,k=block_length_all)  
   }
   if(boot_method_all=="monthly"){
    boot_data = bootstrap_month(boot_data,boot_prop_all)  
   }
  }
  
  
  for(i in 1:n_grad){
    emp_prob_1[[i]][j] <- mean(boot_data[, 1] > median_diag[i, 1] & boot_data[, 2] > median_diag[i, 2])
    
    emp_prob_2[[i]][j] <- mean(boot_data[, 1] > median_diag2[i, 1] & boot_data[, 2] > median_diag2[i, 2])
    
  }
  
}

lb_1 <- sapply(1:n_grad, function(i) quantile(emp_prob_1[[i]], alpha/2))
ub_1 <- sapply(1:n_grad, function(i) quantile(emp_prob_1[[i]], 1 - alpha/2))
med_1 <- sapply(1:n_grad, function(i) quantile(emp_prob_1[[i]], 0.5))

lb_2 <- sapply(1:n_grad, function(i) quantile(emp_prob_2[[i]], alpha/2))
ub_2 <- sapply(1:n_grad, function(i) quantile(emp_prob_2[[i]], 1 - alpha/2))
med_2 <- sapply(1:n_grad, function(i) quantile(emp_prob_2[[i]], 0.5))

par(mfrow=c(1,2),mgp=c(2.2,1,0),mar=c(5,4,4,2)+0.1)

ang_ind = 1:n_grad

plot(ang_ind,med_1,xlab="Angle Index",ylab = "Probability",ylim=range(med_1,lb_1,ub_1,p),type='n',main="HT isoline",cex.lab=1.3, cex.axis=1.4,cex.main=1.7)
polygon(c(rev(ang_ind), ang_ind), c(rev( lb_1), ub_1), col = 'grey80', border = NA)
lines(ang_ind,med_1,type="l",col=1,lwd=2)
lines(ang_ind, ub_1, lty = 'dashed', col = 'blue',lwd=2)
lines(ang_ind, lb_1, lty = 'dashed', col = 'blue',lwd=2)
lines(ang_ind,rep(p,length(ang_ind)),type="l",lwd=3,col=2)
legend("topleft",legend = c("True Probability","Median Estimate","95% Confidence Intervals"),lty=c('solid','solid','dashed'),col=c(2,1,"blue"),cex=1.2,lwd=2,bg="white")

plot(ang_ind,med_2,xlab="Angle Index",ylab = "Probability",ylim=range(med_2,lb_2,ub_2,p),type='n',main="WT isoline",cex.lab=1.3, cex.axis=1.4,cex.main=1.7)
polygon(c(rev(ang_ind), ang_ind), c(rev( lb_2), ub_2), col = 'grey80', border = NA)
lines(ang_ind,med_2,type="l",col=1,lwd=2)
lines(ang_ind, ub_2, lty = 'dashed', col = 'blue',lwd=2)
lines(ang_ind, lb_2, lty = 'dashed', col = 'blue',lwd=2)
lines(ang_ind,rep(p,length(ang_ind)),type="l",lwd=3,col=2)
legend("topleft",legend = c("True Probability","Median Estimate","95% Confidence Intervals"),lty=c('solid','solid','dashed'),col=c(2,1,"blue"),cex=1.2,lwd=2,bg="white")

res = list("ang_ind" = ang_ind,
           "med_x_ht04" = med_x_1, "lb_x_ht04" = lb_x_1, "ub_x_ht04" = ub_x_1,
           "med_y_ht04" = med_y_1, "lb_y_ht04" = lb_y_1, "ub_y_ht04" = ub_y_1,
           "med_x_wt13" = med_x_2, "lb_x_wt13" = lb_x_2, "ub_x_wt13" = ub_x_2,
           "med_y_wt13" = med_y_2, "lb_y_wt13" = lb_y_2, "ub_y_wt13" = ub_y_2,
           "med_ht04" = med_1, "lb_ht04" = lb_1, "ub_ht04" = ub_1,
           "med_wt13" = med_2, "lb_wt13" = lb_2, "ub_wt13" = ub_2) 
return(res) 
}