#' Derives return curves that capture uncertainty
#'
#' Calculates return level curves using two extremal dependence models: (1) Heffernan and Tawn (2004) (hereon in HT04) and (2) Wadsworth and Tawn (2013) (hereon in WT13) as outlined in Murphy-Barltrop et al. (2023).
#'
#' @param data Data frame of raw data detrended if necessary. First column should be of class \code{Date}.
#' @param q Numeric vector of length one specifying quantile level for fitting GPDs and the HT04 and WT13 models.
#' @param rp Numeric vector of length one specifying return period of interest.
#' @param mu Numeric vector of length one specifying the (average) occurrence frequency of events in Data. Default is 365.25, daily data.
#' @param n_sim Numeric vector of length one specifying the number of simulations for HT model. Default is \code{50}.
#' @param n_grad Numeric vector of length one specifying number of number of rays along which to compute points on the curve. Default is \code{50}.
#' @param boot_method Character vector of length one specifying the bootstrap method. Options are \code{"basic"} (default), \code{"block"} or \code{"monthly"}.
#' @param boot_replace Character vector of length one specifying whether bootstrapping is carried out with \code{"T"} or without \code{"F"} replacement. Only required if \code{boot_method = "basic"}. Default is \code{NA}.
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
#' @param alpha Numeric vector of length one specifying the \code{100(1-alpha)\%} confidence interval. Default is \code{0.1}.
#' @param x_lab Character vector specifying the x-axis label. Default is \code{colnames(data)[2]}.
#' @param y_lab Character vector specifying the y-axis label. Default is \code{colnames(data)[3]}.
#' @param x_lim_min Numeric vector of length one specifying x-axis minimum. Default is \code{min(data[,2],na.rm=T)}.
#' @param x_lim_max Numeric vector of length one specifying x-axis maximum. Default is \code{max(data[,2],na.rm=T)+0.3*diff(range(data[,3],na.rm=T))}.
#' @param y_lim_min Numeric vector of length one specifying y-axis minimum. Default is \code{min(data[,3],na.rm=T)}.
#' @param y_lim_max Numeric vector of length one specifying y-axis maximum. Default is \code{max(data[,3],na.rm=T)+0.3*diff(range(data[,3],na.rm=T))}.
#' @param plot Logical; whether to plot return curves. Default is \code{"TRUE"}.
#' @section Details:
#' The HT04 model is fit to two conditional samples. One sample comprises the declustered time series of the first variable paired with concurrent values of the other variable. The second sample is obtained in the same way but with the variables reversed.
#' @return List comprising the median curve based on the HT04 model \code{median_ht04}, and the upper \code{ub_ht04} and lower \code{lb_ht04} bound of its \code{100(1-alpha)\%} confidence interval. Analogous results for the curve based on the WT13 method \code{median_wt13}, \code{ub_wt13} and \code{lb_wt13}. Plots of the median return level curve and associated \code{100(1-alpha)\%} confidence interval for both extremal dependence models.
#' @export
#' @examples
#' #Data starts on first day of 1948
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
#' return_curve_est(data=S22.Detrend.df.extended[,1:3],
#'                  q=0.985,rp=100,mu=365.25,n_sim=100,
#'                  n_grad=50,n_boot=100,boot_method="monthly",
#'                  boot_replace=NA, block_length=NA, boot_prop=0.8,
#'                  decl_method_x="runs", decl_method_y="runs",
#'                  window_length_x=NA,window_length_y=NA,
#'                  u_x=0.95, u_y=0.95,
#'                  sep_crit_x=36, sep_crit_y=36,
#'                  alpha=0.1, x_lab=NA, y_lab=NA)
return_curve_est = function(data,q,rp,mu,n_sim,n_grad,n_boot,boot_method, boot_replace, block_length, boot_prop, decl_method_x, decl_method_y, window_length_x, window_length_y, u_x=NA, u_y=NA, sep_crit_x=NA, sep_crit_y=NA, alpha=0.1, x_lab=colnames(data)[2], y_lab=colnames(data)[2],x_lim_min=min(data[,2],na.rm=T),x_lim_max=max(data[,2],na.rm=T)+0.3*diff(range(data[,2],na.rm=T)),y_lim_min=min(data[,3],na.rm=T),y_lim_max=max(data[,3],na.rm=T)+0.3*diff(range(data[,3],na.rm=T)),plot=T){

 #Convert return period (rp) to probability
 prob = (1/mu) / rp

 #List to output the bootstrapped curves
 boot_curves = list()
 boot_curves2 = list()

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

 #Fit GPDs
 thresh_x = quantile(boot_df[,2],q,na.rm=T)
 gpd_x = evm(data_x[,1],th = thresh_x)
 gpd_x$rate = length(data_x[,1][data_x[,1]>thresh_x])/(length(na.omit(data_x[,1]))/mu)

 thresh_y = quantile(boot_df[,3],q,na.rm=T)
 gpd_y = evm(data_y[,2],th = thresh_y)
 gpd_y$rate = length(data_y[,2][data_y[,2]>thresh_y])/(length(na.omit(data_y[,2]))/mu)

 data_x = na.omit(data_x)
 data_y = na.omit(data_y)

 #Converting to uniform scale
 data_unif_x = cbind(EmpFun(x = boot_df[,2],r = data_x[,1],mod = gpd_x),EmpFun(x = boot_df[,3], r = data_x[,2], mod = gpd_y))
 data_unif_y = cbind(EmpFun(x = boot_df[,2],r = data_y[,1],mod = gpd_x),EmpFun(x = boot_df[,3], r = data_y[,2], mod = gpd_y))

 data_exp_x = apply(data_unif_x, 2, qexp)
 data_exp_y = apply(data_unif_y, 2, qexp)

 #estimated curves exponential margins with HT model
 curve = heff_tawn_curve_exp(data_exp_x = data_exp_x, data_exp_y = data_exp_y, prob = prob, q=q, nsim=n_sim)

 #estimated curves uniform margins
 curve_unif = apply(curve, 2, pexp)

 #estimated curves original margins
 curve_original = cbind(revTransform(x=curve_unif[,1],data = na.omit(boot_df[,2]),qu=q,th=thresh_x,exp(coefficients(gpd_x)[1]),coefficients(gpd_x)[2]),revTransform(x=curve_unif[,2],data = na.omit(boot_df[,3]),qu=q,th=thresh_y,exp(coefficients(gpd_y)[1]),coefficients(gpd_y)[2]))

 boot_curves[[j]] = curve_original

 #estimated curves exponential margins with WT model
 curve2.1 = wads_tawn_curve_exp(data_exp = data_exp_x,prob=prob,q=q)
 curve2.2 = wads_tawn_curve_exp(data_exp = data_exp_y,prob=prob,q=q)

 rad2.1 = apply(curve2.1,1,function(x){sqrt(sum(x^2))})
 rad2.2 = apply(curve2.2,1,function(x){sqrt(sum(x^2))})

 cons_curve = curve2.1

 cons_curve[rad2.1>rad2.2,] = curve2.1[rad2.1>rad2.2,]
 cons_curve[rad2.1<=rad2.2,] = curve2.2[rad2.1<=rad2.2,]

 #estimated curves uniform margins
 curve_unif2 = apply(cons_curve, 2, pexp)

 #estimated curves original margins
 curve_original2 = cbind(revTransform(x=curve_unif2[,1],data = na.omit(boot_df[,2]),qu=q,th=thresh_x,exp(coefficients(gpd_x)[1]),coefficients(gpd_x)[2]),revTransform(x=curve_unif2[,2],data = na.omit(boot_df[,3]),qu=q,th=thresh_y,exp(coefficients(gpd_y)[1]),coefficients(gpd_y)[2]))

 boot_curves2[[j]] = curve_original2

}

angles = ((n_grad:1)/(n_grad+1))*(pi/2)
x0 = min(data[,2], na.rm=T)
y0 = min(data[,3],na.rm=T)

angles_est_points = array(0,dim=c(n_grad,2,n_boot))
angles_est_points2 = array(0,dim=c(n_grad,2,n_boot))

#Finding interceptions of curves with rays
for(k in 1:n_boot){
  curve_angles = atan((boot_curves[[k]][,2]-y0)/(boot_curves[[k]][,1]-x0))

  for(i in 1:n_grad){
    ind = min(which(angles[i] >= curve_angles))
    x1 = boot_curves[[k]][ind,1] - x0
    x2 = boot_curves[[k]][ind-1,1] - x0
    y1 = boot_curves[[k]][ind,2] - y0
    y2 = boot_curves[[k]][ind-1,2] - y0
    p = (x1*tan(angles[i]) - y1) / ( (y2-y1) - (x2-x1)*tan(angles[i]) )
    xhat = x1 + p*(x2-x1) + x0 #change if computing uncertainty
    yhat = y1 + p*(y2-y1) + y0
    angles_est_points[i,,k] = c(xhat,yhat)
  }

  curve_angles2 = atan((boot_curves2[[k]][,2]-y0)/(boot_curves2[[k]][,1]-x0))

  for(i in 1:n_grad){
    ind = min(which(angles[i] >= curve_angles2))
    x1 = boot_curves2[[k]][ind,1] - x0
    x2 = boot_curves2[[k]][ind-1,1] - x0
    y1 = boot_curves2[[k]][ind,2] - y0
    y2 = boot_curves2[[k]][ind-1,2] - y0
    p = (x1*tan(angles[i]) - y1) / ( (y2-y1) - (x2-x1)*tan(angles[i]) )
    xhat = x1 + p*(x2-x1) + x0 #change if computing uncertainty
    yhat = y1 + p*(y2-y1) + y0
    angles_est_points2[i,,k] = c(xhat,yhat)
  }
}

#Calculating median and 90% confidence interval
median = array(0,dim=c(n_grad,2))
lower_bound = array(0,dim=c(n_grad,2))
upper_bound = array(0,dim=c(n_grad,2))

median2 = array(0,dim=c(n_grad,2))
lower_bound2 = array(0,dim=c(n_grad,2))
upper_bound2 = array(0,dim=c(n_grad,2))

for(i in 1:n_grad){
  median[i,] = angles_est_points[i,,order(angles_est_points[i,2,1:n_boot])[n_boot/2]]
  lower_bound[i,] = angles_est_points[i,,order(angles_est_points[i,2,1:n_boot])[n_boot*alpha/2]]
  upper_bound[i,] = angles_est_points[i,,order(angles_est_points[i,2,1:n_boot])[n_boot*(1-alpha/2)]]

  median2[i,] = angles_est_points2[i,,order(angles_est_points2[i,2,1:n_boot])[n_boot/2]]
  lower_bound2[i,] = angles_est_points2[i,,order(angles_est_points2[i,2,1:n_boot])[n_boot*alpha/2]]
  upper_bound2[i,] = angles_est_points2[i,,order(angles_est_points2[i,2,1:n_boot])[n_boot*(1-alpha/2)]]
}

if(plot==T){
par(mfrow=c(2,1))
par(mar=c(4.2,4.5,0.1,0.1))

#plotting results
plot(data[,c(2,3)],pch=16,col="Grey",xlab=x_lab,ylab=y_lab,cex.lab=1.2, cex.axis=1.2,cex.main=1.5,xlim=c(x_lim_min,x_lim_max),ylim=c(y_lim_min,y_lim_max),main="")
lines(median,col="Red3",lwd=2)
lines(upper_bound,lty=2,lwd=2)
lines(lower_bound,lty=2,lwd=2)

plot(data[,c(2,3)],pch=16,col="Grey",xlab=x_lab,ylab=y_lab,cex.lab=1.2, cex.axis=1.2,cex.main=1.5,xlim=c(x_lim_min,x_lim_max),ylim=c(y_lim_min,y_lim_max),main="")
lines(median2,col="green3",lwd=2)
lines(upper_bound2,lty=2,lwd=2)
lines(lower_bound2,lty=2,lwd=2)
}

res = list("median_ht04" = median, "ub_ht04" = upper_bound, "lb_ht04" = lower_bound, "median_wt13" = median2, "ub_wt13"= upper_bound2, "lb_wt13" = lower_bound2)
return(res)
}
