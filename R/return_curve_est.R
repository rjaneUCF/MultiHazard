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
#' @param boot_replace Character vector of length one specifying whether bootstrapping is carried out with \code{"TRUE"} or without \code{"FALSE"} replacement. Only required if \code{boot_method = "basic"}. Default is \code{NA}.
#' @param block_length Numeric vector of length one specifying block length. Only required if \code{boot_method = "block"}. Default is \code{NA}.
#' @param boot_prop Numeric vector of length one specifying the minimum proportion of non-missing values of at least of the variables for a month to be included in the bootstrap. Only required if \code{boot_method = "monthly"}. Default is \code{0.8}.
#' @param n_boot Numeric vector of length one specifying number of bootstrap samples. Default is \code{100}.
#' @param decl_method_x Character vector of length one specifying the declustering method to apply to the first variable. Options are the storm window approach \code{"window"} (default) and the runs method \code{"runs"}.
#' @param decl_method_y Character vector of length one specifying the declustering method to apply to the second variable. Options are the storm window approach \code{"window"} (default) and the runs method \code{"runs"}.
#' @param window_length_x Numeric vector of length one specifying the storm window length to apply during the declustering of the  first variable if \code{decl_method_x = "window"}.
#' @param window_length_y Numeric vector of length one specifying the storm window length to apply during the declustering of the second variable if \code{decl_method_y = "window"}.
#' @param u_x Numeric vector of length one specifying the (quantile) threshold to adopt in the declustering of the first variable if \code{decl_method_x = "runs"}. Default is \code{NA}.
#' @param u_y Numeric vector of length one specifying the (quantile) threshold to adopt in the declustering of the second variable if \code{decl_method_y = "runs"}. Default is \code{NA}.
#' @param sep_crit_x Numeric vector of length one specifying the separation criterion to apply during the declustering of the first variable if \code{decl_method_x = "runs"}. Default is \code{NA}.
#' @param sep_crit_y Numeric vector of length one specifying the separation criterion to apply during the declustering of the second variable if \code{decl_method_y = "runs"}. Default is \code{NA}.
#' @param alpha Numeric vector of length one specifying the \code{100(1-alpha)\%} confidence interval. Default is \code{0.1}.
#' @param most_likely Character vector of length one specifying whether to estimate the relative likelihood of events along the curves. For the ht04 curve probabilites are estimated by simulating from the ht04 model while for the wt13 curve a two-sided conditional sampling (using q as the threshold for both samples) copula theory is adopted. Default is \code{FALSE}.
#' @param n_interp Numeric vector of length one specifying the resolution of the interpolation of the curves Default is \code{1000} thus the curve will be composed of \code{1000} points. The interpolation is only carried out if \code{most-likely = TRUE}.
#' @param n Numeric vector of length one specifying the size of the sample from the fitted joint distributions used to estimate the density along an return curves. Default is \code{10^6}
#' @param n_ensemble Numeric vector of length one specifying the number of possible design events sampled along the two curves. Default is \code{0}
#' @param x_lab Character vector specifying the x-axis label. Default is \code{colnames(data)[2]}.
#' @param y_lab Character vector specifying the y-axis label. Default is \code{colnames(data)[3]}.
#' @param x_lim_min Numeric vector of length one specifying x-axis minimum. Default is \code{min(data[,2],na.rm=TRUE)}.
#' @param x_lim_max Numeric vector of length one specifying x-axis maximum. Default is \code{max(data[,2],na.rm=TRUE)+0.3*diff(range(data[,3],na.rm=TRUE))}.
#' @param y_lim_min Numeric vector of length one specifying y-axis minimum. Default is \code{min(data[,3],na.rm=TRUE)}.
#' @param y_lim_max Numeric vector of length one specifying y-axis maximum. Default is \code{max(data[,3],na.rm=TRUE)+0.3*diff(range(data[,3],na.rm=TRUE))}.
#' @param plot Logical; whether to plot return curves. Default is \code{"TRUE"}.
#' @section Details:
#' The HT04 model is fit to two conditional samples. One sample comprises the declustered time series of the first variable paired with concurrent values of the other variable. The second sample is obtained in the same way but with the variables reversed.
#' @return List comprising the median curve based on the HT04 model \code{median_ht04}, and the upper \code{ub_ht04} and lower \code{lb_ht04} bound of its \code{100(1-alpha)\%} confidence interval. Analogous results for the curve based on the WT13 method \code{median_wt13}, \code{ub_wt13} and \code{lb_wt13}. If \code{plot=TRUE} the median return level curve and associated \code{100(1-alpha)\%} confidence intervals are plotted for both extremal dependence models. If \code{most-likely=TRUE}, the relative probability of events on the two curves is also returned \code{contour_ht04} and \code{contour_wt13} along with the "most-likely"design event \code{most_likely_ht04} and \code{most_likely_wt13}, and an ensemble of possible "design events"sampled along the curve \code{ensemble_ht04} and \code{ensemble_wt13}.
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
return_curve_est = function(data,q,rp,mu,n_sim,n_grad,n_boot,boot_method, boot_replace, block_length, boot_prop, decl_method_x, decl_method_y, window_length_x, window_length_y, u_x=NA, u_y=NA, sep_crit_x=NA, sep_crit_y=NA, alpha=0.1, most_likely=FALSE, n_interp=1000, n=10^6, n_ensemble=0, x_lab=colnames(data)[2], y_lab=colnames(data)[3],x_lim_min=min(data[,2],na.rm=TRUE),x_lim_max=max(data[,2],na.rm=TRUE)+0.3*diff(range(data[,2],na.rm=TRUE)),y_lim_min=min(data[,3],na.rm=TRUE),y_lim_max=max(data[,3],na.rm=TRUE)+0.3*diff(range(data[,3],na.rm=TRUE)),plot=TRUE){

  if (missing(data) || is.null(data)) {
    stop("data is missing.")
  }

  if (!is.data.frame(data)) {
    stop("data must be a data.frame.")
  }

   #Date column validation
  if (!inherits(data[,1], c("Date", "POSIXct", "POSIXt"))) {
    stop("First column of data must be Date/POSIXct format.")
  }


  #Quantile parameter must be numeric and between 0 and 1
  if (missing(q)) {
    stop("q is missing.")
  }

  if (!is.numeric(q) || length(q) != 1) {
    stop("q must be a single numeric value.")
  }

  if (q <= 0 || q >= 1) {
    stop("q must be between 0 and 1 (exclusive).")
  }

  #Return period parameter must be a single numeric value
  if (missing(rp)) {
    stop("rp is missing.")
  }

  if (!is.numeric(rp) || length(rp) != 1 || rp < 1) {
    stop("rp must be a single positive numeric value.")
  }

  #Rate parameter must be a single positive numeric value
  if (missing(mu)) {
    stop("mu (rate parameter) is missing.")
  }

  if (!is.numeric(mu) || length(mu) != 1 || mu <= 0) {
    stop("mu (rate parameter) must be a single positive numeric value.")
  }

  #Simulation parameters
  if (!is.numeric(n_sim) || length(n_sim) != 1 || n_sim <= 0) {
    stop("n_sim must be a positive integer.")
  }

  if (!is.numeric(n_grad) || length(n_grad) != 1 || n_grad <= 0) {
    stop("n_grad must be a positive integer.")
  }

  if (!is.numeric(n_boot) || length(n_boot) != 1 || n_boot <= 0 || n_boot %% 2 != 0) {
    stop("n_boot must be a positive even integer.")
  }


  # Bootstrap method validation
  valid_boot_methods <- c("basic", "block", "monthly")
  if (!boot_method %in% valid_boot_methods) {
    stop(paste("boot_method must be one of:", paste(valid_boot_methods, collapse = ", ")))
  }

  if (boot_method == "block") {
    if (missing(block_length) || !is.numeric(block_length) || block_length <= 0) {
      stop("block_length must be positive when using block bootstrap.")
    }
    if (block_length >= nrow(data)) {
      stop("block_length must be less than data length.")
    }
  }

  if (boot_method == "monthly") {
    if (missing(boot_prop) || !is.numeric(boot_prop) || boot_prop <= 0 || boot_prop > 1) {
      stop("boot_prop must be between 0 and 1 for monthly bootstrap.")
    }
  }

  # Declustering method type validation
  valid_decl_methods <- c("window", "runs")
  if (!decl_method_x %in% valid_decl_methods) {
    stop(paste("decl_method_x must be one of:", paste(valid_decl_methods, collapse = ", ")))
  }

  if (!decl_method_y %in% valid_decl_methods) {
    stop(paste("decl_method_y must be one of:", paste(valid_decl_methods, collapse = ", ")))
  }

  # Declustering (Window method) window length validation
  if (decl_method_x == "window") {
    if (missing(window_length_x)) {
      stop("Error: window_length_x must be specified when using window declustering.")
    }
    if (missing(window_length_x) || !is.numeric(window_length_x) || window_length_x <= 0) {
      stop("Error: window_length_x must be positive when using window declustering.")
    }
  }

  if (decl_method_y == "window") {
    if (missing(window_length_y)) {
      stop("Error: window_length_y must be specified when using window declustering.")
    }
    if (missing(window_length_y) || !is.numeric(window_length_y) || window_length_y <= 0) {
      stop("Error: window_length_y must be positive when using window declustering.")
    }
  }

  # Declustering (runs method) threshold validation
  if (decl_method_x == "runs") {
    if (is.na(u_x) || !is.numeric(u_x)) {
      stop("u_x threshold must be numeric when using runs declustering.")
    }
    if (is.na(sep_crit_x) || !is.numeric(sep_crit_x)) {
      stop("sep_crit_x must be numeric when using runs declustering.")
    }
    if (sep_crit_x <= 0) {
      stop("sep_crit_x must be positive when using runs declustering.")
    }
  }

  if (decl_method_y == "runs") {
    if (is.na(u_y) || !is.numeric(u_y)) {
      stop("u_y threshold must be numeric when using runs declustering.")
    }
    if (is.na(sep_crit_y) || !is.numeric(sep_crit_y)) {
      stop("sep_crit_y must be numeric when using runs declustering.")
    }
    if (sep_crit_y <= 0) {
      stop("sep_crit_y must be positive when using runs declustering.")
    }
  }

  #Alpha paramter validation
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) {
    stop("alpha must be between 0 and 1 (exclusive).")
  }

  #Additional paramter validation
  if (!is.logical(most_likely) || length(most_likely) != 1) {
    stop("most_likely must be TRUE or FALSE.")
  }

  if (!is.numeric(n_interp) || length(n_interp) != 1 || n_interp <= 0) {
    stop("n_interp must be a positive integer.")
  }

  if (!is.numeric(n) || length(n) != 1 || n <= 0) {
    stop("n must be a positive integer.")
  }

  if (!is.numeric(n_ensemble) || length(n_ensemble) != 1 || n_ensemble < 0) {
    stop("n_ensemble must be a non-negative integer.")
  }

  #Data quality checks
  if (any(is.infinite(data[,2]), na.rm = TRUE) || any(is.infinite(data[,3]), na.rm = TRUE)) {
    stop("data contains infinite values.")
  }

  #Declustering (runs method) threshold range validation
  if (decl_method_x == "runs" && !is.na(u_x)) {
    if (u_x < 0 || u_x > 1 ) {
      stop("u_x threshold must be between 0 and 1.")
    }
  }

  if (decl_method_y == "runs" && !is.na(u_y)) {
    if (u_y < 0 || u_y > 1) {
      stop("u_y threshold must be between 0 and 1.")
    }
  }

 #Convert return period (rp) to probability
 prob = (1/mu) / rp

 #List to output the bootstrapped curves
 boot_curves = list()
 boot_curves2 = list()

 #Result vector for most-likely event
 contour_ht04 = NA
 contour_wt13 = NA
 most_likely_ht04 = NA
 most_likely_wt13 = NA
 ensemble_ht04 = NA
 ensemble_wt13 = NA

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
 thresh_x = quantile(boot_df[,2],q,na.rm=TRUE)
 gpd_x = evm(data_x[,1],th = thresh_x)
 gpd_x$rate = length(data_x[,1][data_x[,1]>thresh_x])/(length(na.omit(data_x[,1]))/mu)

 thresh_y = quantile(boot_df[,3],q,na.rm=TRUE)
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
x0 = min(data[,2], na.rm=TRUE)
y0 = min(data[,3],na.rm=TRUE)

angles_est_points = array(0,dim=c(n_grad,2,n_boot))
angles_est_points2 = array(0,dim=c(n_grad,2,n_boot))

#Finding interceptions of curves with rays
for(k in 1:n_boot){
  curve_angles = atan((boot_curves[[k]][,2]-y0)/(boot_curves[[k]][,1]-x0))

  for(i in 1:n_grad){
   match = which(angles[i] >= curve_angles)
   if(length(match)>0){
    ind = min(match)
    x1 = boot_curves[[k]][ind,1] - x0
    x2 = boot_curves[[k]][ind-1,1] - x0
    y1 = boot_curves[[k]][ind,2] - y0
    y2 = boot_curves[[k]][ind-1,2] - y0
    p = (x1*tan(angles[i]) - y1) / ( (y2-y1) - (x2-x1)*tan(angles[i]) )
    xhat = x1 + p*(x2-x1) + x0 #change if computing uncertainty
    yhat = y1 + p*(y2-y1) + y0
    angles_est_points[i,,k] = c(xhat,yhat)
   } else{
    angles_est_points[i,,k] = c(NA,NA)
   }
  }

  curve_angles2 = atan((boot_curves2[[k]][,2]-y0)/(boot_curves2[[k]][,1]-x0))

for(i in 1:n_grad){
 match = which(angles[i] >= curve_angles2)
  if(length(match)>0){
   ind = min(match)
   x1 = boot_curves2[[k]][ind,1] - x0
   x2 = boot_curves2[[k]][ind-1,1] - x0
   y1 = boot_curves2[[k]][ind,2] - y0
   y2 = boot_curves2[[k]][ind-1,2] - y0
   p = (x1*tan(angles[i]) - y1) / ( (y2-y1) - (x2-x1)*tan(angles[i]) )
   xhat = x1 + p*(x2-x1) + x0 #change if computing uncertainty
   yhat = y1 + p*(y2-y1) + y0
   angles_est_points2[i,,k] = c(xhat,yhat)
  } else{
   angles_est_points2[i,,k] = c(NA,NA)
  }
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
  median[i,] = angles_est_points[i,,order(angles_est_points[i,2,1:n_boot])[round(n_boot/2)]]
  lower_bound[i,] = angles_est_points[i,,order(angles_est_points[i,2,1:n_boot])[round(n_boot*alpha/2)]]
  upper_bound[i,] = angles_est_points[i,,order(angles_est_points[i,2,1:n_boot])[round(n_boot*(1-alpha/2))]]

  median2[i,] = angles_est_points2[i,,order(angles_est_points2[i,2,1:n_boot])[round(n_boot/2)]]
  lower_bound2[i,] = angles_est_points2[i,,order(angles_est_points2[i,2,1:n_boot])[round(n_boot*alpha/2)]]
  upper_bound2[i,] = angles_est_points2[i,,order(angles_est_points2[i,2,1:n_boot])[round(n_boot*(1-alpha/2))]]
}

colnames(median) <- c(names(data)[2],names(data)[3])
colnames(median2) <- c(names(data)[2],names(data)[3])
colnames(lower_bound) <- c(names(data)[2],names(data)[3])
colnames(lower_bound2) <- c(names(data)[2],names(data)[3])
colnames(upper_bound) <- c(names(data)[2],names(data)[3])
colnames(upper_bound2) <- c(names(data)[2],names(data)[3])


#(relative) Probabilities implied by the data for the points composing the isoline. Probabilities are scaled to [0,1].
if(most_likely==TRUE){

  #Interpolating the isolines
  x_u = (median[,1] - min(median[,1]))/(max(median[,1])-min(median[,1]))
  y_u = (median[,2] - min(median[,2]))/(max(median[,2])-min(median[,2]))
  x.diff<-c(0,diff(x_u))
  y.diff<-c(0,diff(y_u))
  d<-sqrt(x.diff^2 + y.diff^2)
  #Linearly interpolate the x values with respect to their cumulative distance along the isoline.
  median_approx_x<-approx(x=cumsum(d),y=x_u,xout=seq(0,sum(d),length.out=n_interp))$y * (max(median[,1])-min(median[,1])) + min(median[,1])
  #Linearly interpolate the y values with respect to their cumulative distance along the isoline.
  median_approx_y<-approx(x=cumsum(d),y=y_u,xout=seq(0,sum(d),length.out=n_interp))$y  * (max(median[,2])-min(median[,2])) + min(median[,2])
  median = data.frame(median_approx_x,median_approx_y)

  x_u = (median2[,1] - min(median2[,1]))/(max(median2[,1])-min(median2[,1]))
  y_u = (median2[,2] - min(median2[,2]))/(max(median2[,2])-min(median2[,2]))
  x.diff<-c(0,diff(x_u))
  y.diff<-c(0,diff(y_u))
  d<-sqrt(x.diff^2 + y.diff^2)
  #Linearly interpolate the x values with respect to their cumulative distance along the isoline.
  median2_approx_x<-approx(x=cumsum(d),y=x_u,xout=seq(0,sum(d),length.out=n_interp))$y * (max(median2[,1])-min(median2[,1])) + min(median2[,1])
  #Linearly interpolate the y values with respect to their cumulative distance along the isoline.
  median2_approx_y<-approx(x=cumsum(d),y=y_u,xout=seq(0,sum(d),length.out=n_interp))$y  * (max(median2[,2])-min(median2[,2])) + min(median2[,2])
  median2 = data.frame(median2_approx_x,median2_approx_y)

  colnames(median) <- c(names(data)[2],names(data)[3])
  colnames(median2) <- c(names(data)[2],names(data)[3])

  #Decluster time series

  #Variable 1
  if(decl_method_x == "window"){
    decl = Decluster_SW(Data=data[,1:2], Window_Width = window_length_x)
    x_dec_df = data.frame(as.Date(data[,1]),decl$Declustered)
  }

  if(decl_method_x == "runs"){
    decl = Decluster(Data=data[,2], u=u_x, SepCrit = sep_crit_x, mu = mu)
    x_dec_df = data.frame(as.Date(data[,1]),decl$Declustered)
  }

  colnames(x_dec_df) = colnames(data[,1:2])
  x_dec_df$Date = as.Date(x_dec_df$Date)

  #Variable 2
  if(decl_method_y == "window"){
    decl = Decluster_SW(Data=data[,c(1,3)], Window_Width = window_length_y)
    y_dec_df = data.frame(as.Date(data[,1]),decl$Declustered)
  }

  if(decl_method_y == "runs"){
    decl = Decluster(Data=data[,3], u = u_y, SepCrit = sep_crit_y, mu = mu)
    y_dec_df = data.frame(as.Date(data[,1]),decl$Declustered)
  }

  colnames(y_dec_df) = colnames(data[,c(1,3)])
  y_dec_df$Date = as.Date(y_dec_df$Date)

  #Dataframes
  df = data.frame("Date"=seq.Date(as.Date(min(data$Date,data$Date)),as.Date(max(data$Date,data$Date)),by="day"))
  df_1 = left_join(df,x_dec_df,by="Date")
  decl_df = left_join(df_1,y_dec_df,by="Date")
  colnames(decl_df) = colnames(data)

  #Simulating from ht04 model for ht04 curve
  gpds = Migpd_Fit(Data = decl_df[,2:3],
                   Data_Full =  data[,2:3],
                   mqu = q)

  ht04.sim= HT04(data_Detrend_Dependence_df = data[,2:3],
                 data_Detrend_Declustered_df = decl_df[,2:3],
                 u_Dependence = q,
                 Migpd =gpds,
                 mu=mu,
                 N=n,
                 Margins="laplace",
                 V=10,
                 Maxit=10000)$x.sim

  #copula approach for wt13 method
  #Dataframes for marginal dist.
  data_x = data.frame(decl_df[,2],data[,3])
  data_y = data.frame(data[,2],decl_df[,3])

  colnames(data_x) = colnames(data[,2:3])
  colnames(data_y) = colnames(data[,2:3])


  #Fit GPDs
  thresh_x = quantile(data[,2],q,na.rm=TRUE)
  gpd_x = evm(data_x[,1],th = thresh_x)
  gpd_x$rate = length(data_x[,1][data_x[,1]>thresh_x])/(length(na.omit(data_x[,1]))/mu)
  thresh_y = quantile(data[,3],q,na.rm=TRUE)
  gpd_y = evm(data_y[,2],th = thresh_y)
  gpd_y$rate = length(data_y[,2][data_y[,2]>thresh_y])/(length(na.omit(data_y[,2]))/mu)

  #Fit copula
  cop.con.x = BiCopSelect(pobs(data_x[which(data_x[,1]>thresh_x),1]),pobs(data_x[which(data_x[,1]>thresh_x),2]))
  cop.con.y = BiCopSelect(pobs(data_y[which(data_y[,2]>thresh_y),1]),pobs(data_y[which(data_y[,2]>thresh_y),2]))

  #Simulate from copula
  cop.sim.x.u = BiCopSim(round(n*nrow(data_x[which(data_x[,1]>thresh_x),])/(nrow(data_x[which(data_x[,1]>thresh_x),])+nrow(data_y[which(data_y[,2]>thresh_y),])),0),cop.con.x)
  cop.sim.y.u = BiCopSim(round(n*nrow(data_y[which(data_y[,2]>thresh_y),])/(nrow(data_x[which(data_x[,1]>thresh_x),])+nrow(data_y[which(data_y[,2]>thresh_y),])),0),cop.con.y)

  #Transform to original scale
  cop.sim.x = data.frame(qgpd(cop.sim.x.u[,1],u=thresh_x,exp(coefficients(gpd_x)[1]),coefficients(gpd_x)[2]),
                         revTransform(x=cop.sim.x.u[,2],data = na.omit(data[,3]),qu=q,th=thresh_y,exp(coefficients(gpd_y)[1]),coefficients(gpd_y)[2]))

  cop.sim.y = data.frame(revTransform(x=cop.sim.y.u[,1],data = na.omit(data[,2]),qu=q,th=thresh_x,exp(coefficients(gpd_x)[1]),coefficients(gpd_x)[2]),
                         qgpd(cop.sim.y.u[,2],u=thresh_y,exp(coefficients(gpd_y)[1]),coefficients(gpd_y)[2]))

  colnames(cop.sim.x) = colnames(data[,2:3])
  colnames(cop.sim.y) = colnames(data[,2:3])
  cop.sim = rbind(cop.sim.x,cop.sim.y)

  #Kernel density estimates at points on curve
  prediction_ht04<-kde(x=ht04.sim, eval.points=median)$estimate
  prediction_wt13<-kde(x=cop.sim, eval.points=median2)$estimate

  #Convert probabilities on curve to relative probabilities
  contour_ht04 <- (prediction_ht04-min(prediction_ht04))/(max(prediction_ht04)-min(prediction_ht04))
  contour_wt13 <- (prediction_wt13-min(prediction_wt13))/(max(prediction_wt13)-min(prediction_wt13))

  #Find "most-likely" event
  most_likely_ht04<-data.frame(as.numeric(median[which(prediction_ht04==max(prediction_ht04,na.rm=TRUE)),1]),as.numeric(median[which(prediction_ht04==max(prediction_ht04,na.rm=TRUE)),2]))
  most_likely_wt13<-data.frame(as.numeric(median2[which(prediction_wt13==max(prediction_wt13,na.rm=TRUE)),1]),as.numeric(median2[which(prediction_wt13==max(prediction_wt13,na.rm=TRUE)),2]))

  colnames(most_likely_ht04) <- c(names(data)[2],names(data)[3])
  colnames(most_likely_wt13) <- c(names(data)[2],names(data)[3])

  #Generate an ensemble of "design events"
  if(n_ensemble > 0){
   ensemble_ht04<-median[sample(1:length(prediction_ht04[prediction_ht04>0]),size = n_ensemble, replace = TRUE, prob=prediction_ht04[prediction_ht04>0]),]
   ensemble_wt13<-median2[sample(1:length(prediction_wt13[prediction_wt13>0]),size = n_ensemble, replace = TRUE, prob=prediction_wt13[prediction_wt13>0]),]

   colnames(ensemble_ht04) <- c(names(data)[2],names(data)[3])
   colnames(ensemble_wt13) <- c(names(data)[2],names(data)[3])
  }
}


if(plot==TRUE){
par(mfrow=c(2,1))
par(mar=c(4.2,4.5,0.1,0.1))

#plotting results
plot(data[,c(2,3)],pch=16,col="Grey",xlab=x_lab,ylab=y_lab,cex.lab=1.2, cex.axis=1.2,cex.main=1.5,xlim=c(x_lim_min,x_lim_max),ylim=c(y_lim_min,y_lim_max),main="")
lines(median,col="Red3",lwd=2)
if(most_likely==TRUE){
  points(median,col=rev(heat.colors(150))[20:120][1+100*contour_ht04],lwd=2,pch=16,cex=1.75)
  points(most_likely_ht04,pch=18,cex=1.75,col="dark green")
}
if(n_ensemble>0){
  points(ensemble_ht04,pch=16,cex=1)
}
lines(upper_bound,lty=2,lwd=2)
lines(lower_bound,lty=2,lwd=2)

plot(data[,c(2,3)],pch=16,col="Grey",xlab=x_lab,ylab=y_lab,cex.lab=1.2, cex.axis=1.2,cex.main=1.5,xlim=c(x_lim_min,x_lim_max),ylim=c(y_lim_min,y_lim_max),main="")
lines(median2,col="green3",lwd=2)
if(most_likely==TRUE){
  points(median2,col=rev(heat.colors(150))[20:120][1+100*contour_wt13],lwd=2,pch=16,cex=1.75)
  points(most_likely_wt13,pch=18,cex=1.75,col="dark green")
}
if(n_ensemble>0){
  points(ensemble_wt13,pch=16,cex=1)
}
lines(upper_bound2,lty=2,lwd=2)
lines(lower_bound2,lty=2,lwd=2)
}

res = list("median_ht04" = median, "ub_ht04" = upper_bound, "lb_ht04" = lower_bound, "contour_ht04" = contour_ht04, "most_likely_ht04" = most_likely_ht04, "ensemble_ht04" = ensemble_ht04, "median_wt13" = median2, "ub_wt13"= upper_bound2, "lb_wt13" = lower_bound2, "contour_wt13" = contour_wt13, "most_likely_wt13" = most_likely_wt13, "ensemble_wt13" = ensemble_wt13)
return(res)
}
