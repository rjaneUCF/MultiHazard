#' Inverse PIT GPD
#'
#' Transforms a uniform (0,1) sample to the original scale by invoking the inverse Probability Integral Transform (PIT).
#' Realizations above a high threshold are transformed through a user-specified Generalized Pareto Distribution (GPD) while those below are transformed through the empirical distribution.
#'
#' @param u Vector of the uniform random variates.
#' @param Data Vector of the observations.
#' @param Data_Declust Vector of the declustered observations.
#' @param q Numeric vector of length one, giving the quantile of \code{Data} above which the GPD is fit.
#' @return A vector of \code{u} transformed to the specified GPD.
#' @export
#' @examples
#' #First decluster the rainfall series to find the 500 events
#' #with the highest peaks
#' S13.Rainfall.Declust = Decluster(Data=S13.Detrend.df$Rainfall,
#'                                  SepCrit=24*3, u=0.99667)
#' #Generate some uniform (0,1) random variates
#' unif = runif(100,0,1)
#' #Transform the unifrom variate to the original scale
#' x.sim = inverse_pit_gpd(unif,S13.Detrend.df$Rainfall,S13.Rainfall.Declust$Declsutered,0.95)
#' #Plotting the empirical distribution functions of the sample and observations
#' plot(S13.Detrend.df$Rainfall[order(S13.Detrend.df$Rainfall)],
#'     (1:length(S13.Detrend.df$Rainfall))/length(S13.Detrend.df$Rainfall))
#' points(x.sim[order(x.sim)],1:length(x.sim)/length(x.sim),col=2)
inverse_pit_gpd = function(u,Data,Data_Declust,q){
  thres = quantile(Data,q)
  mod = evm(Data_Declust,th = thres,show=FALSE)
  x = c()
  x[u>q] = qgpd((u[u>q]-q)/(1-q),mod$threshold, sigma = exp(mod$par[1]), xi = mod$par[2])
  x[u<=q] = quantile(Data,u[u<=q])
  return(x)
}
