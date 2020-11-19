#' GPD parameter stability plots
#'
#' Plots showing the stability of the GPD scale and shape parameter estimates across a specified range of thresholds.
#'
#' @param Data Numeric vector containing the declusted data.
#' @param Data_Full Numeric vector containing the non-declustered data.
#' @param u Numeric vector of GPD thresholds; given as a quantiles \code{[0,1]} of \code{Data} vector. Default is \code{0.9} to \code{0.999} in intervals of \code{0.001}.
#' @param Plot Logical; indicating whether to plot diagnostics. Default is \code{FALSE}.
#' @return Plot of the shape and modified scale parameter estimates along with their errors bars over the range of specified thresholds.
#' @seealso \code{\link{Decluster}}
#' @export
#' @examples
#' GPD_Parameter_Stability_Plot(Data = S20.Detrend.Declustered.df$Rainfall,
#'                              Data_Full= na.omit(S20.Detrend.df$Rainfall),
#'                              u=seq(0.9,0.999,0.001))
GPD_Parameter_Stability_Plot<-function(Data,Data_Full,u=0.95,PLOT=FALSE,xlab_hist="Data",y_lab="Data"){

    shape_upper<-numeric(length(u))
    shape.estimate<-numeric(length(u))
    shape_lower<-numeric(length(u))

    scale_upper<-numeric(length(u))
    scale.estimate<-numeric(length(u))
    scale_lower<-numeric(length(u))

    for(i in 1:length(u)){
    gpd<-evm(na.omit(Data), th=quantile(Data_Full,u[i]),penalty = "gaussian",priorParameters = list(c(0, 0), matrix(c(100^2, 0, 0, 0.25), nrow = 2)))

    shape.estimate[i] <- gpd$coefficients[[2]]
    shape_upper[i] <- shape.estimate[i] + 1.96*gpd$se[[2]]
    shape_lower[i] <- shape.estimate[i] -1.96*gpd$se[[2]]

    scale<-sqrt(gpd$se[[1]]^2 - 2 * u[i] * gpd$cov[1,2] + (u[i] * gpd$se[[1]])^2)

    scale.estimate[i] <- gpd$coefficients[[1]]-shape.estimate[i]*u[i]
    scale_upper[i] <- scale.estimate[i] - 1.96 * scale
    scale_lower[i] <- scale.estimate[i] +  1.96 * scale

    }
plot(u,shape.estimate,ylim=c(min(shape_lower,na.rm=T)-0.1,max(shape_upper,na.rm=T)+0.1),pch=16,ylab="Shape")
for(i in 1:length(u)){
segments(u[i],shape.estimate[i],u[i],shape_upper[i])
segments(u[i],shape.estimate[i],u[i],shape_lower[i])
}

plot(u,scale.estimate,ylim=c(min(scale_lower,na.rm=T)-0.1,max(scale_upper,na.rm=T)+0.1),pch=16,ylab="Modified scale")
for(i in 1:length(u)){
  segments(u[i],scale.estimate[i],u[i],scale_upper[i])
  segments(u[i],scale.estimate[i],u[i],scale_lower[i])
}
}

