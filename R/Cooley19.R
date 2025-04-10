#' Derives bivariate isolines using the non-parametric approach of Cooley et al. (2019).
#'
#' The Cooley et al. (2019) method exploits bivariate regular variation and kernel density estimation to generate isolines of bivariate exceedance probabilities. The function utilizes the \code{ks} and \code{texmex} packages, and works for both asymptotic dependence and independence.
#'
#' @param Data Data frame consisting of two columns.
#' @param Migpd An \code{Migpd} object, containing the generalized  Pareto models fitted (independently) to the variables comprising the columns of \code{Data}.
#' @param p.base Numeric vector of length one specifying the exceedance probability of the base isoline. Default is \code{0.01}.
#' @param p.proj Numeric vector of length one specifying the exceedance probability of the projected isoline. Default is \code{0.001}.
#' @param u Numeric vector of length one specifying the quantile at which to estimate the asymptotic nature of the data i.e. chi and chibar. Default is \code{0.95}.
#' @param PLOT Logical; indicating whether to plot the base and projected isolines on the original and transformed scale. Default is \code{FALSE}.
#' @param x_lim_min_T Numeric vector of length one specifying the lower x-axis limit of the transformed scale plot. Default is \code{NA}.
#' @param x_lim_max_T Numeric vector of length one specifying the upper x-axis limit of the transformed scale plot. Default is \code{NA}.
#' @param y_lim_min_T Numeric vector of length one specifying the lower y-axis limit of the transformed scale plot. Default is \code{NA}.
#' @param y_lim_max_T Numeric vector of length one specifying the upper y-axis limit of the transformed scale plot. Default is \code{NA}.
#' @param x_lim_min Numeric vector of length one specifying the lower x-axis limit of the plot on the original scale. Default is \code{NA}.
#' @param x_lim_max Numeric vector of length one specifying the upper x-axis limit of the plot on the original scale. Default is \code{NA}.
#' @param y_lim_min Numeric vector of length one specifying the lower y-axis limit of the plot on the original scale. Default is \code{NA}.
#' @param y_lim_max Numeric vector of length one specifying the lower y-axis limit of the plot on the original scale. Default is \code{NA}.
#' @return List comprising a description of the type of (asymptoptic) dependence \code{Asym}, the values the extremal dependence measures \code{Chi} and \code{n.bar}, exceedance probabilities of the base \code{p.base} and projected \code{p.proj} isolines, as well as the points on the base \code{I.base} and projected \code{I.proj} isolines.
#' @seealso \code{\link{Dataframe_Combine}} \code{\link{Decluster}} \code{\link{GPD_Fit}} \code{\link{Migpd_Fit}}
#' @export
#' @examples
#' S20.GPD<-Migpd_Fit(Data=S20.Detrend.Declustered.df[,-1], mqu =c(0.99,0.99,0.99))
#' Cooley19(Data=na.omit(S20.Detrend.df[,3:4]),Migpd=s.Migpd,
#' p.base=0.01,p.proj=0.001,PLOT=TRUE,x_lim_max_T=500,y_lim_max_T=500)
Cooley19<-function(Data,Migpd,p.base=0.01,p.proj=0.001,u=0.95,PLOT=FALSE,x_lim_min_T=NA, x_lim_max_T=NA,y_lim_min_T=NA,y_lim_max_T=NA,x_lim_min=NA,x_lim_max=NA,y_lim_min=NA,y_lim_max=NA){

  #Check the asymptotic properities
  Chi<-chi(Data,nq=1000)
  CHI<-round(Chi$chi[which(abs(Chi$quantile-u)==min(abs(Chi$quantile-u))),2],2)
  n.bar<-(Chi$chibar[which(abs(Chi$quantile-u)==min(abs(Chi$quantile-u)))]+1)/2

  #Fit Kernel cumulative survival function to Data
  Fhat <- kcde(Data,gridsize=1000,tail.flag="upper.tail")

  #Extract contour
  I.base<-contourLines(x=Fhat$eval.points[[1]], y=Fhat$eval.points[[2]],z=Fhat$estimate, levels=p.base)[[1]]
  I.base <- data.frame(I.base$x,I.base$y)
  colnames(I.base)<-c("x","y")

  #Transform observations and contour to frechet (_F) scale
  x_F<- -1/(log(EmpFun(x=Data[,1],r=Data[,1], mod=Migpd$models[[1]])))
  y_F <- -1/(log(EmpFun(x=Data[,2],r=Data[,2], mod=Migpd$models[[2]])))
  I.base.x_F <- -1/(log(EmpFun(x=Data[,1], r=I.base$x, mod=Migpd$models[[1]])))
  I.base.y_F<-  -1/(log(EmpFun(x=Data[,2], r=I.base$y, mod=Migpd$models[[2]])))
  I.base_F <- data.frame(I.base.x_F,I.base.y_F)
  colnames(I.base_F)<-c("x","y")

  #Projected isoline
  if(CHI>0){
    I.proj.x_F <- (p.base/p.proj)*I.base.x_F
    I.proj.y_F <- (p.base/p.proj)*I.base.y_F
    Asym <- "Asymptotic dependence"
  } else{
    I.proj.x_F <- I.base.x_F*(p.base/p.proj)^(-1/n.bar)
    I.proj.y_F <- I.base.y_F*(p.base/p.proj)^(-1/n.bar)
    Asym <- "Asymptotic independence"
  }
  I.proj_F <- data.frame(I.proj.x_F,I.proj.y_F)
  colnames(I.proj_F)<-c("x","y")

  #Transforming back to original scale
  I.proj.x_U<- exp(-1/I.proj.x_F)
  I.proj.y_U<- exp(-1/I.proj.y_F)

  I.proj.x<-as.numeric(quantile(Data[,1],I.proj.x_U))
  I.proj.x[which(I.proj.x_U>Migpd$models[[1]]$mqu)]<-u2gpd(u=I.proj.x_U[which(I.proj.x_U>Migpd$models[[1]]$mqu)],
                                                           p=Migpd$models[[1]]$rate,
                                                           th=Migpd$models[[1]]$threshold,
                                                           sigma=exp(Migpd$models[[1]]$coefficients[1]),
                                                           xi=Migpd$models[[1]]$coefficients[2])

  I.proj.y<-as.numeric(quantile(Data[,2],I.proj.y_U))
  I.proj.y[which(I.proj.y_U>Migpd$models[[2]]$mqu)]<-u2gpd(u=I.proj.y_U[which(I.proj.y_U>Migpd$models[[2]]$mqu)],
                                                           p=Migpd$models[[2]]$rate,
                                                           th=Migpd$models[[2]]$threshold,
                                                           sigma=exp(Migpd$models[[2]]$coefficients[1]),
                                                           xi=Migpd$models[[2]]$coefficients[2])

  I.proj <- data.frame(I.proj.x,I.proj.y)
  colnames(I.proj)<-c("x","y")

  if(PLOT==TRUE){
    x_lim_min_T<-ifelse(is.na(x_lim_min_T)==T,min(x_F),x_lim_min_T)
    x_lim_max_T<-ifelse(is.na(x_lim_max_T)==T,max(x_F),x_lim_max_T)
    y_lim_min_T<-ifelse(is.na(y_lim_min_T)==T,min(y_F),y_lim_min_T)
    y_lim_max_T<-ifelse(is.na(y_lim_max_T)==T,max(y_F),y_lim_max_T)

    x_lim_min<-ifelse(is.na(x_lim_min)==T,min(Data[,1]),x_lim_min)
    x_lim_max<-ifelse(is.na(x_lim_max)==T,max(Data[,1]),x_lim_max)
    y_lim_min<-ifelse(is.na(y_lim_min)==T,min(Data[,2]),y_lim_min)
    y_lim_max<-ifelse(is.na(y_lim_max)==T,max(Data[,2]),y_lim_max)

    par(mfrow=c(2,2))
    par(mar=c(4.2,4.2,2,2))
    plot(Chi)
    plot(x_F,y_F,xlim=c(x_lim_min_T,x_lim_max_T),ylim=c(y_lim_min_T,y_lim_max_T),xlab=paste("transformed ",names(Data)[1],sep=""),ylab=paste("transformed ",names(Data)[2],sep=""),main="Transformed scale")
    legend("topright",c("Base","Proj"),col=c(3,2),lty=1,bg="transparent",bty='n')
    lines(I.base_F$x,I.base_F$y,col=3,lwd=2)
    lines(I.proj_F$x,I.proj_F$y,col=2,lwd=2)

    plot(Data,xlim=c(x_lim_min,x_lim_max),ylim=c(y_lim_min,y_lim_max),main="Original scale")
    legend("topright",c("Base","Proj"),col=c(3,2),lty=1,bg="transparent",bty='n')
    lines(I.base$x,I.base$y,col=3,lwd=2)
    lines(I.proj$x,I.proj$y,col=2,lwd=2)
  }
  res <- list(Asym = Asym, Chi= CHI, n.bar=n.bar, p.base = p.base, p.proj=p.proj, I.base = I.base, I.proj = I.proj)
  return(res)
}
