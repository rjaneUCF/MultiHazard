#' C and D-vine Copula - Simulation
#'
#' Simulating from specified C- and D-vine copula models. Builds on the \code{CDVineSim} in \code{CDVine}.
#'
#' @param Data Data frame containing \code{n} at least partially concurrent time series. First column may be a \code{"Date"} object. Can be \code{Dataframe_Combine} output.
#' @param Vine_Model An \code{RVineMatrix} object i.e., output of \code{Vine_Copula_Fit} specifying the structure and copula families composing the vine copula.
#' @param Marginals An \code{migpd} object containing the d-independent generalized Pareto models.
#' @param mu (average) Number of events per year. Numeric vector of length one. Default is 365.25, daily data.
#' @param N Number of years worth of extremes to be simulated. Numeric vector of length one. Default 10,000 (years).
#' @return List comprising an integer vector specifying the pair-copula families composing the C- or D-vine copula \code{Vine_family}, its parameters \code{Vine_par} and \code{Vine_par2} and type of regular vine \code{Vine_Type}. In addition, data frames of the simulated observations: \code{u.Sim} on the transformed \code{$[0,1]^n$} and \code{x.Sim} the original scales.
#' @seealso \code{\link{Vine_Copula_Fit}}
#' @export
#' @examples
#' #Fitting vine copula
#' S20.Vine<-Vine_Copula_Fit(Data=S20.Detrend.df)
#' #Simulating from fitted copula
#' S20.Vine.Sim<-Vine_Copula_Sim(Data=S20.Detrend.df,Vine_Model=S20.Vine,
#'                               Marginals=S20.Migpd,N=10)
#' #Plotting observed (black) and simulated (red) values
#' S20.Pairs.Plot.Data<-data.frame(rbind(na.omit(S20.Detrend.df[,-1]),S22.Vine.Sim$x.Sim),
#'                                 c(rep("Observation",nrow(na.omit(S20.Detrend.df))),
#'                                 rep("Simulation",nrow(S20.Vine.Sim$x.Sim))))
#' colnames(S20.Pairs.Plot.Data)<-c(names(S20.Detrend.df)[-1],"Type")
#' pairs(S20.Pairs.Plot.Data[,1:3],
#'       col=ifelse(S20.Pairs.Plot.Data$Type=="Observation","Black","Red"),
#'       upper.panel=NULL)
Vine_Copula_Sim<-function(Data,Vine_Model,Marginals,mu=365.25,N=10000){

  #Number of extreme events
  No.events<-mu*N

  if(class(Data[,1])=="Date" | class(Data[,1])=="factor"){
  #Simulating from copula on the transformed scale
  RMV<- RVineMatrix(Matrix = Vine_Model$Structure, family = Vine_Model$Family,
                    par = Vine_Model$Par, par2 = Vine_Copula$Par2)
  u<-RVineSim(No.events, RMV)
  colnames(u)<-names(Data[2:ncol(Data)])
  #Transforming d-dimensional simulations to the origional scale using the marginal distribution supplied
  x<-matrix(0,nrow=nrow(u),ncol=ncol(u))
  for(i in 1:(ncol(Data)-1)){
   x[,i]<-as.numeric(quantile(na.omit(Data[,(i+1)]),u[,i]))
   x[which(x[,i]>(Marginals$models[i][[1]]$threshold)),i]<-u2gpd(u[which(x[,i]>Marginals$models[i][[1]]$threshold),i], p = Marginals$models[[i]]$rate, th=Marginals$models[i][[1]]$threshold, sigma=exp(Marginals$models[i][[1]]$par[1]),xi=Marginals$models[i][[1]]$par[2])
  }
  x<-data.frame(x)
  colnames(x)<-names(Data[2:ncol(Data)])
  } else {
  #Simulating from copula on the transformed scale
  RMV <- RVineMatrix(Matrix = Vine_Model$Structure, family = Vine_Model$Family,
                     par = Vine_Model$Par, par2 = Vine_Model$Par2)
  u<-RVineSim(No.events,RMV)
  colnames(u)<-names(Data)
  #Transforming d-dimensional simulations to the origional scale using the marginal distribution supplied
  x<-matrix(0,nrow=nrow(u),ncol=ncol(u))
   for(i in 1:(ncol(Data))){
     x[,i]<-as.numeric(quantile(na.omit(Data[,(i)]),u[,i]))
     x[which(x[,i]>(Marginals$models[i][[1]]$threshold)),i]<-u2gpd(u[which(x[,i]>Marginals$models[i][[1]]$threshold),i], p = Marginals$models[[i]]$rate, th=Marginals$models[i][[1]]$threshold, sigma=exp(Marginals$models[i][[1]]$par[1]),xi=Marginals$models[i][[1]]$par[2])
   }
  x<-data.frame(x)
  colnames(x)<-names(Data)
  }
 #Compose list of outputs
 res<-list("Vine_Model" = Vine_Model, "u.Sim"=u,"x.Sim"=x)
 return(res)
}
