#' C and D-vine Copula - Simulation
#'
#' Simulating from specified C- and D-vine copula models. Builds on the \code{CDVineSim} in \code{CDVine}.
#'
#' @param Data Dataframe containing \code{n} at least partially concurrent time series. First column may be a \code{"Date"} object. Can be \code{Dataframe_Combine} output.
#' @param Marginals An \code{migpd} object containing the d-independent generalized Pareto models.
#' @param Vine_family A n*(n-1)/2 integer vector specifying the pair-copula families defining the fitting C- or a D-vine copula models. Can be specified as the \code{Family} agument of a \code{Vine_Copula_Fit} object. See help file of the \code{CDVineSim} function to find out copula represented by integers 0-40.
#' @param Vine_par A n*(n-1)/2 vector of pair-copula parameters.
#' @param Vine_par2 A n*(n-1)/2 vector of second parameters for pair-copula families with two parameters.
#' @param Vine_Type Type of the vine model: \itemize{
#'   \item 1 or "CVine" = C-vine
#'   \item 2 or "DVine" = D-vine
#' } Can be specified as the \code{Type} argument of a \code{Vine_Copula_Fit} object.
#' @param mu (average) Number of events per year. Numeric vector of length one. Default is 365.25, daily data.
#' @param N Number of years worth of extremes to be simulated. Numeric vector of length one. Default 10,000 (years).
#' @return List comprising an integer vector specifing the pair-copula families composing the C- or D-vine copula \code{Vine_family}, its paraeters \code{Vine_par} and \code{Vine_par2} and type of regular vine \code{Vine_Type}. In addition, dataframes of the simulated observations: \code{u.Sim} on the transformed \code{[0,1]^n} and \code{x.Sim} the origional scales.
#' @seealso Detrend_Declustered_Combine CD_Vine_Select migpd.edit
#' @export
#' @examples
#' S22.Vine.Sim<-Vine_Copula_Sim(Data=S22.Detrend.df,Marginals=S22.GPD,Vine_family=S22.Vine$Family, Vine_par=S22.Vine$Par, Vine_par2=S22.Vine$Par2, Vine_Type="DVine",N=10)
#'
#' S22.Pairs.Plot.Data<-data.frame(rbind(na.omit(S22.Detrend.df[,-1]),S22.Vine.Sim$x.Sim),c(rep("Observation",nrow(na.omit(S22.Detrend.df))),rep("Simulation",nrow(S22.Vine.Sim$x.Sim))))
#' colnames(S22.Pairs.Plot.Data)<-c(names(S22.Detrend.df)[-1],"Type")
#' pairs(S22.Pairs.Plot.Data[,1:3],col=ifelse(S22.Pairs.Plot.Data$Type=="Observation","Black","Red"),upper.panel=NULL)

Vine_Copula_Sim<-function(Data,Marginals,Vine_family, Vine_par, Vine_par2, Vine_Type="DVine",mu=365.25,N=10000){

  #Number of extreme events
  No.events<-mu*N

  if(class(Data[,1])=="Date" | class(Data[,1])=="factor"){
  #Simulating from copula on the transformed scale
  u<-CDVineSim(No.events, family=Vine_family, par=Vine_par, par2=Vine_par2, type="CVine")
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
  u<-CDVineSim(No.events, family=Vine_family, par=Vine_par, par2=Vine_par2, type="CVine")
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
 res<-list("Vine_family" = Vine_family, "Vine_par"= Vine_par, "Vine_par2" = Vine_par2, "Vine_Type" = Vine_Type,"u.Sim"=u,"x.Sim"=x)
 return(res)
}
