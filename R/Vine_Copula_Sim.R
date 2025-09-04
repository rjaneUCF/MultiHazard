#' C and D-vine Copula - Simulation
#'
#' Simulating from specified C- and D-vine copula models. Function is a repackaging of the \code{RVineMatrix} and \code{RVineMatrix} functions from the \code{VineCopula} package into a single function.
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
#' #Fitting GPD to independent cluster maxima
#' S20.Migpd<-Migpd_Fit(Data=S20.Detrend.Declustered.df,
#'                      Data_Full=S20.Detrend.df,
#'                      mqu=c(0.975,0.975,0.9676))
#' #Fitting vine copula
#' S20.Vine<-Vine_Copula_Fit(Data=S20.Detrend.df)
#' #Simulating from fitted copula
#' S20.Vine.Sim<-Vine_Copula_Sim(Data=S20.Detrend.df,Vine_Model=S20.Vine,
#'                               Marginals=S20.Migpd,N=10)
#' #Plotting observed (black) and simulated (red) values
#' S20.Pairs.Plot.Data<-data.frame(rbind(na.omit(S20.Detrend.df[,-1]),S20.Vine.Sim$x.Sim),
#'                                 c(rep("Observation",nrow(na.omit(S20.Detrend.df))),
#'                                 rep("Simulation",nrow(S20.Vine.Sim$x.Sim))))
#' colnames(S20.Pairs.Plot.Data)<-c(names(S20.Detrend.df)[-1],"Type")
#' pairs(S20.Pairs.Plot.Data[,1:3],
#'       col=ifelse(S20.Pairs.Plot.Data$Type=="Observation","Black","Red"),
#'       upper.panel=NULL)
Vine_Copula_Sim<-function(Data,Vine_Model,Marginals,mu=365.25,N=10000){

  # Check: Data must be a data frame or matrix
  if (!is.data.frame(Data) && !is.matrix(Data)) {
    stop("Error: Data must be a data frame or matrix.")
  }

  # Check: At least 2 columns (to allow for time + variables or just variables)
  if (ncol(Data) < 2) {
    stop("Error: Data must have at least two columns.")
  }

  # Check: At least 5 rows
  if (nrow(Data) < 5) {
    stop("Error: Data must contain at least 5 rows.")
  }

  # Check: Vine_Model is a list with expected names
  expected_names <- c("Structure", "Family", "Par", "Par2")
  if (!is.list(Vine_Model) || !all(expected_names %in% names(Vine_Model))) {
    stop("Error: Vine_Model must be a list with names 'Structure', 'Family', 'Par', and 'Par2'.")
  }

  # Check: Structure is square matrix
  if (!is.matrix(Vine_Model$Structure) || nrow(Vine_Model$Structure) != ncol(Vine_Model$Structure)) {
    stop("Error: Vine_Model$Structure must be a square matrix.")
  }

  # Check: Family, Par, Par2 are same dimensions as Structure
  dims <- dim(Vine_Model$Structure)
  if (!all(dim(Vine_Model$Family) == dims,
           dim(Vine_Model$Par) == dims,
           dim(Vine_Model$Par2) == dims)) {
    stop("Error: Vine_Model$Family, Par, and Par2 must match dimensions of Structure.")
  }

  # Check: Marginals is a list with $models and $mqu
  if (!is.list(Marginals) || !("models" %in% names(Marginals)) || !("mqu" %in% names(Marginals))) {
    stop("Error: Marginals must be a list containing 'models' and 'mqu'.")
  }

  # Check: Marginals$models length matches number of variables
  n.vars <- if (inherits(Data[,1], c("Date", "POSIXct", "factor", "character"))) ncol(Data) - 1 else ncol(Data)
  if (length(Marginals$models) != n.vars) {
    stop("Error: Marginals$models must have the same length as number of variables in Data.")
  }

  # Check: mu is numeric and positive
  if (!is.numeric(mu) || length(mu) != 1 || mu <= 0) {
    stop("Error: mu must be a single positive number.")
  }

  # Check: N is integer and positive
  if (!is.numeric(N) || length(N) != 1 || N <= 0 || N != as.integer(N)) {
    stop("Error: N must be a single positive integer.")
  }

  #Number of extreme events
  No.events<-round(mu*N,0)

  if(inherits(Data[,1], c("Date", "factor"))){
    #Simulating from copula on the transformed scale
    RMV<- RVineMatrix(Matrix = Vine_Model$Structure, family = Vine_Model$Family,
                      par = Vine_Model$Par, par2 = Vine_Model$Par2)
    u<-RVineSim(No.events, RMV)
    colnames(u)<-names(Data[2:ncol(Data)])
    #Transforming d-dimensional simulations to the original scale using the marginal distribution supplied
    x<-matrix(0,nrow=nrow(u),ncol=ncol(u))
    for(i in 1:(ncol(Data)-1)){
      x[,i]<-as.numeric(quantile(na.omit(Data[,(i+1)]),u[,i]))
      x[which(x[,i]>(Marginals$models[i][[1]]$threshold)),i]<-u2gpd(u[which(x[,i]>Marginals$models[i][[1]]$threshold),i], p = 1-Marginals$mqu, th=Marginals$models[i][[1]]$threshold, sigma=exp(Marginals$models[i][[1]]$par[1]),xi=Marginals$models[i][[1]]$par[2])
    }
    x<-data.frame(x)
    colnames(x)<-names(Data[2:ncol(Data)])
  } else {
    #Simulating from copula on the transformed scale
    RMV <- RVineMatrix(Matrix = Vine_Model$Structure, family = Vine_Model$Family,
                       par = Vine_Model$Par, par2 = Vine_Model$Par2)
    u<-RVineSim(No.events,RMV)
    colnames(u)<-names(Data)
    #Transforming d-dimensional simulations to the original scale using the marginal distribution supplied
    x<-matrix(0,nrow=nrow(u),ncol=ncol(u))
    for(i in 1:(ncol(Data))){
      x[,i]<-as.numeric(quantile(na.omit(Data[,(i)]),u[,i]))
      x[which(x[,i]>(Marginals$models[i][[1]]$threshold)),i]<-u2gpd(u[which(x[,i]>Marginals$models[i][[1]]$threshold),i], p = 1-Marginals$mqu, th=Marginals$models[i][[1]]$threshold, sigma=exp(Marginals$models[i][[1]]$par[1]),xi=Marginals$models[i][[1]]$par[2])
    }
    x<-data.frame(x)
    colnames(x)<-names(Data)
  }
  #Compose list of outputs
  res<-list("Vine_Model" = Vine_Model, "u.Sim"=u,"x.Sim"=x)
  return(res)
}
