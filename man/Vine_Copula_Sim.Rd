% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Vine_Copula_Sim.R
\name{Vine_Copula_Sim}
\alias{Vine_Copula_Sim}
\title{C and D-vine Copula - Simulation}
\usage{
Vine_Copula_Sim(Data, Marginals, Vine_family, Vine_par, Vine_par2,
  Vine_Type = "DVine", mu = 365.25, N = 10000)
}
\arguments{
\item{Data}{Data frame containing \code{n} at least partially concurrent time series. First column may be a \code{"Date"} object. Can be \code{Dataframe_Combine} output.}

\item{Marginals}{An \code{migpd} object containing the d-independent generalized Pareto models.}

\item{Vine_family}{A n*(n-1)/2 integer vector specifying the pair-copula families defining the fitting C- or a D-vine copula models. Can be specified as the \code{Family} argument of a \code{Vine_Copula_Fit} object. See help file of the \code{CDVineSim} function to find the integers representing the different copula families.}

\item{Vine_par}{A n*(n-1)/2 vector of pair-copula parameters.}

\item{Vine_par2}{A n*(n-1)/2 vector of second parameters for pair-copula families with two parameters.}

\item{Vine_Type}{Type of the vine model: \itemize{
  \item 1 or "CVine" = C-vine
  \item 2 or "DVine" = D-vine
} Can be specified as the \code{Type} argument of a \code{Vine_Copula_Fit} object.}

\item{mu}{(average) Number of events per year. Numeric vector of length one. Default is 365.25, daily data.}

\item{N}{Number of years worth of extremes to be simulated. Numeric vector of length one. Default 10,000 (years).}
}
\value{
List comprising an integer vector specifying the pair-copula families composing the C- or D-vine copula \code{Vine_family}, its parameters \code{Vine_par} and \code{Vine_par2} and type of regular vine \code{Vine_Type}. In addition, data frames of the simulated observations: \code{u.Sim} on the transformed \code{$[0,1]^n$} and \code{x.Sim} the original scales.
}
\description{
Simulating from specified C- and D-vine copula models. Builds on the \code{CDVineSim} in \code{CDVine}.
}
\examples{
#Fitting vine copula
S20.Vine<-Vine_Copula_Fit(Data=S20.Detrend.df, FamilySet=NA,
                          Type="DVine", SelCrit="AIC",Indeptest=FALSE,
                          Level=0.05)
#Simulating from fitted copula
S20.Vine.Sim<-Vine_Copula_Sim(Data=S20.Detrend.df,Marginals=S20.Migpd,
                              Vine_family=S20.Vine$Family, Vine_par=S20.Vine$Par,
                              Vine_par2=S20.Vine$Par2, Vine_Type="DVine",N=10)
#Plotting observed (black) and simulated (red) values
S20.Pairs.Plot.Data<-data.frame(rbind(na.omit(S20.Detrend.df[,-1]),S22.Vine.Sim$x.Sim),
                                c(rep("Observation",nrow(na.omit(S20.Detrend.df))),
                                rep("Simulation",nrow(S20.Vine.Sim$x.Sim))))
colnames(S20.Pairs.Plot.Data)<-c(names(S20.Detrend.df)[-1],"Type")
pairs(S20.Pairs.Plot.Data[,1:3],
      col=ifelse(S20.Pairs.Plot.Data$Type=="Observation","Black","Red"),
      upper.panel=NULL)
}
\seealso{
\code{\link{Vine_Copula_Fit}}
}
