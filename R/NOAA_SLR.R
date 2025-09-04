#' NOAA sea-level rise scenarios
#'
#' Time (in years) for a specified amount of sea-level rise (SLR) to occur at Miami Beach according to the five SLR scenarios in NOAA 2017 report titled "Global and Regional Sea Level Rise Scenarios for the United States".
#'
#' @param OsWL_req Numeric vector of SLR required.
#' @param SLR_scen Character vector specifying which of the NOAA (2017) scenarios to consider. Options include \code{High}, Intermediate high \code{Int.High}, \code{Intermediate}, Intermediate low (\code{Int.Low}) and \code{Low}.
#' @param Input_unit Character vector of length one; specifying units of SLR. Default is meters \code{"m"}, other option is feet \code{"ft"}.
#' @param Year.Inital Character vector of length one; specifying the current year.
#' @return List comprising the specified \code{Threshold} as the quantile of the conditioning variable above which declustered excesses are paired with co-occurrences of the other variable, the resulting two-dimensional sample \code{data} and \code{name} of the conditioning variable.
#' @export
#' @examples
#' NOAA_SLR(OsWL_req=seq(0,1,0.01),SLR_scen = c("High","Intermediate","Low"),
#'                    Input_unit="m")
NOAA_SLR<-function(OsWL_req,SLR_scen = c("High","Intermediate","Low"),Input_unit="m",Year.Inital=2020){
  res<-matrix(0,nrow=length(OsWL_req),ncol=length(SLR_scen))
  scenario<-which(colnames(NOAA2017)==SLR_scen)
  for(j in 1:length(scenario)){
    sp.Initial<-spline(NOAA2017$Year,ifelse(Input_unit=="m",1,3.28084)*NOAA2017[,scenario[j]], n = 201)
    sp<-spline(NOAA2017$Year[-which(NOAA2017$Year<Year.Initial)],ifelse(Input_unit=="m",1,3.28084)*(NOAA2017[,scenario[j]][-which(NOAA2017$Year<Year.Initial)]-sp.Initial$y[which(sp.Initial$x==Year.Initial)]), n = 201)
    plot(sp$x,sp$y)
    for(i in 1:length(OsWL_req)){
      res[i,j]<-sp$x[min(which((sp$y-OsWL_req[i])>0))]-Year.Initial
    }
  }
  return(res)
}
