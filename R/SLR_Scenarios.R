#' Sea level rise scenarios
#'
#' Time (in years) for a specified change in sea level according to various sea level projections. Contained within the function are: (1) the three scenarios for Key West in the Southeast Florida Regional Climate Change Compact, (2) those for Miami Beach in "Global and Regional Sea Level Rise Scenarios for the United States" NOAA et al. (2017) and (3) those in the Interagency Sea Level Rise Scenario Tool (NOAA et al. 2022) for Naples and Miami Beach. Users can also input scenarios of their choice.
#'
#' @param SeaLevelRise Numeric vector of length one, specifying the sea level rise required.
#' @param Scenario Character vector of length one, specifying the sea level rise scenarios to be adopted. Options are \code{"Compact"} for those for Key West in the Southeast Florida Regional Climate Change Compact, \code{"NOAA2017"} for those in "Global and Regional Sea Level Rise Scenarios for the United States" at Miami Beach used in Jane et al. (2020), \code{"NOAA2022"} for those for Miami Beach and Naples in the Interagency Sea Level Rise Scenario Tool, or \code{NA} if a set of scenarios are specified by the user (see \code{New_Scenario}).
#' @param Unit Character vector of length one, specifying units of \code{SeaLevelRise}. Options are meters \code{m} and Inches \code{"Inches"}. Default is \code{"m"}.
#' @param Year Numeric vector of length one, specifying the current year. Default is \code{2022}.
#' @param Location Character vector of length one, specifying the location associated with the scenarios. Projections for  \code{"Key West"} (Compact), \code{"Miami Beach"} (NOAA2017 AND NOAA2022) and \code{"Naples"} (NOAA2022) are contained within the package. If a user specified scenarios are employed, set to the name of the site. Default is \code{"Key West"}.
#' @param New_Scenario Dataframe containing sea level rise scenarios. First column must be a year and the scenarios provided in the remaining columns. For the color scale to correlate with the severity of the scenarios they should be listed from most to least severe i.e., the highest SLR scenario should appear in column 2. All entries must be numeric.
#' @return For \code{"Compact"}, \code{"NOAA2017"} and \code{"NOAA2022"} a list length of time for \code{SeaLevelRise} of sea level rise is expected to arise under the \code{High}, \code{Intermediate} and \code{Low}. For user specified scenarios, the time for \code{SeaLevelRise} to occur under each is returned as \code{SLR_Year}. Upper panel: A plot of the scenarios. Scenarios are in bold until the time the SeaLevelRise is reached and are transparent thereafter. Lower panel: A plot showing the number of years before is expected to occur.
#' @examples
#' #Calculate the estimated time required for 0.45m of SLR in Key West according to the scenarios
#' in the Southeast Florida Regional Climate Change Compact
#' SLRScenarios(0.45)
#' #Calculate the estimated time required for 0.8 inches of SLR in Naples according
#' to the scenarios in the 2022 Interagency Sea Level Rise Scenario Tool
#' SLRScenarios(0.45,Scenario="NOAA2022", Unit = "Inches", Location="Naples")
#' #Read in the scenarios for Fort Myers downloaded
#' from https://sealevel.nasa.gov/task-force-scenario-tool/?psmsl_id=1106
#' SeaLevelRise.2022<-read.csv("sl_taskforce_scenarios_psmsl_id_1106_Fort_Myers.csv")
#' #Convert data to the appropriate format for the SLRScenarios function
#' #i.e. first column years, following columns the scenarios most to least extreme,
#' converted from millimeters to meters
#' SeaLevelRise.2022_input<-data.frame(Year=seq(2020,2150,10),
#'                                     "High"=as.numeric(SeaLevelRise.2022[14,-(1:5)])/1000,
#'                                     "Medium"=as.numeric(SeaLevelRise.2022[8,-(1:5)])/1000,
#'                                     "Low"=as.numeric(SeaLevelRise.2022[2,-(1:5)])/1000)
#' #Calculate the estimated time required for 0.8 inches of SLR at Fort Myers
#' SLR_Scenarios(SeaLevelRise=0.8, Scenario="Other", Unit = "m", Year=2022,
#'               Location="Fort Myers", New_Scenario=SeaLevelRise.2022_input)
SLR_Scenarios<-function(SeaLevelRise, Scenario="Compact", Unit = "m", Year=2022, Location="Key West", New_Scenario=NA){

 #Colors
 mypalette<-brewer.pal(7,"Set1")

 #Layout of graphical output
 layout(matrix(c(1,1,2), 3, 1, byrow = TRUE))
 par(mar=c(4.2,6,1,1))

 #Plotting the SLR scenarios
 y_axis<-ifelse(Unit=="m",paste('Projected sea level rise at', Location, '(m)'),past('Projected sea level rise at', Location,'(ft)'))

 if(Scenario=="Compact"){
  #Composing the AR5 scenario
  SeaLevelRise.AR5<-data.frame("Year"=1:4,"Medium"=1:4)
  SeaLevelRise.AR5$Year<-c(1992,2030,2060,2100)
  SeaLevelRise.AR5$Int<-0.0254*c(0,6,14,31)
  #Using splines to interpolate the scenarios
  spline.high<-spline(NOAAetal2012$Year,ifelse(Unit=="m",1,0.00328084)*NOAAetal2012$High, xout=seq(Year,2100,0.25))
  spline.high$y<-spline.high$y-spline.high$y[1]
  spline.int<-spline(USACE2013$Year,ifelse(Unit=="m",1,0.00328084)*USACE2013$High, xout=seq(Year,2100,0.25))
  spline.int$y<-spline.int$y-spline.int$y[1]
  spline.low<-spline(SeaLevelRise.AR5$Year,ifelse(Unit=="m",1,0.00328084)*SeaLevelRise.AR5$Int, xout=seq(Year,2100,0.25))
  spline.low$y<-spline.low$y-spline.low$y[1]
 }
 res<-("x"=spline.high$x,"y"=spline.high$y)
 return(res)
}

