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
 if(Scenario=="NOAA2017"){
   spline.high<-spline(NOAAetal2017$Year,ifelse(Unit=="m",1,0.00328084)*NOAAetal2017$High, xout=seq(Year,2100,0.25))
   spline.high$y<-spline.high$y-spline.high$y[1]
   spline.int<-spline(NOAAetal2017$Year,ifelse(Unit=="m",1,0.00328084)*NOAAetal2017$Intermediate, xout=seq(Year,2100,0.25))
   spline.int$y<-spline.int$y-spline.int$y[1]
   spline.low<-spline(NOAAetal2017$Year,ifelse(Unit=="m",1,0.00328084)*NOAAetal2017$Low, xout=seq(Year,2100,0.25))
   spline.low$y<-spline.low$y-spline.low$y[1]
 }

 if(Scenario=="NOAA2022"){
   NOAAetal2022.Year<-c(2020,2030,2040,2050,2060,2070,2080,2090,2100)
   x<-ifelse(Location=="Miami Beach",1,2)
   spline.high<-spline(NOAAetal2022.Year,ifelse(Unit=="m",1,0.00328084)*NOAAetal2022[c(29,14)[x],6:14]/1000, xout=seq(Year,2100,0.25))
   spline.high$y<-spline.high$y-spline.high$y[1]
   spline.int<-spline(NOAAetal2022.Year,ifelse(Unit=="m",1,0.00328084)*NOAAetal2022[c(23,8)[x],6:14]/1000, xout=seq(Year,2100,0.25))
   spline.int$y<-spline.int$y-spline.int$y[1]
   spline.low<-spline(NOAAetal2022.Year,ifelse(Unit=="m",1,0.00328084)*NOAAetal2022[c(17,2)[x],6:14]/1000, xout=seq(Year,2100,0.25))
   spline.low$y<-spline.low$y-spline.low$y[1]
 }

 if(Scenario=="Other"){
   if((ncol(New_Scenario)-1)==5){
     mypalette<-c("Dark Red",brewer.pal(7,"Set1")[1:3],"Dark Green")
   }
   max.Year<-max(New_Scenario[,1])
   splines<-vector("list",ncol(New_Scenario)-1)
   for(i in 1:(ncol(New_Scenario)-1)){
     splines[[i]]<-spline(New_Scenario[,1],ifelse(Unit=="m",1,0.00328084)*New_Scenario[,i+1], xout=seq(Year,max.Year,0.25))
   }
   max.SLR<-max(unlist(lapply(splines, `[[`, 2)),na.rm=T)
   min.SLR<-min(0,min(unlist(lapply(splines, `[[`, 2)),na.rm=T))
   plot(0,xlab="Year",ylab=y_axis,type='n',xlim=c(Year,max.Year),ylim=c(min.SLR,ifelse(Unit=="m",1,3.28084)*max.SLR*1.25),cex.lab=1.5,cex.axis=1.5)
   for(i in 1:(ncol(New_Scenario)-1)){
     lines(splines[[i]]$x,splines[[i]]$y,col=alpha(mypalette[i],0.2),lwd=5)
     max<-which(abs(splines[[i]]$y-SeaLevelRise)==min(abs(splines[[i]]$y-SeaLevelRise)))
     lines(splines[[i]]$x[1:max],splines[[i]]$y[1:max],col=mypalette[i],lwd=5)
   }
 }

 if(Scenario=="Compact" | Scenario=="NOAA2017" | Scenario=="NOAA2022"){
   plot(0,xlab="Year",ylab=y_axis,type='n',xlim=c(Year,2100),ylim=c(0,ifelse(Unit=="m",1,3.28084)*max(spline.high$y,spline.int$y,spline.low$y)*1.5),cex.lab=1.5,cex.axis=1.5)
   lines(spline.high$x,spline.high$y,col=alpha(mypalette[1],0.2),lwd=5)
   lines(spline.int$x,spline.int$y,col=alpha(mypalette[2],0.2),lwd=5)
   lines(spline.low$x,spline.low$y,col=alpha(mypalette[3],0.2),lwd=5)
   #Scenarios are bold up until the required SeaLevelRise is expected to occur
   max.high<-which(abs(spline.high$y-SeaLevelRise)==min(abs(spline.high$y-SeaLevelRise)))
   lines(spline.high$x[1:max.high],spline.high$y[1:max.high],col=mypalette[1],lwd=5)
   max.int<-which(abs(spline.int$y-SeaLevelRise)==min(abs(spline.int$y-SeaLevelRise)))
   lines(spline.int$x[1:max.int],spline.int$y[1:max.int],col=mypalette[2],lwd=5)
   max.low<-which(abs(spline.low$y-SeaLevelRise)==min(abs(spline.low$y-SeaLevelRise)))
   lines(spline.low$x[1:max.low],spline.low$y[1:max.low],col=mypalette[3],lwd=5)
 }

 #Plotting the estimated time for the required SeaLevelRise to occur
 par(las=1)

 if(Scenario=="Compact"){
   plot(0,xlab="Number of years",ylab="",type='n',xlim=c(Year,2100),ylim=c(0,3),cex.lab=1.5,cex.axis=1.5,yaxt="n",xaxt="n",bty="n")
   axis(1,at=seq(Year,2100,20),seq(0,2100-Year,20),cex.axis=1.5)
   mtext(c("NOAA et al. (2012)","High","USACE 2013","High","IPCC AR5","Medium"),2,-4.15,at=c(2.4,2.2,1.4,1.2,0.4,0.2))
 }

 if(Scenario=="NOAA2017"){
   plot(0,xlab="Number of years",ylab="",type='n',xlim=c(Year,2100),ylim=c(0,3),cex.lab=1.5,cex.axis=1.5,yaxt="n",xaxt="n",bty="n")
   axis(1,at=seq(Year,2100,20),seq(0,2100-Year,20),cex.axis=1.5)
   mtext(c("High","Intermediate","Low"),2,-4.15,at=c(2.25,1.25,0.25))
 }

 if(Scenario=="NOAA2022"){
   plot(0,xlab="Number of years",ylab="",type='n',xlim=c(Year,2100),ylim=c(0,3),cex.lab=1.5,cex.axis=1.5,yaxt="n",xaxt="n",bty="n")
   axis(1,at=seq(Year,2100,20),seq(0,2100-Year,20),cex.axis=1.5)
   mtext(c("High","Intermediate","Low"),2,-4.15,at=c(2.25,1.25,0.25))
 }

 if(Scenario=="Compact" | Scenario=="NOAA2017" | Scenario=="NOAA2022"){
   rect(Year,2,spline.int$x[1:max.high],2.5,col=mypalette[1],border=NA)
   High<-spline.high$x[max.high]
   if(High>2100){
    text(spline.high$x[max.high]+1.5,2.25,"> 80",cex=1.5,font=3)
   } else{
    text(spline.high$x[max.high]+1.5,2.25,paste(High-Year),cex=1.5,font=3)
   }
   rect(Year,1,spline.int$x[1:max.int],1.5,col=mypalette[2],border=NA)
   Intermediate<-spline.int$x[max.int]
   if(Intermediate>2100){
    text(spline.int$x[max.int]+1.5,1.25,"> 80",cex=1.5,font=3)
   } else{
     text(spline.int$x[max.int]+1.5,1.25,paste(Intermediate-Year),cex=1.5,font=3)
   }
   rect(Year,0,spline.low$x[1:max.low],0.5,col=mypalette[3],border=NA)
   Low<-spline.low$x[max.low]
   if(Low>2100){
     text(spline.low$x[max.low]+1.5,0.25,"> 80",cex=1.5,font=3)
   } else{
    text(spline.low$x[max.low]+1.5,0.25,paste(Low-Year),cex=1.5,font=3)
   }
   res<-list("High" = High, "Intermediate" = Intermediate, "Low" = Low)
 }

 if(Scenario=="Other"){
  SLR_Year<-numeric(ncol(New_Scenario)-1)
  plot(0,xlab="Number of years",ylab="",type='n',xlim=c(Year,max.Year),ylim=c(0,ncol(New_Scenario)-1),cex.lab=1.5,cex.axis=1.5,yaxt="n",xaxt="n",bty="n")
  axis(1,at=seq(Year,max.Year,20),seq(0,max.Year-Year,20),cex.axis=1.5)
  mtext(colnames(New_Scenario)[-1],2,ifelse((max.Year-Year)<100,-4.15,-1.1),at=rev(seq(0,ncol(New_Scenario)-2,1))+0.25)
  for(i in 1:(ncol(New_Scenario)-1)){
   j<-rev(1:(ncol(New_Scenario)-1))[i]
   max<-which(abs(splines[[j]]$y-SeaLevelRise)==min(abs(splines[[j]]$y-SeaLevelRise)))
   rect(Year,seq(0,ncol(New_Scenario)-2,1)[i],splines[[i]]$x[1:max],seq(0,ncol(New_Scenario)-2,1)[i]+0.5,col=mypalette[j],border=NA)
   SLR_Year[i]<-round(splines[[i]]$x[max],0)
   if((max.Year-Year)>100){
    text(splines[[i]]$x[max]+2.6,seq(0,ncol(New_Scenario)-2,1)[i]+0.25,paste('>',SLR_Year[i]-Year),cex=1.5,font=3)
   } else{
    text(splines[[i]]$x[max]+1.5,seq(0,ncol(New_Scenario)-2,1)[i]+0.25,paste(SLR_Year[i]-Year),cex=1.5,font=3)
   }
  }
  res<-list("SLR_Year" = SLR_Year)
 }

 return(res)
}


