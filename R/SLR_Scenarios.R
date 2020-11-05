#' Sea level rise scenarios in the Southeast Florida Regional Climate Change Compact:
#'
#' Calculates and plots time required for sea level rise to reach a specified level according to the three scenarios in the Compact.
#'
#' @param data A data frame with \code{n} columns, each comprising a declustered and if necessary detrended time series to be modelled.
#' @param SeaLevelRise Numeric vector of length one, sea level rise required.
#' @return An object of class \code{"migpd"}. There are \code{coef}, \code{print}, \code{plot}, \code{ggplot} and \code{summary} functions available.
#' @export
#' @examples
#' SLRScenarios(0.45)
SLR_Scenarios<-function(SeaLevelRise, Unit = "m"){
mypalette<-brewer.pal(7,"Set1")
layout(matrix(c(1,1,2), 3, 1, byrow = TRUE))
par(mar=c(4.2,6,1,1))

SeaLevelRise.AR5<-data.frame("Year"=1:4,"Medium"=1:4)
SeaLevelRise.AR5$Year<-c(1992,2030,2060,2100)
SeaLevelRise.AR5$Int<-0.0254*c(0,6,14,31)
SeaLevelRise.AR5$Int<-SeaLevelRise.AR5$Int-spline(SeaLevelRise.AR5$Year,SeaLevelRise.AR5$Int,n=109)$y[which(spline(SeaLevelRise.AR5$Year,SeaLevelRise.AR5$Int,n=109)$x==2019)]


y_axis<-ifelse(Unit=="m","Projected sea level rise at Miami Beach (m)","Projected sea level rise at Miami Beach (ft)")
plot(0,xlab="Year",ylab=y_axis,type='n',xlim=c(2020,2100),ylim=c(0,ifelse(Unit=="m",1,3.28084)*2.5),cex.lab=1.5,cex.axis=1.5)
lines(spline(SeaLevelRise.AR5$Year,ifelse(Unit=="m",1,3.28084)*SeaLevelRise.AR5$Int)$x[-which(spline(SeaLevelRise.AR5$Year,ifelse(Unit=="m",1,3.28084)*SeaLevelRise.AR5$Int)$x<2020)],
      spline(SeaLevelRise.AR5$Year,ifelse(Unit=="m",1,3.28084)*SeaLevelRise.AR5$Int)$y[-which(spline(SeaLevelRise.AR5$Year,ifelse(Unit=="m",1,3.28084)*SeaLevelRise.AR5$Int)$x<2020)],col=alpha(mypalette[3],0.2),lwd=5)
lines(spline(SeaLevelRise.USACE2013$Year[-which(SeaLevelRise.USACE2013$Year<2020)],ifelse(Unit=="m",1,3.28084)*SeaLevelRise.USACE2013$High[-which(SeaLevelRise.USACE2013$Year<2020)]),col=alpha(mypalette[2],0.2),lwd=5)
lines(spline(SeaLevelRise.NOAA2012$Year[-which(SeaLevelRise.NOAA2012$Year<2020)],ifelse(Unit=="m",1,3.28084)*SeaLevelRise.NOAA2012$High[-which(SeaLevelRise.NOAA2012$Year<2020)]),col=alpha(mypalette[1],0.2),lwd=5)




lines(spline(SeaLevelRise.AR5$Year,SeaLevelRise.AR5$Int, n = 109)$x[29:(which(abs(spline(SeaLevelRise.AR5$Year,SeaLevelRise.AR5$Int, n = 109)$y-SeaLevelRise)==min(abs(spline(SeaLevelRise.AR5$Year,SeaLevelRise.AR5$Int, n = 109)$y-SeaLevelRise))))],
      (spline(SeaLevelRise.AR5$Year,SeaLevelRise.AR5$Int, n = 109)$y[29:(which(abs(spline(SeaLevelRise.AR5$Year,SeaLevelRise.AR5$Int, n = 109)$y-SeaLevelRise)==min(abs(spline(SeaLevelRise.AR5$Year,SeaLevelRise.AR5$Int, n = 109)$y-SeaLevelRise))))])*ifelse(Unit=="m",1,3.28084),col=mypalette[3],lwd=5)
lines(spline(SeaLevelRise.USACE2013$Year[-which(SeaLevelRise.USACE2013$Year<2020)],ifelse(Unit=="m",1,3.28084)*SeaLevelRise.USACE2013$High[-which(SeaLevelRise.USACE2013$Year<2020)], n = 100)$x[1:(which(abs(spline(SeaLevelRise.USACE2013$Year[-which(SeaLevelRise.USACE2013$Year<2020)],ifelse(Unit=="m",1,3.28084)*SeaLevelRise.USACE2013$High[-which(SeaLevelRise.USACE2013$Year<2020)], n = 100)$y-SeaLevelRise)==min(abs(spline(SeaLevelRise.USACE2013$Year[-which(SeaLevelRise.USACE2013$Year<2020)],ifelse(Unit=="m",1,3.28084)*SeaLevelRise.USACE2013$High[-which(SeaLevelRise.USACE2013$Year<2020)], n = 100)$y-SeaLevelRise))))],
      spline(SeaLevelRise.USACE2013$Year[-which(SeaLevelRise.USACE2013$Year<2020)],ifelse(Unit=="m",1,3.28084)*SeaLevelRise.USACE2013$High[-which(SeaLevelRise.USACE2013$Year<2020)], n = 100)$y[1:(which(abs(spline(SeaLevelRise.USACE2013$Year[-which(SeaLevelRise.USACE2013$Year<2020)],ifelse(Unit=="m",1,3.28084)*SeaLevelRise.USACE2013$High[-which(SeaLevelRise.USACE2013$Year<2020)], n = 100)$y-SeaLevelRise)==min(abs(spline(SeaLevelRise.USACE2013$Year[-which(SeaLevelRise.USACE2013$Year<2020)],ifelse(Unit=="m",1,3.28084)*SeaLevelRise.USACE2013$High[-which(SeaLevelRise.USACE2013$Year<2020)], n = 100)$y-SeaLevelRise))))],col=mypalette[2],lwd=5)
lines(spline(SeaLevelRise.NOAA2012$Year[-which(SeaLevelRise.NOAA2012$Year<2020)],ifelse(Unit=="m",1,3.28084)*SeaLevelRise.NOAA2012$High[-which(SeaLevelRise.NOAA2012$Year<2020)], n = 201)$x[1:(which(abs(spline(SeaLevelRise.NOAA2012$Year[-which(SeaLevelRise.NOAA2012$Year<2020)],ifelse(Unit=="m",1,3.28084)*SeaLevelRise.NOAA2012$High[-which(SeaLevelRise.NOAA2012$Year<2020)], n = 201)$y-SeaLevelRise)==min(abs(spline(SeaLevelRise.NOAA2012$Year[-which(SeaLevelRise.NOAA2012$Year<2020)],ifelse(Unit=="m",1,3.28084)*SeaLevelRise.NOAA2012$High[-which(SeaLevelRise.NOAA2012$Year<2020)], n = 201)$y-SeaLevelRise))))],
      spline(SeaLevelRise.NOAA2012$Year[-which(SeaLevelRise.NOAA2012$Year<2020)],ifelse(Unit=="m",1,3.28084)*SeaLevelRise.NOAA2012$High[-which(SeaLevelRise.NOAA2012$Year<2020)], n = 201)$y[1:(which(abs(spline(SeaLevelRise.NOAA2012$Year[-which(SeaLevelRise.NOAA2012$Year<2020)],ifelse(Unit=="m",1,3.28084)*SeaLevelRise.NOAA2012$High[-which(SeaLevelRise.NOAA2012$Year<2020)], n = 201)$y-SeaLevelRise)==min(abs(spline(SeaLevelRise.NOAA2012$Year[-which(SeaLevelRise.NOAA2012$Year<2020)],ifelse(Unit=="m",1,3.28084)*SeaLevelRise.NOAA2012$High[-which(SeaLevelRise.NOAA2012$Year<2020)], n = 201)$y-SeaLevelRise))))],col=mypalette[1],lwd=5)

par(las=1)
plot(0,xlab="Number of years",ylab="",type='n',xlim=c(2020,2100),ylim=c(0,3),cex.lab=1.5,cex.axis=1.5,yaxt="n",xaxt="n",bty="n")
rect(2020,2,spline(SeaLevelRise.NOAA2012$Year,ifelse(Unit=="m",1,3.28084)*SeaLevelRise.NOAA2012$High, n = 201)$x[(which(abs(spline(SeaLevelRise.NOAA2012$Year,ifelse(Unit=="m",1,3.28084)*SeaLevelRise.NOAA2012$High, n = 201)$y-SeaLevelRise)==min(abs(spline(SeaLevelRise.NOAA2012$Year,ifelse(Unit=="m",1,3.28084)*SeaLevelRise.NOAA2012$High, n = 201)$y-SeaLevelRise))))],2.5,col=mypalette[1],border=NA)
text(spline(SeaLevelRise.NOAA2012$Year,ifelse(Unit=="m",1,3.28084)*SeaLevelRise.NOAA2012$High, n = 201)$x[(which(abs(spline(SeaLevelRise.NOAA2012$Year,ifelse(Unit=="m",1,3.28084)*SeaLevelRise.NOAA2012$High, n = 201)$y-SeaLevelRise)==min(abs(spline(SeaLevelRise.NOAA2012$Year,ifelse(Unit=="m",1,3.28084)*SeaLevelRise.NOAA2012$High, n = 201)$y-SeaLevelRise))))]+1.5,2.25,paste(round(spline(SeaLevelRise.NOAA2012$Year,ifelse(Unit=="m",1,3.28084)*SeaLevelRise.NOAA2012$High, n = 201)$x[(which(abs(spline(SeaLevelRise.NOAA2012$Year,ifelse(Unit=="m",1,3.28084)*SeaLevelRise.NOAA2012$High, n = 201)$y-SeaLevelRise)==min(abs(spline(SeaLevelRise.NOAA2012$Year,ifelse(Unit=="m",1,3.28084)*SeaLevelRise.NOAA2012$High, n = 201)$y-SeaLevelRise))))],0)-2020),cex=1.5,font=3)
High<-spline(SeaLevelRise.NOAA2012$Year,ifelse(Unit=="m",1,3.28084)*SeaLevelRise.NOAA2012$High, n = 201)$x[(which(abs(spline(SeaLevelRise.NOAA2012$Year,ifelse(Unit=="m",1,3.28084)*SeaLevelRise.NOAA2012$High, n = 201)$y-SeaLevelRise)==min(abs(spline(SeaLevelRise.NOAA2012$Year,ifelse(Unit=="m",1,3.28084)*SeaLevelRise.NOAA2012$High, n = 201)$y-SeaLevelRise))))]
rect(2020,1,spline(SeaLevelRise.USACE2013$Year,ifelse(Unit=="m",1,3.28084)*SeaLevelRise.USACE2013$High, n = 201)$x[(which(abs(spline(SeaLevelRise.USACE2013$Year,ifelse(Unit=="m",1,3.28084)*SeaLevelRise.USACE2013$High, n = 201)$y-SeaLevelRise)==min(abs(spline(SeaLevelRise.USACE2013$Year,ifelse(Unit=="m",1,3.28084)*SeaLevelRise.USACE2013$High, n = 201)$y-SeaLevelRise))))],1.5,col=mypalette[2],border=NA)
text(spline(SeaLevelRise.USACE2013$Year,ifelse(Unit=="m",1,3.28084)*SeaLevelRise.USACE2013$High, n = 201)$x[(which(abs(spline(SeaLevelRise.USACE2013$Year,ifelse(Unit=="m",1,3.28084)*SeaLevelRise.USACE2013$High, n = 201)$y-SeaLevelRise)==min(abs(spline(SeaLevelRise.USACE2013$Year,ifelse(Unit=="m",1,3.28084)*SeaLevelRise.USACE2013$High, n = 201)$y-SeaLevelRise))))]+1.5,1.25,paste(round(spline(SeaLevelRise.USACE2013$Year,ifelse(Unit=="m",1,3.28084)*SeaLevelRise.USACE2013$High, n = 201)$x[(which(abs(spline(SeaLevelRise.USACE2013$Year,ifelse(Unit=="m",1,3.28084)*SeaLevelRise.USACE2013$High, n = 201)$y-SeaLevelRise)==min(abs(spline(SeaLevelRise.USACE2013$Year,ifelse(Unit=="m",1,3.28084)*SeaLevelRise.USACE2013$High, n = 201)$y-SeaLevelRise))))],0)-2020),cex=1.5,font=3)
Intermedite<-spline(SeaLevelRise.USACE2013$Year,ifelse(Unit=="m",1,3.28084)*SeaLevelRise.USACE2013$High, n = 201)$x[(which(abs(spline(SeaLevelRise.USACE2013$Year,ifelse(Unit=="m",1,3.28084)*SeaLevelRise.USACE2013$High, n = 201)$y-SeaLevelRise)==min(abs(spline(SeaLevelRise.USACE2013$Year,ifelse(Unit=="m",1,3.28084)*SeaLevelRise.USACE2013$High, n = 201)$y-SeaLevelRise))))]
rect(2020,0,spline(SeaLevelRise.AR5$Year,ifelse(Unit=="m",1,3.28084)*SeaLevelRise.AR5$Int, n = 201)$x[(which(abs(spline(SeaLevelRise.AR5$Year,ifelse(Unit=="m",1,3.28084)*SeaLevelRise.AR5$Int, n = 201)$y-SeaLevelRise)==min(abs(spline(SeaLevelRise.AR5$Year,ifelse(Unit=="m",1,3.28084)*SeaLevelRise.AR5$Int, n = 201)$y-SeaLevelRise))))],0.5,col=mypalette[3],border=NA)
text(spline(SeaLevelRise.AR5$Year,ifelse(Unit=="m",1,3.28084)*SeaLevelRise.AR5$Int, n = 201)$x[(which(abs(spline(SeaLevelRise.AR5$Year,ifelse(Unit=="m",1,3.28084)*SeaLevelRise.AR5$Int, n = 201)$y-SeaLevelRise)==min(abs(spline(SeaLevelRise.AR5$Year,ifelse(Unit=="m",1,3.28084)*SeaLevelRise.AR5$Int, n = 201)$y-SeaLevelRise))))]+1.5,0.25,paste(ifelse((round(spline(SeaLevelRise.AR5$Year,ifelse(Unit=="m",1,3.28084)*SeaLevelRise.AR5$Int, n = 201)$x[(which(abs(spline(SeaLevelRise.AR5$Year,ifelse(Unit=="m",1,3.28084)*SeaLevelRise.AR5$Int, n = 201)$y-SeaLevelRise)==min(abs(spline(SeaLevelRise.AR5$Year,ifelse(Unit=="m",1,3.28084)*SeaLevelRise.AR5$Int, n = 201)$y-SeaLevelRise))))],0)-2020)>79,"> 80",round(spline(SeaLevelRise.AR5$Year,ifelse(Unit=="m",1,3.28084)*SeaLevelRise.AR5$Int, n = 201)$x[(which(abs(spline(SeaLevelRise.AR5$Year,ifelse(Unit=="m",1,3.28084)*SeaLevelRise.AR5$Int, n = 201)$y-SeaLevelRise)==min(abs(spline(SeaLevelRise.AR5$Year,ifelse(Unit=="m",1,3.28084)*SeaLevelRise.AR5$Int, n = 201)$y-SeaLevelRise))))],0)-2020)),cex=1.5,font=3)
Low<-spline(SeaLevelRise.AR5$Year,ifelse(Unit=="m",1,3.28084)*SeaLevelRise.AR5$Int, n = 201)$x[(which(abs(spline(SeaLevelRise.AR5$Year,ifelse(Unit=="m",1,3.28084)*SeaLevelRise.AR5$Int, n = 201)$y-SeaLevelRise)==min(abs(spline(SeaLevelRise.AR5$Year,ifelse(Unit=="m",1,3.28084)*SeaLevelRise.AR5$Int, n = 201)$y-SeaLevelRise))))]
axis(1,at=c(2020,2040,2060,2080,2100),c("0","20","40","60","80"),cex.axis=1.5)
#axis(2,at=c(2.25,1.25,0.25),c("High","Intermediate","Low"),cex.axis=1.5)
mtext(c("NOAA et al (2012)","High","USACE 2013","High","IPCC AR5","Medium"),2,-4.15,at=c(2.4,2.2,1.4,1.2,0.4,0.2))
res<-list("High" = High, "Intermediate" = Intermediate, "Low" = Low)
return(res)
}
