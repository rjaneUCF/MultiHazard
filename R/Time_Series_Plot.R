#' Rainfall and O-sWL time series plots
#'
#' Plots a user specified number of synthetic events where at least O-sWL or rainfall peak exceeds a high threshold. .
#'
#' @param Sample Data frame containing the simulated events. Columns (and their names) required by the function are rainfall peak (Rainfall), O-sWL peak (OsWL), their lag time (Lag), and the ID of the sampled rainfall event (samp).
#' @param Rainfall_Series Data frame with rows comprising time series of rainfall totals associated with cluster maximum of the rainfall series.
#' @param Oswl_Time_Series Data frame with rows comprising the water level curves associated with the simulated events in \code{Sample}.
#' @param Con_Variable Character vector of length one specifying the conditioning variable of the events in \code{Sample}.
#' @param Buffer Numeric vector of length one specifying the extension of the x-axis before and after the rainfall event when \code{Con_Variable == "Rainfall"}. Default is \code{6}.
#' @param Intensity Numeric vector specifying the "intensity" of the O-sWL events in \code{Sample}. Default is \code{NA}.
#' @param Event_ID Numeric vector specifying the events in \code{Sample} to be plot.
#' @param Row Numeric vector of length one specifying the number of rows of subplots in the Figure.
#' @param Col Numeric vector of length one specifying the number of columns of subplots in the Figure. Product of \code{Row} and \code{Col} must be equal to or greater than \code{Event_ID}.
#' @param Mar Numeric vector of length one specifying the margin at the (bottom,left,top,right) of the subplots. Default is \code{c(4.2,4.5,1.5,3.5)}.
#' @return Figure containing a (Row * Col) matrix of subplots each displaying the hyetogaph (grey bars) and water level curve (blue lines) comprising an event.
#' @seealso \code{\link{U_Sample}} \code{\link{WL_Curve}}
#' @export
#' @examples
Time_Series_Plot<-function(Rainfall_Series,Oswl_Time_Series,Sample,Con_Variable,Buffer=6,Intensity=NA,Event_ID=1:16,Row=4,Col=4,Mar=c(4.2,4.5,1.5,3.5)){
 #Figure layout
 par(mfrow=c(Row,Col))
 par(mar=Mar)
 #Loop repeated for each event specified in Event_ID
 for(i in 1:length(Event_ID)){
  j<-Event_ID[i]
  #Hourly rainfall totals sampled using Serinaldi and Kilsby (2013) method
  x_event<-Rainfall_Series[Sample$Samp[j],]
  #Replace peak in sampled event with the simulated peak
  x_event[which(Rainfall_Series[Sample$Samp[j],]==max(Rainfall_Series[Sample$Samp[i],],na.rm=T))] <- S13.rainfall.sample$Xp[i]
  #Find time of peak hourly rainfall
  xp.pos.max<-which(x_event==max(x_event,na.rm=T))[1]
  if(Con_Variable == "Rainfall"){
   #Finding x-axis limits
   x.lower<-min(-Buffer,-xp.pos.max,Sample$Lag[j]-Buffer)
   x.upper<-max(Buffer,Sample$D[j]-xp.pos.max,S13.rainfall.sample$Lag[j]+Buffer)
   #Plotting hyetograph
   plot(0,xlim=c(x.lower-1,x.upper+1),ylim=c(0,max(0.0001,x_event,na.rm=T)),xlab="",ylab="",yaxt='n',type='n')
   for(k in 1:length(x_event)){
    rect(k-xp.pos.max-0.5,
         0,
         k-xp.pos.max+0.5,
         x_event[k],col="dark grey")
   }
   #y-axis tick marks
   axis(side=2, at = pretty(c(0,range(na.omit(x_event)))))
   #y-axis label
   mtext ("Hourly rainfall total [Inches]", side=2, line=2.5,cex=0.6)
   #Plotting water level curve on top of hyetograph
   par(new = TRUE)
   plot(seq(-144,144,1)+Sample$Lag[j],Oswl_Time_Series[j,],xlim=c(x.lower-1,x.upper+1),type='l',ylim=c(0,Sample$OsWL[j]),yaxt='n',xlab="",ylab="",col="blue")
   #2nd y-axis tick marks
   axis(4,labels=seq(0,Sample$OsWL[j],0.5),at=seq(0,Sample$OsWL[j],0.5))
   #Axis labels
   mtext("Time about peak (Hour)", side=1, line=2.25,cex=0.6)
   mtext("Hourly O-sWL [Ft NGVD 29]", side=4, line=2.5,cex=0.6)
   #Vertical blue line denoting time of peak O-sWL
   abline(v=Sample$Lag[j],col="blue",lty=2)
  }
  if(Con_Variable == "OsWL"){
   #Finding x-axis limits
   x.lower<-min(-Buffer,Sample$Lag[j]-xp.pos.max)
   x.upper<-max(Buffer,Sample$Lag[j]+(Sample$D[j]-xp.pos.max))
   #Plotting the hyetograph
   plot(0,xlim=c(x.lower-1,x.upper+1),ylim=c(0,max(0.0001,x_event,na.rm=T)),xlab="",ylab="",yaxt='n',type='n')
   for(k in 1:length(x_event)){
    rect(k-xp.pos.max-0.5+Sample$Lag[j],
         0,
         k-xp.pos.max+0.5+Sample$Lag[j],
         x_event[k],col="dark grey")
   }
   #y-axis ticks
   axis(side=2, at = pretty(c(0,max(na.omit(x_event))),n=5))
   mtext("Hourly rainfall total [Inches]", side=2, line=2.5,cex=0.6)
   #Vertical line denoting time of hourly rainfall peak
   abline(v=Sample$Lag[i],col="Grey",lty=2)
   #Plotting water level curve
   par(new = TRUE)
   plot(seq(-(length(Oswl_Time_Series[1,])-1)/2,(length(Oswl_Time_Series[1,])-1)/2,1),oswl.ts.oswl[i,],xlim=c(x.lower-1,x.upper+1),type='l',ylim=c(0,S13.oswl.sample$OsWL[i]),xaxt='n',yaxt='n',xlab="",ylab="",col="blue")
   #2nd y-axis tick marks
   axis(4,labels=seq(0,Sample$OsWL[i],0.5),at=seq(0,Sample$OsWL[i],0.5))
   #Axis labels
   mtext("Time about peak (Hour)", side=1, line=2.25,cex=0.6)
   mtext("Hourly O-sWL [Ft NGVD 29]", side=4, line=2.5,cex=0.6)
  }
  #Print total rainfall volume in upper left-hand corner of plot
  mtext(paste(round(Sample[i,]$V,2),'"',sep=""),side = 3, adj = 0, cex=0.75)
  #Print "intensity" in upper right-hand corner of plot
  mtext(paste(round(intensity.oswl[i],0)),side = 3, adj = 1, cex=0.75)
 }
}
