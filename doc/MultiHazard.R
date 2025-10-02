## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE, collapse = TRUE,
  comment = "#>")
library(MultiHazard)
library(dplyr)
library(scales)
library(copula)

## ----eval=FALSE---------------------------------------------------------------
# install.packages("remotes")
# remotes::install_github("rjaneUCF/MultiHazard")

## -----------------------------------------------------------------------------
#Viewing first few rows of in the groundwater level records
head(G_3356)

## -----------------------------------------------------------------------------
head(G_3355)

## -----------------------------------------------------------------------------
#Converting Date column to "Date"" object
G_3356$Date<-seq(as.Date("1985-10-23"), as.Date("2019-05-29"), by="day")
G_3355$Date<-seq(as.Date("1985-08-20"), as.Date("2019-06-02"), by="day")
#Converting column containing the readings to a "numeric"" object
G_3356$Value<-as.numeric(as.character(G_3356$Value))
G_3355$Value<-as.numeric(as.character(G_3355$Value))

## -----------------------------------------------------------------------------
#Merge the two dataframes by date
GW_S20<-left_join(G_3356,G_3355,by="Date")
colnames(GW_S20)<-c("Date","G3356","G3355")
#Carrying out imputation
Imp<-Imputation(Data=GW_S20,Variable="G3356",
                x_lab="G-3355 (ft NGVD 29)", y_lab="G-3356 (ft NGVD 29)")

## -----------------------------------------------------------------------------
head(Imp$Data)

## -----------------------------------------------------------------------------
Imp$Model

## -----------------------------------------------------------------------------
G_3356_ValueFilled_NA<-which(is.na(Imp$Data$ValuesFilled)==TRUE)
length(G_3356_ValueFilled_NA) 

## -----------------------------------------------------------------------------
G3356_approx<-approx(seq(1,length(Imp$Data$ValuesFilled),1),Imp$Data$ValuesFilled,
                     xout=seq(1,length(Imp$Data$ValuesFilled),1))
Imp$Data$ValuesFilled[which(is.na(Imp$Data$ValuesFilled)==TRUE)]<-
  G3356_approx$y[which(is.na(Imp$Data$ValuesFilled)==TRUE)]

## -----------------------------------------------------------------------------
#Creating a data frame with the imputed series alongside the corresponding dates 
G_3356_Imp<-data.frame(Imp$Data$Date,Imp$Data$ValuesFilled)
colnames(G_3356_Imp)<-c("Date","ValuesFilled")
#Detrending
G_3356_Detrend<-Detrend(Data=G_3356_Imp,PLOT=TRUE,x_lab="Date",
                        y_lab="Groundwater level (ft NGVD 29)")

## -----------------------------------------------------------------------------
head(G_3356_Detrend)

## -----------------------------------------------------------------------------
S20.Groundwater.Detrend.df<-data.frame(as.Date(GW_S20$Date),G_3356_Detrend)
colnames(S20.Groundwater.Detrend.df)<-c("Date","Groundwater")

## -----------------------------------------------------------------------------
G_3356.Declustered<-Decluster(Data=G_3356_Detrend,u=0.95,SepCrit=3,mu=365.25)

## -----------------------------------------------------------------------------
G_3356_Imp$Detrend<-G_3356_Detrend
plot(as.Date(G_3356_Imp$Date),G_3356_Imp$Detrend,col="Grey",pch=16,
     cex=0.25,xlab="Date",ylab="Groundwater level (ft NGVD 29)")
abline(h=G_3356.Declustered$Threshold,col="Dark Green")
points(as.Date(G_3356_Imp$Date[G_3356.Declustered$EventsMax]),
       G_3356.Declustered$Declustered[G_3356.Declustered$EventsMax],
       col="Red",pch=16,cex=0.5)

## -----------------------------------------------------------------------------
G_3356.Declustered$Threshold

## -----------------------------------------------------------------------------
G_3356.Declustered$Rate

## -----------------------------------------------------------------------------
S20.Groundwater.Detrend.Declustered<-G_3356.Declustered$Declustered

## -----------------------------------------------------------------------------
#Changing names of the data frames
S20.Rainfall.df<-Perrine_df
S20.OsWL.df<-S20_T_MAX_Daily_Completed_Detrend_Declustered[,c(2,4)]
#Converting Date column to "Date"" object
S20.Rainfall.df$Date<-as.Date(S20.Rainfall.df$Date)
S20.OsWL.df$Date<-as.Date(S20.OsWL.df$Date)

## -----------------------------------------------------------------------------
S20.OsWL.Detrend<-Detrend(Data=S20.OsWL.df,Method = "window",PLOT=FALSE,
                          x_lab="Date",y_lab="O-sWL (ft NGVD 29)")

## -----------------------------------------------------------------------------
S20.OsWL.Detrend.df<-data.frame(as.Date(S20.OsWL.df$Date),S20.OsWL.Detrend)
colnames(S20.OsWL.Detrend.df)<-c("Date","OsWL")

## -----------------------------------------------------------------------------
#Declustering rainfall and O-sWL series
S20.Rainfall.Declustered<-Decluster(Data=S20.Rainfall.df$Value,u=0.95,SepCrit=3)$Declustered
S20.OsWL.Detrend.Declustered<-Decluster(Data=S20.OsWL.Detrend,u=0.95,SepCrit=3,mu=365.25)$Declustered

## -----------------------------------------------------------------------------
S20.OsWL.Detrend.Declustered.df<-data.frame(S20.OsWL.df$Date,S20.OsWL.Detrend.Declustered)
colnames(S20.OsWL.Detrend.Declustered.df)<-c("Date","OsWL")
S20.Rainfall.Declustered.df<-data.frame(S20.Rainfall.df$Date,S20.Rainfall.Declustered)
colnames(S20.Rainfall.Declustered.df)<-c("Date","Rainfall")
S20.Groundwater.Detrend.Declustered.df<-data.frame(G_3356$Date,S20.Groundwater.Detrend.Declustered)
colnames(S20.Groundwater.Detrend.Declustered.df)<-c("Date","Groundwater")

## -----------------------------------------------------------------------------
S20.Detrend.df<-Dataframe_Combine(data.1<-S20.Rainfall.df,
                                  data.2<-S20.OsWL.Detrend.df,
                                  data.3<-S20.Groundwater.Detrend.df,
                                  names=c("Rainfall","OsWL","Groundwater"), n=3)
S20.Detrend.Declustered.df<-Dataframe_Combine(data.1<-S20.Rainfall.Declustered.df,
                                              data.2<-S20.OsWL.Detrend.Declustered.df,
                                              data.3<-S20.Groundwater.Detrend.Declustered.df,
                                              names=c("Rainfall","OsWL","Groundwater"), n=3)

## -----------------------------------------------------------------------------
S20.Rainfall.Declustered.SW<-Decluster_SW(Data=S20.Rainfall.df,Window_Width=7)

## -----------------------------------------------------------------------------
plot(S20.Rainfall.df$Date,S20.Rainfall.df$Value,pch=16,cex=0.5,
     xlab="Date",ylab="Total daily rainfall (Inches)")
points(S20.Rainfall.df$Date,S20.Rainfall.Declustered.SW$Declustered,pch=16,col=2,cex=0.5)

## -----------------------------------------------------------------------------
S20.OsWL.Declustered.SW<-Decluster_SW(Data=S20.OsWL.df,Window_Width=3)

## -----------------------------------------------------------------------------
#Declustering
S20.Rainfall.Declustered.S.SW<-Decluster_S_SW(Data=S20.Rainfall.df, 
                                              Window_Width_Sum=7, Window_Width=7)
#First twenty values of the weekly totals
S20.Rainfall.Declustered.S.SW$Totals[1:20]

## -----------------------------------------------------------------------------
#First ten values of the declustered weekly totals
S20.Rainfall.Declustered.S.SW$Declustered[1:20]

## -----------------------------------------------------------------------------
plot(S20.Rainfall.df$Date,S20.Rainfall.Declustered.S.SW$Totals,pch=16,cex=0.5,
     xlab="Date",ylab="Total weekly rainfall (Inches)")
points(S20.Rainfall.df$Date,S20.Rainfall.Declustered.S.SW$Declustered,pch=16,col=2,cex=0.5)

## -----------------------------------------------------------------------------
GPD_Fit(Data=S20.Detrend.Declustered.df$Rainfall,Data_Full=na.omit(S20.Detrend.df$Rainfall),
        u=0.997,PLOT=TRUE,xlab_hist="Rainfall (Inches)",y_lab="Rainfall (Inches)")

## -----------------------------------------------------------------------------
S20.Rainfall.Solari<-GPD_Threshold_Solari(Event=S20.Rainfall.Declustered.SW$Declustered,
                                          Data=S20.Detrend.df$Rainfall)

## -----------------------------------------------------------------------------
S20.Rainfall.Solari$Candidate_Thres

## -----------------------------------------------------------------------------
Rainfall.Thres.Quantile<-ecdf(S20.Detrend.df$Rainfall)(S20.Rainfall.Solari$Candidate_Thres)

## -----------------------------------------------------------------------------
Solari.Sel<-GPD_Threshold_Solari_Sel(Event=S20.Rainfall.Declustered.SW$Declustered,
                                    Data=S20.Detrend.df$Rainfall,
                                    Solari_Output=S20.Rainfall.Solari,
                                    Thres=S20.Rainfall.Solari$Candidate_Thres,
                                    RP_Max=100)

## ----warning=FALSE, message=FALSE, error=FALSE--------------------------------
S20.OsWL.Solari<-GPD_Threshold_Solari(Event=S20.OsWL.Declustered.SW$Declustered,
                                      Data=S20.Detrend.df$OsWL)

## -----------------------------------------------------------------------------
S20.OsWL.Solari$Candidate_Thres

## -----------------------------------------------------------------------------
OsWL.Thres.Quantile<-ecdf(S20.Detrend.df$OsWL)(S20.OsWL.Solari$Candidate_Thres)

## -----------------------------------------------------------------------------
Solari.Sel<-GPD_Threshold_Solari_Sel(Event=S20.OsWL.Declustered.SW$Declustered,
                                     Data=S20.Detrend.df$OsWL,
                                     Solari_Output=S20.OsWL.Solari,
                                     Thres=S20.OsWL.Solari$Candidate_Thres,
                                     RP_Max=100)

## -----------------------------------------------------------------------------
S20.Kendall.Results<-Kendall_Lag(Data=S20.Detrend.df,GAP=0.2)

## -----------------------------------------------------------------------------
S20.Kendall.Results$Value$Rainfall_OsWL

## -----------------------------------------------------------------------------
S20.Kendall.Results$Test$Rainfall_OsWL_Test

## -----------------------------------------------------------------------------
Copula_Threshold_2D(Data_Detrend=S20.Detrend.df[,-c(1,4)],
Data_Declust=S20.Detrend.Declustered.df[,-c(1,4)],
y_lim_min=-0.075, y_lim_max =0.25,
Upper=c(2,9), Lower=c(2,10),GAP=0.15)

## -----------------------------------------------------------------------------
S20.Rainfall<-Con_Sampling_2D(Data_Detrend=S20.Detrend.df[,-c(1,4)],
Data_Declust=S20.Detrend.Declustered.df[,-c(1,4)],
Con_Variable="Rainfall",u = Rainfall.Thres.Quantile)
Diag_Non_Con(Data=S20.Rainfall$Data$OsWL,Omit=c("Gum","RGum"),x_lab="O-sWL (ft NGVD 29)",y_lim_min=0,y_lim_max=1.5)

## -----------------------------------------------------------------------------
Diag_Non_Con_Sel(Data=S20.Rainfall$Data$OsWL,x_lab="O-sWL (ft NGVD 29)",
y_lim_min=0,y_lim_max=1.5,Selected="Logis")

## -----------------------------------------------------------------------------
S20.OsWL<-Con_Sampling_2D(Data_Detrend=S20.Detrend.df[,-c(1,4)],
Data_Declust=S20.Detrend.Declustered.df[,-c(1,4)],
Con_Variable="OsWL",u=OsWL.Thres.Quantile)
S20.OsWL$Data$Rainfall<-S20.OsWL$Data$Rainfall+runif(length(S20.OsWL$Data$Rainfall),0.001,0.01)
Diag_Non_Con_Trunc(Data=S20.OsWL$Data$Rainfall+0.001,x_lab="Rainfall (Inches)",
y_lim_min=0,y_lim_max=2)

## -----------------------------------------------------------------------------
Diag_Non_Con_Trunc_Sel(Data=S20.OsWL$Data$Rainfall+0.001,x_lab="Rainfall (Inches)",
y_lim_min=0,y_lim_max=2,
Selected="BS")

## -----------------------------------------------------------------------------
S20.Copula.Rainfall<-Copula_Threshold_2D(Data_Detrend=S20.Detrend.df[,-c(1,4)],
Data_Declust=S20.Detrend.Declustered.df[,-c(1,4)],
u1=Rainfall.Thres.Quantile,u2=NA,
y_lim_min=0,y_lim_max=0.25, GAP=0.075)$Copula_Family_Var1

## -----------------------------------------------------------------------------
S20.Copula.OsWL<-Copula_Threshold_2D(Data_Detrend=S20.Detrend.df[,-c(1,4)],
Data_Declust=S20.Detrend.Declustered.df[,-c(1,4)],
u1=NA,u2=OsWL.Thres.Quantile,
y_lim_min=0,y_lim_max=0.25,GAP=0.075)$Copula_Family_Var2

## -----------------------------------------------------------------------------
S20.Bivariate<-Design_Event_2D(Data=S20.Detrend.df[,-c(1,4)], 
Data_Con1=S20.Rainfall$Data, 
Data_Con2=S20.OsWL$Data, 
u1=Rainfall.Thres.Quantile, 
u2=OsWL.Thres.Quantile, 
Copula_Family1=S20.Copula.Rainfall,
Copula_Family2=S20.Copula.OsWL, 
Marginal_Dist1="Logis", Marginal_Dist2="BS",
x_lab="Rainfall (Inches)",y_lab="O-sWL (ft NGVD 29)",
RP=100,N=10^7,N_Ensemble=10)

## -----------------------------------------------------------------------------
S20.Bivariate$MostLikelyEvent$`100`

## -----------------------------------------------------------------------------
S20.Bivariate$FullDependence$`100`

## -----------------------------------------------------------------------------
#Fitting the marginal distribution
#See next section for information on the Migpd_Fit function
S20.GPD<-Migpd_Fit(Data=S20.Detrend.Declustered.df[,2:3], Data_Full = S20.Detrend.df[,2:3], 
                   mqu =c(0.99,0.99))
#10-year exceedance probability for daily data
p.10<-(1/365.25)/10
#10-year exceedance probability for daily data
p.100<-(1/365.25)/100
#Calculating the isoline
isoline<-Cooley19(Data=na.omit(S20.Detrend.df[,2:3]),Migpd=S20.GPD,
                  p.base=p.10,p.proj=p.100,PLOT=TRUE,x_lim_max_T=15000,y_lim_max_T=15000)

## -----------------------------------------------------------------------------
#Adding dates to complete final month of combined records
final.month = data.frame(seq(as.Date("2019-02-04"),as.Date("2019-02-28"),by="day"),NA,NA,NA)
colnames(final.month) = c("Date","Rainfall","OsWL","Groundwater")
S22.Detrend.df.extended = rbind(S22.Detrend.df,final.month)

#Derive return curves
curve = return_curve_est(data=S22.Detrend.df.extended[,1:3],
                         q=0.985,rp=100,mu=365.25,n_sim=100,
                         n_grad=50,n_boot=100,boot_method="monthly",
                         boot_replace=NA, block_length=NA, boot_prop=0.8,
                         decl_method_x="runs", decl_method_y="runs",
                         window_length_x=NA,window_length_y=NA,
                         u_x=0.95, u_y=0.95,
                         sep_crit_x=36, sep_crit_y=36,
                         most_likely=T, n_ensemble=10,
                         alpha=0.1, x_lab=NA, y_lab=NA)

## -----------------------------------------------------------------------------
head(curve$median_wt13)

## -----------------------------------------------------------------------------
head(curve$ub_wt13)
head(curve$lb_wt13)

## -----------------------------------------------------------------------------
#"Most-likely" design event
curve$most_likely_wt13

#Ensemble of ten design events 
curve$ensemble_wt13

## ----fig.width=6, fig.height=4------------------------------------------------
#Diagnostic plots for the return curves
curve = return_curve_diag(data=S22.Detrend.df.extended[,1:3],
                          q=0.985,rp=1,mu=365.25,n_sim=100,
                          n_grad=50,n_boot=100,boot_method="monthly",
                          boot_replace=NA, block_length=NA, boot_prop=0.8,
                          decl_method_x="runs", decl_method_y="runs",
                          window_length_x=NA,window_length_y=NA,
                          u_x=0.95, u_y=0.95,
                          sep_crit_x=36, sep_crit_y=36,
                          alpha=0.1,
                          boot_method_all="block", boot_replace_all=NA,
                          block_length_all=14)

## -----------------------------------------------------------------------------
S20.Migpd<-Migpd_Fit(Data=S20.Detrend.Declustered.df[,-1],Data_Full = S20.Detrend.df[,-1],
                     mqu=c(0.975,0.975,0.9676))
summary(S20.Migpd)

## -----------------------------------------------------------------------------
S20.Gaussian<-Standard_Copula_Fit(Data=S20.Detrend.df,Copula_Type="Gaussian")

## -----------------------------------------------------------------------------
S20.Gaussian.Sim<-Standard_Copula_Sim(Data=S20.Detrend.df,Marginals=S20.Migpd,
                                      Copula=S20.Gaussian,N=100)

## -----------------------------------------------------------------------------
S20.Pairs.Plot.Data<-data.frame(rbind(na.omit(S20.Detrend.df[,-1]),S20.Gaussian.Sim$x.Sim),
                                c(rep("Observation",nrow(na.omit(S20.Detrend.df))),
                                  rep("Simulation",nrow(S20.Gaussian.Sim$x.Sim))))
colnames(S20.Pairs.Plot.Data)<-c(names(S20.Detrend.df)[-1],"Type")
pairs(S20.Pairs.Plot.Data[,1:3],
      col=ifelse(S20.Pairs.Plot.Data$Type=="Observation","Black",alpha("Red",0.3)),
      upper.panel=NULL,pch=16)

## -----------------------------------------------------------------------------
Standard_Copula_Sel(Data=S20.Detrend.df)

## -----------------------------------------------------------------------------
S20.Vine<-Vine_Copula_Fit(Data=S20.Detrend.df)

## -----------------------------------------------------------------------------
S20.Vine.Sim<-Vine_Copula_Sim(Data=S20.Detrend.df,Vine_Model=S20.Vine,Marginals=S20.Migpd,N=100)

## -----------------------------------------------------------------------------
S20.Pairs.Plot.Data<-data.frame(rbind(na.omit(S20.Detrend.df[,-1]),S20.Vine.Sim$x.Sim),
                                c(rep("Observation",nrow(na.omit(S20.Detrend.df))),
                                  rep("Simulation",nrow(S20.Vine.Sim$x.Sim))))
colnames(S20.Pairs.Plot.Data)<-c(names(S20.Detrend.df)[-1],"Type")
pairs(S20.Pairs.Plot.Data[,1:3],
      col=ifelse(S20.Pairs.Plot.Data$Type=="Observation","Black",alpha("Red",0.3)),
      upper.panel=NULL,pch=16)

## -----------------------------------------------------------------------------
S20.HT04<-HT04(data_Detrend_Dependence_df=S20.Detrend.df,
               data_Detrend_Declustered_df=S20.Detrend.Declustered.df,
               u_Dependence=0.995,Migpd=S20.Migpd,mu=365.25,N=1000)

## -----------------------------------------------------------------------------
S20.HT04$Model$Rainfall

## -----------------------------------------------------------------------------
S20.HT04$Prop

## -----------------------------------------------------------------------------
S20.HT04.Sim<-S20.HT04$x.sim

## -----------------------------------------------------------------------------
S20.Pairs.Plot.Data<-data.frame(rbind(na.omit(S20.Detrend.df[,-1]),S20.HT04.Sim),
                                c(rep("Observation",nrow(na.omit(S20.Detrend.df))),
                                  rep("Simulation",nrow(S20.HT04.Sim))))
colnames(S20.Pairs.Plot.Data)<-c(names(S20.Detrend.df)[-1],"Type")
pairs(S20.Pairs.Plot.Data[,1:3],
      col=ifelse(S20.Pairs.Plot.Data$Type=="Observation","Black",alpha("Red",0.2)),
      upper.panel=NULL,pch=16)

## -----------------------------------------------------------------------------
#Difference in O-sWL between the most-likely and full dependence events
Diff<-S20.Bivariate$FullDependence$`100`$OsWL-S20.Bivariate$MostLikelyEvent$`100`$OsWL
Diff

## -----------------------------------------------------------------------------
#Time in years for the sea level rise to occur
SLR_Scenarios(SeaLevelRise=Diff,Unit="m")

## -----------------------------------------------------------------------------
#Sea level rise scenarios for Fort Myers
head(sl_taskforce_scenarios_psmsl_id_1106_Fort_Myers)

#Formatting to a data frame that can be interpreted by the tool
SeaLevelRise.2022_input<-data.frame(Year=seq(2020,2150,10),
"High"=as.numeric(sl_taskforce_scenarios_psmsl_id_1106_Fort_Myers[14,-(1:5)])/1000,
"Int_Medium"=as.numeric(sl_taskforce_scenarios_psmsl_id_1106_Fort_Myers[11,-(1:5)])/1000,
"Medium"=as.numeric(sl_taskforce_scenarios_psmsl_id_1106_Fort_Myers[8,-(1:5)])/1000,
"Int_Low"=as.numeric(sl_taskforce_scenarios_psmsl_id_1106_Fort_Myers[5,-(1:5)])/1000,
"Low"=as.numeric(sl_taskforce_scenarios_psmsl_id_1106_Fort_Myers[2,-(1:5)])/1000)

#Finding time in years for 0.8m of sea level rise to occur
SLR_Scenarios(SeaLevelRise=0.8, Scenario="Other", Unit = "m", Year=2022, 
              Location="Fort Myers", New_Scenario=SeaLevelRise.2022_input)

## -----------------------------------------------------------------------------
#Decluster O-sWL series at S-13 using a runs method
S13.OsWL.Declust = Decluster(Data=S13.Detrend.df$OsWL,
                            SepCrit=24*7, u=0.99667)
#Calculate O-sWL of the identified cluster maximum
intensity = Intensity(Data=S13.Detrend.df[,c(1,3)],
                      Cluster_Max=S13.OsWL.Declust$EventsMax,
                      Base_Line=2)

## -----------------------------------------------------------------------------
#Plotting water levels
#Converting Date_Time column to POSIXct class
S13.Detrend.df$Date_Time <- as.POSIXct(S13.Detrend.df$Date_Time, 
                                         format = "%Y-%m-%d %H:%M:%S")
plot(S13.Detrend.df$Date_Time[(S13.OsWL.Declust$EventsMax[1]-48):(S13.OsWL.Declust$EventsMax[1]+48)], 
     S13.Detrend.df$OsWL[(S13.OsWL.Declust$EventsMax[1]-48):(S13.OsWL.Declust$EventsMax[1]+48)],
     xlab="Time", ylab="O-sWL (ft NGVD 29)",type='l',lwd=1.5)

#Adding purple points denoting preceding and following high tides
points(S13.Detrend.df$Date_Time[intensity$Pre.High[1]], 
     S13.Detrend.df$OsWL[intensity$Pre.High[1]],pch=16,cex=1.5,col="Purple")
points(S13.Detrend.df$Date_Time[intensity$Fol.High[1]], 
       S13.Detrend.df$OsWL[intensity$Fol.High[1]],pch=16,cex=1.5,col="Purple")

#Adding orange points denoting preceding and following low tides
points(S13.Detrend.df$Date_Time[intensity$Pre.Low[1]], 
       S13.Detrend.df$OsWL[intensity$Pre.Low[1]],pch=16,cex=1.5,col="Orange")
points(S13.Detrend.df$Date_Time[intensity$Fol.Low[1]], 
       S13.Detrend.df$OsWL[intensity$Fol.Low[1]],pch=16,cex=1.5,col="Orange")

#Recall BaseLine=2
baseline= 2

# Only create polygon showing intensity
above<- S13.Detrend.df$OsWL[intensity$Pre.Low[1]:intensity$Fol.Low[1]] > baseline

if(any(above)) {
  runs <- rle(above)
  ends <- cumsum(runs$lengths)
  starts <- c(1, ends[-length(ends)] + 1)
  
  for(j in which(runs$values)) {
    start_idx <- starts[j]
    end_idx <- ends[j]
    
    x_seg <- S13.Detrend.df$Date_Time[intensity$Pre.Low[1]:intensity$Fol.Low[1]][start_idx:end_idx]
    y_seg <- S13.Detrend.df$OsWL[intensity$Pre.Low[1]:intensity$Fol.Low[1]][start_idx:end_idx]
    
    polygon(x = c(x_seg, rev(x_seg)),
            y = c(y_seg, rep(baseline, length(x_seg))), 
            col = "dark grey", border = NA)
  }
}

## -----------------------------------------------------------------------------
#Four synthetic events all with intensity of 60 units
sim.peaks = c(3.4,4,4.2,5)
sim.intensity = c(60,60,60,60)

#Generating the water level curves
oswl_ts_oswl = WL_Curve(Data = S13.Detrend.df,
                        Cluster_Max = S13.OsWL.Declust$EventsMax,
                        Pre_Low = intensity$Pre.Low,
                        Fol_Low = intensity$Fol.Low,
                        Thres = S13.OsWL.Declust$Threshold, Limit = 45,
                        Peak = sim.peaks,
                        Base_Line=2,
                        Intensity = sim.intensity)

## -----------------------------------------------------------------------------
#Plot the water level curves of the observed peaks
plot(-144:144,
     S13.Detrend.df$OsWL[(S13.OsWL.Declust$EventsMax[1]-144):
                         (S13.OsWL.Declust$EventsMax[1]+144)],
     xlab="Time relative to peak(hour)", ylab="O-sWL (ft NGVD 29)", type='l',ylim=c(-2,6))
for(i in 2:length(S13.OsWL.Declust$EventsMax)){
  lines(-144:144,
        S13.Detrend.df$OsWL[(S13.OsWL.Declust$EventsMax[i]-144):
                            (S13.OsWL.Declust$EventsMax[i]+144)])
}
#Superimpose the curves generated for the four synthetic events
for(i in 1:4){
  lines(-144:144,oswl_ts_oswl$Series[i,],col=2)
}

## -----------------------------------------------------------------------------
#First decluster the rainfall series to find the 500 events
#with the highest peaks
S13.Rainfall.Declust = Decluster(Data=S13.Detrend.df$Rainfall,
                                 SepCrit=24*3, u=0.99667)
#Hourly peaks
peaks = S13.Detrend.df$Rainfall[S13.Rainfall.Declust$EventsMax]
#Set very small rainfall measurements to zero.
#Assumed to be the result of uncertainty in measuring equipment.
S13.Detrend.df$Rainfall[which(S13.Detrend.df$Rainfall<0.01)] = 0
#Find NAs in rainfall series
z = which(is.na(S13.Detrend.df$Rainfall)==T)
#Temporarily set NAs to zero
S13.Detrend.df$Rainfall[z] = 0
#Find times where there is 6-hours of no rainfall
no.rain = rep(NA,length(S13.Detrend.df$Rainfall))
for(i in 6:length(S13.Detrend.df$Rainfall)){
  no.rain[i] = ifelse(sum(S13.Detrend.df$Rainfall[(i-5):i])==0,i,NA)
}
#Remove NAs from results vector as these correspond to times where there is
#rainfall at certain points in the 6 hour period.
no.rain = na.omit(no.rain)
#Reset missing values in the rainfall record back to NA
S13.Detrend.df$Rainfall[z] = NA
#Find the start and end times of the 500 events.
start = rep(NA,length(S13.Rainfall.Declust$EventsMax))
end = rep(NA,length(S13.Rainfall.Declust$EventsMax))
for(i in 1:length(S13.Rainfall.Declust$EventsMax)){
 start[i] = max(no.rain[which(no.rain<S13.Rainfall.Declust$EventsMax[i])])
 end[i] = min(no.rain[which(no.rain>S13.Rainfall.Declust$EventsMax[i])])
}
start = start + 1
end = end - 6
d = end - start + 1 #Duration
#Simulate some peaks by sampling observed peaks with replacement
#I.e., applying the method exactly as in Serinaldi and Kilsby (2013)
sim.peak = sample(peaks,size=500,replace=TRUE)

## -----------------------------------------------------------------------------
sample = U_Sample(Data=S13.Detrend.df$Rainfall,
                  Cluster_Max=S13.Rainfall.Declust$EventsMax,
                  D=d,Start=start,End=end,
                  Xp=sim.peak)

## -----------------------------------------------------------------------------
#Calculating volume and intensity
v<-rep(NA,500)
for(i in 1:length(S13.Rainfall.Declust$EventsMax)){
  v[i] = sum(S13.Detrend.df$Rainfall[(start[i]):(end[i])])
}
I = v/d

#Putting in a data.frame
observations = data.frame(peaks,d,v,I)
colnames(observations) = c("peak","d","v","I")

## -----------------------------------------------------------------------------
#Observations
observations.u = data.frame(pobs(observations))
colnames(observations.u) = c("peak","d","v","I")

#Sample
sample.u = data.frame(pobs(sample))

## ----fig.width=6, fig.height=8------------------------------------------------
#Layout of the plots
par(mfrow=c(4,6))
par(mar=c(4.2,4.2,0.1,0.1))

#Characteristics of observations on original scale
plot(observations$peak,observations$I,pch=16,xlab="",ylab="",cex.axis=1.25)
mtext(expression('X'[p]),side=1,line=3)
mtext("I",side=2,line=2.3)
plot(observations$peak,observations$v,pch=16,xlab="",ylab="",cex.axis=1.25)
mtext(expression('X'[p]),side=1,line=3)
mtext("V",side=2,line=2.3)
plot(observations$peak,observations$d,pch=16,xlab="",ylab="",cex.axis=1.25)
mtext(expression('X'[p]),side=1,line=3)
mtext("D",side=2,line=2.3)
plot(observations$I,observations$v,pch=16,xlab="",ylab="",cex.axis=1.25)
mtext("I",side=1,line=2.5)
mtext("V",side=2,line=2.3)
plot(observations$I,observations$d,pch=16,xlab="",ylab="",cex.axis=1.25)
mtext("I",side=1,line=2.5)
mtext("D",side=2,line=2.3)
plot(observations$v,observations$d,pch=16,xlab="",ylab="",cex.axis=1.25)
mtext("V",side=1,line=2.5)
mtext("D",side=2,line=2.3)

#Characteristics of sample on original scale
plot(sample$Xp,sample$I,pch=16,xlab="",ylab="",cex.axis=1.25,col=2)
mtext(expression('X'[p]),side=1,line=3)
mtext('I',side=2,line=2.3)
plot(sample$Xp,sample$V,pch=16,xlab="",ylab="",cex.axis=1.25,col=2)
mtext(expression('X'[p]),side=1,line=3)
mtext('V',side=2,line=2.3)
plot(sample$Xp,sample$D,pch=16,xlab="",ylab="",cex.axis=1.25,col=2)
mtext(expression('X'[p]),side=1,line=3)
mtext("D",side=2,line=2.3)
plot(sample$I,sample$V,pch=16,xlab="",ylab="",cex.axis=1.25,col=2)
mtext('I',side=1,line=2.5)
mtext('V',side=2,line=2.3)
plot(sample$I,sample$D,pch=16,xlab="",ylab="",cex.axis=1.25,col=2)
mtext('I',side=1,line=2.5)
mtext("D",side=2,line=2.3)
plot(sample$V,sample$D,pch=16,xlab="",ylab="",cex.axis=1.25,col=2)
mtext('V',side=1,line=2.5)
mtext("D",side=2,line=2.3)

#Characteristics of observations on the [0,1] scale
plot(observations.u$peak,observations.u$I,pch=16,xlab="",ylab="",cex.axis=1.25)
mtext(expression('X'[p]),side=1,line=3)
mtext('I',side=2,line=2.3)
plot(observations.u$peak,observations.u$v,pch=16,xlab="",ylab="",cex.axis=1.25)
mtext(expression('X'[p]),side=1,line=3)
mtext('V',side=2,line=2.3)
plot(observations.u$peak,observations.u$d,pch=16,xlab="",ylab="",cex.axis=1.25)
mtext(expression('X'[p]),side=1,line=3)
mtext('D',side=2,line=2.3)
plot(observations.u$I,observations.u$v,pch=16,xlab="",ylab="",cex.axis=1.25)
mtext('I',side=1,line=2.5)
mtext('V',side=2,line=2.3)
plot(observations.u$I,observations.u$d,pch=16,xlab="",ylab="",cex.axis=1.25)
mtext('I',side=1,line=2.5)
mtext('D',side=2,line=2.3)
plot(observations.u$v,observations.u$d,pch=16,xlab="",ylab="",cex.axis=1.25)
mtext('V',side=1,line=2.5)
mtext('D',side=2,line=2.3)

#Characteristics of sample on the [0,1] scale
plot(sample.u$Xp,sample.u$I,pch=16,xlab="",ylab="",cex.axis=1.25,col=2)
mtext(expression('X'[p]),side=1,line=3)
mtext('I',side=2,line=2.3)
plot(sample.u$Xp,sample.u$V,pch=16,xlab="",ylab="",cex.axis=1.25,col=2)
mtext(expression('X'[p]),side=1,line=3)
mtext('V',side=2,line=2.3)
plot(sample.u$Xp,sample.u$D,pch=16,xlab="",ylab="",cex.axis=1.25,col=2)
mtext(expression('X'[p]),side=1,line=3)
mtext("D",side=2,line=2.3)
plot(sample.u$I,sample.u$V,pch=16,xlab="",ylab="",cex.axis=1.25,col=2)
mtext('I',side=1,line=2.5)
mtext('V',side=2,line=2.3)
plot(sample.u$I,sample.u$D,pch=16,xlab="",ylab="",cex.axis=1.25,col=2)
mtext('I',side=1,line=2.5)
mtext("D",side=2,line=2.3)
plot(sample.u$V,sample.u$D,pch=16,xlab="",ylab="",cex.axis=1.25,col=2)
mtext('V',side=1,line=2.5)
mtext("D",side=2,line=2.3)

