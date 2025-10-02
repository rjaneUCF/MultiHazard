#Append S-22 dataframe with HURDAT data
S22.Detrend.df = HURDAT(Data=S22.Detrend.df,lat.loc=25.669610,lon.loc=80.283732,rad=500)

#Inter event rainfall
S22.Detrend.df$Rainfall[which(S22.Detrend.df$Rainfall<0.01)]<-0
z<-which(is.na(S22.Detrend.df$Rainfall)==T)
S22.Detrend.df$Rainfall[z]<-0
no.rain<-rep(NA,length(S22.Detrend.df$Rainfall))
for(i in 6:length(S22.Detrend.df$Rainfall)){
  no.rain[i]<-ifelse(sum(S22.Detrend.df$Rainfall[(i-5):i])==0,i,NA)
}
no.rain<-na.omit(no.rain)

#Find start and end of all rainfall events
start = no.rain[which(diff(no.rain)>1)]
end = no.rain[which(diff(no.rain)>1)+1]

start = start + 1
end = end - 6
d = end - start + 1

S22.Detrend.df$Rainfall[z]<-NA

#Find volume
v<-rep(NA,length(start))
for(i in 1:length(start)){
  v[i] = sum(S22.Detrend.df$Rainfall[(start[i]):(end[i])])
}

event.mid = start+ round(d/2,0)

#Putting volume into a dataframe
S22.Detrend.df.declust = data.frame(S22.Detrend.df[,1],rep(NA,nrow(S22.Detrend.df)),S22.Detrend.df[,c(3,5)])
colnames(S22.Detrend.df.declust)<-c("Date","Rainfall_vol","OsWL","Name")
S22.Detrend.df.declust$Rainfall_vol = rep(0, nrow(S22.Detrend.df.declust))
S22.Detrend.df.declust$Rainfall_vol[event.mid] = v

#Obtaining the final conditional sample
rainfall.thres = 5
events = which(v>rainfall.thres)
wl = rep(NA,length(events))
tc = rep(NA,length(events))
for(j in 1:length(events)){
  if(d[events[j]]<72){
    wl[j] = max(S22.Detrend.df$OsWL[(max(1,event.mid[events[j]]-5)):(min(event.mid[events[j]]+5,nrow(S22.Detrend.df)))])
  }else{
    wl[j] = max(S22.Detrend.df$OsWL[start[events[j]]:end[events[j]]])
  }
  tc[j] = ifelse(any(!is.na(S22.Detrend.df$Name[(max(1,event.mid[events[j]]-5)):(min(event.mid[events[j]]+5,nrow(S22.Detrend.df)))])),"tc","nontc")
}

Sample_df = na.omit(data.frame(v[events],wl, tc, event.mid[events]))
colnames(Sample_df) = c("Rainfall_vol","OsWL","type","x.con")

con.sample.rainfall.tc = list(Threshold = rainfall.thres, Data = Sample_df[which(Sample_df$type=="tc"),1:2], Con_Variable = "rainfall_vol",  x.con = Sample_df$x.con[which(Sample_df$type=="tc")])
con.sample.rainfall.nontc = list(Threshold = rainfall.thres, Data = Sample_df[which(Sample_df$type=="nontc"),1:2], Con_Variable = "rainfall_vol",  x.con = Sample_df$x.con[which(Sample_df$type=="nontc")])

cop.rainfall.tc <- BiCopSelect(pobs(con.sample.rainfall.tc$Data$Rainfall_vol), pobs(con.sample.rainfall.tc$Data$OsWL), familyset = NA, selectioncrit = "AIC",
                               indeptest = FALSE, level = 0.05, weights = NA,
                               rotations = TRUE, se = FALSE, presel = TRUE,
                               method = "mle")$family

cop.rainfall.nontc <- BiCopSelect(pobs(con.sample.rainfall.nontc$Data$Rainfall_vol), pobs(con.sample.rainfall.nontc$Data$OsWL), familyset = NA, selectioncrit = "AIC",
                                  indeptest = FALSE, level = 0.05, weights = NA,
                                  rotations = TRUE, se = FALSE, presel = TRUE,
                                  method = "mle")$family

#Sample conditioned on water level

#Obtaining the final conditional sample
oswl.thres = 3
x.con = which(S22.Detrend.Declustered.df$OsWL>oswl.thres)
rain = rep(NA,length(x.con))
tc = rep(NA,length(x.con))
for(j in 1:length(x.con)){
  rain[j] = ifelse(any(is.na(S22.Detrend.df$Rainfall[(x.con[j]-3):(x.con[j]+3)])==T),NA,sum(S22.Detrend.df$Rainfall[(x.con[j]-3):(x.con[j]+3)]))
  tc[j] = ifelse(any(!is.na(S22.Detrend.df$Name[(x.con[j]-3):(x.con[j]+3)])),"tc","nontc")
}


Sample_df = na.omit(data.frame(rain,S22.Detrend.Declustered.df$OsWL[x.con],tc,x.con))
colnames(Sample_df) = c("Rainfall_vol","OsWL","type","x.con")

con.sample.oswl.tc = list(Threshold = oswl.thres, Data = Sample_df[which(Sample_df$type=="tc"),1:2], Con_Variable = "OsWL",  x.con = Sample_df$x.con[which(Sample_df$type=="tc")])
con.sample.oswl.nontc = list(Threshold = oswl.thres, Data = Sample_df[which(Sample_df$type=="nontc"),1:2], Con_Variable = "OsWL",  x.con = Sample_df$x.con[which(Sample_df$type=="nontc")])

cop.oswl.tc <- BiCopSelect(pobs(con.sample.oswl.tc$Data$Rainfall_vol), pobs(con.sample.oswl.tc$Data$OsWL), familyset = NA, selectioncrit = "AIC",
                           indeptest = FALSE, level = 0.05, weights = NA,
                           rotations = TRUE, se = FALSE, presel = TRUE,
                           method = "mle")$family

cop.oswl.nontc <- BiCopSelect(pobs(con.sample.oswl.nontc$Data$Rainfall_vol), pobs(con.sample.oswl.nontc$Data$OsWL), familyset = NA, selectioncrit = "AIC",
                              indeptest = FALSE, level = 0.05, weights = NA,
                              rotations = TRUE, se = FALSE, presel = TRUE,
                              method = "mle")$family

#Find events in both samples for a given population
z1 = numeric(nrow(con.sample.oswl.tc$Data))
z = numeric(nrow(con.sample.rainfall.tc$Data))
for(i in 1:nrow(con.sample.rainfall.tc$Data)){
  for(j in 1:nrow(con.sample.oswl.tc$Data)){
    z1[j] <- ifelse(con.sample.rainfall.tc$x.con[i] > (con.sample.oswl.tc$x.con[j]-3) & con.sample.rainfall.tc$x.con[i] < (con.sample.oswl.tc$x.con[j]+3),1,0)
  }
  z[i] = ifelse(sum(z1)>0,1,0)
}

n_both_1 = sum(z)

z1 = numeric(nrow(con.sample.oswl.nontc$Data))
z = numeric(nrow(con.sample.rainfall.nontc$Data))
for(i in 1:nrow(con.sample.rainfall.nontc$Data)){
  for(j in 1:nrow(con.sample.oswl.nontc$Data)){
    z1[j] <- ifelse(con.sample.rainfall.nontc$x.con[i] > (con.sample.oswl.nontc$x.con[j]-3) & con.sample.rainfall.nontc$x.con[i] < (con.sample.oswl.nontc$x.con[j]+3),1,0)
  }
  z[i] = ifelse(sum(z1)>0,1,0)
}

n_both_2 = sum(z)

#TC sample
Diag_Non_Con(Data=con.sample.rainfall.tc$Data$OsWL, x_lab="O-sWL (in the sample con. on rainfall volume)",
             y_lim_min = 0, y_lim_max = 2.0)
mar.1  = "Gaus"
con.sample.oswl.tc$Data$Rainfall_vol<-con.sample.oswl.tc$Data$Rainfall_vol+runif(length(con.sample.oswl.tc$Data$Rainfall_vol),0.0001,0.001)
Diag_Non_Con_Trunc(Data=con.sample.oswl.tc$Data$Rainfall_vol,x_lab="Rainfall volume (in the sample con. on O-sWL)",
                   y_lim_min = 0, y_lim_max = 0.5, Omit="Weib")
mar.2 = "Exp"


#Non-TC sample
Diag_Non_Con(Data=con.sample.rainfall.nontc$Data$OsWL, x_lab="O-sWL (in the sample con. on rainfall volume)",
             y_lim_min = 0, y_lim_max = 2.0)
mar.3 = "Gaus"
con.sample.oswl.nontc$Data$Rainfall_vol<-con.sample.oswl.nontc$Data$Rainfall_vol+runif(length(con.sample.oswl.nontc$Data$Rainfall_vol),0.0001,0.001)
Diag_Non_Con_Trunc(Data=con.sample.oswl.nontc$Data$Rainfall_vol, x_lab="Rainfall vol (in the sample con. on O-sWL)",
                   y_lim_min = 0, y_lim_max = 2.0)
mar.4 = "Gam(2)"

#Marginal GPDs for TC sample
rainfall.exe = which(S22.Detrend.df.declust$Rainfall_vol > rainfall.thres)
rainfall.type = rep(NA,length(rainfall.exe))
for(i in 1:length(rainfall.type)){
  rainfall.type[i] = ifelse(any(!is.na(S22.Detrend.df$Name[(max(1,(rainfall.exe[i]-3))):(min(rainfall.exe[i]+3,nrow(S22.Detrend.df)))])),"tc","nontc")
}
declust.rainfall.tc = S22.Detrend.df.declust$Rainfall_vol[rainfall.exe[which(rainfall.type=="tc")]]
gpd.rainfall.tc= GPD_Fit(Data = declust.rainfall.tc, Data_Full = na.omit(S22.Detrend.df.declust$Rainfall_vol), u=NA, Thres =rainfall.thres , mu = 365.25,  Method = "Standard", PLOT=F)
rate.rainfall.tc = length(declust.rainfall.tc)/(length(na.omit(S22.Detrend.df.declust$Rainfall_vol))/(365.25))
gpd.rainfall.tc$Rate = rate.rainfall.tc

x.con = which(S22.Detrend.Declustered.df$OsWL>oswl.thres)
oswl.type  = rep(NA,length(x.con))
for(i in 1:length(oswl.type)){
  oswl.type[i] = ifelse(any(!is.na(S22.Detrend.df$Name[(x.con[i]-10):(x.con[i]+10)])),"tc","nontc")
}
declust.oswl.tc = S22.Detrend.Declustered.df$OsWL[x.con[which(oswl.type=="tc")]] #con.sample.wl.tc$Data$wl
gpd.oswl.tc = GPD_Fit(Data = declust.oswl.tc, Data_Full = na.omit(S22.Detrend.df$OsWL), u=NA, Thres = oswl.thres , mu = 365.25, Method = "Standard", PLOT=F)
rate.oswl.tc = length(declust.oswl.tc)/(length(na.omit(S22.Detrend.df$OsWL))/(365.25))
gpd.oswl.tc$Rate = rate.oswl.tc

#Non-TC sample
declust.rainfall.nontc =  S22.Detrend.df.declust$Rainfall_vol[rainfall.exe[which(rainfall.type=="nontc")]]
gpd.rainfall.nontc= GPD_Fit(Data = declust.rainfall.nontc, Data_Full = na.omit(S22.Detrend.df.declust$Rainfall_vol), u=NA, Thres =rainfall.thres , mu = 365.25, Method = "Standard", PLOT=F)
rate.rainfall.nontc = length(declust.rainfall.nontc)/(length(na.omit(S22.Detrend.df.declust$Rainfall_vol))/(365.25))
gpd.rainfall.nontc$Rate = rate.rainfall.nontc

declust.oswl.nontc =  S22.Detrend.Declustered.df$OsWL[x.con[which(oswl.type=="nontc")]]
gpd.oswl.nontc = GPD_Fit(Data = declust.oswl.nontc, Data_Full = na.omit(S22.Detrend.df$OsWL), u=NA, Thres = oswl.thres , mu = 365.25, Method = "Standard", PLOT=F, min.RI=1)
rate.oswl.nontc = length(declust.oswl.nontc)/(length(na.omit(S22.Detrend.df$OsWL))/(365.25))
gpd.oswl.nontc$Rate = rate.oswl.nontc

#Renaming columns of S22.Detrend.df (so function will not return error message)
colnames(S22.Detrend.df) = c("Date","Rainfall_vol","OsWL","Groundwater")

test_that("Design_Event_2D_Multi_Pop works", {

  result <- Design_Event_2D_Multi_Pop(Data=S22.Detrend.df[,c(1:3)],
                                      Data_Con1=con.sample.rainfall.tc$Data, Data_Con2=con.sample.rainfall.tc$Data,
                                      Data_Con3=con.sample.oswl.tc$Data, Data_Con4=con.sample.oswl.tc$Data,
                                      u1=NA,u2=NA,u3=NA,u4=NA,
                                      N_Both_1=n_both_1, N_Both_2=n_both_2,
                                      Copula_Family1=cop.rainfall.tc, Copula_Family2=cop.rainfall.nontc,
                                      Copula_Family3=cop.oswl.nontc, Copula_Family4=cop.oswl.nontc,
                                      Marginal_Dist1=mar.1, Marginal_Dist2=mar.2,
                                      Marginal_Dist3=mar.3, Marginal_Dist4=mar.4,
                                      GPD1 = gpd.rainfall.tc,
                                      Rate_Con1 = rate.rainfall.tc,
                                      GPD2 = gpd.oswl.tc,
                                      Rate_Con2 = rate.oswl.tc,
                                      GPD3 = gpd.rainfall.nontc,
                                      Rate_Con3 = rate.rainfall.nontc,
                                      GPD4 = gpd.oswl.nontc,
                                      Rate_Con4 = rate.oswl.nontc,
                                      Con1 = "Rainfall_vol", Con2 = "OsWL",
                                      Con3 = "Rainfall_vol", Con4 = "OsWL",
                                      Grid_x_min = 0 ,Grid_x_max = 700, Grid_y_min = -2,
                                      Grid_y_max = 3, Grid_x_interval=0.1, Grid_y_interval=0.01,
                                      x_lab = "Rainfall (in)", y_lab = "Water level (ft NGVD88)",
                                      RP=100,N=10^4,N_Ensemble=10,
                                      Plot_Quantile_Isoline=FALSE,
                                      x_lim_min = 0,
                                      x_lim_max = 100,
                                      y_lim_min = -2,
                                      y_lim_max = 10)

  # Checking type of output
  expect_type(result, 'list')
  expect_named(result, c("FullDependence", "MostLikelyEvent", "Ensemble", "Isoline","Contour",
                         "Quantile_Isoline_1","Quantile_Isoline_2","Threshold_1","Threshold_2"))


  #Checking length of outputs
  expect_true(nrow(result$FullDependence$`100`) == 1 & ncol(result$FullDependence$`100`)==2)
  expect_true(nrow(result$MostLikelyEvent$`100`) == 1 & ncol(result$MostLikelyEvent$`100`)==2)
  expect_equal(nrow(result$Ensemble$`100`),10)
  expect_equal(nrow(result$Isoline$`100`),1915)
  expect_equal(length(result$Contour$`100`),1915)
  expect_equal(length(result$Threshold_1),1)
  expect_equal(length(result$Threshold_2),1)

  #Checking column names of outputs
  expect_equal(colnames(S22.Detrend.df[,2:3]), colnames(result$FullDependence$`100`))
  expect_equal(colnames(S22.Detrend.df[,2:3]), colnames(result$MostLikelyEvent$`100`))
  expect_equal(colnames(S22.Detrend.df[,2:3]), colnames(result$Ensemble$`100`))
  expect_equal(colnames(S22.Detrend.df[,2:3]), colnames(result$Isoline$`100`))
})

test_that("Invalid inputs", {

  expect_error(
    Design_Event_2D_Multi_Pop(Data=S22.Detrend.df[,c(1:3)],
                              Data_Con1=con.sample.rainfall.tc$Data, Data_Con2=con.sample.rainfall.tc$Data,
                              Data_Con3=con.sample.oswl.tc$Data, Data_Con4=con.sample.oswl.tc$Data,
                              u1=NA,u2=NA,u3=NA,u4=NA,
                              N_Both_1=n_both_1, N_Both_2=n_both_2,
                              Copula_Family1=105, Copula_Family2=cop.rainfall.nontc,
                              Copula_Family3=cop.oswl.nontc, Copula_Family4=cop.oswl.nontc,
                              Marginal_Dist1=mar.1, Marginal_Dist2=mar.2,
                              Marginal_Dist3=mar.3, Marginal_Dist4=mar.4,
                              GPD1 = gpd.rainfall.tc,
                              Rate_Con1 = rate.rainfall.tc,
                              GPD2 = gpd.oswl.tc,
                              Rate_Con2 = rate.oswl.tc,
                              GPD3 = gpd.rainfall.nontc,
                              Rate_Con3 = rate.rainfall.nontc,
                              GPD4 = gpd.oswl.nontc,
                              Rate_Con4 = rate.oswl.nontc,
                              Con1 = "Rainfall_vol", Con2 = "OsWL",
                              Con3 = "Rainfall_vol", Con4 = "OsWL",
                              Grid_x_min = 0 ,Grid_x_max = 700, Grid_y_min = -2,
                              Grid_y_max = 3, Grid_x_interval=0.1, Grid_y_interval=0.01,
                              x_lab = "Rainfall (in)", y_lab = "Water level (ft NGVD88)",
                              RP=100,N=10^4,N_Ensemble=10,
                              Plot_Quantile_Isoline=FALSE,
                              x_lim_min = 0,
                              x_lim_max = 100,
                              y_lim_min = -2,
                              y_lim_max = 10),
    "Invalid Copula_Family1.")

  expect_error(
    Design_Event_2D_Multi_Pop(Data=S22.Detrend.df[,c(1:3)],
                              Data_Con1=con.sample.rainfall.tc$Data, Data_Con2=con.sample.rainfall.tc$Data,
                              Data_Con3=con.sample.oswl.tc$Data, Data_Con4=con.sample.oswl.tc$Data,
                              u1=NA,u2=NA,u3=NA,u4=NA,
                              N_Both_1=n_both_1,
                              Copula_Family1=cop.rainfall.tc, Copula_Family2=cop.rainfall.nontc,
                              Copula_Family3=cop.oswl.nontc, Copula_Family4=cop.oswl.nontc,
                              Marginal_Dist1=mar.1, Marginal_Dist2=mar.2,
                              Marginal_Dist3=mar.3, Marginal_Dist4=mar.4,
                              GPD1 = gpd.rainfall.tc,
                              Rate_Con1 = rate.rainfall.tc,
                              GPD2 = gpd.oswl.tc,
                              Rate_Con2 = rate.oswl.tc,
                              GPD3 = gpd.rainfall.nontc,
                              Rate_Con3 = rate.rainfall.nontc,
                              GPD4 = gpd.oswl.nontc,
                              Rate_Con4 = rate.oswl.nontc,
                              Con1 = "Rainfall_vol", Con2 = "OsWL",
                              Con3 = "Rainfall_vol", Con4 = "OsWL",
                              Grid_x_min = 0 ,Grid_x_max = 700, Grid_y_min = -2,
                              Grid_y_max = 10, Grid_x_interval=0.1, Grid_y_interval=0.01,
                              x_lab = "Rainfall (in)", y_lab = "Water level (ft NGVD88)",
                              RP=100,N=10^4,N_Ensemble=10,
                              Plot_Quantile_Isoline=FALSE,
                              x_lim_min = 0,
                              x_lim_max = 100,
                              y_lim_min = -2,
                              y_lim_max = 10),
    "N_Both_2 parameter is required and cannot be missing.")
})
