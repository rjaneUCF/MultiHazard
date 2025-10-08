#a
library('MultiHazard')
library('dplyr')
library('gamlss')
library('gamlss.mx')

rain_example_df = read.csv('F:/OneDrive - University of Central Florida/Documents/Shiny/Data/Miami_Airport_Rainfall_df.csv')[,2:3]
tw_example_df = read.csv('F:/OneDrive - University of Central Florida/Documents/Shiny/Data/S22_Tailwater_df.csv')[,2:3]

#Detrend
rainfall_det_df  = rain_example_df
rainfall_det_df$Date = as.Date(rainfall_det_df$Date)
colnames(rainfall_det_df)= c("Date","rain")

tw_det = Detrend(Data=tw_example_df,Method="Linear",Window_Width= NA, End_Length = NA, PLOT=FALSE)
tw_det_df  = data.frame(tw_example_df[,1],tw_det)
colnames(tw_det_df)= c("Date","tw")
tw_det_df$Date = as.Date(tw_det_df$Date)

#Decluster series
decl = Decluster(Data=rainfall_det_df[,2], u = 0.95, SepCrit = 3, mu = 365.25)
rain_dec_df = data.frame(as.Date(rainfall_det_df[,1]),decl$Declustered)
colnames(rain_dec_df) = c("Date","rain")
rain_dec_df$Date = as.Date(rain_dec_df$Date)

decl = Decluster(Data=tw_det_df[,2], u = 0.95, SepCrit = 3, mu = 365.25)
tw_dec_df = data.frame(as.Date(tw_det_df[,1]),decl$Declustered)
colnames(tw_dec_df) = c("Date","tw")
tw_dec_df$Date = as.Date(tw_dec_df$Date)

#Dataframes
df = data.frame("Date"=seq.Date(as.Date(min(rainfall_det_df$Date,tw_det_df$Date)),as.Date(max(rainfall_det_df$Date,tw_det_df$Date)),by="day"))
df_1 = left_join(df,rainfall_det_df,by="Date")
data_detrend_df = left_join(df_1,tw_det_df,by="Date")
colnames(data_detrend_df) = c("Date","rain","tw")

df = data.frame("Date"=seq.Date(as.Date(min(rainfall_det_df$Date,tw_det_df$Date)),as.Date(max(rainfall_det_df$Date,tw_det_df$Date)),by="day"))
df_1 = left_join(df,rain_dec_df,by="Date")
data_decl_df = left_join(df_1,tw_dec_df,by="Date")
colnames(data_decl_df) = c("Date","rain","tw")

con.sample.rain = Con_Sampling_2D(Data_Detrend=data_detrend_df, Data_Declust=data_decl_df, 
                                  Con_Variable="rain", u = 0.985)

con.sample.tw = Con_Sampling_2D(Data_Detrend=data_detrend_df, Data_Declust=data_decl_df, 
                                Con_Variable="tw", u = 0.985)

cop.rain = Copula_Threshold_2D(Data_Detrend=data_detrend_df,Data_Declust=data_decl_df,u1 = 0.985)$Copula_Family_Var1
cop.tw = Copula_Threshold_2D(Data_Detrend=data_detrend_df,Data_Declust=data_decl_df,u2 = 0.985)$Copula_Family_Var2

  
#The Diag_Non_Con_Trunc() function assesses the fit of multiple non-extreme value marginal distributions for bounded  variables.
Diag_Non_Con(Data=con.sample.rain$Data$tw, x_lab="Tail water level (in the sample con. on Rainfall)",
             y_lim_min = 0, y_lim_max = 3.0, Omit = c("RGum","Gum"))
#Gaussian best fitting
mar.1  = "Logis"

#The Diag_Non_Con_Trunc() function assesses the fit of multiple non-extreme value marginal distributions for bounded  variables.
con.sample.tw$Data$rain<-con.sample.tw$Data$rain+runif(length(con.sample.tw$Data$rain),0.0001,0.001)
Diag_Non_Con_Trunc(Data=con.sample.tw$Data$rain, x_lab="Rainfall (in the sample con. on Surge)",
                   y_lim_min = 0, y_lim_max = 4.0)
mar.2 = "Gam(2)"

#Marginal GPDs
declust.rain = data_decl_df$rain[data_decl_df$rain>con.sample.rain$Threshold]
gpd.rain= GPD_Fit(Data = declust.rain, Data_Full = na.omit(data_detrend_df$rain), u=NA, Thres =con.sample.rain$Threshold , Method = "Standard")
rate.rain = length(declust.rain)/(length(na.omit(data_detrend_df$rain))/(365.25))
gpd.rain$Rate = rate.rain

declust.tw = data_decl_df$tw[which(data_decl_df$tw>=con.sample.tw$Threshold)]
gpd.tw = GPD_Fit(Data = declust.tw, Data_Full = na.omit(data_detrend_df$tw), u=NA, Thres = con.sample.tw$Threshold , Method = "Standard")
rate.tw = length(declust.tw)/(length(na.omit(data_detrend_df$tw))/(365.25))
gpd.tw$Rate = rate.tw

Design.event<-Design_Event_2D(Data=data_detrend_df[,2:3],
                              Data_Con1= con.sample.rain$Data,Data_Con2=con.sample.tw$Data,
                              Thres1=con.sample.rain$Threshold, Thres2=con.sample.tw$Threshold,
                              Copula_Family1=cop.rain,Copula_Family2=cop.tw,
                              Marginal_Dist1=mar.1, Marginal_Dist2=mar.2,
                              GPD1 = gpd.rain,
                              Rate_Con1 = rate.rain,
                              GPD2 = gpd.tw,
                              Rate_Con2 = rate.tw,
                              Con1 = "rain", Con2 = "tw", 
                              x_lab = "Rainfall (in)", y_lab = "Water level (ft NGVD88)",
                              mu = 365.25,
                              x_lim_min = 0,
                              x_lim_max = 20,
                              y_lim_min = -2,
                              y_lim_max = 8,
                              RP=100,
                              N=10^1)


Design.event.grid<-Design_Event_2D_Grid(Data=data_detrend_df[,2:3],
                                        Data_Con1= con.sample.rain$Data,Data_Con2=con.sample.tw$Data,
                                        Thres1=con.sample.rain$Threshold, Thres2=con.sample.tw$Threshold,
                                        Copula_Family1=cop.rain,Copula_Family2=cop.tw,
                                        Marginal_Dist1=mar.1, Marginal_Dist2=mar.2,
                                        GPD1 = gpd.rain,
                                        Rate_Con1 = rate.rain,
                                        GPD2 = gpd.tw,
                                        Rate_Con2 = rate.tw,
                                        Grid_x_min = 0,
                              Grid_x_max = 20,
                              Grid_y_min = 0.5,
                              Grid_y_max = 8,
                              Grid_x_interval = 0.01,
                              Grid_y_interval = 0.01,
                              Con1 = "rain", Con2 = "tw", 
                              x_lab = "Rainfall (in)", y_lab = "Water level (ft NGVD88)",
                              mu = 365.25,
                              x_lim_min = 0,
                              x_lim_max = 20,
                              y_lim_min = -2,
                              y_lim_max = 8,
                              RP=100,
                              N=10^6)
  

##Plots
png("F:/OneDrive - University of Central Florida/Documents/JOSS/Figure_1_draft.png",width=9,height=9,res=300,units='in')
par(mfrow=c(2,2))
par(mar=c(4.2,4.5,0.5,0.5))
plot(0,type='n',xlim=c(15,500),ylim=c(0.22,2.5),xaxt='n',yaxt='n',xlab="",ylab="")
axis(1,seq(0,500,100),labels=NA)
mtext(seq(0,500,100),side=1,at=seq(0,500,100),line=0.7)
axis(2,seq(0.5,2.5,0.5),labels=NA)
mtext(seq(0.5,2.5,0.5),side=2,at=seq(0.5,2.5,0.5),line=0.7)

text(13.5,2.52,'(a)')
mtext(side=1,'Rainfall (mm)',line=2.35)
mtext(side=2, 'Water level (m NGVD29)',line=2.5)

points(25.4*data_detrend_df$rain,0.3048*data_detrend_df$tw,pch=16,col="Grey")
points(25.4*con.sample.rain$Data$rain,0.3048*con.sample.rain$Data$tw,pch=1,col="Blue")
points(25.4*con.sample.tw$Data$rain,0.3048*con.sample.tw$Data$tw,pch=4,col="Red")
points(25.4*Design.event$Quantile_Isoline_1$`100`$rain[which(Design.event$Quantile_Isoline_1$`100`$rain<12.5)],
      0.3048*Design.event$Quantile_Isoline_1$`100`$tw[which(Design.event$Quantile_Isoline_1$`100`$rain<12.5)],col="Blue",pch=16)
points(25.4*Design.event.grid$Isoline$`100`$rain[which(Design.event.grid$Isoline$`100`$rain>12.5)],
       0.3048*Design.event.grid$Isoline$`100`$tw[which(Design.event.grid$Isoline$`100`$rain>12.5)],col="Blue",pch=16)
points(25.4*Design.event$Quantile_Isoline_2$`100`$rain,0.3048*Design.event$Quantile_Isoline_2$`100`$tw,col="Red",pch=16)
box()

#b
plot(0,type='n',xlim=c(15,500),ylim=c(0.2,2.5),xaxt='n',yaxt='n',xlab="",ylab="")

text(13.5,2.52,'(b)')
axis(1,seq(0,500,100),labels=NA)
mtext(seq(0,500,100),side=1,at=seq(0,500,100),line=0.7)
axis(2,seq(0.5,2.5,0.5),labels=NA)
mtext(seq(0.5,2.5,0.5),side=2,at=seq(0.5,2.5,0.5),line=0.7)

mtext(side=1,'Rainfall (mm)',line=2.35)
mtext(side=2, 'Water level (m NGVD29)',line=2.5)
points(25.4*data_detrend_df$rain,0.3048*data_detrend_df$tw,pch=16,col="Grey")
points(25.4*con.sample.rain$Data$rain,0.3048*con.sample.rain$Data$tw,pch=1,col="Blue")
points(25.4*con.sample.tw$Data$rain,0.3048*con.sample.tw$Data$tw,pch=4,col="Red")
points(25.4*Design.event.grid$Isoline$`100`$rain,0.3048*Design.event.grid$Isoline$`100`$tw, col=rev(heat.colors(150))[20:120][1+100*Design.event.grid$Contour$`100`],pch=16)
points(25.4*Design.event.grid$MostLikelyEvent$`100`$rain,0.3048*Design.event.grid$MostLikelyEvent$`100`$tw,pch=18,cex=1.2)
box()


#cd
#Load necessary packages
library(dplyr)
library(MultiHazard)
library(copula)
library(VineCopula)
library(gamlss)
library(gamlss.mx)
library(lubridate)
library(texmex)

# Functions ---------------------------------------------------------------
Laplace_inverse = function(u){ #inverse function of standard Laplace cdf for HT model
  x = c()
  x[u<=0.5] = log(2*u[u<=0.5])
  x[u>0.5] = -log(2*(1-u[u>0.5]))
  return(x)
}

Laplace_cdf = function(x){ #standard Laplace cdf function for HT model
  u = c()
  u[x<0] = exp(x[x<0])/2
  u[x>=0] = 1-exp(-x[x>=0])/2
  return(u)
}

HeffTawnNegLL=function(X,Y,par){ #negative log likelihood for estimation of HT model parameters
  #X is conditioning variable, Y is other variable
  alpha=par[1]
  beta=par[2]
  sig=par[3]
  mu=par[4]
  if(alpha < -1 || alpha > 1 || beta > 1 || beta < -5 ){return(1e10)}
  if(sig<=0){return(1e10)}
  negloglik <- -sum(dnorm(Y,alpha*X+mu*((X)^beta),sig*((X)^beta),log=T))
  if(is.finite(negloglik)){
    return(negloglik)
  } else {
    return(1e10)
  }
}

heff_tawn_root = function(y,prob,Ypar,YZ,nsim){ #function for uniroot, aiming to find the point on the return curve where y=x
  q = 1 - Laplace_cdf(y)
  sample_y = rexp(nsim) + y
  sample_z = sample(YZ,nsim,replace = T)
  sample_x = (sample_y^Ypar[2])*sample_z + Ypar[1]*sample_y
  x = quantile(sample_x,(1-prob/q)) 
  return(x-y)
}

heff_tawn_curve_exp = function(data_exp_x,data_exp_y,prob,q,nsim){ #function for estimating return curve for data on standard exponential margins using HT model
  #prob is curve survival probability (small)
  #q is quantile cdf probability for fitting HT model
  #nsim is number of simulation for HT model
  
  #transform data_exp to Laplace margins
  
  data_u_x = apply(data_exp_x,2,pexp)
  data_x = apply(data_u_x,2,Laplace_inverse)
  
  data_u_y = apply(data_exp_y,2,pexp)
  data_y = apply(data_u_y,2,Laplace_inverse)
  
  #fitting the HT model, conditioning on one variable
  
  ux = quantile(data_x[,1],q)
  uy = quantile(data_y[,2],q)
  dataX = data_x[data_x[,1]>ux,]
  dataY = data_y[data_y[,2]>uy,]
  
  #conditioning on large values of Y
  Yopt = optim(fn=HeffTawnNegLL,X=dataY[,2],Y=dataY[,1],par=rep(1/2,4),control = list(maxit=10000)) #opposite way to make things confusing
  Ypar = Yopt$par
  YZ = (dataY[,1] - Ypar[1]*dataY[,2])/(dataY[,2]^Ypar[2])
  
  #conditioning on large values of X
  Xopt = optim(fn=HeffTawnNegLL,X=dataX[,1],Y=dataX[,2],par=rep(1/2,4),control = list(maxit=10000))
  Xpar = Xopt$par
  XZ = (dataX[,2] - Xpar[1]*dataX[,1])/(dataX[,1]^Xpar[2])
  
  #Finding points on curve where y=x
  test = tryCatch(uniroot(heff_tawn_root,interval = c(Laplace_inverse(q),Laplace_inverse(1-prob-10^(-10))),prob=prob,Ypar=Ypar,YZ=YZ,nsim=10000),error = function(e){1})
  test2 = tryCatch(uniroot(heff_tawn_root,interval = c(Laplace_inverse(q),Laplace_inverse(1-prob-10^(-10))),prob=prob,Ypar=Xpar,YZ=XZ,nsim=10000),error = function(e){1})
  
  if(is.list(test) & is.list(test2)){ #if roots can be found for both conditioning variables
    
    y = test$root #extract values
    x = test2$root
    
    val = mean(c(x,y))#take the average of these values and use this as the crossover point from region 1 to region 2
    qvalsY = Laplace_cdf(val) #cdf probability value
    qvalsX = Laplace_cdf(val)
    
    qvalsY = exp(seq(log(1-prob-10^(-10)),log(qvalsY),length.out=100)) #create sequence of exceedance probabilities for Y variable in region 1
    y = Laplace_inverse(qvalsY) #find corresponding values of Y
    x = rep(NA,length.out=100) 
    for(i in 1:length(y)){ #find corresponding values of X variable using fitted model
      sample_y = rexp(nsim) + y[i]
      sample_z = sample(YZ,nsim,replace = T)
      sample_x = (sample_y^Ypar[2])*sample_z + Ypar[1]*sample_y
      x[i] = quantile(sample_x,(1-prob/(1-qvalsY[i]))) 
    }
    
  } else { #if rootfinder gives an error, we must try to find the exceedance probability manually or just consider y values above the q-th quantile
    
    qvalsY = exp(seq(log(1 - prob - 10^(-10)),log(min(1-sqrt(prob),q)),length.out=1000)) #exceedance probabilities
    y = Laplace_inverse(qvalsY) #correspoding values of y
    x = rep(NA,length.out = 1000) 
    i = 1
    sample_y = rexp(nsim) + y[i]
    sample_z = sample(YZ,nsim,replace = T)
    sample_x = (sample_y^Ypar[2])*sample_z + Ypar[1]*sample_y
    x[i] = quantile(sample_x,(1-prob/(1-qvalsY[i])))
    
    while(y[i]>x[i] & i<1000){ #finding the crossover point where we go from region 1 to region 2 (if it exists)
      #using model to obtain estimates of x values on curve
      i = i+1
      sample_y = rexp(nsim) + y[i]
      sample_z = sample(YZ,nsim,replace = T)
      sample_x = (sample_y^Ypar[2])*sample_z + Ypar[1]*sample_y
      x[i] = quantile(sample_x,(1-prob/(1-qvalsY[i]))) 
    }
    
    qvalsX = Laplace_cdf(x[i])
    qvalsY = exp(seq(log(qvalsY[1]),log(qvalsY[i-1]),length.out=100)) #selecting probabilities up to the crossover point
    y = Laplace_inverse(qvalsY) #corresponding values of y
    x = rep(NA,length.out=100) 
    for(i in 1:length(y)){#using model to obtain estimates of x values on curve
      sample_y = rexp(nsim) + y[i]
      sample_z = sample(YZ,nsim,replace = T)
      sample_x = (sample_y^Ypar[2])*sample_z + Ypar[1]*sample_y
      x[i] = quantile(sample_x,(1-prob/(1-qvalsY[i]))) 
    }
  }
  
  qvalsX = exp(seq(log(max(qvalsX,q)),log(1 - prob - 10^(-10)),length.out=100)) #selecting exceedance probabilities of X variable such that curve point is in region 2
  x = c(x, Laplace_inverse(qvalsX))
  y = c(y,rep(NA,length.out = length(qvalsX))) 
  
  for(i in (length(qvalsX)+1):(2*length(qvalsX))){#using model to obtain estimates of y values on curve
    sample_x = rexp(nsim) + x[i]
    sample_z = sample(XZ,nsim,replace = T)
    sample_y = (sample_x^Xpar[2])*sample_z + Xpar[1]*sample_x
    y[i] = quantile(sample_y,(1-prob/(1-qvalsX[i-length(qvalsX)])),na.rm = T) 
  }
  
  #transforming the curve back to standard exponential margins
  curve = cbind(x,y)
  curve = apply(curve,2,Laplace_cdf)
  curve = apply(curve,2,qexp)
  
  #imposing property 3.1
  end = qexp(1-prob)
  curve = rbind(c(0,end),curve,c(end,0))
  
  #imposing property 3.4 
  for(i in 2:(length(curve[,1])-1)){
    if(curve[i+1,1]<curve[i,1]){
      curve[i+1,1] = curve[i,1]
    }
    if(curve[i+1,2]>curve[i,2]){
      curve[i+1,2] = curve[i,2]
    }
  }
  return(curve)#returns matrix containing curve estimate
}

EmpFun <- function(x, r, mod) {
  
  X<-x
  #Empirical distribution function
  F <-ecdf(x)
  Femp<-F(r)
  
  #Empirical distribution function
  ox <- order(x)
  names(x) <- 1:length(x)
  x <- sort(x)
  run <- rle(x)
  p <- cumsum(run$lengths)/length(x)  #divisor
  p <- rep(p, run$lengths)
  p <- p[order(as.character(names(x)))]
  x <- x[order(as.character(names(x)))]
  #Femp <- p
  
  #dpareto
  th <-mod$threshold
  sigma <- exp(mod$coefficients[1])
  xi <- mod$coefficients[2]
  q <- F(th)
  Para <- 1-(1-q) * pgpd(r,u=th,sigma=exp(mod$coefficients[1]),xi=mod$coefficients[2],lower.tail = FALSE)
  
  #Pareto distribution above threshold and emprical distribution below
  res <- Femp
  res[r>th] <- Para[r>th]
  
  #Ensuring results are in the order of x
  res <- as.numeric(res)
  res
}

u2gpd <-
  function( u , p = 1, th=0, sigma, xi){
    ( ( ( 1 - u )/p )^(-xi) -1 ) * sigma/xi + th
  }


revTransform <-
  function (x, data, qu, th = 0, sigma = 1, xi = 0, method = "mixture") {
    if (!is.element(method, c("mixture", "empirical")))
      stop("method should be 'mixture' or 'empirical'")
    data<-na.omit(data)
    qu<-ecdf(na.omit(data))(th)
    n <- length(data)
    probs <- (1:n)/(n + 1)
    px <- sapply(x, function(x, p) p[abs(x - p) == min(abs(x - p))][1], p = probs) # take 1st item in case of ties
    px <- as.integer(round(px * (1 + n)))
    res <- sort(data)[px]
    if (method == "mixture") {
      res[res > th] <- u2gpd(x[res > th], p=1-qu, th = th, sigma = sigma, xi = xi)
    }
    res[order(x)] <- sort(res)
    res
  }

##Pre-processing data

rain_example_df = read.csv('F:/OneDrive - University of Central Florida/Documents/Shiny/Data/Miami_Airport_Rainfall_df.csv')[,2:3]
tw_example_df = read.csv('F:/OneDrive - University of Central Florida/Documents/Shiny/Data/S22_Tailwater_df.csv')[,2:3]

#Convert to metric units
#rain_example_df$Value = 25.4*rain_example_df$Value
#tw_example_df$ValueFilled = 0.3048*tw_example_df$ValueFilled

#Detrend
rainfall_det_df  = rain_example_df
rainfall_det_df$Date = as.Date(rainfall_det_df$Date)
colnames(rainfall_det_df)= c("Date","rain")

tw_det = Detrend(Data=tw_example_df,Method="Linear",Window_Width= NA, End_Length = NA, PLOT=FALSE)
tw_det_df  = data.frame(tw_example_df[,1],tw_det)
colnames(tw_det_df)= c("Date","tw")
tw_det_df$Date = as.Date(tw_det_df$Date)

#Decluster series
decl = Decluster(Data=rainfall_det_df[,2], u = 0.95, SepCrit = 3, mu = 365.25)
rain_dec_df = data.frame(as.Date(rainfall_det_df[,1]),decl$Declustered)
colnames(rain_dec_df) = c("Date","rain")
rain_dec_df$Date = as.Date(rain_dec_df$Date)

decl = Decluster(Data=tw_det_df[,2], u = 0.95, SepCrit = 3, mu = 365.25)
tw_dec_df = data.frame(as.Date(tw_det_df[,1]),decl$Declustered)
colnames(tw_dec_df) = c("Date","tw")
tw_dec_df$Date = as.Date(tw_dec_df$Date)

#Dataframes
df = data.frame("Date"=seq.Date(as.Date(min(rainfall_det_df$Date,tw_det_df$Date)),as.Date(max(rainfall_det_df$Date,tw_det_df$Date)),by="day"))
df_1 = left_join(df,rainfall_det_df,by="Date")
data_detrend_df = left_join(df_1,tw_det_df,by="Date")
colnames(data_detrend_df) = c("Date","rain","tw")

df = data.frame("Date"=seq.Date(as.Date(min(rainfall_det_df$Date,tw_det_df$Date)),as.Date(max(rainfall_det_df$Date,tw_det_df$Date)),by="day"))
df_1 = left_join(df,rain_dec_df,by="Date")
data_decl_df = left_join(df_1,tw_dec_df,by="Date")
colnames(data_decl_df) = c("Date","rain","tw")

prob = (1/365.25) / 100  #10^(-3) #curve probabilities
q = 0.985 #quantile level for fitting GPDs and estimating extremal dependence structures
nsim = 100000 #number of simulations for HT model

#gas data
data_x = data_detrend_df[,-1] #data.frame(data_decl_df$rain,data_detrend_df$tw)
data_y = data_detrend_df[,-1]  #data.frame(data_detrend_df$rain,data_decl_df$tw)

colnames(data_x) = c("rain","tw")
colnames(data_y) = c("rain","tw")

data_x = na.omit(data_x)
data_y = na.omit(data_y)

#Fit GPDs
thresh_x = quantile(data_detrend_df$rain,q,na.rm=T)
gpd_x = evm(data_decl_df$rain,th=thresh_x, penalty = "gaussian", 
            priorParameters = list(c(0, 0), matrix(c(100^2, 0, 0, 0.25), nrow = 2)))
gpd_x$rate =  length(data_decl_df$rain[which(data_decl_df$rain>=con.sample.rain$Threshold)])/(length(na.omit(data_detrend_df$rain))/(365.25))

thresh_y = quantile(data_detrend_df$tw,q,na.rm=T)
gpd_y = evm(data_decl_df$tw,th=thresh_y, penalty = "gaussian", 
            priorParameters = list(c(0, 0), matrix(c(100^2, 0, 0, 0.25), nrow = 2)))
gpd_y$rate = length(data_decl_df$tw[which(data_decl_df$tw>=con.sample.tw$Threshold)])/(length(na.omit(data_detrend_df$tw))/(365.25))

#Converting to uniform scale
data_unif_x = cbind(EmpFun(x = na.omit(data_detrend_df$rain), r = data_x[,1], mod = gpd_x),EmpFun(x = na.omit(data_detrend_df$tw), r = data_x[,2],mod = gpd_y))
data_unif_y = cbind(EmpFun(x = na.omit(data_detrend_df$rain), r = data_y[,1], mod = gpd_x),EmpFun(x = na.omit(data_detrend_df$tw), r = data_y[,2],mod = gpd_y))

data_exp_x = apply(data_unif_x, 2, qexp)
data_exp_y = apply(data_unif_y, 2, qexp)

#estimated curves exponential margins
curve = heff_tawn_curve_exp(data_exp_x = data_exp_x, data_exp_y = data_exp_y, prob = prob, q=q, nsim=nsim)

#estimated curves uniform margins
curve_unif = apply(curve, 2, pexp)

#estimated curves original margins
curve_original = cbind(revTransform(x=curve_unif[,1],data = na.omit(data_detrend_df$rain),qu=q,th=thresh_x,exp(coefficients(gpd_x)[1]),coefficients(gpd_x)[2]),revTransform(x=curve_unif[,2],data = na.omit(data_detrend_df$tw),qu=q,th=thresh_y,exp(coefficients(gpd_y)[1]),coefficients(gpd_y)[2]))

#Monthly block bootstrap

#
extra = data.frame(seq(as.Date("2019-01-27"),as.Date("2019-01-31"),by="day"),rep("NA",5),rep("NA",5))
colnames(extra) = c("Date","rain","tw")

data_detrend_df = rbind(data_detrend_df,extra)

m.y = data.frame(month(data_detrend_df$Date),year(data_detrend_df$Date))
colnames(m.y) = c("month","year")

m.y.complete = data.frame(rep(1:12,length(min(m.y$year):max(m.y$year))),rep(min(m.y$year):max(m.y$year),each=12))

#indicator variable denoting which month has rainfall, tw or both data types avaliable
I = rep(NA, nrow(m.y.complete))

for(i in 1:nrow(m.y.complete)){
  d = data_detrend_df[which(month(data_detrend_df$Date) == m.y.complete[i,1] & year(data_detrend_df$Date) == m.y.complete[i,2]),]
  r = length(which(is.na(d$rain) == FALSE & is.na(d$tw))) / nrow(d)
  tw = length(which(is.na(d$rain) & is.na(d$tw)== FALSE)) / nrow(d)
  b = length(which(is.na(d$rain)  == FALSE & is.na(d$tw)  == FALSE)) / nrow(d)
  
  I[i] = ifelse(b>0.8,"B",ifelse(r>0.8,"R",ifelse(tw>0.8,"tw",NA)))
  
}

#Combinin information into a matrix
m.y.complete = data.frame(rep(1:12,length(min(m.y$year):max(m.y$year))),rep(min(m.y$year):max(m.y$year),each=12),I)
colnames(m.y.complete) = c("month","year","indicator")

z1 = rep(NA, length(seq(1948,2024,4)))
for(i in 1:length(seq(1948,2024,4))){ 
  z1[i] = ifelse(length(which(m.y.complete$month==2 & m.y.complete$year==seq(1948,2024,4)[i]))>0,which(m.y.complete$month==2 & m.y.complete$year==seq(1948,2024,4)[i]),NA)
}

m.y.complete$month[z1] = 13

boot_curves = list()
for(j in 1:100){
  
  #Forbootstrap, do not pick months where no data is avaliable
  z = which(is.na(m.y.complete$indicator))
  
  #Bootstrap by month
  boot.m.y = rep(NA, nrow(m.y.complete))
  
  for(i in 1:nrow(m.y.complete[-z,])){
    boot.m.y[i] = sample(which(m.y.complete$month ==  m.y.complete$month[i] &  (m.y.complete$indicator == m.y.complete$indicator[i] | m.y.complete$indicator == "B")),1)
  }
  
  #Bootsrap sample of months and years
  m.y.boot = m.y.complete[boot.m.y,1:2]
  
  #Composing bootstrap sample of rainfall and tw
  
  boot_df = data_detrend_df
  
  for(i in 1:nrow(m.y.complete[-z,])){
    boot_df[which(month(data_detrend_df$Date) == m.y.complete[i,1] & year(data_detrend_df$Date) == m.y.complete[i,2]),2:3] = data_detrend_df[which(month(data_detrend_df$Date) == m.y.boot[i,1] & year(data_detrend_df$Date) == m.y.boot[i,2]),2:3]
    print(i)
  }
  
  boot_df$rain = as.numeric(boot_df$rain)
  boot_df$tw = as.numeric(boot_df$tw) 
  
  #Decluster series
  decl = Decluster(Data=boot_df[,2], u = 0.95, SepCrit = 3, mu = 365.25)
  boot_rain_dec_df = data.frame(as.Date(boot_df[,1]),decl$Declustered)
  colnames(boot_rain_dec_df) = c("Date","rain")
  boot_rain_dec_df$Date = as.Date(boot_rain_dec_df$Date)
  
  decl = Decluster(Data=boot_df[,3], u = 0.95, SepCrit = 3, mu = 365.25)
  boot_tw_dec_df = data.frame(as.Date(boot_df[,1]),decl$Declustered)
  colnames(boot_tw_dec_df) = c("Date","tw")
  boot_tw_dec_df$Date = as.Date(boot_tw_dec_df$Date)
  
  #Dataframes
  df = data.frame("Date"=seq.Date(as.Date(min(boot_df$Date,boot_df$Date)),as.Date(max(boot_df$Date,boot_df$Date)),by="day"))
  df_1 = left_join(df,boot_rain_dec_df,by="Date")
  boot_decl_df = left_join(df_1,boot_tw_dec_df,by="Date")
  colnames(boot_decl_df) = c("Date","rain","tw")
  
  #gas data
  data_x = boot_df[,-1] #data.frame(data_decl_df$rain,data_detrend_df$tw)
  data_y = boot_df[,-1]  #data.frame(data_detrend_df$rain,data_decl_df$tw)
  
  colnames(data_x) = c("rain","tw")
  colnames(data_y) = c("rain","tw")
  
  data_x = na.omit(data_x)
  data_y = na.omit(data_y)
  
  #Fit GPDs
  thresh_x = quantile(boot_df$rain,q,na.rm=T)
  gpd_x = evm(boot_decl_df$rain,th=thresh_x, penalty = "gaussian", 
              priorParameters = list(c(0, 0), matrix(c(100^2, 0, 0, 0.25), nrow = 2)))
  gpd_x$rate =  length(boot_decl_df$rain[which(boot_decl_df$rain>=con.sample.rain$Threshold)])/(length(na.omit(boot_df$rain))/(365.25))
  
  thresh_y = quantile(boot_df$tw,q,na.rm=T)
  gpd_y = evm(boot_decl_df$tw,th=thresh_y, penalty = "gaussian", 
              priorParameters = list(c(0, 0), matrix(c(100^2, 0, 0, 0.25), nrow = 2)))
  gpd_y$rate = length(boot_decl_df$tw[which(boot_decl_df$tw>=con.sample.tw$Threshold)])/(length(na.omit(boot_df$tw))/(365.25))
  
  
  #Converting to uniform scale
  data_unif_x = cbind(EmpFun(x = boot_df$rain,r = data_x[,1],mod = gpd_x),EmpFun(x = boot_df$tw, r = data_x[,2], mod = gpd_y))
  data_unif_y = cbind(EmpFun(x = boot_df$rain,r = data_y[,1],mod = gpd_x),EmpFun(x = boot_df$tw, r = data_y[,2], mod = gpd_y))
  
  data_exp_x = apply(data_unif_x, 2, qexp)
  data_exp_y = apply(data_unif_y, 2, qexp)
  
  #estimated curves exponential margins
  curve = heff_tawn_curve_exp(data_exp_x = data_exp_x, data_exp_y = data_exp_y, prob = prob, q=q, nsim=nsim)
  
  #estimated curves uniform margins
  curve_unif = apply(curve, 2, pexp)
  
  #estimated curves original margins
  curve_original = cbind(revTransform(x=curve_unif[,1],data = na.omit(boot_df$rain),qu=q,th=thresh_x,exp(coefficients(gpd_x)[1]),coefficients(gpd_x)[2]),revTransform(x=curve_unif[,2],data = na.omit(boot_df$tw),qu=q,th=thresh_y,exp(coefficients(gpd_y)[1]),coefficients(gpd_y)[2]))
  
  boot_curves[[j]] = curve_original
}

data_detrend_df$rain = as.numeric(data_detrend_df$rain)
data_detrend_df$tw = as.numeric(data_detrend_df$tw)

plot(data_detrend_df[,c(2,3)],pch=16,col="grey",main="",xlab="",ylab="",xaxt='n',yaxt='n',xlim=c(15/25.4,500/25.4),ylim=c(0.18/0.3048,2.5/0.3048))
axis(1,seq(0,500,100)/25.4,labels=NA)
mtext(seq(0,500,100),side=1,at=seq(0,500,100)/25.4,line=0.7)
axis(2,seq(0.5,2.5,0.5)/0.3048,labels=NA)
mtext(seq(0.5,2.5,0.5),side=2,at=seq(0.5,2.5,0.5)/0.3048,line=0.7)

mtext(side=1,'Rainfall (mm)',line=2.35)
mtext(side=2, 'Water level (m NGVD29)',line=2.5)

abline(h=25.4*min(data_detrend_df$rain),lwd=3,col=2)
abline(v=0.3048*min(data_detrend_df$tw),lwd=3,col=2)

no_grad = 25 #number of gradients
angles = ((no_grad:1)/(no_grad+1))*(pi/2)
x0 = min(data_detrend_df$rain, na.rm=T)
y0 = min(data_detrend_df$tw,na.rm=T)
for(i in 1:no_grad){
  abline(a = y0-x0*tan(angles[i]),b=tan(angles[i]),lwd=2)
}
text(485/24.5,0.16/0.3048,'(c)')
box()

angles_est_points = array(0,dim=c(no_grad,2,100))

for(k in 1:100){
  curve_angles = atan((boot_curves[[k]][,2]-y0)/(boot_curves[[k]][,1]-x0))
  
  for(i in 1:no_grad){
    ind = min(which(angles[i] >= curve_angles))
    x1 = boot_curves[[k]][ind,1] - x0
    x2 = boot_curves[[k]][ind-1,1] - x0
    y1 = boot_curves[[k]][ind,2] - y0
    y2 = boot_curves[[k]][ind-1,2] - y0
    p = (x1*tan(angles[i]) - y1) / ( (y2-y1) - (x2-x1)*tan(angles[i]) )
    xhat = x1 + p*(x2-x1) + x0 #change if computing uncertainty
    yhat = y1 + p*(y2-y1) + y0
    angles_est_points[i,,k] = c(xhat,yhat)
  }
}


median = array(0,dim=c(no_grad,2))
lower_bound = array(0,dim=c(no_grad,2))
upper_bound = array(0,dim=c(no_grad,2))

for(i in 1:no_grad){
  median[i,] = angles_est_points[i,,order(angles_est_points[i,2,1:100])[50]]
  lower_bound[i,] = angles_est_points[i,,order(angles_est_points[i,2,1:100])[5]]
  upper_bound[i,] = angles_est_points[i,,order(angles_est_points[i,2,1:100])[95]]
}


plot(0,type='n',xlim=c(15,500),ylim=c(0.2,2.5),xaxt='n',yaxt='n',xlab="",ylab="")
points(25.4*data_detrend_df[,2],0.3048*data_detrend_df[,3],pch=16,col="grey")
axis(1,seq(0,500,100),labels=NA)
mtext(seq(0,500,100),side=1,at=seq(0,500,100),line=0.7)
axis(2,seq(0,2.5,0.5),labels=NA)
mtext(seq(0.5,2.5,0.5),side=2,at=seq(0.5,2.5,0.5),line=0.7)

mtext(side=1,'Rainfall (mm)',line=2.35)
mtext(side=2, 'Water level (m NGVD29)',line=2.5)

lines(25.4*median[,1],0.3048*median[,2],col="blue",lwd=2)
lines(25.4*upper_bound[,1],0.3048*upper_bound[,2],lty=2,lwd=2)
lines(25.4*lower_bound[,1],0.3048*lower_bound[,2],lty=2,lwd=2)

text(502,0.18,'(d)')
box()
dev.off()