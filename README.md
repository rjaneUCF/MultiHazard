# MultiHazard

```{r setup, include=FALSE}
library("MultiHazard")
library("dplyr")
library("scales")
library("statmod")
library("CDVine")
```
## 1. Introduction
The `MultiHazard` package provides tools for stationary multivariate statistical modeling such as of the joint distribution of MULTIple co-occurring HAZARDs. The package contains functions for pre-processing data including imputing missing values, detrending and declustering time series (Section 2) as well as analyzing pairwise correlations over a range of lags (Section 3). Functionality is also built in to conditionally sample a bivariate dataset (given one of the variables is above a predetermined threshold), selecting the best fitting amongst an array of parametric (extreme and non-extreme, truncated and non-truncated) marginal distributions or copulas, which allow the dependence structure between the set of variables to be modelled independenly of marginal distributions (Section 3). Estimation of joint probability contours using the method of overlaying (conditional) contours given in Bender et al. (2016) and subsequently for a given return period extracting design events assuming full dependence, as well as the "most likely" or an ensemble of possible design events once accounting for dependence is possible (Section 4). The package also provides the capability of fitting and simulating synthetic records from three higher dimensional approaches - standard (elliptic/Archimedean) copulas, Pair Copula Constructions (PCCs) and the conditional threshold exceedance approach of Heffernan and Tawn (2004) (Section 5).  

## 2. Pre-processing
### Imputation

Well G_3356 represents the groundwater level at Site S20, however, it contains missing values. Lets impute missing values in the record at Well G_3356 using recordings at nearby Well G_3355.  Firstly, reading in the two time series.
```{r}
#Viewing first few rows of in the groundwater level records
head(G_3356)
head(G_3355)
#Converting Date column to "Date"" object
G_3356$Date<-seq(as.Date("1985-10-23"), as.Date("2019-05-29"), by="day")
G_3355$Date<-seq(as.Date("1985-08-20"), as.Date("2019-06-02"), by="day")
#Converting column containing the readings to a "numeric"" object
G_3356$Value<-as.numeric(as.character(G_3356$Value))
G_3355$Value<-as.numeric(as.character(G_3355$Value))
```
Warning message confirms there are NAs in the record at Well G_3356. Before carrying out the imputation the two data frames need to be merged.
```{r}
#Merge the two dataframes by date
GW_S20<-left_join(G_3356,G_3355,by="Date")
colnames(GW_S20)<-c("Date","G3356","G3355")
#Carrying out imputation
Imp<-Imputation(Data=GW_S20,Variable="G3356",x_lab="G-3355 (ft NGVD 29)", y_lab="G-3356 (ft NGVD 29)")
```
The completed record is given in the `ValuesFilled` column of the data frame outputted as the `Data` object while the linear regression model including its coefficient of determinant are given by the `model` output argument. 
```{r}
head(Imp$Data)
Imp$Model
```

Are any values still NA?
  
```{r}
G_3356_ValueFilled_NA<-which(is.na(Imp$Data$ValuesFilled)==TRUE)
length(G_3356_ValueFilled_NA) 
```
Linear interpolating the three remaining NAs.

```{r}
G3356_approx<-approx(seq(1,length(Imp$Data$ValuesFilled),1),Imp$Data$ValuesFilled,
                     xout=seq(1,length(Imp$Data$ValuesFilled),1))
Imp$Data$ValuesFilled[which(is.na(Imp$Data$ValuesFilled)==TRUE)]<-
  G3356_approx$y[which(is.na(Imp$Data$ValuesFilled)==TRUE)]
```

### Detrending

In the analysis completed O-sWL (Ocean-side Water Level) and groundwater level series are subsequently detrended. The Detrend() function uses either a linear fit covering the entire data (Method=`linear`) or moving averge window (Method=`window`) of a specified length (`Window_Width`) to remove trends from a time series. The residuals are added to the final `End_Length` observations.
The default  Detrend() parameters specify a moving average (Method=`window`) three month window (Window_Width=`89`), to remove any seasonality from the time series. The deafult is then to add the residuals to the average of the final five years of observations (End_Length=`1826`) to bring the record to the present day level, accounting for the Perigian tide in the case of O-sWL. The mean of the observations over the first three months were subtracted from the values during this period before the present day (5 year) average was added. The following R code detrends the record at Well G_3356. Note the function requires a Date object and the completed series.
```{r}
#Cresaring a data from with the imputed series alongside the corresponding dates 
G_3356_Imp<-data.frame(Imp$Data$Date,Imp$Data$ValuesFilled)
colnames(G_3356_Imp)<-c("Date","ValuesFilled")
#Detrending
G_3356_Detrend<-Detrend(Data=G_3356_Imp,PLOT=TRUE,x_lab="Date",
                        y_lab="Groundwater level (ft NGVD 29)")
```
Output of the Detrend_3Month() function is simply the detrended time series.
```{r}
head(G_3356_Detrend)
```
Creating a data frame containing the detrended groundwater series at site S_20 i.e. G_3356_Detrend and their corresponding dates
```{r}
S20.Groundwater.Detrend.df<-data.frame(as.Date(GW_S20$Date),G_3356_Detrend)
colnames(S20.Groundwater.Detrend.df)<-c("Date","Groundwater")
```

### Declustering

The Decluster() function declusters a time series using a threshold u specified as a quantile of the completed series and  separation criterion SepCrit to ensure independent events. If mu=365.25 then SepCrit denotes the minimum number of days readings must remain below the threshold before a new event is defined.

```{r Multi}
G_3356.Declustered<-Decluster(Data=G_3356_Detrend,u=0.95,SepCrit=3,mu=365.25)
```

Plot showing the completed, detrended record at Well G-3356 (grey circles) along with cluster maxima (red circles) identified using a 95% threshold (green line) and three day separation criterion.

```{r}
G_3356_Imp$Detrend<-G_3356_Detrend
plot(as.Date(G_3356_Imp$Date),G_3356_Imp$Detrend,col="Grey",pch=16,
     cex=0.25,xlab="Date",ylab="Groundwater level (ft NGVD 29)")
abline(h=G_3356.Declustered$Threshold,col="Dark Green")
points(as.Date(G_3356_Imp$Date[G_3356.Declustered$EventsMax]),
       G_3356.Declustered$Declustered[G_3356.Declustered$EventsMax],col="Red",pch=16,cex=0.5)
```

Other outputs from the Decluster() function include the threshold on the original scale  
```{r}
G_3356.Declustered$Threshold
```
and the number of events per year
```{r}
G_3356.Declustered$EventsPerYear
```
In preparation for later work, lets assign the detrended and declustered groundwater series at site S20 a name.
```{r}
S20.Groundwater.Detrend.Declustered<-G_3356.Declustered$Declustered
```


Reading in the other rainfall and O-sWL series at Site S20
```{r}
S20.Rainfall.df<-Perrine_df
S20.Rainfall.df$Date<-as.Date(S20.Rainfall.df$Date)
S20.Rainfall.df$Na.Removed<-S20.Rainfall.df$Value
S20.Rainfall.df$Na.Removed[which(is.na(S20.Rainfall.df$Value)==TRUE)]<-0
S20.OsWL.df<-S20_T_MAX_Daily_Completed_Detrend_Declustered[,c(2,4)]
S20.OsWL.df$Date<-as.Date(S20.OsWL.df$Date)
```
Detrending and declustering the rainfall and O-sWL series at Site S20
```{r}
S20.OsWL.Detrend<-Detrend(Data=S20.OsWL.df,Method = "window",PLOT=FALSE,
                          x_lab="Date",y_lab="O-sWL (ft NGVD 29)")
```
Creating a dataframe with the date alongside the detrended OsWL series
```{r}
S20.OsWL.Detrend.df<-data.frame(as.Date(S20.OsWL.df$Date),S20.OsWL.Detrend)
colnames(S20.OsWL.Detrend.df)<-c("Date","OsWL")
```
Declustering rainfall and O-sWL series at site S20,
```{r, warning=FALSE}
#Setting missing values to zero 
S20.Rainfall.Declustered<-Decluster(Data=S20.Rainfall.df$Na.Removed,u=0.95,SepCrit=3)$Declustered
S20.Rainfall.Declustered[which(is.na(S20.Rainfall.df$Value)==TRUE)]<-NA
#S20.O-sWL does not contain any missing values and so can be declustered without
S20.OsWL.Detrend.Declustered<-Decluster(Data=S20.OsWL.Detrend,u=0.95,SepCrit=3,mu=365.25)$Declustered
```
Creating data frames with the date alongside declustered series
```{r}
S20.OsWL.Detrend.Declustered.df<-data.frame(S20.OsWL.df$Date,S20.OsWL.Detrend.Declustered)
colnames(S20.OsWL.Detrend.Declustered.df)<-c("Date","OsWL")
S20.Rainfall.Declustered.df<-data.frame(S20.Rainfall.df$Date,S20.Rainfall.Declustered)
colnames(S20.Rainfall.Declustered.df)<-c("Date","Rainfall")
S20.Groundwater.Detrend.Declustered.df<-data.frame(G_3356$Date,S20.Groundwater.Detrend.Declustered)
colnames(S20.Groundwater.Detrend.Declustered.df)<-c("Date","OsWL")
```
Use the Dataframe_Combine function to create data frames containing all observations of the original, detrended if necessary, and declustered time series. On dates where not all variables are observed, missing values are assigned NA.  
```{r}
S20.Detrend.df<-Dataframe_Combine(data.1<-S20.Rainfall.df,
                                  data.2<-S20.OsWL.Detrend.df,
                                  data.3<-S20.Groundwater.Detrend.df,
                                  names=c("Rainfall","OsWL","Groundwater"))
S20.Detrend.Declustered.df<-Dataframe_Combine(data.1<-S20.Rainfall.Declustered.df,
                                              data.2<-S20.OsWL.Detrend.Declustered.df,
                                              data.3<-S20.Groundwater.Detrend.Declustered.df,
                                              names=c("Rainfall","OsWL","Groundwater"))
```
#### Fit GPD

The GPD_Fit() function fits a generalized Pareto distribution (GPD) to observations above a threshold u, specified as a quantile of the completed time series.  To fit the distribution the   GPD_Fit() function requires the declustered series as its Data argument and the entire completed series, detrended if necessary, as its Data.Full argument. The completed series is required to calcuate the value on the original scale corresponding to u. If PLOT=="TRUE" then diagnostic plots are produced to allow an assessment of the fit.

```{r}
GPD_Fit(Data=S20.Detrend.Declustered.df$Rainfall,Data_Full=na.omit(S20.Detrend.df$Rainfall),
        u=0.997,PLOT="TRUE",xlab_hist="O-sWL (ft NGVD 29)",y_lab="O-sWL (ft NGVD 29)")
```

## 3. Correlation analysis
We can use the `Kendall_Lag` function to view the Kendall's rank correlations coefficient $\tau$ between the time seres over a range of lags
```{r}
S20.Kendall.Results<-Kendall_Lag(Data=S20.Detrend.df,GAP=0.2)
```
Lets pull out the Kendall correlation coefficient values between rainfall and O-sWL for lags of $-5,...,0,..,5$ applied to the latter quantity 
```{r}
S20.Kendall.Results$Value$Rainfall_OsWL
```
and the corresponding p-values testing the null hypothesis $\tau =0$
```{r}
S20.Kendall.Results$Test$Rainfall_OsWL_Test
```

## 4. Bivariate Analysis
In the report the 2D analysis considers the two forcings currently accounted for in structural design assessments undertaken by SFWMD: rainfall and O-sWL. The 2D analysis commences with the well-established two-sided conditional sampling approach, where excesses of a conditioning variable are paired with co-occurring values of another variable to create two samples. For each sample the marginals (one extreme, one non-extreme) and joint distribution are then modeled. 

The two (conditional) joint distributions are modeled independently of the marginals by using a copula. The Copula_Threshold_2D() function explores the sensitivity of the best fitting copula, in terms of Akaike Information Criterion (AIC), to allow the practitioner to make an informed choice with regards to threshold selection. It undertakes the conditional sampling described above and reports the best fitting bivariate copula. The procedure is carried out for a single or range of thresholds specified by the `Thres` argument and the procedure is automatically repeated with the variables switched.


```{r, fig.height=6, fig.align="center", warning=FALSE}
Copula_Threshold_2D(Data_Detrend=S20.Detrend.df[,-c(1,4)],
                    Data_Declust=S20.Detrend.Declustered.df[,-c(1,4)],
                    y_lim_min=-0.075, y_lim_max =0.25,
                    Upper=c(2,9), Lower=c(2,10),GAP=0.15)
```
The Diag_Non_Con() function is designed to aid in the selection of the appropriate (non-extreme) unbounded marginal distribution for the non-conditioned variable.  
```{r, fig.height=6, fig.align="center", warning=FALSE}
S20.Rainfall<-Con_Sampling_2D(Data_Detrend=S20.Detrend.df[,-c(1,4)],
                              Data_Declust=S20.Detrend.Declustered.df[,-c(1,4)],
                              Con_Variable="Rainfall",Thres=0.97)
Diag_Non_Con(Data=S20.Rainfall$Data$OsWL,x_lab="O-sWL (ft NGVD)",y_lim_min=0,y_lim_max=1.5)
```
The Diag_Non_Con_Sel() function, is similar to the Diag_Non_Con() command, but only plots the probability density function and cumulative distribution function of a (single) selected univariate distribution in order to more clearly demonstrate the goodness of fit of a particular distribution. The options are the Gaussian (`Gaus`) and logistic (`Logis`) distributions. 
```{r, fig.height=6, fig.align="center", warning=FALSE}
Diag_Non_Con_Sel(Data=S20.Rainfall$Data$OsWL,x_lab="O-sWL (ft NGVD)",
                 y_lim_min=0,y_lim_max=1.5,Selected="Logis")
```
A generalized Pareto distribution is fitted to the marginal distribution of the conditioning variable i.e. the declustered excesses identified using Con_Sampling_2D(). 
Use `GPD_Fit` to model the OsWL excesses in the `S22.Rainfall$Data` sample
```{r, fig.height=6, fig.align="center", warning=FALSE}
GPD_Fit(Data=S20.Rainfall$Data$Rainfall,Data_Full=S20.Rainfall$Data$Rainfall,u=0,
        PLOT="TRUE",xlab_hist="Rainfall (Inches)",y_lab="Rainfall (Inches)")
```
The process of selecting a conditional sample and fitting marginal distributions is repeated but instead conditioning on O-sWL. The non-conditional variable in this case is (total daily) rainfall, which has a lower bound at zero, and thus requires a suitably truncated distribution. The `Diag_Non_Con_Trunc` fits a selection of truncated distributions to a vector of data. The `Diag_Non_Con_Sel_Trunc` function is analogous to the 'Diag_Non_Con_Sel' function, avalibale distributions are the  Birnbaum-Saunders (`BS`), exponential (`Exp`), gamma (`Gam`), inverse Gaussian (`InvG`), lognormal (`LogN`), Tweedie (`Twe`) and Weibull (`Weib`). 
```{r, fig.height=6, fig.align="center", warning=FALSE}
S20.OsWL<-Con_Sampling_2D(Data_Detrend=S20.Detrend.df[,-c(1,4)],
                          Data_Declust=S20.Detrend.Declustered.df[,-c(1,4)],
                          Con_Variable="OsWL",Thres=0.98)
Diag_Non_Con_Trunc(Data=S20.OsWL$Data$Rainfall+0.001,x_lab="Rainfall (Inches)",
                   y_lim_min=0,y_lim_max=2)
Diag_Non_Con_Trunc_Sel(Data=S20.OsWL$Data$Rainfall+0.001,x_lab="Rainfall (Inches)",
                       y_lim_min=0,y_lim_max=2,
                       Selected="BS")
```
Using `GPD_Fit` to assess the fit of the GPD to the OsWL excesses
```{r, fig.height=6, fig.align="center", warning=FALSE}
GPD_Fit(Data=S20.OsWL$Data$OsWL,Data_Full=S20.OsWL$Data$OsWL,u=0,
        PLOT="TRUE",xlab_hist="OsWL (ft NGVD)",y_lab="OsWL (ft NGVD)")
```
The `Design_Event_2D()` function finds the isoline associated with a particular return period, by overlaying the two corresponding isolines from the joint distributions fitted to the conditional samples using the method in Bender et al. (2016). `Design_Event_2D()` requires the copulas families chosen to model the dependence structure in the two conditional samples as input.
```{r, fig.height=6, fig.align="center", results="hide", warning=FALSE}
S20.Copula.Rainfall<-Copula_Threshold_2D(Data_Detrend=S20.Detrend.df[,-c(1,4)],
                                         Data_Declust=S20.Detrend.Declustered.df[,-c(1,4)],
                                         Thres =0.97,y_lim_min=-0.075,y_lim_max=0.25,
                                         Upper=c(2,9),Lower=c(2,10),GAP=0.15)$Copula_Family_Var1
S20.Copula.OsWL<-Copula_Threshold_2D(Data_Detrend=S20.Detrend.df[,-c(1,4)],
                                     Data_Declust=S20.Detrend.Declustered.df[,-c(1,4)],
                                     Thres =0.97,y_lim_min=-0.075, y_lim_max =0.25,
                                     Upper=c(2,9),Lower=c(2,10),GAP=0.15)$Copula_Family_Var2
```
As input the function requires 
*  Data = Original (detrended) rainfall and O-sWL series
*  Data_Con1/Data_Con2 = two conditionally sampled datasets, 
*  Thres1/Thres2 = two thresholds associated with the conditionally sampled datasets 
*  Copula_Family1/Copula_Family2 two families of the two fitted copulas
*  Marginal_Dist1/Marginal_Dist2 Selected non-extreme margnal distributions
*  HazScen = Hazard scenario (AND/OR)
*  RP = Return Period of interest
*  N = size of the sample from the fitted joint distributions used to estimate the density along the isoline of interest
*  N_Ensemble = size of the ensemble of events sampled along the isoline of interest 
```{r, fig.height=6, fig.align="center", warning=FALSE}
S20.Bivariate<-Design_Event_2D(Data=S20.Detrend.df[,-c(1,4)], Data_Con1=S20.Rainfall$Data, 
                               Data_Con2=S20.OsWL$Data, Thres1=0.97, Thres2=0.97, 
                               Copula_Family1=S20.Copula.Rainfall, Copula_Family2=S20.Copula.OsWL, 
                               Marginal_Dist1="Logis", Marginal_Dist2="Twe",
                               x_lab="Rainfall (mm)",y_lab="O-sWL (m NGVD 29)",
                               RP=100,N=10,N_Ensemble=10)
```
Design event according to the "Most likely" event approach (diamond in the plot)
```{r}
S20.Bivariate$MostLikelyEvent
```
Design event under the assuption of full dependence (Triangle in the plot)
```{r}
S20.Bivariate$FullDependence
```

## 5. Trivariate analysis
In the report three higher dimensional (>3) approaches are implemented to model the joint distribution of rainfall, O-sWL and groundwater level, they are:
  
* Standard (trivariate) copula
* Pair Copula Construction
* Heffernan and Tawn (2004)

In the package, each approach has a `_Fit` and `_Sim` function. The latter requires a `MIGPD` object as its `Marginals` input argument, in order for the simulations on $[0,1]3$ to be transformed back to the original scale. The `Migpd_Fit` command fits independent GPDs to the data in each row of a dataframe (excluding the first column if it is a "Date" object) creating a `MIGPD` object. 

```{r}
S20.Migpd<-Migpd_Fit(Data=S20.Detrend.Declustered.df[,-1],mqu=c(0.975,0.975,0.9676))
summary(S20.Migpd)
```
Standard (trivariate) copula are the most conceptually simple of the copula based models, using a single parametric multivariate probability distribution as the copula. The Standard_Copula_Fit() function fits elliptic (specified by `Gaussian` or `tcop`) or Archimedean (specified by `Gumbel`,`Clayton` or `Frank`) copula to a trivariate dataset. Lets first fit a Gaussian copula
```{r}
S20.Gaussian<-Standard_Copula_Fit(Data=S20.Detrend.df,Copula_Type="Gaussian")
```
From which the Standard_Copula_Sim() function can be used to simulate a synthetic record of N years
```{r}
S20.Gaussian.Sim<-Standard_Copula_Sim(Data=S20.Detrend.df,Marginals=S20.Migpd,
                                      Copula=S20.Gaussian,N=100)
```
Plotting the observationed and simulated values
```{r, fig.height=6, fig.align="center"}
S20.Pairs.Plot.Data<-data.frame(rbind(na.omit(S20.Detrend.df[,-1]),S20.Gaussian.Sim$x.Sim),
                                c(rep("Observation",nrow(na.omit(S20.Detrend.df))),
                                  rep("Simulation",nrow(S20.Gaussian.Sim$x.Sim))))
colnames(S20.Pairs.Plot.Data)<-c(names(S20.Detrend.df)[-1],"Type")
pairs(S20.Pairs.Plot.Data[,1:3],
      col=ifelse(S20.Pairs.Plot.Data$Type=="Observation","Black",alpha("Red",0.3)),
      upper.panel=NULL,pch=16)
```
The Standard_Copula_Sel() function can be used to deduce the best fitting in terms of AIC 
```{r}
Standard_Copula_Sel(Data=S20.Detrend.df)
```
Standard trivariate copulas lack flexibility to model joint distributions where heterogeneous dependencies exist between the variable pairs. Pair copula constructions construct multivariate distribution using a cascade of bivariate copulas (some of which are conditional). As the dimensionality of the problem increases the number of mathematically equally valid decompositions quickly becomes large.  Bedford and Cooke (2001,2002) introduced the regular vine, a graphical model which helps to organize the possible decompositions. The Canonical (C-) and D- vine are two commonly utilized sub-categories of regular vines, in the trivariate case a vine copula is simultaneously a C- and D-vine. Lets fit a regular vine copula model
```{r}
S20.Vine<-Vine_Copula_Fit(Data=S20.Detrend.df)
```
From which the Vine_Copula_Sim() function can be used to simulate a synthetic record of N years
```{r}
S20.Vine.Sim<-Vine_Copula_Sim(Data=S20.Detrend.df,Marginals=S20.Migpd,                              Vine_family=S20.Vine$Family,Vine_par=S20.Vine$Par,Vine_par2=S20.Vine$Par2,N=100)
```
Plotting the observationed and simulated values
```{r, fig.height=5, fig.align="center"}
S20.Pairs.Plot.Data<-data.frame(rbind(na.omit(S20.Detrend.df[,-1]),S20.Vine.Sim$x.Sim),
                                c(rep("Observation",nrow(na.omit(S20.Detrend.df))),
                                  rep("Simulation",nrow(S20.Vine.Sim$x.Sim))))
colnames(S20.Pairs.Plot.Data)<-c(names(S20.Detrend.df)[-1],"Type")
pairs(S20.Pairs.Plot.Data[,1:3],
      col=ifelse(S20.Pairs.Plot.Data$Type=="Observation","Black",alpha("Red",0.3)),
      upper.panel=NULL,pch=16)
```
Finally, lets implement the Heffernan and Tawn (2004) approach, where a non-linear regression model is fitted to the (joint) observations where a (conditioning) variable is above a specified threshold. The regression model typically adopted is  
**Y**_<sub>-i</sub>=**a**Y_<sub>i</sub>+Y<sub>i</sub><sup>**b**</sup>**Z** for Y<sub>i</sub>>v <br />
where **Y** is a set of variables transdformed to a common scale, **Y**<sub>-i</sub> is the set of variables excluding Y<sub>-i</sub>, **a** and **b** are vectors of regression parameters and **Z** is a vector of residuals. The dependence structure, when a specified variable is extreme is thus captured by the regression parameters and the joint residuals. The procedure is repeated conditioning on each variable in turn to build up of the joint distribution when at least one variable is in an extreme state. The `HT04` command fits and simulates N years worth of simulations from the model.
```{r}
S20.HT04<-HT04(data_Detrend_Dependence_df=S20.Detrend.df,
               data_Detrend_Declustered_df=S20.Detrend.Declustered.df,
               u_Dependence=0.995,Migpd=S20.Migpd,mu=365.25,N=1000)
```
Output of the function includes the three conditional `Models`, proportion of occasions where each variable is most extreme given at least one variable is extreme `prop`as well as, the simulations on the transfromed scale `u.Sim` (gumbel by default) and original scale `x.Sim`. Lets view the fitted model when conditioning on rainfall
```{r}
S20.HT04$Model$Rainfall
S20.HT04$Prop
```
and the which the proporiton of the occasions in the original sample that rainfall is the most extreme of the drivers given that at least one driver is extreme. 

The HT04 approach uses rejection sampling to generate synthetic records. The first step involves sampling a variable, conditioned to exceed the `u_Dependence` threshold. A joint residual associated with the corresponding regression is independently sampled and other variables estimated using the fitted regression parameters. If the variable conditioned to be extreme in step one is not the most extreme the sample is rejected. The process is repeated until the relative proportion of simulated events where each variable is a maximum, conditional on being above the threshold, is consistent with the empirical distribution. Labelling the simulations `S20.HT04.Sim` 
```{r}
S20.HT04.Sim<-S20.HT04$x.sim
```
and now plotting the simulations from the HT04 model
```{r, fig.height=5, fig.align="center"}
S20.Pairs.Plot.Data<-data.frame(rbind(na.omit(S20.Detrend.df[,-1]),S20.HT04.Sim),
                                c(rep("Observation",nrow(na.omit(S20.Detrend.df))),
                                  rep("Simulation",nrow(S20.HT04.Sim))))
colnames(S20.Pairs.Plot.Data)<-c(names(S20.Detrend.df)[-1],"Type")
pairs(S20.Pairs.Plot.Data[,1:3],
      col=ifelse(S20.Pairs.Plot.Data$Type=="Observation","Black",alpha("Red",0.2)),
      upper.panel=NULL,pch=16)
```





