
# MultiHazard <img src="https://user-images.githubusercontent.com/15319503/195926656-9d3d37b5-86ab-4d4b-9e6d-3c70d5399c73.png" align="right" height="200"/>

The `MultiHazard` package provides tools for stationary multivariate
statistical modeling, for example, to estimate joint occurrence
probabilities of MULTIple co-occurring HAZARDs. The package contains
functions for pre-processing data, including imputing missing values,
detrending and declustering time series as well as analyzing pairwise
correlations over a range of lags. Functionality is built in to
implement the conditional sampling - copula theory approach described in
[Jane et al. (2020)](https://doi.org/10.5194/nhess-20-2681-2020)
including the automated threshold selection method from [Solari et
al. (2017)](https://doi.org/10.1002/2016WR019426). There is a function
that calculates joint probability contours using the method of
overlaying conditional contours given in [Bender et
al. (2016)](https://doi.org/10.1080/02626667.2015.1052816) and extracts
design events such as the “most likely” event or an ensemble of possible
design events. The package also includes methods from [Murphy-Barltrop
et al. (2023)](https://doi.org/10.1002/env.2797) and [Murphy-Barltrop et
al. (2024)](https://doi.org/10.1007/s10687-024-00490-4) for deriving
isolines using the [Heffernan and Tawn
(2004)](https://doi.org/10.1111/j.1467-9868.2004.02050.x) [HT04](#ht04)
and [Wadsworth and Tawn (2013)](https://doi.org/10.3150/12-BEJ471)
\[WT13\] models, together with a novel bootstrap procedure for
quantifying sampling uncertainty in the isolines. Three higher
dimensional approaches — standard (elliptic/Archimedean) copulas, Pair
Copula Constructions (PCCs) and a conditional threshold exceedance
approach (HT04) — are coded. Finally, the package can be implemented to
derive temporally coherent extreme events comprising a hyetograph and
water level curve for simulated peak rainfall and peak sea level events,
as outlined in (Report).  
\### Citation:

> Jane, R., Wahl, T., Peña, F., Obeysekera, J., Murphy-Barltrop, C.,
> Ali, J., Maduwantha, P., Li, H., and Malagón Santos, V. (under review)
> MultiHazard: Copula-based Joint Probability Analysis in R. Journal of
> Open Source Software. \[under revision\]

<br>

### Applications of package:

> Li, H., Jane, R. A., Eilander, D., Enríquez, A. R., Haer, T., and
> Ward, P. J. (2025). Assessing the spatial correlation of potential
> compound flooding in the United States, EGUsphere \[preprint\],
> <https://doi.org/10.5194/egusphere-2025-2993>.

> Amorim, R., Villarini, G., Kim, H., Jane, R., and Wahl, T. (2025). A
> Practitioner’s approach to process-driven modeling of compound
> rainfall and storm surge extremes for coastal Texas, J. Hydrol. Eng.,
> 30(5), 04025025, <https://doi.org/10.1061/JHYEFF.HEENG-648>.

> Maduwantha, P., Wahl, T., Santamaria-Aguilar, S., Jane, R., Booth, J.
> F., Kim, H., and Villarini, G. (2024). A multivariate statistical
> framework for mixed storm types in compound flood analysis, Nat.
> Hazards Earth Syst. Sci., 24, 4091–4107,
> <https://doi.org/10.5194/nhess-24-4091-2024>.

> Nasr, A. A., Wahl, T., Rashid, M. M., Jane, R., Camus, P. and Haigh,
> I. D. (2023). Temporal changes in dependence between compound coastal
> and inland flooding drivers around the contiguous United States
> coastline, Weather Clim. Extrem., 41, 100594,
> <https://doi.org/10.1016/j.wace.2023.100594>.

> Kim, H., Villarini, G., Jane, R., Wahl, T., Misra, S., and Michalek,
> A. (2023). On the generation of high‐resolution probabilistic design
> events capturing the joint occurrence of rainfall and storm surge in
> coastal basins, Int. J. Climatol, 43(2), 761-771,
> <https://doi.org/10.1002/joc.7825>.

> Kim, T., Villarini, G., Kim, H., Jane, R., and Wahl, T. (2023). On the
> compounding of nitrate loads and discharge, J. Environ. Qual., 52,
> 706–717. <https://doi.org/10.1002/jeq2.20458>.

> Peña, F., Obeysekera, J., Jane R., Nardi, F., Maran, C., Cadogan, A.,
> de Groen, F., and Melesse, A. (2023). Investigating compound flooding
> in a low elevation coastal karst environment using multivariate
> statistical and 2D hydrodynamic modeling, Weather Clim. Extrem., 39,
> 100534. <https://doi.org/10.1016/j.wace.2022.100534>.

> Jane, R., Cadavid, L., Obeysekera, J., and Wahl, T. (2020).
> Multivariate statistical modelling of the drivers of compound flood
> events in South Florida, Nat. Hazards Earth Syst. Sci., 20, 2681–2699,
> <https://doi.org/10.5194/nhess-20-2681-2020>.

<br>

## Community guidelines

Contributions to the `MultiHazard` package are welcome! Please feel free
to submit issues or pull requests on GitHub.

## Installation

Install the latest version of this package by entering the following in
R:

``` r
install.packages("remotes")
remotes::install_github("rjaneUCF/MultiHazard")
```

## 1. Introduction

The `MultiHazard` package provides tools for stationary multivariate
statistical modeling, for example, to estimate the joint distribution of
MULTIple co-occurring HAZARDs. This document is designed to explain and
demonstrate the functions contained within the package. Section 1 looks
at the functions concerned with pre-processing the data including
imputing missing values. Section 2 illustrates the functions for
detrending and declustering time series while Section 3 introduces a
function that analyzes pairwise correlations over a range of lags.
Section 4 shows how the conditional sampling - copula theory approach in
[Jane et al. (2020)](https://doi.org/10.5194/nhess-20-2681-2020) can be
implemented including the automated threshold selection method in
[Solari et al. (2017)](https://doi.org/10.1002/2016WR019426). Functions
for selecting the best fitting among an array of (non-extreme, truncated
and non-truncated) parametric marginal distributions, and copulas to
model the dependence structure are demonstrated in this section. Section
4 also contains an explanation of the function that derives the joint
probability contours according to the method of overlaying (conditional)
contours given in [Bender et
al. (2016)](https://doi.org/10.1080/02626667.2015.1052816), and extracts
design events such as the “most likely” event or an ensemble of possible
design events. Section 4 also introduces the functions that generate
isolines using the methods from [Murphy-Barltrop et
al. (2023)](https://doi.org/10.1002/env.2797) and [Murphy-Barltrop et
al. (2024)](https://doi.org/10.1007/s10687-024-00490-4), and implements
a novel bootstrap procedure for quantifying sampling uncertainty in the
isolines. Section 5 introduces the functions for fitting and simulating
synthetic events from three higher-dimensional approaches - standard
(elliptic/Archimedean) copulas, Pair Copula Constructions (PCCs) and the
conditional threshold exceedance approach of HT04. Section 6 describes a
function that calculates the time for a user-specified height of sea
level rise to occur under various scenarios. Lastly, Section 7 shows the
simulation of temporally coherent extreme rainfall and ocean water level
events.

## 2. Pre-processing

### Imputation

Well G_3356 represents the groundwater level at Site S20, however, it
contains missing values. Let’s impute missing values in the record at
Well G_3356 using recordings at nearby Well G_3355. First, reading in
the two time series.

``` r
#Viewing first few rows of in the groundwater level records
head(G_3356)
```

    ##         Date Value
    ## 1 1985-10-23  2.46
    ## 2 1985-10-24  2.47
    ## 3 1985-10-25  2.41
    ## 4 1985-10-26  2.37
    ## 5 1985-10-27  2.63
    ## 6 1985-10-28  2.54

``` r
head(G_3355)
```

    ##         Date Value
    ## 1 1985-08-20  2.53
    ## 2 1985-08-21  2.50
    ## 3 1985-08-22  2.46
    ## 4 1985-08-23  2.43
    ## 5 1985-08-24  2.40
    ## 6 1985-08-25  2.37

``` r
#Converting Date column to "Date"" object
G_3356$Date<-seq(as.Date("1985-10-23"), as.Date("2019-05-29"), by="day")
G_3355$Date<-seq(as.Date("1985-08-20"), as.Date("2019-06-02"), by="day")
#Converting column containing the readings to a "numeric"" object
G_3356$Value<-as.numeric(as.character(G_3356$Value))
G_3355$Value<-as.numeric(as.character(G_3355$Value))
```

Warning message confirms there are NAs in the record at Well G_3356.
Before carrying out the imputation the two data frames need to be
merged.

``` r
#Merge the two dataframes by date
GW_S20<-left_join(G_3356,G_3355,by="Date")
colnames(GW_S20)<-c("Date","G3356","G3355")
#Carrying out imputation
Imp<-Imputation(Data=GW_S20,Variable="G3356",
                x_lab="G-3355 (ft NGVD 29)", y_lab="G-3356 (ft NGVD 29)")
```

![](README_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

The completed record is given in the `ValuesFilled` column of the data
frame outputted as the `Data` object while the linear regression model
including its coefficient of determinant are given by the `model` output
argument.

``` r
head(Imp$Data)
```

    ##         Date G3356 G3355 ValuesFilled
    ## 1 1985-10-23  2.46  2.87         2.46
    ## 2 1985-10-24  2.47  2.85         2.47
    ## 3 1985-10-25  2.41  2.82         2.41
    ## 4 1985-10-26  2.37  2.79         2.37
    ## 5 1985-10-27  2.63  2.96         2.63
    ## 6 1985-10-28  2.54  2.96         2.54

``` r
Imp$Model
```

    ## 
    ## Call:
    ## lm(formula = data[, variable] ~ data[, Other.variable])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -1.60190 -0.16654  0.00525  0.16858  1.73824 
    ## 
    ## Coefficients:
    ##                        Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)            0.275858   0.012366   22.31   <2e-16 ***
    ## data[, Other.variable] 0.700724   0.004459  157.15   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2995 on 11825 degrees of freedom
    ##   (445 observations deleted due to missingness)
    ## Multiple R-squared:  0.6762, Adjusted R-squared:  0.6762 
    ## F-statistic: 2.47e+04 on 1 and 11825 DF,  p-value: < 2.2e-16

Are any values still NA?

``` r
G_3356_ValueFilled_NA<-which(is.na(Imp$Data$ValuesFilled)==TRUE)
length(G_3356_ValueFilled_NA) 
```

    ## [1] 3

Linear interpolating the three remaining NAs.

``` r
G3356_approx<-approx(seq(1,length(Imp$Data$ValuesFilled),1),Imp$Data$ValuesFilled,
                     xout=seq(1,length(Imp$Data$ValuesFilled),1))
Imp$Data$ValuesFilled[which(is.na(Imp$Data$ValuesFilled)==TRUE)]<-
  G3356_approx$y[which(is.na(Imp$Data$ValuesFilled)==TRUE)]
```

### Detrending

In the analysis completed O-sWL (Ocean-side Water Level) and groundwater
level series are subsequently detrended. The Detrend() function uses
either a linear fit covering the entire data (`Method`=`linear`) or
moving average window (`Method`=`window`) of a specified length
(`Window_Width`) to remove trends from a time series. The residuals are
added to the final `End_Length` observations. The default `Detrend()`
parameters specify a moving average (`Method`=`window`) three month
window (`Window_Width`=`89`), to remove any seasonality from the time
series. The default is then to add the residuals to the average of the
final five years of observations (`End_Length`=`1826`) to bring the
record to the present day level, accounting for the Perigean tide in the
case of O-sWL. The mean of the observations over the first three months
were subtracted from the values during this period before the present
day (5-year) average was added. The following R code detrends the record
at Well G_3356. Note the function requires a Date object and the
completed series.

``` r
#Creating a data frame with the imputed series alongside the corresponding dates 
G_3356_Imp<-data.frame(Imp$Data$Date,Imp$Data$ValuesFilled)
colnames(G_3356_Imp)<-c("Date","ValuesFilled")
#Detrending
G_3356_Detrend<-Detrend(Data=G_3356_Imp,PLOT=TRUE,x_lab="Date",
                        y_lab="Groundwater level (ft NGVD 29)")
```

![](README_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

Output of the function is simply the detrended time series.

``` r
head(G_3356_Detrend)
```

    ## [1] 2.394700 2.411588 2.360033 2.327588 2.595588 2.520255

Creating a data frame containing the detrended groundwater series at
site S_20 i.e. G_3356_Detrend and their corresponding dates

``` r
S20.Groundwater.Detrend.df<-data.frame(as.Date(GW_S20$Date),G_3356_Detrend)
colnames(S20.Groundwater.Detrend.df)<-c("Date","Groundwater")
```

### Declustering

The `Decluster()` function declusters a time series using a threshold u
specified as a quantile of the completed series and separation criterion
`SepCrit` to ensure independent events. If `mu`=`365.25` then `SepCrit`
denotes the minimum number of days readings must remain below the
threshold before a new event is defined.

``` r
G_3356.Declustered<-Decluster(Data=G_3356_Detrend,u=0.95,SepCrit=3,mu=365.25)
```

Plot showing the completed, detrended record at Well G-3356 (grey
circles) along with cluster maxima (red circles) identified using a 95%
threshold (green line) and three day separation criterion.

``` r
G_3356_Imp$Detrend<-G_3356_Detrend
plot(as.Date(G_3356_Imp$Date),G_3356_Imp$Detrend,col="Grey",pch=16,
     cex=0.25,xlab="Date",ylab="Groundwater level (ft NGVD 29)")
abline(h=G_3356.Declustered$Threshold,col="Dark Green")
points(as.Date(G_3356_Imp$Date[G_3356.Declustered$EventsMax]),
       G_3356.Declustered$Declustered[G_3356.Declustered$EventsMax],
       col="Red",pch=16,cex=0.5)
```

![](README_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

Other outputs from the `Decluster()` function include the threshold on
the original scale

``` r
G_3356.Declustered$Threshold
```

    ## [1] 2.751744

and the number of events per year

``` r
G_3356.Declustered$Rate
```

    ## [1] 7.857399

In preparation for later work, lets assign the detrended and declustered
groundwater series at site S20 a name.

``` r
S20.Groundwater.Detrend.Declustered<-G_3356.Declustered$Declustered
```

Reading in the other rainfall and O-sWL series at Site S20

``` r
#Changing names of the data frames
S20.Rainfall.df<-Perrine_df
S20.OsWL.df<-S20_T_MAX_Daily_Completed_Detrend_Declustered[,c(2,4)]
#Converting Date column to "Date"" object
S20.Rainfall.df$Date<-as.Date(S20.Rainfall.df$Date)
S20.OsWL.df$Date<-as.Date(S20.OsWL.df$Date)
```

Detrending and declustering the rainfall and O-sWL series at Site S20

``` r
S20.OsWL.Detrend<-Detrend(Data=S20.OsWL.df,Method = "window",PLOT=FALSE,
                          x_lab="Date",y_lab="O-sWL (ft NGVD 29)")
```

Creating a dataframe with the date alongside the detrended OsWL series

``` r
S20.OsWL.Detrend.df<-data.frame(as.Date(S20.OsWL.df$Date),S20.OsWL.Detrend)
colnames(S20.OsWL.Detrend.df)<-c("Date","OsWL")
```

Declustering rainfall and O-sWL series at site S20,

``` r
#Declustering rainfall and O-sWL series
S20.Rainfall.Declustered<-Decluster(Data=S20.Rainfall.df$Value,u=0.95,SepCrit=3)$Declustered
S20.OsWL.Detrend.Declustered<-Decluster(Data=S20.OsWL.Detrend,u=0.95,SepCrit=3,mu=365.25)$Declustered
```

Creating data frames with the date alongside declustered series

``` r
S20.OsWL.Detrend.Declustered.df<-data.frame(S20.OsWL.df$Date,S20.OsWL.Detrend.Declustered)
colnames(S20.OsWL.Detrend.Declustered.df)<-c("Date","OsWL")
S20.Rainfall.Declustered.df<-data.frame(S20.Rainfall.df$Date,S20.Rainfall.Declustered)
colnames(S20.Rainfall.Declustered.df)<-c("Date","Rainfall")
S20.Groundwater.Detrend.Declustered.df<-data.frame(G_3356$Date,S20.Groundwater.Detrend.Declustered)
colnames(S20.Groundwater.Detrend.Declustered.df)<-c("Date","Groundwater")
```

Use the `Dataframe_Combine()` function to create data frames containing
all observations of the original, detrended if necessary, and
declustered time series. On dates where not all variables are observed,
missing values are assigned NA.

``` r
S20.Detrend.df<-Dataframe_Combine(data.1<-S20.Rainfall.df,
                                  data.2<-S20.OsWL.Detrend.df,
                                  data.3<-S20.Groundwater.Detrend.df,
                                  names=c("Rainfall","OsWL","Groundwater"))
S20.Detrend.Declustered.df<-Dataframe_Combine(data.1<-S20.Rainfall.Declustered.df,
                                              data.2<-S20.OsWL.Detrend.Declustered.df,
                                              data.3<-S20.Groundwater.Detrend.Declustered.df,
                                              names=c("Rainfall","OsWL","Groundwater"))
```

The package contains two other declustering functions. The
`Decluster_SW()` function declusters a time series via a storm window
approach. A moving window of length (`Window_Width`) is moved over the
time series, if the maximum value is located at the center of the window
then the value is considered a peak and retained, otherwise it is set
equal to NA. For a seven day window at S20:

``` r
S20.Rainfall.Declustered.SW<-Decluster_SW(Data=S20.Rainfall.df,Window_Width=7)
```

Plotting the original and detrended series:

``` r
plot(S20.Rainfall.df$Date,S20.Rainfall.df$Value,pch=16,cex=0.5,
     xlab="Date",ylab="Total daily rainfall (Inches)")
points(S20.Rainfall.df$Date,S20.Rainfall.Declustered.SW$Declustered,pch=16,col=2,cex=0.5)
```

![](README_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

Repeating the analysis for the O-sWL with a 3-day window.

``` r
S20.OsWL.Declustered.SW<-Decluster_SW(Data=S20.OsWL.df,Window_Width=3)
```

The `Decluster_S_SW()` function declusters a summed time series via a
storm window approach. First a moving window of width
(`Window_Width_Sum`) travels across the data and each time the values
are summed. As with the `Decluster_SW()` function a moving window of
length `(Window_Width`) is then moved over the time series, if the
maximum value in a window is located at its center then the value
considered a peak and retained, otherwise it is set equal to NA. To
decluster weekly precipitation totals using a seven day storm window at
S20:

``` r
#Declustering
S20.Rainfall.Declustered.S.SW<-Decluster_S_SW(Data=S20.Rainfall.df, 
                                              Window_Width_Sum=7, Window_Width=7)
#First twenty values of the weekly totals
S20.Rainfall.Declustered.S.SW$Totals[1:20]
```

    ##  [1]   NA   NA   NA 0.23 0.19 0.10 1.56 1.94 2.04 2.04 2.04 2.91 3.02 2.75 3.15
    ## [16] 3.05 3.05 3.05 2.18 2.03

``` r
#First ten values of the declustered weekly totals
S20.Rainfall.Declustered.S.SW$Declustered[1:20]
```

    ##  [1]   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA 3.15
    ## [16]   NA   NA   NA   NA   NA

Plotting the original and detrended series:

``` r
plot(S20.Rainfall.df$Date,S20.Rainfall.Declustered.S.SW$Totals,pch=16,cex=0.5,
     xlab="Date",ylab="Total weekly rainfall (Inches)")
points(S20.Rainfall.df$Date,S20.Rainfall.Declustered.S.SW$Declustered,pch=16,col=2,cex=0.5)
```

![](README_files/figure-gfm/unnamed-chunk-29-1.png)<!-- -->

#### Fit GPD

The `GPD_Fit()` function fits a generalized Pareto distribution (GPD) to
observations above a threshold u, specified as a quantile of the
completed time series. To fit the distribution the `GPD_Fit()` function
requires the declustered series as its Data argument and the entire
completed series, detrended if necessary, as its `Data_Full` argument.
The completed series is required to calculate the value on the original
scale corresponding to u. If `PLOT=TRUE` then diagnostic plots are
produced to allow an assessment of the fit.

``` r
GPD_Fit(Data=S20.Detrend.Declustered.df$Rainfall,Data_Full=na.omit(S20.Detrend.df$Rainfall),
        u=0.997,PLOT=TRUE,xlab_hist="Rainfall (Inches)",y_lab="Rainfall (Inches)")
```

![](README_files/figure-gfm/unnamed-chunk-30-1.png)<!-- -->

    ## $Threshold
    ## [1] 3.5465
    ## 
    ## $Rate
    ## [1] 1.021938
    ## 
    ## $sigma
    ## [1] 1.421037
    ## 
    ## $xi
    ## [1] 0.1696554
    ## 
    ## $sigma.SE
    ## [1] 0.2365755
    ## 
    ## $xi.SE
    ## [1] 0.1781031

#### Solari (2017) automated threshold selection

[Solari et al. (2017)](https://doi.org/10.1002/2016WR019426) proposes a
methodology for automatic threshold estimation, based on an
EDF-statistic and a goodness of fit test to test the null hypothesis
that exceedances of a high threshold come from a GPD distribution.

EDF-statistics measure the distance between the empirical distribution
$F_n$ obtained from the sample and the parametric distribution $F(x)$.
The Anderson Darling $A^2$ statistic is an EDF-statistic, which assigns
more weight to the tails of the data than similar measures. [Sinclair et
al. (1990)](https://doi.org/10.1080/03610929008830405) proposed the
right-tail weighted Anderson Darling statistic $A_R^2$ which allocates
more weight to the upper tail and less to the lower tail of the
distribution than $A^2$ and is given by:
$${A}_{R}^{2}= -\frac{n}{2} -  \sum\_{i=1}^{n} \left\[\left(2-\frac{(2i-1)}{n}\right)log(1-z\_{i})+2z\_{i}\right\] $$
where $z=F(x)$ and $n$ is the sample size. The approach in [Solari et
al. (2017)](https://doi.org/10.1002/2016WR019426) is implemented as
follows:

1.  A time series is declustered using the storm window approach to
    identify independent peaks.
2.  Candidate thresholds are defined by ordering the peaks and removing
    any repeated values. A GPD is fit to all the peaks above each
    candidate threshold. The right-tail weighted Anderson-Darling
    statistic $A_R^2$ and its corresponding p-value are calculated for
    each threshold.
3.  The threshold that minimizes one minus the p-value is then selected.

The `GPD_Threshold_Solari()` function carries out these steps.

``` r
S20.Rainfall.Solari<-GPD_Threshold_Solari(Event=S20.Rainfall.Declustered.SW$Declustered,
                                          Data=S20.Detrend.df$Rainfall)
```

![](README_files/figure-gfm/unnamed-chunk-31-1.png)<!-- -->

The optimum threshold according to the Solari approach is

``` r
S20.Rainfall.Solari$Candidate_Thres
```

    ## [1] 2.6

``` r
Rainfall.Thres.Quantile<-ecdf(S20.Detrend.df$Rainfall)(S20.Rainfall.Solari$Candidate_Thres)
```

The `GPD_Threshold_Solari_Sel()` allows the goodness-of-fit at a
particular threshold (`Thres`) to be investigated in more detail. Let’s
study the fit of the threshold selected by the method.

``` r
Solari.Sel<-GPD_Threshold_Solari_Sel(Event=S20.Rainfall.Declustered.SW$Declustered,
                                    Data=S20.Detrend.df$Rainfall,
                                    Solari_Output=S20.Rainfall.Solari,
                                    Thres=S20.Rainfall.Solari$Candidate_Thres,
                                    RP_Max=100)
```

![](README_files/figure-gfm/unnamed-chunk-34-1.png)<!-- -->

Repeating the automated threshold selection procedure for O-sWL.

``` r
S20.OsWL.Solari<-GPD_Threshold_Solari(Event=S20.OsWL.Declustered.SW$Declustered,
                                      Data=S20.Detrend.df$OsWL)
```

    ## Error in solve.default(family$info(o)) : 
    ##   system is computationally singular: reciprocal condition number = 5.53436e-17

    ## Error in diag(o$cov) : invalid 'nrow' value (too large or NA)

    ## Error in solve.default(family$info(o)) : 
    ##   system is computationally singular: reciprocal condition number = 4.03099e-17

    ## Error in diag(o$cov) : invalid 'nrow' value (too large or NA)

    ## Error in solve.default(family$info(o)) : 
    ##   system is computationally singular: reciprocal condition number = 8.36628e-17

    ## Error in diag(o$cov) : invalid 'nrow' value (too large or NA)

    ## Error in solve.default(family$info(o)) : 
    ##   system is computationally singular: reciprocal condition number = 7.23864e-17

    ## Error in diag(o$cov) : invalid 'nrow' value (too large or NA)

    ## Error in solve.default(family$info(o)) : 
    ##   Lapack routine dgesv: system is exactly singular: U[2,2] = 0

    ## Error in diag(o$cov) : invalid 'nrow' value (too large or NA)

    ## Error in solve.default(family$info(o)) : 
    ##   system is computationally singular: reciprocal condition number = 8.26056e-17

    ## Error in diag(o$cov) : invalid 'nrow' value (too large or NA)

![](README_files/figure-gfm/unnamed-chunk-35-1.png)<!-- -->

``` r
S20.OsWL.Solari$Candidate_Thres
```

    ## [1] 2.032

``` r
OsWL.Thres.Quantile<-ecdf(S20.Detrend.df$OsWL)(S20.OsWL.Solari$Candidate_Thres)
```

and checking the fit of the GPD at the selected threshold.

``` r
Solari.Sel<-GPD_Threshold_Solari_Sel(Event=S20.OsWL.Declustered.SW$Declustered,
                                     Data=S20.Detrend.df$OsWL,
                                     Solari_Output=S20.OsWL.Solari,
                                     Thres=S20.OsWL.Solari$Candidate_Thres,
                                     RP_Max=100)
```

![](README_files/figure-gfm/unnamed-chunk-38-1.png)<!-- -->

## 3. Correlation analysis

We can use the `Kendall_Lag()` function to view the Kendall’s rank
correlations coefficient $\tau$ between the time series over a range of
lags

``` r
S20.Kendall.Results<-Kendall_Lag(Data=S20.Detrend.df,GAP=0.2)
```

![](README_files/figure-gfm/unnamed-chunk-39-1.png)<!-- -->

Let’s pull out the Kendall correlation coefficient values between
rainfall and O-sWL for lags of $−5, ..., 0, ..,5$ applied to the latter
quantity

``` r
S20.Kendall.Results$Value$Rainfall_OsWL
```

    ##  [1]  0.046483308  0.051860955  0.051392298  0.051311970  0.054097316
    ##  [6]  0.058316831  0.061388245  0.035305812  0.004206059 -0.014356749
    ## [11] -0.025993095 -0.030431776 -0.029481162

and the corresponding p-values testing the null hypothesis $\tau=0$

``` r
S20.Kendall.Results$Test$Rainfall_OsWL_Test
```

    ##  [1] 5.819748e-13 1.014698e-15 8.196547e-16 1.030482e-15 5.160274e-17
    ##  [6] 1.887733e-19 1.987221e-20 2.232548e-07 4.033739e-01 7.669577e-02
    ## [11] 1.186248e-03 6.352872e-05 3.864403e-05

## 4. Bivariate Analysis

#### Two-sided conditional sampling - copula theory method

In the report the 2D analysis considers the two forcings currently
accounted for in structural design assessments undertaken by SFWMD:
rainfall and O-sWL. The 2D analysis commences with the well-established
two-sided conditional sampling approach, where excesses of a
conditioning variable are paired with co-occurring values of another
variable to create two samples. For each sample the marginals (one
extreme, one non-extreme) and joint distribution are then modeled.

The two (conditional) joint distributions are modeled independently of
the marginals by using a copula. The `Copula_Threshold_2D()` function
explores the sensitivity of the best fitting copula, in terms of Akaike
Information Criterion (AIC), to allow the practitioner to make an
informed choice with regards to threshold selection. It undertakes the
conditional sampling described above and reports the best fitting
bivariate copula. The procedure is carried out for a single or range of
thresholds specified by the `u` argument and the procedure is
automatically repeated with the variables switched.

``` r
Copula_Threshold_2D(Data_Detrend=S20.Detrend.df[,-c(1,4)],
Data_Declust=S20.Detrend.Declustered.df[,-c(1,4)],
y_lim_min=-0.075, y_lim_max =0.25,
Upper=c(2,9), Lower=c(2,10),GAP=0.15)
```

![](README_files/figure-gfm/unnamed-chunk-42-1.png)<!-- -->

    ## $Kendalls_Tau1
    ##  [1] 0.05627631 0.05803451 0.05900376 0.08072261 0.10731477 0.14151449
    ##  [7] 0.14427232 0.14762199 0.13101587 0.05056147
    ## 
    ## $p_value_Var1
    ##  [1] 7.698074e-03 9.251524e-03 1.424015e-02 1.974145e-03 1.559739e-04
    ##  [6] 1.171141e-05 4.256557e-05 2.055448e-04 6.000832e-03 4.413023e-01
    ## 
    ## $N_Var1
    ##  [1] 1008  904  776  661  559  432  363  286  200  107
    ## 
    ## $Copula_Family_Var1
    ##  [1]  6  6  6  6  6  6  6 13  6  6
    ## 
    ## $Kendalls_Tau2
    ##  [1] 0.1113049 0.1359921 0.1377104 0.1561184 0.1579352 0.1359861 0.1183870
    ##  [8] 0.1463056 0.1482198 0.1904729
    ## 
    ## $p_value_Var2
    ##  [1] 2.720713e-05 2.549673e-06 1.199322e-05 1.143046e-05 1.879301e-04
    ##  [6] 1.335664e-02 5.199742e-02 3.145804e-02 6.891926e-02 8.220838e-02
    ## 
    ## $N_Var2
    ##  [1] 760 639 535 416 290 168 136 110  77  43
    ## 
    ## $Copula_Family_Var2
    ##  [1]   6   6   6   6 204 204 204 204   1 204

The `Diag_Non_Con()` function is designed to aid in the selection of the
appropriate (non-extreme) unbounded marginal distribution for the
non-conditioned variable.

``` r
S20.Rainfall<-Con_Sampling_2D(Data_Detrend=S20.Detrend.df[,-c(1,4)],
Data_Declust=S20.Detrend.Declustered.df[,-c(1,4)],
Con_Variable="Rainfall",u = Rainfall.Thres.Quantile)
Diag_Non_Con(Data=S20.Rainfall$Data$OsWL,Omit=c("Gum","RGum"),x_lab="O-sWL (ft NGVD 29)",y_lim_min=0,y_lim_max=1.5)
```

![](README_files/figure-gfm/unnamed-chunk-43-1.png)<!-- -->

    ## $AIC
    ##   Distribution      AIC
    ## 1         Gaus 81.27871
    ## 2          Gum       NA
    ## 3         Lapl 92.82204
    ## 4        Logis 84.19031
    ## 5         RGum       NA
    ## 
    ## $Best_fit
    ## [1] "Gaus"

The `Diag_Non_Con_Sel()` function, is similar to the `Diag_Non_Con()`
command, but only plots the probability density function and cumulative
distribution function of a (single) selected univariate distribution in
order to more clearly demonstrate the goodness of fit of a particular
distribution. The options are the Gaussian (`Gaus`), Gumbel (`Gum`),
Laplace (`Lapl`), logistic (`Logis`) and reverse Gumbel (`RGum`)
distributions.

``` r
Diag_Non_Con_Sel(Data=S20.Rainfall$Data$OsWL,x_lab="O-sWL (ft NGVD 29)",
y_lim_min=0,y_lim_max=1.5,Selected="Logis")
```

![](README_files/figure-gfm/unnamed-chunk-44-1.png)<!-- -->

A generalized Pareto distribution is fitted to the marginal distribution
of the conditioning variable i.e. the declustered excesses identified
using `Con_Sampling_2D()`.

The process of selecting a conditional sample and fitting marginal
distributions is repeated but instead conditioning on O-sWL. The
non-conditional variable in this case is (total daily) rainfall, which
has a lower bound at zero, and thus requires a suitably truncated
distribution. The `Diag_Non_Con_Trunc` fits a selection of truncated
distributions to a vector of data. The `Diag_Non_Con_Sel_Trunc` function
is analogous to the `Diag_Non_Con_Sel` function, available distributions
are the Birnbaum-Saunders (`BS`), exponential (`Exp`), gamma (`Gam(2)`),
inverse Gaussian (`InvG`), lognormal (`LogN`), Tweedie (`Twe`) and
Weibull (`Weib`). If the gamlss and gamlss.mx packages are loaded then
the three-parameter gamma (`Gam(3)`), two-parameter mixed gamma
(`GamMix(2)`) and three-parameter mixed gamma (`GamMix(3)`)
distributions are also tested.

``` r
S20.OsWL<-Con_Sampling_2D(Data_Detrend=S20.Detrend.df[,-c(1,4)],
Data_Declust=S20.Detrend.Declustered.df[,-c(1,4)],
Con_Variable="OsWL",u=OsWL.Thres.Quantile)
S20.OsWL$Data$Rainfall<-S20.OsWL$Data$Rainfall+runif(length(S20.OsWL$Data$Rainfall),0.001,0.01)
Diag_Non_Con_Trunc(Data=S20.OsWL$Data$Rainfall+0.001,x_lab="Rainfall (Inches)",
y_lim_min=0,y_lim_max=2)
```

![](README_files/figure-gfm/unnamed-chunk-45-1.png)<!-- -->

    ## $AIC
    ##   Distribution         AIC
    ## 1           BS   0.3278007
    ## 2          Exp 109.1167201
    ## 3       Gam(2)  30.8279174
    ## 4       Gam(3)  31.0401538
    ## 5        LNorm  26.4125980
    ## 6        TNorm 216.4250444
    ## 7          Twe  25.6933297
    ## 8         Weib  28.8704164
    ## 
    ## $Best_fit
    ## [1] "BS"

``` r
Diag_Non_Con_Trunc_Sel(Data=S20.OsWL$Data$Rainfall+0.001,x_lab="Rainfall (Inches)",
y_lim_min=0,y_lim_max=2,
Selected="BS")
```

![](README_files/figure-gfm/unnamed-chunk-46-1.png)<!-- -->

The `Design_Event_2D()` function finds the isoline associated with a
particular return period, by overlaying the two corresponding isolines
from the joint distributions fitted to the conditional samples using the
method in Bender et al. (2016). `Design_Event_2D()` requires the copulas
families chosen to model the dependence structure in the two conditional
samples as input.

``` r
S20.Copula.Rainfall<-Copula_Threshold_2D(Data_Detrend=S20.Detrend.df[,-c(1,4)],
Data_Declust=S20.Detrend.Declustered.df[,-c(1,4)],
u1=Rainfall.Thres.Quantile,u2=NA,
y_lim_min=0,y_lim_max=0.25, GAP=0.075)$Copula_Family_Var1
```

![](README_files/figure-gfm/unnamed-chunk-47-1.png)<!-- -->

``` r
S20.Copula.OsWL<-Copula_Threshold_2D(Data_Detrend=S20.Detrend.df[,-c(1,4)],
Data_Declust=S20.Detrend.Declustered.df[,-c(1,4)],
u1=NA,u2=OsWL.Thres.Quantile,
y_lim_min=0,y_lim_max=0.25,GAP=0.075)$Copula_Family_Var2
```

![](README_files/figure-gfm/unnamed-chunk-48-1.png)<!-- -->

As input the function requires

- `Data` = Original (detrended) rainfall and O-sWL series
- `Data_Con1`/`Data_Con2` = two conditionally sampled data sets,
- `u1`/`u2` or `Thres1`/`Thres2` = two thresholds associated with the
  conditionally sampled data sets
- `Copula_Family1`/`Copula_Family2` two families of the two fitted
  copulas
- `Marginal_Dist1`/`Marginal_Dist2`Selected non-extreme marginal
  distributions
- `RP` = Return Period of interest
- `N` = size of the sample from the fitted joint distributions used to
  estimate the density along the isoline of interest
- `N_Ensemble` = size of the ensemble of events sampled along the
  isoline of interest

``` r
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
```

![](README_files/figure-gfm/unnamed-chunk-49-1.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-49-2.png)<!-- -->

Design event according to the “Most likely” event approach (diamond in
the plot)

``` r
S20.Bivariate$MostLikelyEvent$`100`
```

    ##    Rainfall     OsWL
    ## 1 0.1800065 3.173927

Design event under the assumption of full dependence (triangle in the
plot)

``` r
S20.Bivariate$FullDependence$`100`
```

    ##   Rainfall OsWL
    ## 1    14.56 3.31

#### Cooley (2019) projection method

[Cooley et al. (2019)](https://doi.org/10.1007/s10687-019-00348-0) puts
forward a non-parametric approach for constructing the isoline
associated with exceedance probability $p$. The approach centers around
constructing a base isoline with a larger exceedance probability
$p_{base} < p$ and projecting it to more extreme levels. $p_{base}$
should be small enough to be representative of the extremal dependence
but large enough for sufficient data to be involved in the estimation
procedure.

The approach begins by approximating the joint survival function via a
kernel density estimator from which the base isoline is derived. For the
marginal distributions, a GPD is fit above a sufficiently high threshold
to allow extrapolation into the tails and the empirical distribution is
used below the threshold. Unless the joint distribution of the two
variables is regularly varying, a marginal transformation is required
for the projection. The two marginals are thus transformed to Frechet
scales. For asymptotic dependence, on the transformed scale the isoline
with exceedance probability $p$ can be obtained as
$l_{T}(p)=s^{−1}l_T(p_{base})$ where $\frac{p_{base}}{p} > 1$. For the
case of asymptotic independence,
$l_{T}(p)=s^{\frac{1}{\eta}}l_{T}(p_{base})$, where $\eta$ is the tail
dependence coefficient. Applying the inverse Frechet transformation
gives the isoline on the original scale.

Let’s estimate the 100-year (p=0.01) rainfall-OsWL isoline at S20 using
the 10-year isoline as the base isoline.

``` r
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
```

![](README_files/figure-gfm/unnamed-chunk-52-1.png)<!-- -->

#### Cooley (2019) projection method

When the dependence structure between two variables differs
significantly across conditional samples — as captured by their
respective copulas — it can produce an inflection point in certain
regions of the probability space when overlaying the partial isolines.
While these effects are mathematically sound, they represent
methodological artifacts rather than natural phenomena, potentially
limiting the physical interpretability of results in these regions of
the joint distribution. [Murphy-Barltrop et
al. (2023)](https://doi.org/10.1002/env.2797) proposed new estimation
methods for bivariate isolines that avoid fitting copulas while
accounting for both asymptotic independence (AI) and asymptotic
dependence (AD). Their techniques exploit existing bivariate extreme
value models — specifically the [Heffernan and Tawn
(2004)](https://doi.org/10.1111/j.1467-9868.2004.02050.x) [HT04](#ht04)
and [Wadsworth and Tawn (2013)](https://doi.org/10.3150/12-BEJ471)
\[WT13\] approaches. An additional benefit of these new approaches is
that it is possible to derive confidence intervals that account for the
sampling uncertainty.

When employing the HT04 model to construct the isoline with exceedance
probability $p$, the first step is to convert the variables $(X,Y)$ to
the Laplace scale $(X_L,Y_L)$. The HT04 model is a conditional
exceedance model and is therefore fit twice, conditioning on both $X_L$
and $Y_L$ separately, thus allowing us to estimate the curve in
different regions. In particular, we consider the regions defined by
$R_{Y_L}$ where $y_L>x_L$ and $R_{X_L}$ where \$ y_L x_L\$. For region
$R_{Y_L}$, we start by selecting a high threshold $u_{Y_L}$ such that
$P(Y_L>u_{Y_L})>p$ and fit the HT04 model to observations where
$y_L>u_{Y_L}$. Next, a decreasing series of thresholds is defined in the
interval $(u_{Y_L},F^{-1}{_Y_L}(1-p))$ where $u_{Y_L}$ is the minimal
quantile for which the fitted model is valid and $F^{-1}{_Y_L}(1-p)$ is
the limit that values can attain on this curve. For a given quantile in
the interval $y*$, use the model to simulate from the conditional
distribution $X_L | Y_L>y*$ and estimate $x*$ the $(1-p/q)$th quantile.
Since $q=P(Y_L>y*)$ and
$P(X_L>x* ,Y_L>y*) = P(X_L>x* | Y_L>y*)P(Y_L>y*) = p/q \times q = p$,
the point $(x*,y*)$ lies on the isoline. The process is repeated until a
value $x**$ where $x**\leq y*$ is obtained or we exhaust all values in
the interval. A very similar procedure is implemented for region
$R_{X_L}$, this time selecting $u_{X_L}$, a high quantile of $X_L$ for
quantiles in the interval $(x^{'}**,F^{-1}{_X_L}(1-p))$ where
$(x^{'}**=x**)$ if it exists and $u_{X_L}$ otherwise.

The WT13 model is defined in terms of
$T_w := \min\left\{\frac{X}{w},\frac{Y}{1-w}\right\}$ where $(X,Y)$ are
a pair of variables with standard exponential marginal distributions,
for any ray in $[0,1]$,
$$\Pr(T_w > t) = L(t|w) \exp(-\lambda(w)t), \quad \lambda(w) \geq \max(w,1-w)$$
where $L(\cdot|w)$ is slowly varying for each ray $w \in [0,1]$ and
$\lambda(w)$ is termed the angular dependence function. To find the
isoline with exceedance probability $p < p^*$, a set $W$ of equally
spaced rays on $[0,1]$ is defined and for each ray estimating the
angular dependence function via the Hill estimator $\hat{\lambda}(w)$
using observations above the $95\%$ quantile of the variable $T_w$. For
any large $u$, WT13 states that
$$\Pr(T_w > t + u | T_w > u) \approx \exp(-t\hat{\lambda}(w))$$ for any
$w \in [0,1]$ and $t > 0$. If $u$ is the $(1-p^*)^{th}$ quantile of
$T_w$ i.e., $\Pr(T_w > u) = p^*$, then
$$p=\Pr(T_w > t + u) = \Pr(T_w > t + u | T_w > u)\Pr(T_w > u) = p^*\exp(-t\hat{\lambda}(w))$$.
Re-arranging for $t$ gives
$t = -\frac{1}{\hat{\lambda}(w)}\log\left(\frac{p}{p^*}\right)$. An
estimate for the return curve with exceedance probability $p$ is
obtained by letting $(x,y) = (w(t+u),(1-w)(t+u))$.

The confidence intervals for isolines are estimated using a bootstrap
procedure following the methodology outlined in Murphy-Barltrop et
al. (2023). In the procedure a set of rays are defined at equally-spaced
angles within the interval $(0,\pi/2)$. Each ray intersects the isoline
exactly once, and since the angle is fixed, the L2 norm represents the
radial distance to the intersection point. Next, the observed data set
is bootstrapped a large number of times (e.g., 1000 iterations) to
generate multiple samples of the same size as the original data set. For
each bootstrap sample, the L2 norm is calculated for each ray’s
intersection with the fitted isoline. The confidence intervals are then
derived by computing the relevant percentiles (e.g., 2.5th and 97.5th
percentiles for 95% confidence intervals) of these L2 norms across all
bootstrap iterations.

The `return_curve_est` function derives isolines for a given `rp`. The
quantiles `q` of the GPDs of the HT04 and WT13 models must be specified
along with the average occurrence frequency of the events in the data
`mu`. The methods for declustering the time series and the associated
parameters are also required. The bootstrapping procedure for estimate
the sampling uncertainty can be carried out using a basic
(`boot_method = "basic"`), block (`boot_method = "block"`) or monthly
(`boot_method = "month"`) bootstrap. The latter two are recommend where
there is a temporal dependence in the data. For the basic bootstrap
whether the sampling is carried out with `boot_replace = T` or without
`boot_replace = F` replacement must be specified while the block
bootstrap require `"block_length"`. The number of number of rays along
which to compute points on the curve for each sample `n_grad` is another
input. For the HT04 model the number of simulations `n_sim` is required.

``` r
#Adding dates to complete final month of combined records
final.month = data.frame(seq(as.Date("2019-02-04"),as.Date("2019-02-28"),by="day"),NA,NA,NA)
colnames(final.month) = c("Date","Rainfall","OsWL","Groundwater")
S22.Detrend.df.extended = rbind(S22.Detrend.df,final.month)

#Derive return curves
return_curve_est(data=S22.Detrend.df.extended[,1:3],
                 q=0.985,rp=100,mu=365.25,n_sim=100,
                 n_grad=50,n_boot=100,boot_method="monthly",
                 boot_replace=NA, block_length=NA, boot_prop=0.8,
                 decl_method_x="runs", decl_method_y="runs",
                 window_length_x=NA,window_length_y=NA,
                 u_x=0.95, u_y=0.95,
                 sep_crit_x=36, sep_crit_y=36,
                 alpha=0.1, x_lab=NA, y_lab=NA)
```

![](README_files/figure-gfm/unnamed-chunk-53-1.png)<!-- -->

    ## $median_ht04
    ##         Rainfall     OsWL
    ##  [1,]  0.2725040 9.865885
    ##  [2,]  0.5455255 9.865881
    ##  [3,]  0.7922189 9.570535
    ##  [4,]  1.0317140 9.353040
    ##  [5,]  1.2296552 8.942701
    ##  [6,]  1.3849222 8.429811
    ##  [7,]  1.5758608 8.216731
    ##  [8,]  1.7572968 8.008114
    ##  [9,]  1.9451421 7.857596
    ## [10,]  2.0931769 7.600905
    ## [11,]  2.2765065 7.481387
    ## [12,]  2.4204634 7.269067
    ## [13,]  2.6087465 7.184538
    ## [14,]  2.7941337 7.094371
    ## [15,]  2.9404456 6.926342
    ## [16,]  3.1504169 6.887953
    ## [17,]  3.3646027 6.848794
    ## [18,]  3.4862962 6.651693
    ## [19,]  3.6343396 6.505900
    ## [20,]  3.8278617 6.428591
    ## [21,]  4.0756698 6.418190
    ## [22,]  4.2273612 6.274489
    ## [23,]  4.4803610 6.250631
    ## [24,]  4.6467729 6.118398
    ## [25,]  4.9250194 6.100226
    ## [26,]  5.1553572 6.020100
    ## [27,]  5.3680588 5.914763
    ## [28,]  5.6948940 5.900217
    ## [29,]  6.0098046 5.857203
    ## [30,]  6.2587821 5.747544
    ## [31,]  6.5697446 5.671757
    ## [32,]  6.9092458 5.599363
    ## [33,]  7.1688477 5.459894
    ## [34,]  7.5377731 5.373066
    ## [35,]  7.7488042 5.182151
    ## [36,]  8.1486326 5.078668
    ## [37,]  8.4708745 4.918352
    ## [38,]  9.0040015 4.832198
    ## [39,]  9.5249821 4.711129
    ## [40,] 10.0154773 4.550449
    ## [41,] 10.4497319 4.345430
    ## [42,] 10.6132674 4.040866
    ## [43,] 10.8847916 3.758766
    ## [44,] 11.2794552 3.491370
    ## [45,] 11.7976285 3.226490
    ## [46,] 12.7835612 3.005507
    ## [47,] 14.4394320 2.809120
    ## [48,] 15.1039166 2.420715
    ## [49,] 15.3741795 1.969378
    ## [50,] 15.6015263 1.501809
    ## 
    ## $ub_ht04
    ##         Rainfall      OsWL
    ##  [1,]  0.8501701 28.615392
    ##  [2,]  1.4992689 25.329177
    ##  [3,]  2.0863271 23.536186
    ##  [4,]  2.3642352 20.114207
    ##  [5,]  2.8030110 19.078428
    ##  [6,]  3.1641276 17.947719
    ##  [7,]  3.6030472 17.473147
    ##  [8,]  3.9821019 16.853903
    ##  [9,]  4.2891425 16.095903
    ## [10,]  4.3929281 14.830034
    ## [11,]  4.6204176 14.132924
    ## [12,]  5.0278851 13.999599
    ## [13,]  5.4366962 13.865833
    ## [14,]  5.5010300 12.977994
    ## [15,]  5.9030018 12.875957
    ## [16,]  6.2783324 12.712872
    ## [17,]  6.5059788 12.289817
    ## [18,]  6.7615388 11.941391
    ## [19,]  7.1337748 11.787077
    ## [20,]  7.5867378 11.738598
    ## [21,]  8.0555633 11.688421
    ## [22,]  8.2643075 11.291217
    ## [23,]  8.5329455 10.980827
    ## [24,]  8.5738818 10.426231
    ## [25,]  8.8805880 10.179541
    ## [26,]  8.8934847  9.644830
    ## [27,]  9.2248345  9.430679
    ## [28,]  9.3157598  9.002385
    ## [29,]  9.3157598  8.517495
    ## [30,]  9.5338083  8.220729
    ## [31,] 10.1491445  8.205562
    ## [32,] 10.6703821  8.091583
    ## [33,] 10.9616251  7.808282
    ## [34,] 11.2172658  7.497422
    ## [35,] 11.3803246  7.132236
    ## [36,] 11.8769562  6.935152
    ## [37,] 12.0277241  6.554762
    ## [38,] 12.8104580  6.443333
    ## [39,] 12.8560704  6.001599
    ## [40,] 13.1140467  5.642343
    ## [41,] 13.6772725  5.372185
    ## [42,] 14.1980717  5.060831
    ## [43,] 14.6079921  4.695188
    ## [44,] 14.9709441  4.299818
    ## [45,] 15.6484186  3.946327
    ## [46,] 16.1776618  3.532369
    ## [47,] 17.1554323  3.145434
    ## [48,] 20.2899258  2.901269
    ## [49,] 20.9244288  2.311706
    ## [50,] 21.7228761  1.690406
    ## 
    ## $lb_ht04
    ##         Rainfall     OsWL
    ##  [1,]  0.1115031 4.640226
    ##  [2,]  0.2232181 4.640226
    ##  [3,]  0.3341600 4.627292
    ##  [4,]  0.4433195 4.601288
    ##  [5,]  0.5393527 4.495699
    ##  [6,]  0.6471972 4.483330
    ##  [7,]  0.7500000 4.445736
    ##  [8,]  0.8493669 4.398200
    ##  [9,]  0.9535623 4.372555
    ## [10,]  1.0635139 4.364222
    ## [11,]  1.1701847 4.341874
    ## [12,]  1.2571509 4.266211
    ## [13,]  1.3710114 4.260273
    ## [14,]  1.4874482 4.254201
    ## [15,]  1.6067786 4.247977
    ## [16,]  1.7217833 4.227499
    ## [17,]  1.8492957 4.224205
    ## [18,]  1.9811506 4.220799
    ## [19,]  2.1178374 4.217268
    ## [20,]  2.2432055 4.190013
    ## [21,]  2.3861845 4.180950
    ## [22,]  2.5415551 4.179531
    ## [23,]  2.7046752 4.178041
    ## [24,]  2.8674884 4.166616
    ## [25,]  3.0452195 4.161618
    ## [26,]  3.2270361 4.150274
    ## [27,]  3.4296662 4.147685
    ## [28,]  3.6461007 4.144919
    ## [29,]  3.8782559 4.141952
    ## [30,]  4.1284014 4.138755
    ## [31,]  4.3612459 4.108393
    ## [32,]  4.6035764 4.071571
    ## [33,]  4.8377701 4.016551
    ## [34,]  5.1033017 3.967524
    ## [35,]  5.4525263 3.949076
    ## [36,]  5.7660829 3.892299
    ## [37,]  5.9842037 3.774302
    ## [38,]  6.3400725 3.704653
    ## [39,]  6.9067564 3.696823
    ## [40,]  7.2367588 3.571266
    ## [41,]  7.6485902 3.454323
    ## [42,]  8.2967624 3.381763
    ## [43,]  8.6260027 3.190658
    ## [44,]  9.2496310 3.046831
    ## [45,]  9.5855036 2.812972
    ## [46,]  9.6476315 2.518720
    ## [47,]  9.6946646 2.221591
    ## [48,]  9.7417957 1.923841
    ## [49,]  9.7891150 1.624903
    ## [50,] 10.8681232 1.355974
    ## 
    ## $contour_ht04
    ## [1] NA
    ## 
    ## $most_likely_ht04
    ## [1] NA
    ## 
    ## $ensemble_ht04
    ## [1] NA
    ## 
    ## $median_wt13
    ##         Rainfall     OsWL
    ##  [1,]  0.2174872 8.080187
    ##  [2,]  0.4353877 8.080187
    ##  [3,]  0.6169207 7.678766
    ##  [4,]  0.7583855 7.145696
    ##  [5,]  0.8801972 6.691455
    ##  [6,]  1.0567857 6.674435
    ##  [7,]  1.1894129 6.452156
    ##  [8,]  1.3483722 6.382237
    ##  [9,]  1.4803901 6.224162
    ## [10,]  1.6054017 6.067614
    ## [11,]  1.7419147 5.964326
    ## [12,]  1.8512270 5.799699
    ## [13,]  1.9948415 5.734130
    ## [14,]  2.1495497 5.693324
    ## [15,]  2.3103665 5.660972
    ## [16,]  2.4733301 5.627057
    ## [17,]  2.6447220 5.601924
    ## [18,]  2.8281630 5.588771
    ## [19,]  2.9950823 5.541164
    ## [20,]  3.1388596 5.455267
    ## [21,]  3.2760584 5.359333
    ## [22,]  3.4470886 5.304841
    ## [23,]  3.6003309 5.223455
    ## [24,]  3.7771600 5.164478
    ## [25,]  3.9572244 5.102154
    ## [26,]  4.1121721 5.008560
    ## [27,]  4.2742793 4.917651
    ## [28,]  4.4602781 4.842463
    ## [29,]  4.6837884 4.790162
    ## [30,]  4.8901536 4.714004
    ## [31,]  5.1209287 4.646161
    ## [32,]  5.3824202 4.587652
    ## [33,]  5.5342156 4.447772
    ## [34,]  5.6863434 4.304143
    ## [35,]  5.8629621 4.169475
    ## [36,]  6.2874911 4.151930
    ## [37,]  6.6014917 4.058299
    ## [38,]  6.8669290 3.927653
    ## [39,]  7.1707615 3.799099
    ## [40,]  7.5461536 3.680293
    ## [41,]  7.9927741 3.563815
    ## [42,]  8.5968172 3.467136
    ## [43,]  9.3217665 3.365649
    ## [44,] 10.1797219 3.250524
    ## [45,] 11.1148438 3.098855
    ## [46,] 12.0719021 2.895037
    ## [47,] 12.7943696 2.605417
    ## [48,] 13.6735249 2.288170
    ## [49,] 15.1475945 1.955403
    ## [50,] 15.5711021 1.500872
    ## 
    ## $ub_wt13
    ##         Rainfall      OsWL
    ##  [1,]  0.5841520 19.981153
    ##  [2,]  0.8011109 14.009755
    ##  [3,]  1.1348452 13.268062
    ##  [4,]  1.2774295 11.337389
    ##  [5,]  1.5188744 10.805883
    ##  [6,]  1.7698301 10.488886
    ##  [7,]  2.0320734 10.299862
    ##  [8,]  2.3035441 10.179984
    ##  [9,]  2.5443861  9.963721
    ## [10,]  2.7899565  9.791190
    ## [11,]  3.0651150  9.719296
    ## [12,]  3.3503064  9.669269
    ## [13,]  3.4969932  9.283104
    ## [14,]  3.5881752  8.820274
    ## [15,]  3.8175458  8.687797
    ## [16,]  4.0648409  8.590827
    ## [17,]  4.3014521  8.471465
    ## [18,]  4.4841564  8.263292
    ## [19,]  4.6605715  8.054639
    ## [20,]  4.8089732  7.814566
    ## [21,]  4.9790610  7.614473
    ## [22,]  5.2568562  7.553846
    ## [23,]  5.4889440  7.427854
    ## [24,]  5.7903955  7.372892
    ## [25,]  6.0738102  7.284955
    ## [26,]  6.2581071  7.089397
    ## [27,]  6.4943324  6.941497
    ## [28,]  6.8130966  6.858235
    ## [29,]  7.1573154  6.780601
    ## [30,]  7.5305511  6.707939
    ## [31,]  7.9580057  6.654487
    ## [32,]  8.1747724  6.437932
    ## [33,]  8.4993048  6.283678
    ## [34,]  9.0434065  6.242344
    ## [35,]  9.6339064  6.194429
    ## [36,] 10.1953515  6.097813
    ## [37,] 10.5766647  5.887169
    ## [38,] 10.7638374  5.577073
    ## [39,] 10.8255413  5.214968
    ## [40,] 11.2820719  4.996780
    ## [41,] 11.4755750  4.671774
    ## [42,] 11.6304732  4.330286
    ## [43,] 11.8876721  4.011001
    ## [44,] 12.4592543  3.749750
    ## [45,] 13.1536583  3.479976
    ## [46,] 13.8964438  3.178258
    ## [47,] 15.4699538  2.936726
    ## [48,] 18.1153016  2.699761
    ## [49,] 19.7846771  2.241408
    ## [50,] 21.7630021  1.691642
    ## 
    ## $lb_wt13
    ##         Rainfall     OsWL
    ##  [1,]  0.1107750 4.616592
    ##  [2,]  0.2217604 4.616592
    ##  [3,]  0.3331685 4.616592
    ##  [4,]  0.4452145 4.616592
    ##  [5,]  0.5529042 4.582999
    ##  [6,]  0.6531681 4.515271
    ##  [7,]  0.7461882 4.428331
    ##  [8,]  0.8529508 4.412449
    ##  [9,]  0.9624612 4.403831
    ## [10,]  1.0550979 4.337767
    ## [11,]  1.1486620 4.280797
    ## [12,]  1.2493609 4.246103
    ## [13,]  1.3467666 4.202993
    ## [14,]  1.4426755 4.156884
    ## [15,]  1.5309737 4.095741
    ## [16,]  1.6322580 4.060782
    ## [17,]  1.7301412 4.017824
    ## [18,]  1.8498303 4.008710
    ## [19,]  1.9784527 4.006916
    ## [20,]  2.0962331 3.982392
    ## [21,]  2.2146435 3.953793
    ## [22,]  2.3349206 3.922746
    ## [23,]  2.4459100 3.876009
    ## [24,]  2.5610344 3.830452
    ## [25,]  2.6733916 3.778157
    ## [26,]  2.8134238 3.749209
    ## [27,]  2.9407300 3.701961
    ## [28,]  3.0696393 3.651037
    ## [29,]  3.1874477 3.586061
    ## [30,]  3.3483644 3.549698
    ## [31,]  3.5221369 3.514400
    ## [32,]  3.7180301 3.484787
    ## [33,]  3.9182322 3.447198
    ## [34,]  4.1364642 3.409320
    ## [35,]  4.3712941 3.368466
    ## [36,]  4.6601292 3.341600
    ## [37,]  4.9834539 3.313884
    ## [38,]  5.3193828 3.272633
    ## [39,]  5.7343904 3.242646
    ## [40,]  6.2239891 3.214380
    ## [41,]  6.7081508 3.155147
    ## [42,]  7.2342770 3.079460
    ## [43,]  7.9783234 3.027760
    ## [44,]  8.5830397 2.900846
    ## [45,]  8.9075630 2.686243
    ## [46,]  9.4409716 2.486641
    ## [47,]  9.5271464 2.200847
    ## [48,]  9.5989228 1.910602
    ## [49,]  9.7266879 1.621053
    ## [50,] 10.8440608 1.355233
    ## 
    ## $contour_wt13
    ## [1] NA
    ## 
    ## $most_likely_wt13
    ## [1] NA
    ## 
    ## $ensemble_wt13
    ## [1] NA

The `return_curve_diag()` function calculates the empirical probability
of observing data within the survival regions defined by a subset of
points on the return curve. If the curve is a good fit, the empirical
probabilities should closely match the probabilities associated with the
return level curve. The procedure which is introduced in
[Murphy-Barltrop et al. (2023)](https://doi.org/10.1002/env.2797) uses
bootstrap resampling of the original data set to obtain confidence
intervals for the empirical estimates. Since observations are required
in the survival regions to estimate the empirical probabilities, it is
recommended that this function be run for shorter return periods that
may usually be considered for design e.g., for a 1-year return period
rather than a 50-year return period. The inputs are almost identical to
those of the `return_curve_est` function. The only additional arguments
are `boot_method_all` which details the bootstrapping procedure - basic
or block - to use when estimating the distribution of empirical
(survival) probabilities from the original data set (without any
declustering), `boot_replace_all` which specifies whether to sample with
replacement in the case of a basic bootstrap and `block_length_all`
which specifies the block length for the block bootstrap.

``` r
#Diagnostic plots for the return curves
return_curve_diag(data=S22.Detrend.df.extended[,1:3],
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
```

![](README_files/figure-gfm/unnamed-chunk-54-1.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-54-2.png)<!-- -->

    ## $ang_ind
    ##  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
    ## [26] 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50
    ## 
    ## $med_x_ht04
    ##          50%          50%          50%          50%          50%          50% 
    ## 0.0006084191 0.0006084191 0.0006544559 0.0006544559 0.0008180699 0.0008180699 
    ##          50%          50%          50%          50%          50%          50% 
    ## 0.0008180699 0.0007953076 0.0009059629 0.0010795842 0.0009434740 0.0011120695 
    ##          50%          50%          50%          50%          50%          50% 
    ## 0.0010729798 0.0009409752 0.0009409752 0.0009156780 0.0009156780 0.0010345279 
    ##          50%          50%          50%          50%          50%          50% 
    ## 0.0011898565 0.0013030996 0.0013452227 0.0013218726 0.0013218726 0.0013563565 
    ##          50%          50%          50%          50%          50%          50% 
    ## 0.0013563565 0.0017153997 0.0020197312 0.0021537329 0.0022487314 0.0024821999 
    ##          50%          50%          50%          50%          50%          50% 
    ## 0.0026623731 0.0025714720 0.0026684807 0.0028866649 0.0025950124 0.0024626504 
    ##          50%          50%          50%          50%          50%          50% 
    ## 0.0022923336 0.0024368624 0.0023565938 0.0025847505 0.0022677330 0.0020197312 
    ##          50%          50%          50%          50%          50%          50% 
    ## 0.0019710547 0.0017870465 0.0014678985 0.0014011072 0.0014011072 0.0014011072 
    ##          50%          50% 
    ## 0.0014011072 0.0014011072 
    ## 
    ## $lb_x_ht04
    ##           5%           5%           5%           5%           5%           5% 
    ## 0.0001559758 0.0001559758 0.0001605507 0.0001605507 0.0003227813 0.0003227813 
    ##           5%           5%           5%           5%           5%           5% 
    ## 0.0003227813 0.0003095776 0.0003112069 0.0004451470 0.0003059909 0.0004341306 
    ##           5%           5%           5%           5%           5%           5% 
    ## 0.0003095776 0.0003095776 0.0003095776 0.0003095776 0.0003095776 0.0004396828 
    ##           5%           5%           5%           5%           5%           5% 
    ## 0.0006323351 0.0007550595 0.0007550595 0.0006189964 0.0006189964 0.0007620331 
    ##           5%           5%           5%           5%           5%           5% 
    ## 0.0007620331 0.0009485027 0.0011949167 0.0013185367 0.0014616903 0.0016116659 
    ##           5%           5%           5%           5%           5%           5% 
    ## 0.0019062535 0.0018396393 0.0018641614 0.0019693507 0.0016735975 0.0016090197 
    ##           5%           5%           5%           5%           5%           5% 
    ## 0.0013806019 0.0016626796 0.0015605128 0.0016923991 0.0014179549 0.0012596441 
    ##           5%           5%           5%           5%           5%           5% 
    ## 0.0011304462 0.0011705413 0.0009225815 0.0008075203 0.0008075203 0.0008075203 
    ##           5%           5% 
    ## 0.0008075203 0.0008075203 
    ## 
    ## $ub_x_ht04
    ##         95%         95%         95%         95%         95%         95% 
    ## 0.001053636 0.001053636 0.001218433 0.001218433 0.001521155 0.001521155 
    ##         95%         95%         95%         95%         95%         95% 
    ## 0.001521155 0.001459837 0.001582203 0.001784878 0.001594407 0.001821059 
    ##         95%         95%         95%         95%         95%         95% 
    ## 0.001642715 0.001582341 0.001582341 0.001506158 0.001506158 0.001650799 
    ##         95%         95%         95%         95%         95%         95% 
    ## 0.001858359 0.002100453 0.002100453 0.002100453 0.002100453 0.002174913 
    ##         95%         95%         95%         95%         95%         95% 
    ## 0.002174913 0.002500250 0.002876487 0.002999889 0.003150180 0.003326434 
    ##         95%         95%         95%         95%         95%         95% 
    ## 0.003517416 0.003505186 0.003634996 0.003639340 0.003406239 0.003327032 
    ##         95%         95%         95%         95%         95%         95% 
    ## 0.003184183 0.003270572 0.003172529 0.003364172 0.002919595 0.002673306 
    ##         95%         95%         95%         95%         95%         95% 
    ## 0.002636455 0.002504713 0.002158624 0.002088294 0.002088294 0.002088294 
    ##         95%         95% 
    ## 0.002088294 0.002088294 
    ## 
    ## $med_y_ht04
    ##          50%          50%          50%          50%          50%          50% 
    ## 0.0007510036 0.0008459802 0.0009781058 0.0010830149 0.0014138831 0.0014138831 
    ##          50%          50%          50%          50%          50%          50% 
    ## 0.0014138831 0.0013269131 0.0014194491 0.0015849995 0.0013977765 0.0013962731 
    ##          50%          50%          50%          50%          50%          50% 
    ## 0.0012812347 0.0011867551 0.0011867551 0.0010877857 0.0010877857 0.0010861449 
    ##          50%          50%          50%          50%          50%          50% 
    ## 0.0011829890 0.0012019283 0.0012726731 0.0011793777 0.0011991061 0.0011816523 
    ##          50%          50%          50%          50%          50%          50% 
    ## 0.0011816523 0.0011993022 0.0011781735 0.0012084666 0.0013183272 0.0014793738 
    ##          50%          50%          50%          50%          50%          50% 
    ## 0.0015085568 0.0014581685 0.0014867612 0.0016514369 0.0015973623 0.0014306939 
    ##          50%          50%          50%          50%          50%          50% 
    ## 0.0012948527 0.0012976480 0.0012368601 0.0013193349 0.0011335074 0.0009607686 
    ##          50%          50%          50%          50%          50%          50% 
    ## 0.0009568012 0.0007667033 0.0007486656 0.0006517866 0.0006517866 0.0006517866 
    ##          50%          50% 
    ## 0.0006517866 0.0006517866 
    ## 
    ## $lb_y_ht04
    ##           5%           5%           5%           5%           5%           5% 
    ## 0.0004173214 0.0004359424 0.0006284637 0.0007320108 0.0008573863 0.0008573863 
    ##           5%           5%           5%           5%           5%           5% 
    ## 0.0008573863 0.0007611562 0.0008385809 0.0009808704 0.0007674073 0.0008390004 
    ##           5%           5%           5%           5%           5%           5% 
    ## 0.0007322699 0.0006327980 0.0006327980 0.0006266845 0.0006266845 0.0005463466 
    ##           5%           5%           5%           5%           5%           5% 
    ## 0.0006285697 0.0006371878 0.0006526641 0.0006328313 0.0006371878 0.0006283161 
    ##           5%           5%           5%           5%           5%           5% 
    ## 0.0006283161 0.0006459428 0.0006377206 0.0006476511 0.0007513505 0.0008596804 
    ##           5%           5%           5%           5%           5%           5% 
    ## 0.0009048082 0.0008545877 0.0009463742 0.0009773935 0.0008699622 0.0008386623 
    ##           5%           5%           5%           5%           5%           5% 
    ## 0.0007509012 0.0007509012 0.0006465261 0.0007378331 0.0006316929 0.0005268623 
    ##           5%           5%           5%           5%           5%           5% 
    ## 0.0005268623 0.0003332827 0.0003168284 0.0002112189 0.0002112189 0.0002112189 
    ##           5%           5% 
    ## 0.0002112189 0.0002112189 
    ## 
    ## $ub_y_ht04
    ##         95%         95%         95%         95%         95%         95% 
    ## 0.001191033 0.001296899 0.001508626 0.001695223 0.002047025 0.002047025 
    ##         95%         95%         95%         95%         95%         95% 
    ## 0.002047025 0.001923055 0.002031342 0.002259620 0.002095779 0.001947403 
    ##         95%         95%         95%         95%         95%         95% 
    ## 0.001921990 0.001848758 0.001848758 0.001815213 0.001815213 0.001811306 
    ##         95%         95%         95%         95%         95%         95% 
    ## 0.001834436 0.001918100 0.002017329 0.001821955 0.001939500 0.001855796 
    ##         95%         95%         95%         95%         95%         95% 
    ## 0.001855796 0.001924429 0.001830253 0.001852121 0.002167001 0.002275162 
    ##         95%         95%         95%         95%         95%         95% 
    ## 0.002141939 0.002072522 0.002179313 0.002195886 0.002287798 0.002006521 
    ##         95%         95%         95%         95%         95%         95% 
    ## 0.001946262 0.001841247 0.001740272 0.001950301 0.001702998 0.001428483 
    ##         95%         95%         95%         95%         95%         95% 
    ## 0.001428483 0.001209324 0.001221095 0.001113035 0.001113035 0.001113035 
    ##         95%         95% 
    ## 0.001113035 0.001113035 
    ## 
    ## $med_x_wt13
    ##          50%          50%          50%          50%          50%          50% 
    ## 0.0006544559 0.0014290799 0.0015274177 0.0015974441 0.0022452121 0.0023227809 
    ##          50%          50%          50%          50%          50%          50% 
    ## 0.0023010219 0.0021840882 0.0020823680 0.0022573072 0.0023978756 0.0020772110 
    ##          50%          50%          50%          50%          50%          50% 
    ## 0.0020790025 0.0018191232 0.0018356726 0.0018356726 0.0019692554 0.0021822160 
    ##          50%          50%          50%          50%          50%          50% 
    ## 0.0021822160 0.0023335444 0.0023694854 0.0025628942 0.0026338407 0.0028158739 
    ##          50%          50%          50%          50%          50%          50% 
    ## 0.0028642887 0.0028755930 0.0031983204 0.0030432543 0.0035010552 0.0037125044 
    ##          50%          50%          50%          50%          50%          50% 
    ## 0.0037114015 0.0037935261 0.0037599425 0.0037665167 0.0037107533 0.0034463021 
    ##          50%          50%          50%          50%          50%          50% 
    ## 0.0032534528 0.0031564616 0.0029915384 0.0029821995 0.0028945000 0.0026938420 
    ##          50%          50%          50%          50%          50%          50% 
    ## 0.0024282896 0.0020320329 0.0018400171 0.0014011072 0.0014011072 0.0014011072 
    ##          50%          50% 
    ## 0.0014011072 0.0014011072 
    ## 
    ## $lb_x_wt13
    ##           5%           5%           5%           5%           5%           5% 
    ## 0.0001605507 0.0005995816 0.0007492852 0.0007788278 0.0011859237 0.0013997909 
    ##           5%           5%           5%           5%           5%           5% 
    ## 0.0012930950 0.0011919523 0.0010731418 0.0012309048 0.0012830644 0.0010646805 
    ##           5%           5%           5%           5%           5%           5% 
    ## 0.0011721206 0.0009686772 0.0010259266 0.0010259266 0.0011986649 0.0012929680 
    ##           5%           5%           5%           5%           5%           5% 
    ## 0.0012929680 0.0013247680 0.0016163232 0.0017923751 0.0017543220 0.0019420311 
    ##           5%           5%           5%           5%           5%           5% 
    ## 0.0020558114 0.0019955584 0.0023361025 0.0022664588 0.0026329181 0.0027657115 
    ##           5%           5%           5%           5%           5%           5% 
    ## 0.0027839893 0.0028148745 0.0028058161 0.0028044857 0.0027649580 0.0025895684 
    ##           5%           5%           5%           5%           5%           5% 
    ## 0.0022664588 0.0022492888 0.0020507249 0.0020457853 0.0019674228 0.0017524166 
    ##           5%           5%           5%           5%           5%           5% 
    ## 0.0015414023 0.0013604571 0.0012862137 0.0008075203 0.0008075203 0.0008075203 
    ##           5%           5% 
    ## 0.0008075203 0.0008075203 
    ## 
    ## $ub_x_wt13
    ##         95%         95%         95%         95%         95%         95% 
    ## 0.001218433 0.002270751 0.002215939 0.002484984 0.003170058 0.003167238 
    ##         95%         95%         95%         95%         95%         95% 
    ## 0.003108398 0.003079801 0.002848633 0.003126007 0.003386497 0.002972094 
    ##         95%         95%         95%         95%         95%         95% 
    ## 0.003043104 0.002789921 0.002570931 0.002570931 0.002738524 0.002897331 
    ##         95%         95%         95%         95%         95%         95% 
    ## 0.002897331 0.003010713 0.003243099 0.003451704 0.003501642 0.003674309 
    ##         95%         95%         95%         95%         95%         95% 
    ## 0.003801415 0.003838336 0.004299089 0.004057981 0.004432688 0.004606472 
    ##         95%         95%         95%         95%         95%         95% 
    ## 0.004666239 0.004687848 0.004859386 0.004794224 0.004452368 0.004191582 
    ##         95%         95%         95%         95%         95%         95% 
    ## 0.004024828 0.003974477 0.003720799 0.003890983 0.003722593 0.003477667 
    ##         95%         95%         95%         95%         95%         95% 
    ## 0.003080480 0.002754735 0.002670038 0.002088294 0.002088294 0.002088294 
    ##         95%         95% 
    ## 0.002088294 0.002088294 
    ## 
    ## $med_y_wt13
    ##          50%          50%          50%          50%          50%          50% 
    ## 0.0009781058 0.0017091325 0.0018987601 0.0019774319 0.0022023782 0.0022717445 
    ##          50%          50%          50%          50%          50%          50% 
    ## 0.0022505470 0.0023304368 0.0022235572 0.0022930824 0.0024445985 0.0022935611 
    ##          50%          50%          50%          50%          50%          50% 
    ## 0.0021542286 0.0018746988 0.0018079635 0.0018780969 0.0016944057 0.0016944057 
    ##          50%          50%          50%          50%          50%          50% 
    ## 0.0015745554 0.0015663560 0.0017668795 0.0019020520 0.0020200242 0.0020272077 
    ##          50%          50%          50%          50%          50%          50% 
    ## 0.0022030669 0.0021550811 0.0023482213 0.0022797614 0.0024107446 0.0024658388 
    ##          50%          50%          50%          50%          50%          50% 
    ## 0.0025276864 0.0024173918 0.0024316423 0.0024763315 0.0021663809 0.0019765111 
    ##          50%          50%          50%          50%          50%          50% 
    ## 0.0019315385 0.0018263864 0.0016727312 0.0017211724 0.0017084961 0.0014439809 
    ##          50%          50%          50%          50%          50%          50% 
    ## 0.0012230268 0.0009607686 0.0008590604 0.0006517866 0.0006517866 0.0006517866 
    ##          50%          50% 
    ## 0.0006517866 0.0006517866 
    ## 
    ## $lb_y_wt13
    ##           5%           5%           5%           5%           5%           5% 
    ## 0.0006284637 0.0011747384 0.0011684295 0.0012054607 0.0015096481 0.0016177423 
    ##           5%           5%           5%           5%           5%           5% 
    ## 0.0014833392 0.0013985426 0.0013661645 0.0013768209 0.0015833658 0.0014904481 
    ##           5%           5%           5%           5%           5%           5% 
    ## 0.0014914155 0.0010851347 0.0010933807 0.0011549803 0.0010629177 0.0010629177 
    ##           5%           5%           5%           5%           5%           5% 
    ## 0.0009565809 0.0009569938 0.0010759940 0.0012732986 0.0012816382 0.0012895206 
    ##           5%           5%           5%           5%           5%           5% 
    ## 0.0014934649 0.0014980853 0.0014917135 0.0014917135 0.0014727884 0.0014734271 
    ##           5%           5%           5%           5%           5%           5% 
    ## 0.0015896788 0.0015734057 0.0016100527 0.0017567323 0.0014020472 0.0014057991 
    ##           5%           5%           5%           5%           5%           5% 
    ## 0.0012760237 0.0012465747 0.0010399796 0.0011524195 0.0010780516 0.0008429796 
    ##           5%           5%           5%           5%           5%           5% 
    ## 0.0007306372 0.0005407639 0.0004294186 0.0002112189 0.0002112189 0.0002112189 
    ##           5%           5% 
    ## 0.0002112189 0.0002112189 
    ## 
    ## $ub_y_wt13
    ##         95%         95%         95%         95%         95%         95% 
    ## 0.001508626 0.002369864 0.002558398 0.002594220 0.002766746 0.002926860 
    ##         95%         95%         95%         95%         95%         95% 
    ## 0.002906091 0.003099516 0.002992636 0.003009645 0.003251391 0.003178344 
    ##         95%         95%         95%         95%         95%         95% 
    ## 0.003033735 0.002716324 0.002767102 0.002666670 0.002493802 0.002493802 
    ##         95%         95%         95%         95%         95%         95% 
    ## 0.002395005 0.002286509 0.002677348 0.002811581 0.002911085 0.002885364 
    ##         95%         95%         95%         95%         95%         95% 
    ## 0.003044405 0.002944144 0.003336478 0.003212588 0.003227902 0.003335437 
    ##         95%         95%         95%         95%         95%         95% 
    ## 0.003567953 0.003381854 0.003235036 0.003245798 0.003048707 0.002795783 
    ##         95%         95%         95%         95%         95%         95% 
    ## 0.002609490 0.002495132 0.002302311 0.002354904 0.002423208 0.002011937 
    ##         95%         95%         95%         95%         95%         95% 
    ## 0.001844490 0.001429143 0.001318773 0.001113035 0.001113035 0.001113035 
    ##         95%         95% 
    ## 0.001113035 0.001113035 
    ## 
    ## $med_ht04
    ##          50%          50%          50%          50%          50%          50% 
    ## 0.0007002801 0.0007414731 0.0009062449 0.0009886307 0.0011534025 0.0011534025 
    ##          50%          50%          50%          50%          50%          50% 
    ## 0.0013181743 0.0014005602 0.0014005602 0.0015653320 0.0015653320 0.0015653320 
    ##          50%          50%          50%          50%          50%          50% 
    ## 0.0014829461 0.0014829461 0.0014005602 0.0012357884 0.0011534025 0.0012357884 
    ##          50%          50%          50%          50%          50%          50% 
    ## 0.0014829461 0.0015653320 0.0018124897 0.0016477179 0.0018124897 0.0016477179 
    ##          50%          50%          50%          50%          50%          50% 
    ## 0.0015653320 0.0018124897 0.0020184544 0.0023068051 0.0023891910 0.0027187346 
    ##          50%          50%          50%          50%          50%          50% 
    ## 0.0025951557 0.0023068051 0.0023891910 0.0023068051 0.0019772615 0.0018124897 
    ##          50%          50%          50%          50%          50%          50% 
    ## 0.0016477179 0.0016477179 0.0014829461 0.0014829461 0.0013181743 0.0011534025 
    ##          50%          50%          50%          50%          50%          50% 
    ## 0.0010710166 0.0009062449 0.0008238590 0.0007002801 0.0007002801 0.0007002801 
    ##          50%          50% 
    ## 0.0007002801 0.0007002801 
    ## 
    ## $lb_ht04
    ##           5%           5%           5%           5%           5%           5% 
    ## 0.0003295436 0.0004860768 0.0004943154 0.0006549679 0.0008238590 0.0008238590 
    ##           5%           5%           5%           5%           5%           5% 
    ## 0.0009062449 0.0009021256 0.0009845115 0.0010710166 0.0009062449 0.0009845115 
    ##           5%           5%           5%           5%           5%           5% 
    ## 0.0009021256 0.0008238590 0.0008197397 0.0007414731 0.0007373538 0.0007414731 
    ##           5%           5%           5%           5%           5%           5% 
    ## 0.0009886307 0.0010668973 0.0012316691 0.0011534025 0.0012357884 0.0010710166 
    ##           5%           5%           5%           5%           5%           5% 
    ## 0.0010710166 0.0014005602 0.0014829461 0.0017301038 0.0018907563 0.0019772615 
    ##           5%           5%           5%           5%           5%           5% 
    ## 0.0018948756 0.0017301038 0.0017301038 0.0017301038 0.0014005602 0.0013140550 
    ##           5%           5%           5%           5%           5%           5% 
    ## 0.0010668973 0.0010710166 0.0009845115 0.0009062449 0.0008238590 0.0006590872 
    ##           5%           5%           5%           5%           5%           5% 
    ## 0.0006549679 0.0005767013 0.0005725820 0.0004119295 0.0004119295 0.0004119295 
    ##           5%           5% 
    ## 0.0004119295 0.0004119295 
    ## 
    ## $ub_ht04
    ##         95%         95%         95%         95%         95%         95% 
    ## 0.001075136 0.001153403 0.001400560 0.001404680 0.001730104 0.001730104 
    ##         95%         95%         95%         95%         95%         95% 
    ## 0.001894876 0.002059647 0.002142033 0.002310924 0.002306805 0.002232658 
    ##         95%         95%         95%         95%         95%         95% 
    ## 0.002150272 0.002224419 0.002063767 0.001894876 0.001812490 0.001894876 
    ##         95%         95%         95%         95%         95%         95% 
    ## 0.002224419 0.002306805 0.002553963 0.002389191 0.002553963 0.002471577 
    ##         95%         95%         95%         95%         95%         95% 
    ## 0.002306805 0.002640468 0.002801120 0.003052397 0.003138903 0.003550832 
    ##         95%         95%         95%         95%         95%         95% 
    ## 0.003381941 0.003130664 0.003134783 0.003052397 0.002640468 0.002475696 
    ##         95%         95%         95%         95%         95%         95% 
    ## 0.002306805 0.002389191 0.002146153 0.001981381 0.001734223 0.001569451 
    ##         95%         95%         95%         95%         95%         95% 
    ## 0.001487065 0.001482946 0.001318174 0.001153403 0.001153403 0.001153403 
    ##         95%         95% 
    ## 0.001153403 0.001153403 
    ## 
    ## $med_wt13
    ##          50%          50%          50%          50%          50%          50% 
    ## 0.0009062449 0.0018124897 0.0022244192 0.0022656121 0.0028011204 0.0028835063 
    ##          50%          50%          50%          50%          50%          50% 
    ## 0.0030482781 0.0029658922 0.0028835063 0.0032130499 0.0036249794 0.0037073653 
    ##          50%          50%          50%          50%          50%          50% 
    ## 0.0037897512 0.0034602076 0.0037073653 0.0034602076 0.0036249794 0.0035425935 
    ##          50%          50%          50%          50%          50%          50% 
    ## 0.0034602076 0.0036661724 0.0037897512 0.0039545230 0.0039545230 0.0038309441 
    ##          50%          50%          50%          50%          50%          50% 
    ## 0.0040369089 0.0037897512 0.0037073653 0.0036249794 0.0038721371 0.0040369089 
    ##          50%          50%          50%          50%          50%          50% 
    ## 0.0039133300 0.0037073653 0.0035425935 0.0033778217 0.0031306640 0.0026363487 
    ##          50%          50%          50%          50%          50%          50% 
    ## 0.0023891910 0.0023068051 0.0018948756 0.0019772615 0.0017712968 0.0014829461 
    ##          50%          50%          50%          50%          50%          50% 
    ## 0.0014005602 0.0011534025 0.0009886307 0.0007002801 0.0007002801 0.0007002801 
    ##          50%          50% 
    ## 0.0007002801 0.0007002801 
    ## 
    ## $lb_wt13
    ##           5%           5%           5%           5%           5%           5% 
    ## 0.0004943154 0.0012357884 0.0017218652 0.0017259845 0.0021379140 0.0022202999 
    ##           5%           5%           5%           5%           5%           5% 
    ## 0.0023068051 0.0023850717 0.0021420333 0.0024674576 0.0026363487 0.0027187346 
    ##           5%           5%           5%           5%           5%           5% 
    ## 0.0027187346 0.0024715769 0.0026363487 0.0027146153 0.0028835063 0.0028793870 
    ##           5%           5%           5%           5%           5%           5% 
    ## 0.0028011204 0.0028835063 0.0030441588 0.0032954358 0.0031265447 0.0030482781 
    ##           5%           5%           5%           5%           5%           5% 
    ## 0.0030482781 0.0029617729 0.0028835063 0.0028835063 0.0030482781 0.0032089306 
    ##           5%           5%           5%           5%           5%           5% 
    ## 0.0028835063 0.0027146153 0.0027970012 0.0026363487 0.0023850717 0.0018907563 
    ##           5%           5%           5%           5%           5%           5% 
    ## 0.0016477179 0.0016477179 0.0013181743 0.0013181743 0.0011534025 0.0010668973 
    ##           5%           5%           5%           5%           5%           5% 
    ## 0.0009021256 0.0006590872 0.0006590872 0.0004119295 0.0004119295 0.0004119295 
    ##           5%           5% 
    ## 0.0004119295 0.0004119295 
    ## 
    ## $ub_wt13
    ##         95%         95%         95%         95%         95%         95% 
    ## 0.001400560 0.002718735 0.002887626 0.003048278 0.003464327 0.003789751 
    ##         95%         95%         95%         95%         95%         95% 
    ## 0.003872137 0.003872137 0.003629099 0.004119295 0.004448838 0.004613610 
    ##         95%         95%         95%         95%         95%         95% 
    ## 0.004695996 0.004370572 0.004704235 0.004535344 0.004778382 0.004613610 
    ##         95%         95%         95%         95%         95%         95% 
    ## 0.004535344 0.004539463 0.004700115 0.004947273 0.004860768 0.004864887 
    ##         95%         95%         95%         95%         95%         95% 
    ## 0.005272697 0.004790740 0.004786621 0.004695996 0.004943154 0.005029659 
    ##         95%         95%         95%         95%         95%         95% 
    ## 0.004860768 0.004695996 0.004452958 0.004288186 0.003958642 0.003460208 
    ##         95%         95%         95%         95%         95%         95% 
    ## 0.003213050 0.003048278 0.002718735 0.002640468 0.002471577 0.002063767 
    ##         95%         95%         95%         95%         95%         95% 
    ## 0.001820728 0.001647718 0.001565332 0.001153403 0.001153403 0.001153403 
    ##         95%         95% 
    ## 0.001153403 0.001153403

![](README_files/figure-gfm/unnamed-chunk-54-3.png)<!-- -->

## 5. Trivariate analysis

The package contains three higher dimensional $(>3)$ approaches are
implemented to model the joint distribution of rainfall, O-sWL and
groundwater level. They are:

- Standard (trivariate) copula
- Pair Copula Construction
- Heffernan and Tawn (2004)

#### Standard (trivariate) copula

In the package, each approach has a `_Fit` and `_Sim` function. The
latter requires a `MIGPD` object as its `Marginals` input argument, in
order for the simulations on $\[0,1\]^3$ to be transformed back to the
original scale. The `Migpd_Fit` command fits independent GPDs to the
data in each row of a dataframe (excluding the first column if it is a
“Date” object) creating a `MIGPD` object.

``` r
S20.Migpd<-Migpd_Fit(Data=S20.Detrend.Declustered.df[,-1],Data_Full = S20.Detrend.df[,-1],
                     mqu=c(0.975,0.975,0.9676))
summary(S20.Migpd)
```

    ## 
    ## A collection of 3 generalized Pareto models.
    ## All models converged.
    ## Penalty to the likelihood: gaussian
    ## 
    ## Summary of models:
    ##                   Rainfall      OsWL Groundwater
    ## Threshold        1.6000000 1.9385406   2.8599327
    ## P(X < threshold) 0.9750000 0.9750000   0.9676000
    ## sigma            0.9040271 0.1574806   0.3083846
    ## xi               0.1742220 0.2309118  -0.3441421
    ## Upper end point        Inf       Inf   3.7560295

Standard (trivariate) copula are the most conceptually simple of the
copula based models, using a single parametric multivariate probability
distribution as the copula. The `Standard_Copula_Fit()` function fits
elliptic (specified by `Gaussian` or `tcop`) or Archimedean (specified
by `Gumbel`,`Clayton` or `Frank`) copula to a trivariate dataset. Let
first fit a Gaussian copula

``` r
S20.Gaussian<-Standard_Copula_Fit(Data=S20.Detrend.df,Copula_Type="Gaussian")
```

From which the `Standard_Copula_Sim()` function can be used to simulate
a synthetic record of N years

``` r
S20.Gaussian.Sim<-Standard_Copula_Sim(Data=S20.Detrend.df,Marginals=S20.Migpd,
                                      Copula=S20.Gaussian,N=100)
```

Plotting the observed and simulated values

``` r
S20.Pairs.Plot.Data<-data.frame(rbind(na.omit(S20.Detrend.df[,-1]),S20.Gaussian.Sim$x.Sim),
                                c(rep("Observation",nrow(na.omit(S20.Detrend.df))),
                                  rep("Simulation",nrow(S20.Gaussian.Sim$x.Sim))))
colnames(S20.Pairs.Plot.Data)<-c(names(S20.Detrend.df)[-1],"Type")
pairs(S20.Pairs.Plot.Data[,1:3],
      col=ifelse(S20.Pairs.Plot.Data$Type=="Observation","Black",alpha("Red",0.3)),
      upper.panel=NULL,pch=16)
```

![](README_files/figure-gfm/unnamed-chunk-58-1.png)<!-- -->

The `Standard_Copula_Sel()` function can be used to deduce the best
fitting in terms of AIC

``` r
Standard_Copula_Sel(Data=S20.Detrend.df)
```

    ##     Copula        AIC
    ## 1 Gaussian -1389.1438
    ## 2    t-cop -1434.7850
    ## 3   Gumbel  -876.7537
    ## 4  Clayton  -565.7595
    ## 5    Frank  -854.4284

#### Pair Copula Construction

Standard trivariate copulas lack flexibility to model joint
distributions where heterogeneous dependencies exist between the
variable pairs. Pair copula constructions construct multivariate
distribution using a cascade of bivariate copulas (some of which are
conditional). As the dimensionality of the problem increases the number
of mathematically equally valid decompositions quickly becomes large.
Bedford and Cooke (2001,2002) introduced the regular vine, a graphical
model which helps to organize the possible decompositions. The Canonical
(C-) and D- vine are two commonly utilized sub-categories of regular
vines, in the trivariate case a vine copula is simultaneously a C- and
D-vine. Let’s fit a regular vine copula model

``` r
S20.Vine<-Vine_Copula_Fit(Data=S20.Detrend.df)
```

From which the `Vine_Copula_Sim()` function can be used to simulate a
synthetic record of N years

``` r
S20.Vine.Sim<-Vine_Copula_Sim(Data=S20.Detrend.df,Vine_Model=S20.Vine,Marginals=S20.Migpd,N=100)
```

Plotting the observed and simulated values

``` r
S20.Pairs.Plot.Data<-data.frame(rbind(na.omit(S20.Detrend.df[,-1]),S20.Vine.Sim$x.Sim),
                                c(rep("Observation",nrow(na.omit(S20.Detrend.df))),
                                  rep("Simulation",nrow(S20.Vine.Sim$x.Sim))))
colnames(S20.Pairs.Plot.Data)<-c(names(S20.Detrend.df)[-1],"Type")
pairs(S20.Pairs.Plot.Data[,1:3],
      col=ifelse(S20.Pairs.Plot.Data$Type=="Observation","Black",alpha("Red",0.3)),
      upper.panel=NULL,pch=16)
```

![](README_files/figure-gfm/unnamed-chunk-62-1.png)<!-- -->

#### HT04

Finally, let us implement the Heffernan and Tawn (2004) approach, in
which a non-linear regression model is fitted to joint observations
where a conditioning variable exceeds a specified threshold. The
regression model typically adopted is:
$\textbf{Y}_{-i} = \textbf{a}Y_i + Y_i^{\textbf{b}}\textbf{Z} \hspace{1cm} \text{for} \hspace{1cm} Y_i > v$
where:

$\textbf{Y}$ is a set of variables transformed to a common scale
$\textbf{Y}_{-i}$ is the set of variables excluding $Y_i$ $\textbf{a}$
and $\textbf{b}$ are vectors of regression parameters $\textbf{Z}$ is a
vector of residuals

The dependence structure when a specified variable is extreme is thus
captured by the regression parameters and the joint residuals. The
procedure is repeated, conditioning on each variable in turn, to build
up the joint distribution when at least one variable is in an extreme
state. The HT04 command fits the model and simulates N years’ worth of
data from the fitted model.

``` r
S20.HT04<-HT04(data_Detrend_Dependence_df=S20.Detrend.df,
               data_Detrend_Declustered_df=S20.Detrend.Declustered.df,
               u_Dependence=0.995,Migpd=S20.Migpd,mu=365.25,N=1000)
```

Output of the function includes the three conditional `Models`,
proportion of occasions where each variable is most extreme given at
least one variable is extreme `prop`as well as, the simulations on the
transformed scale `u.Sim` (Gumbel by default) and original scale
`x.Sim`. Let’s view the fitted model when conditioning on rainfall

``` r
S20.HT04$Model$Rainfall
```

    ## mexDependence(x = Migpd, which = colnames(data_Detrend_Dependence_df)[i], 
    ##     dqu = u_Dependence, margins = "laplace", constrain = FALSE, 
    ##     v = V, maxit = Maxit)
    ## 
    ## 
    ## Marginal models:
    ## 
    ## Dependence model:
    ## 
    ## Conditioning on Rainfall variable.
    ## Thresholding quantiles for transformed data: dqu = 0.995
    ## Using laplace margins for dependence estimation.
    ## Log-likelihood = -108.5807 -84.53349 
    ## 
    ## Dependence structure parameter estimates:
    ##     OsWL Groundwater
    ## a 1.0000      0.3295
    ## b 0.6948     -1.5240

``` r
S20.HT04$Prop
```

    ## [1] 0.3552632 0.3157895 0.3289474

and the which the proporiton of the occasions in the original sample
that rainfall is the most extreme of the drivers given that at least one
driver is extreme.

The HT04 approach uses rejection sampling to generate synthetic records.
The first step involves sampling a variable, conditioned to exceed the
`u_Dependence` threshold. A joint residual associated with the
corresponding regression is independently sampled and other variables
estimated using the fitted regression parameters. If the variable
conditioned to be extreme in step one is not the most extreme the sample
is rejected. The process is repeated until the relative proportion of
simulated events where each variable is a maximum, conditional on being
above the threshold, is consistent with the empirical distribution.
Labeling the simulations `S20.HT04.Sim`

``` r
S20.HT04.Sim<-S20.HT04$x.sim
```

and now plotting the simulations from the HT04 model

``` r
S20.Pairs.Plot.Data<-data.frame(rbind(na.omit(S20.Detrend.df[,-1]),S20.HT04.Sim),
                                c(rep("Observation",nrow(na.omit(S20.Detrend.df))),
                                  rep("Simulation",nrow(S20.HT04.Sim))))
colnames(S20.Pairs.Plot.Data)<-c(names(S20.Detrend.df)[-1],"Type")
pairs(S20.Pairs.Plot.Data[,1:3],
      col=ifelse(S20.Pairs.Plot.Data$Type=="Observation","Black",alpha("Red",0.2)),
      upper.panel=NULL,pch=16)
```

![](README_files/figure-gfm/unnamed-chunk-67-1.png)<!-- -->

## 6. Sea Level Rise

The `SLR_Scenarios` function estimates the time required for a
user-specified amount of sea level rise ( `SeaLevelRise`) to occur under
various sea level rise scenarios. The default scenarios are for Key West
from the Southeast Florida Regional Climate Change Compact (2019). Let’s
calculate how long before the O-sWL in the 100-year “most-likely” design
event (see section 4) equals that of the corresponding design event
derived under full dependence.

``` r
#Difference in O-sWL between the most-likely and full dependence events
Diff<-S20.Bivariate$FullDependence$`100`$OsWL-S20.Bivariate$MostLikelyEvent$`100`$OsWL
Diff
```

    ## [1] 0.1360731

``` r
#Time in years for the sea level rise to occur
SLR_Scenarios(SeaLevelRise=Diff,Unit="m")
```

![](README_files/figure-gfm/unnamed-chunk-69-1.png)<!-- -->

    ## $High
    ## [1] 2032
    ## 
    ## $Intermediate
    ## [1] 2034
    ## 
    ## $Low
    ## [1] 2046

Scenarios from the Interagency Sea Level Rise Scenario Tool (2022) for
Miami Beach and Naples can be utilized by changing the `Scenario` and
`Location` arguments. Alternatively, a user can input other sea level
rise scenarios into the function. For example, below we use the
scenarios from the same tool but for Fort Myers.

``` r
SeaLevelRise.2022<-read.csv("F:/OneDrive - University of Central Florida/Documents/sl_taskforce_scenarios_psmsl_id_1106_Fort_Myers.csv")
SeaLevelRise.2022_input<-data.frame(Year=seq(2020,2150,10),
"High"=as.numeric(SeaLevelRise.2022[14,-(1:5)])/1000,
"Int_Medium"=as.numeric(SeaLevelRise.2022[11,-(1:5)])/1000,
"Medium"=as.numeric(SeaLevelRise.2022[8,-(1:5)])/1000,
"Int_Low"=as.numeric(SeaLevelRise.2022[5,-(1:5)])/1000,
"Low"=as.numeric(SeaLevelRise.2022[2,-(1:5)])/1000)
SLR_Scenarios(SeaLevelRise=0.8, Scenario="Other", Unit = "m", Year=2022, 
              Location="Fort Myers", New_Scenario=SeaLevelRise.2022_input)
```

![](README_files/figure-gfm/unnamed-chunk-70-1.png)<!-- -->

    ## $SLR_Year
    ## [1] 2150 2120 2085 2070 2064

## 7. Simulating temporally varying events

By sampling peaks from a multivariate statistical model e.g., a vine
copula that includes the relative lag between the rainfall and water
level peaks as well as their magnitudes temporally-varying synthetic
events can be simulated.

#### Water level curves

The intensity of a storm surge event is defined in [Wahl et
al. 2011](https://doi.org/10.5194/nhess-11-2925-2011) as the area of the
water level curve above a specified base line level from the first low
water level of the preceding tide and last low tide of the following
tide. The `intensity()` function calculates this metric with the default
baseline level set as the mean of the water level time series
(`Base_Line`=`"Mean "`). Let’s calculate the “intensity” of cluster
maxima in the hourly O-sWL time series at control structure S-13.

``` r
#Decluster O-sWL series at S-13 using a runs method
S13.OsWL.Declust = Decluster(Data=S13.Detrend.df$OsWL,
                            SepCrit=24*7, u=0.99667)
#Calculate O-sWL of the identified cluster maximum
intensity = Intensity(Data=S13.Detrend.df[,c(1,3)],
                      Cluster_Max=S13.OsWL.Declust$EventsMax,
                      Base_Line=2)
```

Plotting a subset of the events:

``` r
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
```

![](README_files/figure-gfm/unnamed-chunk-72-1.png)<!-- -->

The `WL_Curve()` function generates a water level curve conditioned on a
water level event peak (`Peak`) and “intensity” (`Intensity`). The
function works by re-scaling observed events curves between four time
units before and after the peak such that the re-scaled peak coincides
with conditioned peak. It then computes the intensity of the re-scaled
events and selects the event with the highest intensity that’s less than
conditioned intensity. Intensity units are then added to ensure the
intensity of the simulated curve matches the conditioned intensity. The
Intensity units, which manifest as increases in water level, are
distributed in proportion with the decreases of water levels from the
peak OsWL i.e., more units are added at lower water levels. No units are
added about the peak. For conditioned event peaks below a user-specified
threshold (`Thres`) an observed curve with an intensity less than some
user-specified limit (`Limit`) is randomly sampled.

``` r
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
```

Superimposing the four simulated water level curves on the observed
curves:

``` r
#Plot the water level curves of the observed peaks
plot(1:289,
     S13.Detrend.df$OsWL[(S13.OsWL.Declust$EventsMax[1]-144):
                         (S13.OsWL.Declust$EventsMax[1]+144)],
     type='l',ylim=c(-2,6))
for(i in 2:length(S13.OsWL.Declust$EventsMax)){
  lines(1:289,
        S13.Detrend.df$OsWL[(S13.OsWL.Declust$EventsMax[i]-144):
                            (S13.OsWL.Declust$EventsMax[i]+144)])
}
#Superimpose the curves generated for the four synthetic events
for(i in 1:4){
  lines(1:289,oswl_ts_oswl$Series[i,],col=2)
}
```

![](README_files/figure-gfm/unnamed-chunk-74-1.png)<!-- -->

#### Hyetographs

[Serinaldi and Kilsby](https://doi.org/10.1002/wrcr.20221) studied the
pairwise relationships between key hyetograph characteristics: maximum
value $Xp$, volume $V$, duration $D$, and average intensity $I$ at $282$
rain gauge stations in central eastern Europe with daily rainfall
records and three 5-minutes record gauges in Italy. The authors found
that at both temporal scales the only stochastic relationship was
between volume and duration, and that simple bootstrap procedures could
be used to generate events that preserve the pairwise relationships
among the characteristics. The `U_Sample` function implements a
bootstrap procedure referred to as the U-boot algorithm in [Serinaldi
and Kilsby (2013)](https://doi.org/10.1002/wrcr.20221) to generate
hyetograhs. For a simulated peak `Xp`, a duration is independently
sampled and a set of non-peaks are sampled at random from one of the
events with the same duration. To implement the method exactly as in
[Serinaldi and Kilsby (2013)](https://doi.org/10.1002/wrcr.20221), set
`Xp` equal to a sample (taken with replacement) of the observed cluster
maximum (peaks). An example application of the function is given below
for rainfall recorded at control structure S-13. First the hourly time
series is declustered to find the characteristics of the rainfall events
containing the $500$ highest peaks:

``` r
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
```

Now the bootstrapping procedure can be carried out:

``` r
sample = U_Sample(Data=S13.Detrend.df$Rainfall,
                  Cluster_Max=S13.Rainfall.Declust$EventsMax,
                  D=d,Start=start,End=end,
                  Xp=sim.peak)
```

Let’s check whether the characteristic of the bootstrapped events match
those of the observed events. First, calculating the volume and
intensity of the $500$ observed events with the highest hourly peaks.

``` r
#Calculating volume and intensity
v<-rep(NA,500)
for(i in 1:length(S13.Rainfall.Declust$EventsMax)){
  v[i] = sum(S13.Detrend.df$Rainfall[(start[i]):(end[i])])
}
I = v/d

#Putting in a data.frame
observations = data.frame(peaks,d,v,I)
colnames(observations) = c("peak","d","v","I")
```

To aid the comparison information relating to the marginal distributions
of the characteristics is removed by transforming the data to the
$[0,1]^2$ scale.

``` r
#Observations
observations.u = data.frame(pobs(observations))
colnames(observations.u) = c("peak","d","v","I")

#Sample
sample.u = data.frame(pobs(sample))
```

Scatter plots showing the pairwise relationships of the characteristics:

``` r
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
```

![](README_files/figure-gfm/unnamed-chunk-79-1.png)<!-- -->

The figure demonstrates that the simulated events preserves the key
characteristics of the observed events.
