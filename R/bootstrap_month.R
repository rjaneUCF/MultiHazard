#' Implements a monthly bootstrap 
#' 
#' Months in which at least one variable exceeds the user-specified minimum proportion of non-missing values are sampled with replacement. February of leap years are treated as a 13th month. 
#'
#' @param data Data frame of raw data detrended if necessary. First column should be of the \code{Date}.
#' @param block_prop Numeric vector of length one specifying the minimum proportion of non-missing values of at least one of the variables for a month to be included in the bootstrap. Default is \code{0.8}.
#' @return Dataframe containing a bootstrap undertaken with replacement that accounts for monthly-scale seasonality.
#' @export
#' @examples
#' #Let's assess the sampling variability in kendall's tau
#' #correlation coefficient between rainfall and OsWL at S-22.
#' 
#' #Data starts on first day of 1948
#' head(S22.Detrend.Declustered.df)
#' 
#' #Dataframe ends on 1948-02-03 
#' tail(S22.Detrend.Declustered.df)
#'
#' #Adding dates to complete final month of combined records
#' final.month = data.frame(seq(as.Date("2019-02-04"),as.Date("2019-02-28"),by="day"),NA,NA,NA)
#' colnames(final.month) = c("Date","Rainfall","OsWL","Groundwater")
#' S22.Detrend.Declustered.df = rbind(S22.Detrend.Declustered.df,final.month)
#'
#' #Generate 100 monthly bootstrap samples of rainfall and OsWL
#' cor = rep(NA,100)
#' for(i in 1:100){
#'  boot_df = bootstrap_month(S22.Detrend.df[,c(1:3)], boot_prop=0.8)
#'  boot_df = na.omit(boot_df)
#'  cor[i] = cor(boot_df$Rainfall, boot_df$OsWL, method="kendall")
#' }
#'
#' #Compare means of bootstrap samples with the mean of the observed data 
#' hist(cor)
#' df = na.omit(S22.Detrend.df[,1:3])
#' abline(v=cor(df$Rainfall,df$OsWL, method="kendall"),col=2,lwd=2)
bootstrap_month = function(data,boot_prop=0.8){

#Finding month and years represented in the observational records
m.y = data.frame(month(data$Date),year(data$Date))
colnames(m.y) = c("month","year")
  
#Continuous records of alld months and years
m.y.complete = data.frame(rep(1:12,length(min(m.y$year):max(m.y$year))),rep(min(m.y$year):max(m.y$year),each=12))
  
#Indicator variable denoting which month has rainfall, tail water or both data types available
I = rep(NA, nrow(m.y.complete))
  
for(i in 1:nrow(m.y.complete)){
 d = data[which(month(data$Date) == m.y.complete[i,1] & year(data$Date) == m.y.complete[i,2]),]
 x = length(which(!is.na(d[,2])  & is.na(d[,3]))) / nrow(d)
 y = length(which(is.na(d[,2]) & !is.na(d[,3]))) / nrow(d)
 b = length(which(!is.na(d[,2]) & !is.na(d[,3]))) / nrow(d)
    
 I[i] = ifelse(b>boot_prop,"B",ifelse(x>boot_prop,"x",ifelse(y>boot_prop,"y",NA)))
   
}
  
#Combining information into a matrix
m.y.complete = data.frame(rep(1:12,length(min(m.y$year):max(m.y$year))),rep(min(m.y$year):max(m.y$year),each=12),I)
colnames(m.y.complete) = c("month","year","indicator")

#Treating February in leap years as a `13th month` 
leap.years = seq(1880,2052,4)
s = leap.years[min(which(leap.years>=min(m.y$year)))]
e = leap.years[max(which(leap.years<=max(m.y$year)))]
leap.years = seq(s,e,4)

#Identifying leap years
z1 = rep(NA, length(leap.years))

for(i in 1:length(leap.years)){ 
 z1[i] = ifelse(length(which(m.y.complete$month==2 & m.y.complete$year==leap.years[i]))>0,which(m.y.complete$month==2 & m.y.complete$year==leap.years[i]),NA)
}
  
m.y.complete$month[z1] = 13

#For bootstrap, do not pick months where no data is available
z = which(!is.na(m.y.complete$indicator))

#Bootstrap by month
boot.m.y = rep(NA, nrow(m.y.complete))

for(i in 1:nrow(m.y.complete[z,])){
  boot.m.y[z[i]] = sample(which(m.y.complete$month ==  m.y.complete$month[z[i]] &  m.y.complete$indicator == m.y.complete$indicator[z[i]]),1)
}

#Bootstrap sample of months and years
m.y.boot = m.y.complete[boot.m.y,1:2]

#Convert '13' back to '02'
m.y.boot$month[which(m.y.boot$month==13)] = 2
m.y.complete$month[which(m.y.complete$month==13)] = 2

#Composing the bootstrap sample 
boot_df = data.frame(data[,1],rep(NA,nrow(data)),rep(NA,nrow(data)))
colnames(boot_df) = colnames(data)

for(i in 1:nrow(m.y.complete[z,])){
  boot_df[which(month(data$Date) == m.y.complete[z[i],1] & year(data$Date) == m.y.complete[z[i],2]),] = data[which(month(data$Date) == m.y.boot[z[i],1] & year(data$Date) == m.y.boot[z[i],2]),]
}

boot_df

}


