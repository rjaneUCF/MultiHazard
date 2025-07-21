#First decluster the rainfall series to find the 500 events with the highest peaks
S13.Rainfall.Declust = Decluster(Data=S13.Detrend.df$Rainfall,
                                 SepCrit=24*3, u=0.99667)
#Set very small rainfall measurements to zero.
#Assumed to be the result of uncertainty in measuring equipment.
S13_Rainfall$Rainfall[which(S13_Rainfall$Rainfall<0.01)] = 0
#Find NAs in rainfall series
z = which(is.na(S13_Rainfall$Rainfall)==T)
#Temporarily set NAs to zero
S13_Rainfall$Rainfall[z] = 0
#Find times where there is 6-hours of no rainfall
no.rain = rep(NA,length(S13_Rainfall$Rainfall))
for(i in 6:length(S13_Rainfall$Rainfall)){
   no.rain[i] = ifelse(sum(S13_Rainfall$Rainfall[(i-5):i])==0,i,NA)
 }
#Remove NAs from results vector as these correspond to times where there is
#rainfall at certain points in the 6 hour period.
no.rain = na.omit(no.rain)
#Reset missing values in the rainfall record back to NA
S13_Rainfall$Rainfall[z] = NA
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
sim.peak = sample(S13.Rainfall.Declust$EventsMax,size=500,replace=TRUE)
#Derive the hyetographs


test_that("U_Sample works", {

 result = U_Sample(Data=S13_Rainfall$Rainfall,
                   Cluster_Max=S13.Rainfall.Declust$EventsMax,
                   D=d,Start=start,End=end,
                   Xp=sim.peak)

 # Checking type of output
 expect_type(result, 'list')
 expect_named(result, c("Xp","D","Samp","V","Vn","I","In","Start","End"))
 expect_type(result$Xp, "double")
 expect_type(result$D, "integer")
 expect_type(result$Samp, "integer")
 expect_equal(nrow(result), length(sim.peak))

 #  Check for no unexpected NA values in key columns
 expect_false(all(is.na(result$Samp)))
 expect_false(all(is.na(result$V)))

 # Check logical constraints
 expect_true(all(result$Start <= result$End, na.rm = TRUE))
 expect_true(all(result$D > 0, na.rm = TRUE))

})



test_that("Invalid inputs", {

  expect_error(
    U_Sample(Data=S13_Rainfall$Rainfall,
             Cluster_Max=S13.Rainfall.Declust$EventsMax,
             D=d[-1],Start=start,End=end,
             Xp=sim.peak),
    "Cluster_Max, D, Start, and End must have the same length.")

  expect_error(
    U_Sample(Data=S13_Rainfall$Rainfall,
             Cluster_Max=S13.Rainfall.Declust$EventsMax,
             D=d,Start=start,End=end),
    "Xp must contain at least one value.")
})
