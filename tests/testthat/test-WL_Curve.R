
S13.OsWL.Declust = Decluster(Data=S13.Detrend.df$OsWL,
                              SepCrit=24*7, u=0.99667)
#Use O-sWL intensity function to obtain index of preceding and following low water levels
intensity = OsWL_Intensity(Data=S13.Detrend.df,Cluster_Max=S13.OsWL.Declust$EventsMax)

#Four synthetic events
sim.peaks = c(3.4,4,4.2,5)
sim.intensity = c(38,48,120,140)


test_that("WL_Curve function works", {
  
  result <- WL_Curve(Data = S13.Detrend.df,
                     Cluster_Max = S13.OsWL.Declust$EventsMax,
                     Pre_Low = intensity$Pre.Low,
                     Fol_Low = intensity$Fol.Low,
                     Thres = S13.OsWL.Declust$Threshold, Limit = 45,
                     Peak = sim.peaks,
                     Intensity = sim.intensity)
  
  # Checking type of output
  expect_type(result, 'list')
  expect_named(result, c("Series","Intensity","Event"))
  
  #Checking length of outputs
  expect_equal(nrow(result$Series),length(sim.intensity))
  expect_equal(ncol(result$Series),2*144+1)
  expect_equal(length(result$Intensity),length(sim.intensity))
  expect_equal(length(result$Event),length(sim.intensity))
  
  #Checking column names of outputs
  expect_true(all(result$Event > 0 & result$Event < length(S13.OsWL.Declust$EventsMax)))
  expect_gte(all(as.numeric(result$Intensity)), 0)
})


#Checking invalid inputs

test_that("Intensity invalid inputs", {
  
  S13.Detrend.df.new.name = S13.Detrend.df
  colnames(S13.Detrend.df.new.name) = c("Date","Rainfall", "Water_level")
  expect_error(WL_Curve(Data = S13.Detrend.df.new.name,
                        Cluster_Max = S13.OsWL.Declust$EventsMax,
                        Pre_Low = intensity$Pre.Low,
                        Fol_Low = intensity$Fol.Low,
                        Thres = S13.OsWL.Declust$Threshold, Limit = 45,
                        Peak = sim.peaks,
                        Intensity = sim.intensity),
               "Data must contain 'OsWL' column.")
  
  expect_error(WL_Curve(Data = S13.Detrend.df,
                        Cluster_Max = S13.OsWL.Declust$EventsMax,
                        Pre_Low = intensity$Pre.Low,
                        Fol_Low = intensity$Fol.Low,
                        Thres = S13.OsWL.Declust$Threshold, Limit = 45,
                        Peak = sim.peaks[-1],
                        Intensity = sim.intensity),
               "Peak and Intensity vectors must be of equal length.")
  
  expect_error(WL_Curve(Data = S13.Detrend.df,
                        Cluster_Max = S13.OsWL.Declust$EventsMax,
                        Pre_Low = intensity$Pre.Low,
                        Fol_Low = intensity$Fol.Low,
                        Thres = S13.OsWL.Declust$Threshold, Limit = 45,
                        Peak = sim.peaks,
                        Intensity = sim.intensity,
                        Base_Line=c(6,7,8)),
               "Base_Line must be a single value")
  
  dummy_index <- rep(1, length(S13.OsWL.Declust$EventsMax))  # all invalid
  expect_error(WL_Curve(Data = S13.Detrend.df,
                        Cluster_Max = dummy_index,
                        Pre_Low = dummy_index,
                        Fol_Low = dummy_index,
                        Thres = S13.OsWL.Declust$Threshold, Limit = 45,
                        Peak = sim.peaks,
                        Intensity = sim.intensity),
               "No Cluster_Max values are within valid bounds")
})
