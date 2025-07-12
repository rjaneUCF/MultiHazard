
S13.OsWL.Declust = Decluster(Data=S13.Detrend.df$OsWL,
                            SepCrit=24*7, u=0.99667)

test_that("Intensity function works", {

  result <- Intensity(Data=S13.Detrend.df,Cluster_Max=S13.OsWL.Declust$EventsMax)

  # Checking type of output
  expect_type(result, 'data.frame')
  expect_named(result, c("Pre.High","Fol.High","Pre.Low","Fol.Low","Intensity"))

  #Checking length of outputs
  expect_equal(length(result$Pre.High),length(S13.OsWL.Declust$EventsMax))
  expect_equal(length(result$Fol.High),length(S13.OsWL.Declust$EventsMax))
  expect_equal(length(result$Pre.Low),length(S13.OsWL.Declust$EventsMax))
  expect_equal(length(result$Fol.Low),length(S13.OsWL.Declust$EventsMax))

  #Checking value of outputs
  expect_true(all(result$Pre.High >= 0))
  expect_true(all(result$Fol.High >= 0))
  expect_true(all(result$Pre.Low >= 0))
  expect_true(all(result$Fol.Low >= 0))
  expect_true(all(result$Intensity >= 0))
})


#Checking invalid inputs

test_that("Intensity invalid inputs", {

  expect_error(Intensity(Data=S13.Detrend.df,Cluster_Max=S13.OsWL.Declust$EventsMax),
               "Data must have at least one row and one column")

  expect_error(Intensity(Data=S13.Detrend.df,Cluster_Max=S13.OsWL.Declust$EventsMax),
               "Cluster_Max must only contain positive integer values that do not exceed length of time series")

  expect_error(Intensity(Data=S13.Detrend.df,Cluster_Max=S13.OsWL.Declust$EventsMax),
               "Base_Line must be a single value")
})
