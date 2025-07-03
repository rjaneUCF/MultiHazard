#Test that functions work as intended

test_that("Conditional Sampling 2Ds functions basic functionality", {

  result<-Con_Sampling_2D(Data_Detrend=S20.Detrend.df[,-c(1,4)],
                                Data_Declust=S20.Detrend.Declustered.df[,-c(1,4)],
                                Con_Variable="Rainfall",u=0.97)
  
  # Check return type
  expect_type(result, "list")
  
  #Check names of output
  expect_named(result, c("Threshold", "Data", "Rate", "Con_Variable", "x.con"))
  
  # Check length of outputs
  expect_length(result$Threshold, 1)
  expect_length(result$Rate, 1)
  expect_length(result$Con_Variable, 1)
  expect_equal(nrow(result$Data), length(x.con))
  
})
