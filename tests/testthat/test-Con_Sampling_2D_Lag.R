#Test that functions work as intended

test_that("Conditional Sampling 2Ds functions basic functionality", {

  result<-Con_Sampling_2D_Lag(Data_Detrend=S20.Detrend.df[,-c(1,4)],
                              Data_Declust=S20.Detrend.Declustered.df[,-c(1,4)],
                              Con_Variable="Rainfall",u=0.97,
                              Lag_Backward = 3, Lag_Forward = 3)

  # Check return type
  expect_type(result, "list")

  #Check names of output
  expect_named(result, c("Threshold", "Data", "Con_Variable", "x.con", "x.noncon"))

  # Check length of outputs
  expect_length(result$Threshold, 1)
  expect_length(result$Con_Variable, 1)
  expect_equal(result$Con_Variable, "Rainfall")
  expect_equal(nrow(result$Data), length(result$x.con))
  expect_equal(nrow(result$Data), length(result$x.noncon))
  expect_lt(max(result$x.con - result$x.noncon), 4)
  expect_gt(min(result$x.con - result$x.noncon), -4)
})

test_that("Con_Sampling_2D function with different(quantile) threshold values", {

  result1 <- Con_Sampling_2D_Lag(Data_Detrend=S20.Detrend.df[,-c(1,4)],
                             Data_Declust=S20.Detrend.Declustered.df[,-c(1,4)],
                             Con_Variable="Rainfall",u=0.97)
  result2 <- Con_Sampling_2D_Lag(Data_Detrend=S20.Detrend.df[,-c(1,4)],
                                 Data_Declust=S20.Detrend.Declustered.df[,-c(1,4)],
                                 Con_Variable="Rainfall",u=0.9)
  expect_true(result1$Threshold >= result2$Threshold)
  expect_true(length(result1$x.con) <= length(result2$x.con))
  expect_true(nrow(result1$Data) <= nrow(result2$Data))
})


#Test reproducibility

test_that("Function is deterministic", {

  result1 <- Con_Sampling_2D_Lag(Data_Detrend=S20.Detrend.df[,-c(1,4)],
                                 Data_Declust=S20.Detrend.Declustered.df[,-c(1,4)],
                                 Con_Variable="Rainfall",u=0.97)
  result2 <- Con_Sampling_2D_Lag(Data_Detrend=S20.Detrend.df[,-c(1,4)],
                                 Data_Declust=S20.Detrend.Declustered.df[,-c(1,4)],
                                 Con_Variable="Rainfall",u=0.97)
  expect_identical(result1, result2)
})

# Test that invalid inputs gives errors
test_that("Invalid inputs produce errors", {
  expect_error(Con_Sampling_2D_Lag(Data_Detrend=S20.Detrend.df[,-c(1,4)],
                                   Data_Declust=S20.Detrend.Declustered.df[,-c(1,4)],
                                   Con_Variable="Rainfall",u=-0.1))

  expect_error(Con_Sampling_2D_Lag(Data_Detrend=S20.Detrend.df[,-c(1,4)],
                                   Data_Declust=S20.Detrend.Declustered.df[,-c(1,4)],
                                   Con_Variable="Rainfall",u="0.97"))

  expect_error(Con_Sampling_2D_Lag(Data_Detrend=S20.Detrend.df[,-c(1,4)],
                                   Data_Declust=S20.Detrend.Declustered.df[,-c(1,4)],
                                   Con_Variable="Rainfall",u=1.5,
                                   Lag_Backward = 3, Lag_Forward = 3))

  expect_error(Con_Sampling_2D_Lag(Data_Detrend=S20.Detrend.df[,4],
                                   Data_Declust=S20.Detrend.Declustered.df[,-c(1,4)],
                                   Con_Variable="Rainfall",u=0.97,
                                   , Lag_Backward = 3, Lag_Forward = 3))

  expect_error(Con_Sampling_2D_Lag(Data_Detrend=S20.Detrend.df[,-c(1,4)],
                                   Data_Declust=S20.Detrend.Declustered.df[,-c(1,4)],
                                   Con_Variable="Invalid",u=0.97,
                                   Lag_Backward = 3, Lag_Forward = 3))

  expect_error(Con_Sampling_2D_Lag(Data_Detrend=S20.Detrend.df[,-c(1,4)],
                                   Data_Declust=S20.Detrend.Declustered.df[,-c(1,4)],
                                   Con_Variable="Invalid",u=0.97,
                                   Lag_Backward = 3, Lag_Forward = 3))

  expect_error(Con_Sampling_2D_Lag(Data_Detrend=S20.Detrend.df[,-c(1,4)],
                                   Data_Declust=S20.Detrend.Declustered.df[,-c(1,4)],
                                   Con_Variable="Invalid",u=0.97,
                                   Lag_Backward = 3, Lag_Forward = "Invalid"))
})
