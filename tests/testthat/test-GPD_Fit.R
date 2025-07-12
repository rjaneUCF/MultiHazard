
test_that("GPD_Fit basic functionality works", {

  result <- GPD_Fit(Data=S20.Detrend.Declustered.df$Rainfall,Data_Full=na.omit(S20.Detrend.df$Rainfall),
                    u=0.997,PLOT=FALSE, xlab_hist="Rainfall (Inches)",y_lab="Rainfall (Inches)")

  expect_type(result, 'list')
  expect_named(result,  c("Threshold","Rate","sigma","xi","sigma.SE","xi.SE"))

  expect_true(is.numeric(result$Threshold))
  expect_true(is.numeric(result$Rate))
  expect_true(is.numeric(result$sigma))
  expect_true(is.numeric(result$xi))
  expect_true(is.numeric(result$sigma.SE))
  expect_true(is.numeric(result$xi.SE))

  expect_gte(result$Rate,  0)
})


test_that("Test for invalid inputs", {

  expect_error(GPD_Fit(Data=S20.Detrend.Declustered.df$Rainfall,Data_Full="Not numeric",
                       u=0.997,PLOT=FALSE,xlab_hist="Rainfall (Inches)",y_lab="Rainfall (Inches)"),
               "Data_Full must be a numeric vector.")

  expect_error(GPD_Fit(Data=S20.Detrend.Declustered.df$Rainfall,Data_Full=na.omit(S20.Detrend.df$Rainfall),
                      u=1.5,PLOT=FALSE,xlab_hist="Rainfall (Inches)",y_lab="Rainfall (Inches)"),
              "u must be a numeric value strictly between 0 and 1.")

  expect_error(GPD_Fit(Data=S20.Detrend.Declustered.df$Rainfall,Data_Full=na.omit(S20.Detrend.df$Rainfall),
                       u=0.997,Method = "Unknown",PLOT=FALSE,xlab_hist="Rainfall (Inches)",y_lab="Rainfall (Inches)"),
               "Method must be Standard or Solari.")

  expect_error(GPD_Fit(Data=S20.Detrend.Declustered.df$Rainfall,Data_Full=na.omit(S20.Detrend.df$Rainfall),
                       u=0.997,PLOT=FALSE,xlab_hist=6,y_lab="Rainfall (Inches)"),
               "xlab_hist and y_lab must be character strings.")
})


