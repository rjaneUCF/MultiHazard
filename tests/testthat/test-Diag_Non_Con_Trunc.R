#Test that functions work as intended

S20.OsWL<-Con_Sampling_2D(Data_Detrend=S20.Detrend.df[,-c(1,4)],
                          Data_Declust=S20.Detrend.Declustered.df[,-c(1,4)],
                          Con_Variable="OsWL",Thres=0.97)

S20.OsWL$Data$Rainfall <- S20.OsWL$Data$Rainfall + runif(length(S20.OsWL$Data$Rainfall),0.001,0.01)

test_that("Diag_Non_Con_Trunc functions basic functionality", {

  result <- suppressWarnings(Diag_Non_Con_Trunc(Data=S20.OsWL$Data$Rainfall ,x_lab="Rainfall (Inches)",
                               y_lim_min=0,y_lim_max=1.5))

  # Check return type
  expect_type(result, "list")

  #Check names of output
  expect_named(result, c("AIC", "Best_fit"))
  expect_named(result$AIC, c("Distribution", "AIC"))

  # Check distribution names
  expected_dists <- c("BS","Exp","Gam(2)","Gam(3)","GamMix(2)","GamMix(3)","LNorm","TNorm","Twe","Weib")
  expect_equal(result$AIC$Distribution, expected_dists)

  # Check length of outputs
  expect_equal(nrow(result$AIC), 10)
  expect_equal(length(result$Best_fit), 1)

  # AIC values should be numeric where not NA
  non_na_aic <- result$AIC$AIC[!is.na(result$AIC$AIC)]
  expect_true(all(is.numeric(non_na_aic)))
  expect_true(all(is.finite(non_na_aic)))

  # Best fit should correspond to minimum AIC
  min_aic_dist <- result$AIC$Distribution[which.min(result$AIC$AIC)]
  expect_equal(result$Best_fit, min_aic_dist)

})


#Test reproducibility

test_that("Function is deterministic", {

  result1 <- Diag_Non_Con_Trunc(Data=S20.OsWL$Data$Rainfall ,x_lab="Rainfall (Inches)",
                                 y_lim_min=0,y_lim_max=1.5)
  result2 <- Diag_Non_Con_Trunc(Data=S20.OsWL$Data$Rainfall ,x_lab="Rainfall (Inches)",
                                y_lim_min=0,y_lim_max=1.5)
  expect_identical(result1$AIC, result2$AIC)
  expect_identical(result1$Best_fit, result2$Best_fit)
})


# Test NA handling
test_that("Handles NA values", {

  data_with_na <- c(S20.OsWL$Data$Rainfall , NA, NA, NA)

  expect_warning(Diag_Non_Con_Trunc(Data = data_with_na,x_lab="Rainfall (Inches)",
                              y_lim_min=0,y_lim_max=1.5),
                 "Removed 3 NA values from Data")


  # Test all NA data - should error
  expect_error(Diag_Non_Con_Trunc(Data = rep(NA, 20), x_lab = "Rainfall (Inches)"),
               "Data must be numeric, got: logical")
})

# Test that invalid inputs gives errors
test_that("Invalid inputs produce errors", {

  expect_error(Diag_Non_Con_Trunc(Data=numeric(0),x_lab="Rainfall (Inches)",
                                  y_lim_min=0,y_lim_max=-1.5),
               "Data is empty")

  expect_error(Diag_Non_Con_Trunc(Data=S20.OsWL$Data$Rainfall [1:5],x_lab="Rainfall (Inches)",
                                  y_lim_min=0,y_lim_max=1.5),
               "Data must have at least 10 non-missing observations, got: 5")

  expect_error(Diag_Non_Con_Trunc(Data=S20.OsWL$Data$Rainfall ,x_lab="Rainfall (Inches)",
                                  y_lim_min=0,y_lim_max=-1.5),
               "y_lim_min must be less than y_lim_max, got: y_lim_min = 0, y_lim_max = -1.5")

  expect_error(Diag_Non_Con_Trunc(Data=S20.OsWL$Data$Rainfall ,x_lab="Rainfall (Inches)",
                                  Omit= "Gaussian", y_lim_min=0,y_lim_max=1.5),
               "Invalid distribution names in Omit: Gaussian")

  data_with_inf <- c(S20.OsWL$Data$Rainfall , Inf, -Inf)
  expect_warning(Diag_Non_Con_Trunc(Data = data_with_inf, x_lab = "Rainfall (Inches)"),
                 "Data contains 2 infinite values")

  expect_error(Diag_Non_Con_Trunc(Data = S20.OsWL$Data$Rainfall , x_lab = "Rainfall (Inches)",
                                  Omit = c("BS","Exp","Gam(2)","Gam(3)","GamMix(2)","GamMix(3)","LNorm","TNorm","Twe","Weib")),
               "Cannot omit all distributions")
})

test_that("Test the plots work", {
  # Test that function runs without fatal errors
  expect_no_error({
    result <- suppressWarnings(
      Diag_Non_Con_Trunc(Data = S20.OsWL$Data$Rainfall, x_lab = "Rainfall (Inches)")
    )
  })
})
  
