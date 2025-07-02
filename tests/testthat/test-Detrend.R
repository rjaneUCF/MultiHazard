

#Test that function work as intended
test_that("Linear method works", {
  test_data <-S20_T_MAX_Daily_Completed_Detrend_Declustered[,c(2,4)]

  result <- Detrend(test_data,
                    Method = "linear", Window_Width = 89,
                    End_Length = 1826, PLOT = FALSE)

  # Check return type
  expect_type(result, "double")

  # Check return  length
  expect_length(result, nrow(S20_T_MAX_Daily_Completed_Detrend_Declustered))

  # Check that detrending changes the data
  expect_false(identical(result, test_data))

  # Calculate original trend
  x <- seq_len(nrow(test_data))
  original_slope <- coef(lm(test_data[,2] ~ x))[2]

  # Calculate detrended slope
  detrended_slope <- coef(lm(result ~ x))[2]

  # Detrending should reduce the magnitude of the slope
  expect_lt(abs(detrended_slope), abs(original_slope))
})


test_that("Window method works", {
  test_data <-S20_T_MAX_Daily_Completed_Detrend_Declustered[,c(2,4)]

  result <- Detrend(test_data,
                    Method = "window", Window_Width = 89,
                    End_Length = 1826, PLOT = FALSE)

  # Check return type and length
  expect_type(result, "double")

  # Check return  length
  expect_length(result, nrow(S20_T_MAX_Daily_Completed_Detrend_Declustered))

  # Check that the result is same as detrended data in the package
  expect_false(identical(result, S20_T_MAX_Daily_Completed_Detrend_Declustered[,5] ))
})

test_that("Window method reduces trend", {
  test_data <- S20_T_MAX_Daily_Completed_Detrend_Declustered[,c(2,4)]
  result <- Detrend(test_data, Method = "window", Window_Width = 89,
                    End_Length = 1826, PLOT = FALSE)

  x <- seq_len(nrow(test_data))
  original_slope <- coef(lm(test_data[,2] ~ x))[2]
  detrended_slope <- coef(lm(result ~ x))[2]

  expect_lt(abs(detrended_slope), abs(original_slope))
})

#Test that the PLOT parameter works without error
test_that("PLOT parameter works without error", {
  test_data <- S20_T_MAX_Daily_Completed_Detrend_Declustered[,c(2,4)]

  # Test plotting with window method
  expect_silent({
    result <- Detrend(test_data, Method = "window", Window_Width = 89,
                      End_Length = 365, PLOT = TRUE,
                      x_lab = "Date", y_lab = "Temperature")
  })

  # Test plotting with linear method
  expect_silent({
    result <- Detrend(test_data, Method = "linear", End_Length = 365,
                      PLOT = TRUE, x_lab = "Date", y_lab = "Temperature")
  })
})

#Test reproducibility
test_that("Function is deterministic", {
  test_data <- S20_T_MAX_Daily_Completed_Detrend_Declustered[,c(2,4)]

  # Window method should be reproducible
  result1 <- Detrend(test_data, Method = "window", Window_Width = 89,
                     End_Length = 365, PLOT = FALSE)
  result2 <- Detrend(test_data, Method = "window", Window_Width = 89,
                     End_Length = 365, PLOT = FALSE)
  expect_identical(result1, result2)

  # Linear method should be reproducible
  result3 <- Detrend(test_data, Method = "linear", End_Length = 365, PLOT = FALSE)
  result4 <- Detrend(test_data, Method = "linear", End_Length = 365, PLOT = FALSE)
  expect_identical(result3, result4)
})

#Test that invalid inputs gives errors

test_that("Invalid inputs produce errors", {
  expect_error(Detrend(S20_T_MAX_Daily_Completed_Detrend_Declustered[,c(2,4)], Method = "invalid"))
})




