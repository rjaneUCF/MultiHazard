
#Test that functions work as intended

test_that("Declsuter functions basic functionality", {
  test_data <-S20_T_MAX_Daily_Completed_Detrend_Declustered[,c(2,5)]

  result <- Decluster(test_data[,2], u = 0.98, SepCrit = 3, mu = 365.25)

  # Check return type
  expect_type(result, "list")

  #Check names of output
  expect_named(result, c("Threshold", "Rate", "EventsMax", "Detrended", "Declustered"))

  # Check length of outputs
  expect_length(result$Threshold, 1)
  expect_length(result$Rate, 1)
  expect_length(result$Detrended, nrow(test_data))
  expect_length(result$Declustered, nrow(test_data))

  # Check that the result is same as detrended data in the package
  expect_false(identical(result, S20_T_MAX_Daily_Completed_Detrend_Declustered[,6]))
})

# Test declustering logic
test_that("Declustering actually works", {
  test_vector <- S20_T_MAX_Daily_Completed_Detrend_Declustered[,5]
  result <- Decluster(test_vector)

  # Should have more NAs in declustered than original
  original_nas <- sum(is.na(test_vector))
  declustered_nas <- sum(is.na(result$Declustered))
  expect_true(declustered_nas >= original_nas)
})

# Test parameter calculation is valid
test_that("Decluster function parameter validation", {

  # Test 8: Parameter validation
  test_data <- S20_T_MAX_Daily_Completed_Detrend_Declustered[,c(2,5)]

  # Test with different u values
  result1 <- Decluster(test_data[,2], u = 0.5)
  result2 <- Decluster(test_data[,2], u = 0.9)

  expect_true(result1$Threshold <= result2$Threshold)

  # Test with different SepCrit values
  result3 <- Decluster(test_data[,2], SepCrit = 1)
  result4 <- Decluster(test_data[,2], SepCrit = 5)

  expect_is(result3, "list")
  expect_is(result4, "list")
})

test_that("Decluster function with different mu values", {

  # Test 9: Different mu values (annual rate calculation)
  test_data <- S20_T_MAX_Daily_Completed_Detrend_Declustered[,c(2,5)]

  result1 <- Decluster(test_data[,2], mu = 365.25)
  result2 <- Decluster(test_data[,2], mu = 1)

  expect_true(result1$Rate != result2$Rate)
  expect_true(result1$Rate > result2$Rate) # (annual rate > daily rate)
})

test_that("Threshold calculation is correct", {
  test_vector <- S20_T_MAX_Daily_Completed_Detrend_Declustered[,5]
  result <- Decluster(test_vector, u = 0.95)
  expected_threshold <- quantile(test_vector, 0.95, na.rm = TRUE)
  expect_equal(result$Threshold, as.numeric(expected_threshold))
})

#Test reproducibility

test_that("Function is deterministic", {
  test_data <- S20_T_MAX_Daily_Completed_Detrend_Declustered[,c(2,5)]

  # Window method should be reproducible
  result1 <- Decluster(test_data[,2], u = 0.98, SepCrit = 3, mu = 365.25)
  result2 <- Decluster(test_data[,2], u = 0.98, SepCrit = 3, mu = 365.25)
  expect_identical(result1, result2)
})


# Test NA handling
test_that("Handles NA values", {
  test_data_with_na <- S20_T_MAX_Daily_Completed_Detrend_Declustered[,c(2,3)]
  result <- Decluster(test_data_with_na[,2], u = 0.95)
  expect_type(result, "list")
  expect_true(any(is.na(result$Detrended)) || any(is.na(result$Declustered)))
})

# Test that invalid inputs gives errors
test_that("Invalid inputs produce errors", {
  test_vector <- S20_T_MAX_Daily_Completed_Detrend_Declustered[,5]
  expect_error(Decluster(test_vector, u = 1.5))      # u > 1
  expect_error(Decluster(test_vector, u = -0.1))     # u < 0
  expect_error(Decluster(letters[1:10]))             # character data
})

# Test data integrity
test_that("Data integrity preserved", {
  test_vector <- S20_T_MAX_Daily_Completed_Detrend_Declustered[,5]
  result <- Decluster(test_vector, u = 0.95)
  expect_equal(length(result$Declustered), length(test_vector))
})
