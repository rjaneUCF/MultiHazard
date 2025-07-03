#Test that functions work as intended

test_that("Annual maximum functions basic functionality", {

  result <- Annual_Max(Data_Detrend = S20_T_MAX_Daily_Completed_Detrend_Declustered[,c(2,5)])

  # Check return type
  expect_type(result, "list")

  # Check names of output
  expect_named(result, c("Event", "AM"))

  # Verify outputs are numeric
  expect_type(result$Event, "double")
  expect_type(result$AM, "double")

  # Check length of outputs
  expect_equal(length(result$Event), length(result$AM))
  expect_less_than(length(result$Event), nrow(S20_T_MAX_Daily_Completed_Detrend_Declustered[,c(2,5)])/365.25 +10)
  expect_greater_than(length(result$Event), nrow(S20_T_MAX_Daily_Completed_Detrend_Declustered[,c(2,5)])/365.25 -10)

})


#Test function with different percent complete values
test_that("Annual maximum function with different Complete_Prop values", {

  result1 <- Annual_Max(Data_Detrend = S20_T_MAX_Daily_Completed_Detrend_Declustered[,c(2,5)], Complete_Prop = 0.9)
  result2 <- Annual_Max(Data_Detrend = S20_T_MAX_Daily_Completed_Detrend_Declustered[,c(2,5)], Complete_Prop = 0.5)

  expect_true(length(result1$Event) <= length(result2$Event))
  expect_true(length(result1$AM) <= length(result2$AM))
})


# Test that invalid inputs gives errors
test_that("Invalid inputs produce errors", {

  expect_error(Annual_Max(Data_Detrend = S20_T_MAX_Daily_Completed_Detrend_Declustered[,c(2,5)],
                          Complete_Prop = -1))
  expect_error(Annual_Max(Data_Detrend = S20_T_MAX_Daily_Completed_Detrend_Declustered[,c(2,5)],
                          Complete_Prop = 1.15))
  expect_error(Annual_Max(Data_Detrend = S20_T_MAX_Daily_Completed_Detrend_Declustered[,5]))
  expect_error(Annual_Max(Data_Detrend = letters[1:10]))
})

#Test reproducibility

test_that("Function is deterministic", {

  #Find annual max
  result1 <- Annual_Max(Data_Detrend = S20_T_MAX_Daily_Completed_Detrend_Declustered[,c(2,5)])
  result2 <- Annual_Max(Data_Detrend = S20_T_MAX_Daily_Completed_Detrend_Declustered[,c(2,5)])

  #See whether outputs are the same
  expect_identical(result1$Event, result2$Event)
  expect_identical(result1$AM, result2$AM)
})
