#Test that functions work as intended

test_that("Annual maximum functions basic functionality", {

  test_data = data.frame(S20_T_MAX_Daily_Completed_Detrend_Declustered$Date,S20_T_MAX_Daily_Completed_Detrend_Declustered$Detrend)
  print(summary(S20_T_MAX_Daily_Completed_Detrend_Declustered$Date))
  test_data[,1] <- as.Date(as.character(test_data[,1]), format = "%m/%d/%Y")
  print(summary(test_data))

  result <- Annual_Max(Data_Detrend = test_data)

  # Check return type
  expect_type(result, "list")

  # Check names of output
  expect_named(result, c("Event", "AM"))

  # Verify outputs are numeric
  expect_type(result$Event, "double")
  expect_type(result$AM, "double")

  # Check length of outputs
  expect_equal(length(result$Event), length(result$AM))
  expect_lt(length(result$Event), nrow(S20_T_MAX_Daily_Completed_Detrend_Declustered[,c(2,5)])/365.25 +10)
  expect_gt(length(result$Event), nrow(S20_T_MAX_Daily_Completed_Detrend_Declustered[,c(2,5)])/365.25 -10)

})


#Test function with different percent complete values
test_that("Annual maximum function with different Complete_Prop values", {

  test_data = data.frame(S20_T_MAX_Daily_Completed_Detrend_Declustered$Date,S20_T_MAX_Daily_Completed_Detrend_Declustered$Detrend)
  test_data[,1] <- as.Date(as.character(test_data[,1]), format = "%m/%d/%Y")

  result1 <- Annual_Max(Data_Detrend = test_data, Complete_Prop = 0.9)
  result2 <- Annual_Max(Data_Detrend = test_data, Complete_Prop = 0.5)

  expect_true(length(result1$Event) <= length(result2$Event))
  expect_true(length(result1$AM) <= length(result2$AM))
})


# Test that invalid inputs gives errors
test_that("Invalid inputs produce errors", {

  test_data = data.frame(S20_T_MAX_Daily_Completed_Detrend_Declustered$Date,S20_T_MAX_Daily_Completed_Detrend_Declustered$Detrend)
  test_data[,1] <- as.Date(as.character(test_data[,1]), format = "%m/%d/%Y")

  expect_error(Annual_Max(Data_Detrend = test_data, Complete_Prop = -1))
  expect_error(Annual_Max(Data_Detrend = test_data, Complete_Prop = 1.15))
  expect_error(Annual_Max(Data_Detrend = test_data[,2]))
  expect_error(Annual_Max(Data_Detrend = letters[1:10]))
})

#Test reproducibility

test_that("Function is deterministic", {

  test_data = data.frame(S20_T_MAX_Daily_Completed_Detrend_Declustered$Date,S20_T_MAX_Daily_Completed_Detrend_Declustered$Detrend)
  test_data[,1] <- as.Date(as.character(test_data[,1]), format = "%m/%d/%Y")

  #Find annual max
  result1 <- Annual_Max(Data_Detrend = test_data)
  result2 <- Annual_Max(Data_Detrend = test_data)

  #See whether outputs are the same
  expect_identical(result1$Event, result2$Event)
  expect_identical(result1$AM, result2$AM)
})
