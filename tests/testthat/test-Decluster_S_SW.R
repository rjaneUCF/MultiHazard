
#Test that functions work as intended
print(summary(S13_Rainfall))
test_that("Decluster_SW functions basic functionality", {
 
  result <- Decluster_S_SW(Data = S13_Rainfall, Window_Width_Sum_24=24, Window_Width = 7*24)
  
  # Check return type
  expect_type(result, "list")
  
  #Check names of output
  expect_named(result, c("Detrend", "Totals", "Declustered", "EventID", "Event_start", "Event_End"))
  
  # Check length of outputs
  expect_equal(length(result$Detrend), nrow(test_data))
  expect_equal(length(result$Totals),  nrow(test_data))
  expect_equal(length(result$Declustered),  nrow(test_data))
  
})

# Test declustering logic
test_that("Declustering actually works", {
  test_data =S13_Rainfall
  #test_data[,1] <- as.Date(test_data[,1], format = "%Y-%m-%d")
  
  result <- Decluster_S_SW(Data = test_data, Window_Width_Sum_24=24, Window_Width = 7*24)
  
  # Should have more NAs in declustered than original
  original_nas <- sum(is.na(test_data[,2]))
  declustered_nas <- sum(is.na(result$Declustered))
  expect_true(declustered_nas >= original_nas)
  
  # Filter out boundary events that can't be tested
  testable_events <- result$EventID[result$EventID > 1 & result$EventID < nrow(test_data)]
  
  if (length(testable_events) > 0) {
    # Test that events are >= neighbors (allows for equal values)
    left_neighbors <- test_data[testable_events - 1, 2]
    right_neighbors <- test_data[testable_events + 1, 2]
    event_values <- test_data[testable_events, 2]
    
    # Remove any comparisons involving NA values
    valid_comparisons <- !is.na(event_values) & !is.na(left_neighbors) & !is.na(right_neighbors)
    
    if (sum(valid_comparisons) > 0) {
      expect_true(all(event_values[valid_comparisons] >= left_neighbors[valid_comparisons]))
      expect_true(all(event_values[valid_comparisons] >= right_neighbors[valid_comparisons]))
    }
  }
})


# Test parameter calculation is valid
test_that("Decluster function parameter validation", {
  test_data =S13_Rainfall
  # Test 8: Parameter validation
  #test_data = data.frame(S20_T_MAX_Daily_Completed_Detrend_Declustered$Date,S20_T_MAX_Daily_Completed_Detrend_Declustered$Detrend)
  #test_data[,1] <- as.Date(test_data[,1], format = "%Y-%m-%d")
  
  # Test with different Window_Width values
  result1 <- Decluster_S_SW(test_data, Window_Width_Sum_24=24, Window_Width = 3*24)
  result2 <- Decluster_S_SW(test_data, Window_Width_Sum_24=24, Window_Width = 7*24)
  
  expect_true(length(result1$EventID) >= length(result2$EventID))
})


#Test reproducibility

test_that("Function is deterministic", {
  #test_data = data.frame(S20_T_MAX_Daily_Completed_Detrend_Declustered$Date,S20_T_MAX_Daily_Completed_Detrend_Declustered$Detrend)
  #test_data[,1] <- as.Date(test_data[,1], format = "%Y-%m-%d")
  test_data =S13_Rainfall
  
  result1 <- Decluster_S_SW(Data = test_data, Window_Width_Sum_24=24, Window_Width = 7*24)
  result2 <- Decluster_S_SW(Data = test_data, Window_Width_Sum_24=24, Window_Width = 7*24)
  expect_identical(result1, result2)
})


# Test NA handling
test_that("Handles NA values", {
  test_data_with_na = data.frame(S20_T_MAX_Daily_Completed_Detrend_Declustered$Date,S20_T_MAX_Daily_Completed_Detrend_Declustered$Detrend)
  test_data_with_na[,1] <- as.Date(test_data_with_na[,1], format = "%Y-%m-%d")
  result <- Decluster_S_SW(test_data_with_na, Window_Width_Sum_24=24, Window_Width = 7*24)
  expect_type(result, "list")
  expect_true(any(is.na(result$Detrend)) || any(is.na(result$Declustered)))
})

# Test that invalid inputs gives errors
test_that("Invalid inputs produce errors", {
  #test_data = data.frame(S20_T_MAX_Daily_Completed_Detrend_Declustered$Date,S20_T_MAX_Daily_Completed_Detrend_Declustered$Detrend)
  #test_data[,1] <- as.Date(test_data[,1], format = "%Y-%m-%d")
  
  test_data =S13_Rainfall
  expect_error(Decluster_S_SW(Data = test_data[,1],Window_Width_Sum_24=24, Window_Width = 7*24),
               "Data must be a data frame")
  
  expect_error(Decluster_S_SW(), "Data parameter is required")
  
  expect_error(Decluster_S_SW(Data = test_data), "Window_Width parameter is required")
  
  expect_error(Decluster_S_SW(Data = "not_a_df", Window_Width_Sum_24=24, Window_Width = 7*24), "Data must be a data frame")
  
  expect_error(Decluster_S_SW(Data = test_data, Window_Width_Sum_24=24, Window_Width = "Invalid"),
               "Window_Width must be numeric")
  
  expect_error(Decluster_S_SW(Data = test_data, Window_Width_Sum_24=24, Window_Width = 35.5),
               "Window_Width must be a single integer value")
  
  expect_error(Decluster_S_SW(Data = test_data, Window_Width_Sum_24=24, Window_Width= -2),
               "Window_Width must be positive")
  
  test_data[,2] <- rep("2",nrow(test_data))                       
  expect_error(Decluster_S_SW(Data = test_data, Window_Width_Sum_24=24, Window_Width = 7*24),
               "Second column of Data must be numeric")
})

