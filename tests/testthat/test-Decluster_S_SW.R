
#Test that functions work as intended
S13_Rainfall_Test = S13_Rainfall[227795:254098,]

test_that("Decluster_SW functions basic functionality", {

  result <- Decluster_S_SW(Data = S13_Rainfall_Test, Window_Width_Sum=24, Window_Width = 7*24)

  # Check return type
  expect_type(result, "list")

  #Check names of output
  expect_named(result, c("Detrend", "Totals", "Declustered", "EventID", "Event_Start", "Event_End"))

  # Check length of outputs
  expect_equal(length(result$Detrend), nrow(S13_Rainfall_Test))
  expect_equal(length(result$Totals),  nrow(S13_Rainfall_Test))
  expect_equal(length(result$Declustered),  nrow(S13_Rainfall_Test))

})

# Test declustering logic
test_that("Declustering actually works", {

  result <- Decluster_S_SW(Data = S13_Rainfall_Test, Window_Width_Sum=24, Window_Width = 7*24)

  # Should have more NAs in declustered than original
  original_nas <- sum(is.na(S13_Rainfall_Test[,2]))
  declustered_nas <- sum(is.na(result$Declustered))
  expect_true(declustered_nas >= original_nas)

  # Filter out boundary events that can't be tested
  testable_events <- result$EventID[result$EventID > 1 & result$EventID < length(result$Totals)]

  if (length(testable_events) > 0) {
    # Get the summed data (this is what was used for event detection)
    summed_data <- result$Totals

    # Get neighbors and event values as VECTORS (not arrays)
    left_neighbors <- summed_data[testable_events - 1]
    right_neighbors <- summed_data[testable_events + 1]
    event_values <- summed_data[testable_events]

    # Remove any comparisons involving NA values
    valid_comparisons <- !is.na(event_values) & !is.na(left_neighbors) & !is.na(right_neighbors)

    if (sum(valid_comparisons) > 0) {
      expect_true(all(event_values[valid_comparisons] >= left_neighbors[valid_comparisons]),
                  info = "Events should be >= left neighbors")
      expect_true(all(event_values[valid_comparisons] >= right_neighbors[valid_comparisons]),
                  info = "Events should be >= right neighbors")
    }
  }
})


# Test parameter calculation is valid
test_that("Decluster function parameter validation", {

  # Test with different Window_Width values
  result1 <- Decluster_S_SW(S13_Rainfall_Test, Window_Width_Sum=24, Window_Width = 3*24)
  result2 <- Decluster_S_SW(S13_Rainfall_Test, Window_Width_Sum=24, Window_Width = 7*24)

  expect_true(length(result1$EventID) >= length(result2$EventID))
})


#Test reproducibility

test_that("Function is deterministic", {

  result1 <- Decluster_S_SW(Data = S13_Rainfall_Test, Window_Width_Sum=24, Window_Width = 7*24)
  result2 <- Decluster_S_SW(Data = S13_Rainfall_Test, Window_Width_Sum=24, Window_Width = 7*24)
  expect_identical(result1, result2)
})


# Test NA handling
test_that("Handles NA values", {
  test_data_with_na = S13_Rainfall_Test
  result <- Decluster_S_SW(test_data_with_na, Window_Width_Sum=24, Window_Width = 7*24)
  expect_type(result, "list")
  expect_true(any(is.na(result$Detrend)) || any(is.na(result$Declustered)))
})

# Test that invalid inputs gives errors
test_that("Invalid inputs produce errors", {

  expect_error(Decluster_S_SW(Data = S13_Rainfall_Test[,1],Window_Width_Sum=24, Window_Width = 7*24),
               "Data must be a data frame")

  expect_error(Decluster_S_SW(), "Data parameter is required")

  expect_error(Decluster_S_SW(Data = S13_Rainfall_Test), "Window_Width parameter is required")

  expect_error(Decluster_S_SW(Data = "not_a_df", Window_Width_Sum=24, Window_Width = 7*24), "Data must be a data frame")

  expect_error(Decluster_S_SW(Data = S13_Rainfall_Test, Window_Width_Sum=24, Window_Width = "Invalid"),
               "Window_Width must be numeric")

  expect_error(Decluster_S_SW(Data = S13_Rainfall_Test, Window_Width_Sum=24, Window_Width = 35.5),
               "Window_Width must be a single integer value")

  expect_error(Decluster_S_SW(Data = S13_Rainfall_Test, Window_Width_Sum=24, Window_Width= -2),
               "Window_Width must be positive")

  S13_Rainfall_modified = S13_Rainfall_Test
  S13_Rainfall_modified[,2] <- rep("2",nrow(S13_Rainfall_Test))
  expect_error(Decluster_S_SW(Data = S13_Rainfall_modified, Window_Width_Sum=24, Window_Width = 7*24),
               "Second column of Data must be numeric")
})

