#Setting up test data

setup_test_data <-function(){
 GW_S20<-data.frame(G_3356,G_3355[65:12336,-1])
 colnames(GW_S20)<-c("Date","G3356","G3355")
 return(GW_S20)
}

#Test that function work as intended
test_that("Linear method works", {

  GW_S20 <- setup_test_data()

  result <- Imputation(Data=GW_S20,Variable="G3356",
                       x_lab="Groundwater level (ft NGVD 29)",
                       y_lab="Groundwater level (ft NGVD 29)")

  # Check return type
  expect_type(result, "list")

  # Check return names
  expect_named(result, c("Data", "Model"))

  # Check return length
  expect_equal(nrow(result$Data), nrow(GW_S20))

  # Check that imputation changes the data
  expect_false(identical(result$Data$ValuesFilled, GW_S20$G3356))

  # Check new column exists
  expect_true("ValuesFilled" %in% names(result$Data))

  # Verify imputed values are numeric
  expect_type(result$Data$ValuesFilled, "double")

  # Check that imputation reduces number of NAs
  expect_lt(sum(is.na(result$Data$ValuesFilled)), sum(is.na(GW_S20$G3356)))
})

# Test model quality
test_that("Linear model is reasonable", {
  GW_S20 <- setup_test_data()

  result <- Imputation(Data=GW_S20, Variable="G3356",
                       x_lab="GW Level", y_lab="GW Level")

  # Check that model summary exists and has expected structure
  expect_s3_class(result$Model, "summary.lm")
  expect_true("coefficients" %in% names(result$Model))
  expect_true("r.squared" %in% names(result$Model))

  # Basic sanity check - R-squared should be between 0 and 1
  expect_gte(result$Model$r.squared, 0)
  expect_lte(result$Model$r.squared, 1)
})


#Test reproducibility
test_that("Function is deterministic", {
  GW_S20 <- setup_test_data()

  # Window method should be reproducible
  result1 <- Imputation(Data=GW_S20,Variable="G3356",
                        x_lab="Groundwater level (ft NGVD 29)",
                        y_lab="Groundwater level (ft NGVD 29)")
  result2 <- Imputation(Data=GW_S20,Variable="G3356",
                        x_lab="Groundwater level (ft NGVD 29)",
                        y_lab="Groundwater level (ft NGVD 29)")
  expect_identical(result1$Data$ValuesFilled, result2$Data$ValuesFilled)

})

# Test that only missing values are imputed
test_that("Only missing values are imputed", {
  GW_S20 <- setup_test_data()

  result <- Imputation(Data=GW_S20, Variable="G3356",
                       x_lab="GW Level", y_lab="GW Level")

  # Non-missing values should remain unchanged
  non_missing_idx <- !is.na(GW_S20$G3356)
  expect_equal(result$Data$ValuesFilled[non_missing_idx],
               GW_S20$G3356[non_missing_idx])
})

#Test that invalid inputs gives errors

test_that("Invalid inputs produce errors", {
  GW_S20 <- setup_test_data()
  expect_error(Imputation(Data=GW_S20,Variable="G666"))
})


# Test with minimal data
test_that("Handles edge cases", {
  # Create minimal test data
  minimal_data <- data.frame(
    Date = as.Date(c("2020-01-01", "2020-01-02", "2020-01-03")),
    var1 = c(1, NA, 3),
    var2 = c(2, 4, 6)
  )

  result <- Imputation(Data=minimal_data, Variable="var1",
                       x_lab="X", y_lab="Y")

  expect_type(result, "list")
  expect_false(is.na(result$Data$ValuesFilled[2]))  # Should impute the NA
})

