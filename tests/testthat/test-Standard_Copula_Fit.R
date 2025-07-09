#Test that functions work as intended

test_that("Gaussian copula works", {
  result <- Standard_Copula_Fit(Data=S20.Detrend.df, Copula_Type="Gaussian")
  expect_type(result, 'S4')
  expect_true(inherits(result, "normalCopula"))

  # Test parameter validity
  params <- getTheta(result)
  expect_true(all(params >= -1 & params <= 1))  # Correlation bounds
  expect_equal(length(params), 3)
})

test_that("t-copula works", {
  result <- Standard_Copula_Fit(Data=S20.Detrend.df, Copula_Type="tcopula")
  expect_type(result, 'S4')
  expect_true(inherits(result, "tCopula"))

  # Test parameter validity
  params <- getTheta(result)
  expect_true(all(params[1:3] >= -1 & params[1:3] <= 1))
  expect_true(params[4] >0)
  expect_equal(length(params), 4)

})

test_that("Gumbel copula works", {
  result <- Standard_Copula_Fit(Data=S20.Detrend.df, Copula_Type="Gumbel")
  expect_type(result, 'S4')
  expect_true(inherits(result, "gumbelCopula"))
  params <- getTheta(result)
  expect_true(all(params >= 1))
})

test_that("Clayton copula works", {
  result <- Standard_Copula_Fit(Data=S20.Detrend.df, Copula_Type="Clayton")
  expect_type(result, 'S4')
  expect_true(inherits(result, "claytonCopula"))
  params <- getTheta(result)
  expect_equal(length(params), 1)
  expect_true(all(params > 0))
})

test_that("Frank copula works", {
  result <- Standard_Copula_Fit(Data=S20.Detrend.df, Copula_Type="Frank")
  expect_type(result, 'S4')
  expect_true(inherits(result, "frankCopula"))
  params <- getTheta(result)
  expect_true(is.finite(params))
  expect_true(is.numeric(params))
  expect_equal(length(params), 1)
})



# Test that invalid inputs gives errors
test_that("Invalid inputs produce errors", {

  expect_error(Standard_Copula_Fit(Data = "", Copula_Type="Gumbel"),
               "Error in Data[, 1] : incorrect number of dimensions")

  expect_error(Standard_Copula_Fit(Data=S20.Detrend.df[1:5,], Copula_Type="Gumbel"),
               "Error: Insufficient non-missing data (need at least 10 complete observations)")

  expect_error(Standard_Copula_Fit(Data=S20.Detrend.df[1:5,], Copula_Type="Not a Copula"),
               "Error: Copula_Type must be one of: Gaussian, tcopula, Gumbel, Clayton, Frank")

  expect_error(Standard_Copula_Fit(Data=6, Copula_Type = "Gaussian"),
               "Error: Data must be a data.frame or matrix")

})

# Test default parameter behavior
test_that("Default Copula_Type works", {
  result <- Standard_Copula_Fit(Data = S20.Detrend.df)  # Should default to Gaussian
  expect_true(inherits(result, "normalCopula"))
})

# Test with different first column types
test_that("Different first column types handled correctly", {
  S20.Detrend.df.modified <- data.frame(
    Category = factor(rep(c("A", "B"), length.out = nrow(S20.Detrend.df))),
    S20.Detrend.df[,2:ncol(S20.Detrend.df)]
  )
  result <- Standard_Copula_Fit(Data = S20.Detrend.df.modified, Copula_Type = "Gaussian")
  expect_equal(result@dimension, ncol(S20.Detrend.df) - 1)
})
