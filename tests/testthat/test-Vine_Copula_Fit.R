

test_that("Vine_Copula_Fit works", {


  result <- Vine_Copula_Fit(Data=S20.Detrend.df)

  # Checking type of output
  expect_type(result, 'list')
  expect_named(result,c("Structure", "Family", "Par", "Par2"))

  # Checking length of outputs
  expect_equal(dim(result$Structure),c(3,3))
  expect_equal(dim(result$Family),c(3,3))
  expect_equal(dim(result$Par),c(3,3))
  expect_equal(dim(result$Par2),c(3,3))

  # Checking type of output
  expect_true(all(is.numeric(result$Structure)))
  expect_true(all(is.numeric(result$Family)))
  expect_true(all(is.numeric(result$Par)))
  expect_true(all(is.numeric(result$Par2)))
})

test_that("Invalid input", {

  expect_error(Vine_Copula_Fit(Data = 5),
               "Error: Data must be a data frame or matrix.")

  expect_error(Vine_Copula_Fit(Data=S20.Detrend.df[1:4,]),
               "Error: Data must contain at least 5 rows.")
})


test_that("Function is deterministic", {

  set.seed(123)
  result1 <- Vine_Copula_Fit(Data=S20.Detrend.df)
  set.seed(123)
  result2 <- Vine_Copula_Fit(Data=S20.Detrend.df)

  expect_identical(result1,result2)
})
