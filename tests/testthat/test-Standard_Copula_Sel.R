

test_that("Basic function works correctly", {
    result <- Standard_Copula_Sel(Data=S20.Detrend.df)
    
    # Check return type
    expect_type(result, "data.frame")
    
    #Check names of output
    expect_equal(colnames(result), c("Copula","AIC"))
   
    # Check distribution names
    expected_cop <- c("Gaussian","t-cop","Gumbel","Clayton","Frank")
    expect_equal(result$Copula, expected_cop)
    
    # Check length of outputs
    expect_equal(nrow(result), 5)        
    
    # AIC values should be numeric where not NA
    non_na_aic <- result$AIC[!is.na(result$AIC)]
    expect_true(all(is.numeric(non_na_aic)))
    expect_true(all(is.finite(non_na_aic)))
    
})


test_that("Invalid inputs produce errors", {
  
  expect_error(Standard_Copula_Sel(),
               "Error: Data is required")
  
  expect_error(Standard_Copula_Sel(Data=1:5),
               "Error: Data must be a data.frame or matrix")
  
  expect_error(Standard_Copula_Sel(Data=S20.Detrend.df[1:5,]),
               "Error: Data must have at least 10 rows")
  
})