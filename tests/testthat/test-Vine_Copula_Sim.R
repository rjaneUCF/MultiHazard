
S20.Migpd<-Migpd_Fit(Data=S20.Detrend.Declustered.df,Data_Full=S20.Detrend.df,mqu=c(0.975,0.975,0.9676))
S20.Vine<-Vine_Copula_Fit(Data=S20.Detrend.df)

test_that("Vine_Copula_Sim works", {

  MU = 365.25
  N_sim = 10
  result <- Vine_Copula_Sim(Data=S20.Detrend.df,Vine_Model=S20.Vine,
                            Marginals=S20.Migpd,N=10,mu=MU)

  # Checking type of output
  expect_type(result, 'list')
  expect_named(result,c("Vine_Model", "u.Sim", "x.Sim"))

  # Checking length of outputs
  expect_equal(nrow(result$u.Sim),round(MU*N_sim,0))
  expect_equal(nrow(result$x.Sim),round(MU*N_sim,0))
  expect_equal(dim(result$u.Sim), dim(result$x.Sim))
  expect_true(all(result$u.Sim>0 & result$u.Sim<1))
  expect_true(all(sapply(result$x.Sim, is.numeric)))

  #Checking column names of outputs
  expect_equal(colnames(S20.Detrend.df)[-1], colnames(result$u.Sim))
  expect_equal(colnames(S20.Detrend.df)[-1], colnames(result$x.Sim))

})

test_that("Check inputs are valid", {

  expect_error(Vine_Copula_Sim(Data=5,Vine_Model=S20.Vine,
                               Marginals=S20.Migpd,N=10),
               "Error: Data must be a data frame or matrix.")

  expect_error(Vine_Copula_Sim(Data=S20.Detrend.df[1:5,],Vine_Model=S20.Vine,
                               Marginals=S20.Migpd,N=10),
               "Error: Data must contain at least 5 rows.")

  invalid.model <- S20.Vine
  invalid_model$Par2 <- NULL
  expect_error(Vine_Copula_Sim(Data=S20.Detrend.df, Vine_Model=invalid_model, Marginals=S20.Migpd, N=10),
               "Vine_Model must be a list with names 'Structure', 'Family', 'Par', and 'Par2'"
  )

})
