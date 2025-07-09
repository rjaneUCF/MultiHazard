#Test that functions work as intended
S20.Migpd<-Migpd_Fit(Data=S20.Detrend.Declustered.df[,-1],mqu=c(0.975,0.975,0.9676))
S20.Gaussian<-Standard_Copula_Fit(Data=S20.Detrend.df,Copula_Type="Gaussian")

test_that("Standard_Copula_Sim works", {

  N_sim = 100
  MU = 365.25
  result <- Standard_Copula_Sim(Data=S20.Detrend.df,Marginals=S20.Migpd,Copula=S20.Gaussian,
                                mu= MU,N = N_sim)

  # Checking type of output
  expect_type(result, 'list')
  expect_named(result, c("u.Sim","x.Sim"))

  #Checking length of outputs
  expect_equal(nrow(result$u.Sim), round(MU*N_sim,0))
  expect_equal(nrow(result$x.Sim), round(MU*N_sim,0))
  expect_true(all(result$u.Sim>0 & result$u.Sim<1))

  #Checking column names of outputs
  expect_equal(colnames(S20.Detrend.df)[-1], colnames(result$u.Sim))
  expect_equal(colnames(S20.Detrend.df)[-1], colnames(result$x.Sim))
})

#Checking invalid inputs

test_that("Standard_Copula_Sim invalid inputs", {

  expect_error(Standard_Copula_Sim(Data=S20.Detrend.df,Marginals=12,Copula=S20.Gaussian,
                                   mu=365.25,N=100),
               "Marginals must be a list"
  )
})
