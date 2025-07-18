S20.Migpd<-Migpd_Fit(Data=S20.Detrend.Declustered.df,Data_Full=S20.Detrend.df,mqu=c(0.975,0.975,0.9676))


test_that("HT04 works", {

  result <- HT04(data_Detrend_Dependence_df=S20.Detrend.df,
                 data_Detrend_Declustered_df=S20.Detrend.Declustered.df,
                 u_Dependence=0.995,Migpd=S20.Migpd,mu=365.25,N=1000)

  # Checking type of output
  expect_type(result, 'list')
  expect_named(result, c("Model", "Prop", "z", "u.sim","x.sim"))


  #Checking length of outputs
  expect_equal(length(result$Model),3)
  expect_equal(length(result$Prop),3)
  expect_equal(length(result$z),3)
  expect_equal(nrow(result$u.sim),round((1-0.995)*365.25*1000,0))
  expect_equal(nrow(result$x.sim),round((1-0.995)*365.25*1000,0))

  #Checking column names of outputs
  expect_equal(paste(colnames(S20.Detrend.df[,2:4]),".trans",sep=""), colnames(result$u.sim))
  expect_equal(colnames(S20.Detrend.df[,2:4]), colnames(result$x.sim))
  expect_equal(colnames(S20.Detrend.df[,2:4]), names(result$z))
})

test_that("Invalid inputs", {

  expect_error(
    HT04(data_Detrend_Dependence_df=S20.Detrend.df,
         data_Detrend_Declustered_df=S20.Detrend.Declustered.df[,1:3],
         u_Dependence=0.995,Migpd=S20.Migpd,mu=365.25,N=1000),
    "Data frames must have the same number of numeric columns after removing date/factor columns.")

  S20.Migpd.modified = S20.Migpd$models
  S20.Migpd.modified$models = S20.Migpd$models[-1]
  expect_error(
    HT04(data_Detrend_Dependence_df=S20.Detrend.df,
         data_Detrend_Declustered_df=S20.Detrend.Declustered.df,
         u_Dependence=0.995,Migpd=S20.Migpd.modified,mu=365.25,N=1000),
    "Number of models.*must match.*number of data columns")
})
