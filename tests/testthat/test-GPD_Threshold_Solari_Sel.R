Rainfall_Declust_SW<-Decluster_SW(Data=S22.Detrend.df[24473:25568,c(1:2)],Window_Width=7)

#Finding an appropriate threshold for the declustered series
S22_OsWL_Solari<-suppressWarnings(GPD_Threshold_Solari(Event=Rainfall_Declust_SW$Declustered,
                                     Data=na.omit(S22.Detrend.df[24473:25568,2])))


test_that("GPD_Threshold_Solari_Sel works", {

  result <- suppressWarnings(GPD_Threshold_Solari_Sel(Event=Rainfall_Declust_SW$Declustered,
                                     Data=S22.Detrend.df[24473:25568,2],
                                     Solari_Output=S22_OsWL_Solari,
                                     Thres=S22_OsWL_Solari$Candidate_Thres))

  # Checking type of output
  expect_type(result, 'list')
  expect_named(result, c("Estimate","CI_Lower","CI_Upper","BOOT"))


  #Checking length of outputs
  expect_equal(length(result$Estimate),32)
  expect_equal(length(result$CI_Lower),32)
  expect_equal(length(result$CI_Upper),32)
  expect_true(nrow(result$BOOT)==10000 & ncol(result$BOOT)==32)
})

test_that("Invalid inputs", {

  expect_error(
    GPD_Threshold_Solari_Sel(Event=Rainfall_Declust_SW$Declustered,
                             Data=S22.Detrend.df[24473:25568,2],
                             Solari_Output = NULL,
                             Thres=S22_OsWL_Solari$Candidate_Thres),
    "Solari_Output is missing.")

  expect_error(
    GPD_Threshold_Solari_Sel(Event=Rainfall_Declust_SW$Declustered,
                             Data=S22.Detrend.df[24473:25568,2],
                             Solari_Output=S22_OsWL_Solari,
                             Thres=S22_OsWL_Solari$Candidate_Thres,RP_Max=-50),
    "RP_Max must be a single positive numeric value.")
})
