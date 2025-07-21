
Rainfall_Declust_SW<-suppressWarnings(Decluster_SW(Data=S22.Detrend.df[,c(1:2)],Window_Width=7))


test_that("GPD_Threshold_Solari works", {


  result <-  GPD_Threshold_Solari(Event=Rainfall_Declust_SW$Declustered,
                                  Data=S22.Detrend.df[,2])

  # Checking type of output
  expect_type(result, 'list')
  expect_named(result, c("Thres_Candidate", "Thres_Candidate_Quantile",
                         "GPD_MLE", "CI_Upper","CI_Lower",
                         "AR2", "AR2_pValue",
                         "Candidate_Thres"))


  #Checking length of outputs
  expect_length(result$Thres_Candidate,248)
  expect_length(result$Thres_Candidate_Quantile,962)
  expect_true(nrow(result$GPD_MLE)==248 & ncol(result$GPD_MLE)==11,)
  expect_true(nrow(result$CI_Upper)==248 & ncol(result$CI_Upper)==11,)
  expect_true(nrow(result$CI_Lower)==248 & ncol(result$CI_Lower)==11,)
  expect_length(result$AR2,248)
  expect_length(result$AR2_pValue,248)
  expect_length(result$Candidate_Thres,1)

  #Checking column names of outputs
  expect_equal(colnames(result$GPD_MLE), c("xi","sigma","u","MRLP","mod_sigma","rate","10","50","100","500","1000"))

  #Checking values of outputs
  expect_true(all(result$Thres_Candidate>min(Rainfall_Declust_SW$Declustered,na.rm=T)))
  expect_true(all(result$Thres_Candidate_Quantile>0 & result$Thres_Candidate_Quantile<=1))
  expect_true(all(result$GPD_MLE[,6]>0))
  expect_true(all(result$GPD_MLE[,8]>result$GPD_MLE[,7]))
})

test_that("Invalid inputs", {

  expect_error(GPD_Threshold_Solari(Event=Rainfall_Declust_SW$Declustered,
                                    Data=rep("Not numeric",5)),
              "Data must be a numeric vector.")


  expect_error(GPD_Threshold_Solari(Event=Rainfall_Declust_SW$Declustered,
                                    Data=S22.Detrend.df[,2],
                                    Alpha=5),
               "Alpha")
})

