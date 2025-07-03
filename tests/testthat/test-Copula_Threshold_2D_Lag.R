#Test that function work as intended
test_that("Copula_Threshold_2D_Lag method works", {
  
  result <- Copula_Threshold_2D_Lag(Data_Detrend=S20.Detrend.df[,-c(1,4)],
                                    Data_Declust=S20.Detrend.Declustered.df[,-c(1,4)],
                                    y_lim_min=-0.075, y_lim_max =0.25,
                                    Upper=c(6,8), Lower=c(6,8),GAP=0.1, PLOT=FALSE)
  
  # Check return type
  expect_type(result, "list")
  
  #Check names of output
  expect_named(result, c("Kendalls_Tau1", "p_value_Var1", "N_Var1", "Copula_Family_Var1",
                         "Kendalls_Tau2", "p_value_Var2", "N_Var2", "Copula_Family_Var2"))
  
  # Check length of outputs
  expect_equal(length(result$Kendalls_Tau1), length(seq(0.9,0.99,0.01)))
  expect_equal(length(result$Kendalls_Tau2), length(seq(0.9,0.99,0.01)))
  expect_equal(length(result$p_value_Var1), length(seq(0.9,0.99,0.01)))
  expect_equal(length(result$p_value_Var2), length(seq(0.9,0.99,0.01)))
  expect_equal(length(result$N_Var1), length(seq(0.9,0.99,0.01)))
  expect_equal(length(result$N_Var2), length(seq(0.9,0.99,0.01)))
  expect_true(all(result$Copula_Family_Var1 %in% c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 13, 14, 16, 17, 18, 19, 20, 23, 24, 26, 27, 28, 29, 30, 33, 34, 36, 37, 38, 39, 40, 104, 114, 124, 134, 204, 214, 224, 234)))
  expect_true(all(result$Copula_Family_Var2 %in% c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 13, 14, 16, 17, 18, 19, 20, 23, 24, 26, 27, 28, 29, 30, 33, 34, 36, 37, 38, 39, 40, 104, 114, 124, 134, 204, 214, 224, 234)))
  expect_true(all(result$Kendalls_Tau1 >-1 & result$Kendalls_Tau1<=1))
  expect_true(all(result$Kendalls_Tau2 >=-1 & result$Kendalls_Tau2 <=1))
  expect_true(all(result$p_value_Var1 >0 & result$p_value_Var1 <=1))
  expect_true(all(result$p_value_Var2 >0 & result$p_value_Var2 <=1))
})

#Test reproducibility
test_that("Function is deterministic", {
  
  result1 <- Copula_Threshold_2D_Lag(Data_Detrend=S20.Detrend.df[,-c(1,4)],
                                     Data_Declust=S20.Detrend.Declustered.df[,-c(1,4)],
                                     y_lim_min=-0.075, y_lim_max =0.25,
                                     Upper=c(6,8), Lower=c(6,8),GAP=0.1, PLOT=FALSE)
  result2 <- Copula_Threshold_2D_Lag(Data_Detrend=S20.Detrend.df[,-c(1,4)],
                                     Data_Declust=S20.Detrend.Declustered.df[,-c(1,4)],
                                     y_lim_min=-0.075, y_lim_max =0.25,
                                     Upper=c(6,8), Lower=c(6,8),GAP=0.1, PLOT=FALSE)
  expect_identical(result1, result2)
  
})



#Test that invalid inputs gives errors

test_that("Invalid inputs produce errors", {
  expect_error(Copula_Threshold_2D_Lag(Data_Detrend=S20.Detrend.df[,-c(1,4)],
                                       Data_Declust=S20.Detrend.Declustered.df[,-c(1,4)],
                                       y_lim_min=-0.075, y_lim_max =0.25,
                                       u1= seq(-1,1,0.1),u2=seq(0.5,0.99,0.01),
                                       Upper=c(6,8), Lower=c(6,8),GAP=0.1, PLOT=FALSE))
  
  expect_error(Copula_Threshold_2D_Lag(Data_Detrend=S20.Detrend.df[,-c(1,4)],
                                       Data_Declust=S20.Detrend.Declustered.df[,-c(1,4)],
                                       y_lim_min=-0.075, y_lim_max =0.25,
                                       u1= seq(0.5,0.99,0.01),u2="Invalid",
                                       Upper=c(6,8), Lower=c(6,8),GAP=0.1, PLOT=FALSE))
  
  expect_error(Copula_Threshold_2D_Lag(Data_Detrend=S20.Detrend.df[,-c(1,4)],
                                       Data_Declust=S20.Detrend.Declustered.df[,4],
                                       y_lim_min=-0.075, y_lim_max =0.25,
                                       u1= seq(-1,1,0.1),u2=seq(0.5,0.99,0.01),
                                       Upper=c(6,8), Lower=c(6,8),GAP=0.1, PLOT=FALSE))
  
  expect_error(Copula_Threshold_2D_Lag(Data_Detrend=data.frame(),
                                       Data_Declust=S20.Detrend.Declustered.df[,-c(1,4)],
                                       y_lim_min=-0.075, y_lim_max =0.25,
                                       u1= seq(-1,1,0.1),u2=seq(0.5,0.99,0.01),
                                       Upper=c(6,8), Lower=c(6,8),GAP=0.1, PLOT=FALSE))
  
  expect_error(
    Copula_Threshold_2D_Lag(Data_Detrend=S20.Detrend.df[,-c(1,4)],
                            Data_Declust=S20.Detrend.Declustered.df[,4],
                            y_lim_min=-0.075, y_lim_max =0.25,
                            u1= seq(0.5,0.99,0.01),u2=seq(0.5,0.99,0.01),
                            Upper=c(6,8), Lower=c(6,8),
                            Lag_Backward_Var1="Invalid", Lag_Forward_Var1=3,
                            Lag_Backward_Var2=3, Lag_Forward_Var2=3,
                            GAP=0.1, PLOT=FALSE),
    "Lag_Backward_Var1 must be numeric" )
  
  expect_error(Copula_Threshold_2D_Lag(Data_Detrend=S20.Detrend.df[,-c(1,4)],
                                       Data_Declust=S20.Detrend.Declustered.df[,4],
                                       y_lim_min=-0.075, y_lim_max =0.25,
                                       u1= seq(0.5,0.99,0.01),u2=seq(0.5,0.99,0.01),
                                       Upper=c(6,8), Lower=c(6,8),
                                       Lag_Backward_Var1=3, Lag_Forward_Var1="Invalid",
                                       Lag_Backward_Var2=3, Lag_Forward_Var2=3,
                                       GAP=0.1, PLOT=FALSE), "Lag_Forward_Var1 must be numeric")
  
  expect_error(Copula_Threshold_2D_Lag(Data_Detrend=S20.Detrend.df[,-c(1,4)],
                                       Data_Declust=S20.Detrend.Declustered.df[,4],
                                       y_lim_min=-0.075, y_lim_max =0.25,
                                       u1= seq(0.5,0.99,0.01),u2=seq(0.5,0.99,0.01),
                                       Upper=c(6,8), Lower=c(6,8),
                                       Lag_Backward_Var1=3, Lag_Forward_Var1=3,
                                       Lag_Backward_Var2="Invalid", Lag_Forward_Var2=3,
                                       GAP=0.1, PLOT=FALSE), "Lag_Backward_Var2 must be numeric")
  
  expect_error(Copula_Threshold_2D_Lag(Data_Detrend=S20.Detrend.df[,-c(1,4)],
                                       Data_Declust=S20.Detrend.Declustered.df[,4],
                                       y_lim_min=-0.075, y_lim_max =0.25,
                                       u1= seq(0.5,0.99,0.01),u2=seq(0.5,0.99,0.01),
                                       Upper=c(6,8), Lower=c(6,8),
                                       Lag_Backward_Var1=3, Lag_Forward_Var1=3,
                                       Lag_Backward_Var2=3, Lag_Forward_Var2="Invalid",
                                       GAP=0.1, PLOT=FALSE), "Lag_Forward_Var2 must be numeric")
  
  expect_error(Copula_Threshold_2D_Lag(Data_Detrend=S20.Detrend.df[,-c(1,4)],
                                       Data_Declust=S20.Detrend.Declustered.df[,4],
                                       y_lim_min=-0.075, y_lim_max =-0.25,
                                       u1= seq(0.5,0.99,0.01),u2=seq(0.5,0.99,0.01),
                                       Upper=c(6,8), Lower=c(6,8),
                                       Lag_Backward_Var1=3, Lag_Forward_Var1=3,
                                       Lag_Backward_Var2=3, Lag_Forward_Var2=3,
                                       GAP=0.1, PLOT=FALSE), "y_lim_min must be less than y_lim_max")
})