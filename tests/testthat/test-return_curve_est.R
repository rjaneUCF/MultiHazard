
final.month = data.frame(seq(as.Date("2019-02-04"),as.Date("2019-02-28"),by="day"),NA,NA,NA)
colnames(final.month) = c("Date","Rainfall","OsWL","Groundwater")
S22.Detrend.df.extended = rbind(S22.Detrend.df,final.month)

test_that("return_curve_est works", {


  result <-  return_curve_est(data=S22.Detrend.df.extended[,1:3],
                              q=0.985,rp=100,mu=365.25,n_sim=100,
                              n_grad=50,n_boot=100,boot_method="monthly",
                              boot_replace=NA, block_length=NA, boot_prop=0.8,
                              decl_method_x="runs", decl_method_y="runs",
                              window_length_x=NA,window_length_y=NA,
                              u_x=0.95, u_y=0.95,
                              sep_crit_x=36, sep_crit_y=36,
                              alpha=0.1, x_lab=NA, y_lab=NA,most_likely=TRUE)

  # Checking type of output
  expect_type(result, 'list')
  expect_named(result, c("median_ht04", "ub_ht04", "lb_ht04", "contour_ht04",
                         "most_likely_ht04", "ensemble_ht04", "median_wt13",
                         "ub_wt13", "lb_wt13", "contour_wt13",
                         "most_likely_wt13", "ensemble_wt13"))


  #Checking length of outputs
  expect_equal(nrow(result$median_ht04),50)
  expect_equal(nrow(result$ub_ht04),50)
  expect_equal(nrow(result$lb_ht04),50)
  expect_equal(nrow(result$contour_ht04),50)
  expect_equal(result$ensemble_ht04,NA)
  expect_true(nrow(result$most_likely_ht04)==1 & ncol(result$most_likely_ht04)==2)
  expect_true(nrow(result$most_likely_wt13)==1 & ncol(result$most_likely_wt13)==2)

  #Checking column names of outputs
  expect_equal(colnames(S22.Detrend.df.extended)[-1], colnames(result$median_ht04))
  expect_equal(colnames(S22.Detrend.df.extended)[-1], colnames(result$ub_ht04))
  expect_equal(colnames(S22.Detrend.df.extended)[-1], colnames(result$lb_ht04))
  expect_equal(colnames(S22.Detrend.df.extended)[-1], colnames(result$contour_ht04))
})


test_that("Invalid inputs", {

            expect_error(
                  return_curve_est(data=S22.Detrend.df.extended[,2:3],
                                   q=0.985,rp=100,mu=365.25,n_sim=100,
                                   n_grad=50,n_boot=100,boot_method="monthly",
                                   boot_replace=NA, block_length=NA, boot_prop=0.8,
                                   decl_method_x="runs", decl_method_y="runs",
                                   window_length_x=NA,window_length_y=NA,
                                   u_x=0.95, u_y=0.95,
                                   sep_crit_x=36, sep_crit_y=36,
                                   alpha=0.1, x_lab=NA, y_lab=NA,most_likely=FALSE),
       "Error: First column of data must be Date/POSIXct format")


  expect_error(
    return_curve_est(data=S22.Detrend.df.extended[,1:3],
                     q=0.985,rp=100,mu=365.25,n_sim=100,
                     n_grad=50,n_boot=100,boot_method="monthly",
                     boot_replace=NA, block_length=NA, boot_prop=0.8,
                     decl_method_x="runs", decl_method_y="runs",
                     window_length_x=NA,window_length_y=NA,
                     u_x=1.5, u_y=0.95,
                     sep_crit_x=36, sep_crit_y=36,
                     alpha=0.1, x_lab=NA, y_lab=NA,most_likely=FALSE),
    "Error: u_x threshold must be between 0 and 1")
})
