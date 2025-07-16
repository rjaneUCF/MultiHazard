

final.month = data.frame(seq(as.Date("2019-02-04"),as.Date("2019-02-28"),by="day"),NA,NA,NA)
colnames(final.month) = c("Date","Rainfall","OsWL","Groundwater")
S22.Detrend.df.extended = rbind(S22.Detrend.df,final.month)

test_that("return_curve_est_diag works", {


  result <-  return_curve_diag(data=S22.Detrend.df.extended[,1:3],
                               q=0.985,rp=1,mu=365.25,n_sim=100,
                               n_grad=50,n_boot=100,boot_method="monthly",
                               boot_replace=NA, block_length=NA, boot_prop=0.8,
                               decl_method_x="runs", decl_method_y="runs",
                               window_length_x=NA,window_length_y=NA,
                               u_x=0.95, u_y=0.95,
                               sep_crit_x=36, sep_crit_y=36,
                               alpha=0.1,
                               boot_method_all="block", boot_replace_all=NA,
                               block_length_all=14)

  # Checking type of output
  expect_type(result, 'list')
  expect_named(result, c("ang_ind",
                         "med_x_ht04", "lb_x_ht04", "ub_x_ht04",
                         "med_y_ht04", "lb_y_ht04", "ub_y_ht04",
                         "med_x_wt13", "lb_x_wt13", "ub_x_wt13",
                         "med_y_wt13", "lb_y_wt13", "ub_y_wt13",
                         "med_ht04", "lb_ht04", "ub_ht04",
                         "med_wt13", "lb_wt13", "ub_wt13"))


  #Checking length of outputs
  expect_equal(length(result$ang_ind),50)
  expect_equal(length(result$med_x_ht04),50)
  expect_equal(length(result$lb_x_ht04),50)
  expect_equal(length(result$ub_x_ht04),50)
  expect_equal(length(result$med_y_ht04),50)
  expect_equal(length(result$lb_y_ht04),50)
  expect_equal(length(result$ub_y_ht04),50)

  #Checking column names of outputs
  expect_equal(colnames(S22.Detrend.df.extended)[-1], names(result$median_ht04))
  expect_equal(colnames(S22.Detrend.df.extended)[-1], names(result$ub_ht04))
  expect_equal(colnames(S22.Detrend.df.extended)[-1], names(result$lb_ht04))
})


test_that("return_curve_est_diag invalid inputs", {

  expect_error(
      return_curve_diag(data=S22.Detrend.df.extended[,1:3],
                               q=0.985,rp=1,mu=365.25,n_sim=100,
                               n_grad=50,n_boot=100,boot_method="monthly",
                               boot_replace=NA, block_length=NA, boot_prop=0.8,
                               decl_method_x="runs", decl_method_y="runs",
                               window_length_x=NA,window_length_y=NA,
                               u_x=1.5, u_y=0.95,
                               sep_crit_x=36, sep_crit_y=36,
                               alpha=0.1,
                               boot_method_all="block", boot_replace_all=NA,
                               block_length_all=14),
             "u_x threshold must be between 0 and 1.")

  expect_error(
    return_curve_diag(data=S22.Detrend.df.extended[,1:3],
                      q=0.985,mu=365.25,n_sim=100,
                      n_grad=50,n_boot=100,boot_method="monthly",
                      boot_replace=NA, block_length=NA, boot_prop=0.8,
                      decl_method_x="runs", decl_method_y="runs",
                      window_length_x=NA,window_length_y=NA,
                      u_x=0.95, u_y=0.95,
                      sep_crit_x=36, sep_crit_y=36,
                      alpha=0.1,
                      boot_method_all="block", boot_replace_all=NA,
                      block_length_all=14),
               " rp is missing.")

})

