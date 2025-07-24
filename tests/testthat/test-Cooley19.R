S20.GPD<-Migpd_Fit(Data=S20.Detrend.Declustered.df[,-c(1,4)], Data_Full=S20.Detrend.df[,-c(1,4)],  mqu =c(0.99,0.99))

test_that("Cooley19 works", {

  result = Cooley19(Data=na.omit(S20.Detrend.df[,2:3]),Migpd=S20.GPD,
                    p.base=0.01,p.proj=0.001,PLOT=TRUE,x_lim_max_T=500,y_lim_max_T=500)


  # Checking type of output
  expect_type(result, 'list')
  expect_named(result, c("Asym", "Chi", "n.bar", "p.base", "p.proj", "I.base", "I.proj"))
  expect_type(result$Asym, "character")
  expect_type(result$Chi, "numeric")
  expect_type(result$n.bar, "numeric")
  expect_type(result$p.base, "numeric")
  expect_type(result$p.proj, "numeric")
  expect_type(result$I.base, "data.frame")
  expect_type(result$I.proj, "data.frame")

  # Check that key outputs are reasonable
  expect_true(result$p.base == 0.01)
  expect_true(result$p.proj == 0.001)
  expect_true(result$Asym %in% c("Asymptotic dependence", "Asymptotic independence"))

  # Check for reasonable data frame sizes (more flexible than exact length)
  expect_true(nrow(result$I.base) > 0)
  expect_true(nrow(result$I.proj) > 0)
  expect_true(ncol(result$I.base) >= 2)
  expect_true(ncol(result$I.proj) >= 2)

  # Check for no unexpected NA values in key outputs
  expect_false(all(is.na(result$Chi)))
  expect_false(all(is.na(result$n.bar)))
  expect_false(is.na(result$Asym))

})

test_that("Invalid inputs", {

  expect_error(
    Cooley19(Data=na.omit(S20.Detrend.df[1:8,2:3]),Migpd=S20.GPD,
             p.base=0.01,p.proj=0.001,PLOT=TRUE,x_lim_max_T=500,y_lim_max_T=500),
    "Data must have at least 10 observations for reliable analysis.")

  expect_error(
    Cooley19(Data=na.omit(S20.Detrend.df[,2:3]),Migpd=S20.GPD,
             p.base=0.01,p.proj=5,PLOT=TRUE,x_lim_max_T=500,y_lim_max_T=500),
    "p.proj must be a single numeric value between 0 and 1 \\(exclusive\\).")

  expect_error(
    Cooley19(Data=na.omit(S20.Detrend.df[,2:3]),
             Migpd=S20.GPD,
             p.base=0.001,
             p.proj=0.01,
             PLOT=TRUE,
             x_lim_max_T=500,
             y_lim_max_T=500),
    "p.proj must be smaller than p.base for extrapolation.")
})
