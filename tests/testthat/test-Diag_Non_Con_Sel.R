#Test that functions work as intended

S20.Rainfall<-Con_Sampling_2D(Data_Detrend=S20.Detrend.df[,-c(1,4)],
                              Data_Declust=S20.Detrend.Declustered.df[,-c(1,4)],
                              Con_Variable="Rainfall",Thres=0.97)



# Test basic functionality for both distributions
test_that("Basic functionality works for all distributions", {

  # Test Gaussian
  expect_no_error(Diag_Non_Con_Sel(Data = S20.Rainfall$Data$OsWL,
                                 x_lab = "O-sWL (ft NGVD 29)",
                                 y_lim_min = 0, y_lim_max = 1.5,
                                 Selected = "Gaus"))

  # Test Gumbel
  expect_no_error(Diag_Non_Con_Sel(Data = S20.Rainfall$Data$OsWL,
                                 x_lab = "O-sWL (ft NGVD 29)",
                                 y_lim_min = 0, y_lim_max = 1.5,
                                 Selected = "Gum"))

  # Test Laplace
  expect_no_error(Diag_Non_Con_Sel(Data = S20.Rainfall$Data$OsWL,
                                 x_lab = "O-sWL (ft NGVD 29)",
                                 y_lim_min = 0, y_lim_max = 1.5,
                                 Selected = "Lapl"))
  # Test Logistic
  expect_no_error(Diag_Non_Con_Sel(Data = S20.Rainfall$Data$OsWL,
                                 x_lab = "O-sWL (ft NGVD 29)",
                                 y_lim_min = 0, y_lim_max = 1.5,
                                 Selected = "Logis"))

  # Test Revere Gumbel
  expect_no_error(Diag_Non_Con_Sel(Data = S20.Rainfall$Data$OsWL,
                                 x_lab = "O-sWL (ft NGVD 29)",
                                 y_lim_min = 0, y_lim_max = 1.5,
                                 Selected = "RGum"))
})

# Test NA handling
test_that("Handles NA values", {

  #data_with_na <- c(S20.Rainfall$Data$OsWL, NA, NA, NA)

  #expect_warning(Diag_Non_Con_Sel(Data = data_with_na,x_lab="O-sWL (ft NGVD 29)",
  #                            y_lim_min=0,y_lim_max=1.5, Omit = c("Gum","RGum"), Selected="Logis"),
  #               "Removed 3 NA values from Data",
  #               all = FALSE)


  # Test all NA data - should error
  expect_error(Diag_Non_Con_Sel(Data = rep(NA, 20), x_lab = "O-sWL (ft NGVD 29)", Selected = "Logis"),
               "Data must be numeric, got: logical")
})

# Test that invalid inputs gives errors
test_that("Invalid inputs produce errors", {

  expect_error(Diag_Non_Con_Sel(Data=numeric(0),x_lab="O-sWL (ft NGVD 29)",
                                y_lim_min=0,y_lim_max=1.5, Selected= "Logis"),
               "Data is empty")

  expect_error(Diag_Non_Con_Sel(Data=S20.Rainfall$Data$OsWL[1:5],x_lab="O-sWL (ft NGVD 29)",
                                y_lim_min=0,y_lim_max=1.5, Selected= "Logis"),
               "Data must have at least 10 non-missing observations, got: 5")

  expect_error(Diag_Non_Con_Sel(Data=S20.Rainfall$Data$OsWL,x_lab="O-sWL (ft NGVD 29)",
                                y_lim_min=0,y_lim_max=-1.5, Selected= "Logis"),
               "y_lim_min must be less than y_lim_max, got: y_lim_min = 0, y_lim_max = -1.5")

  expect_error(Diag_Non_Con_Sel(Data=S20.Rainfall$Data$OsWL,x_lab="O-sWL (ft NGVD 29)",
                                y_lim_min=0,y_lim_max=1.5, Selected= "Exponential"),
               "Invalid distribution names in Selected: Exponential. Valid options are: Gaus, Gum, Lapl, Logis, RGum")

  #data_with_inf <- c(S20.Rainfall$Data$OsWL, Inf, -Inf)
  #expect_warning(Diag_Non_Con_Sel(Data = data_with_inf, x_lab = "O-sWL (ft NGVD 29)", Selected= "Logis"),
  #               "Data contains 2 infinite values. Removing them.",
  #                fixed = TRUE,
  #               all = FALSE)

  expect_error(Diag_Non_Con_Sel(Data=S20.Rainfall$Data$OsWL,x_lab="O-sWL (ft NGVD 29)",
                                Selected= c("Gaus", "Logis"), y_lim_min=0,y_lim_max=1.5),
               "Cannot select more than one distribution.")

    expect_error(Diag_Non_Con_Sel(Data = S20.Rainfall$Data$OsWL,
                                  x_lab = "O-sWL (ft NGVD 29)"),
                 "argument \"Selected\" is missing")

    expect_error(Diag_Non_Con_Sel(Data=S20.Rainfall$Data$OsWL,x_lab="O-sWL (ft NGVD 29)",
                              Omit= "Exponential", y_lim_min=0,y_lim_max=1.5, Selected="Logis"),
                 "Invalid distribution names in Omit: Exponential. Valid options are: Gaus, Gum, Lapl, Logis, RGum")

    expect_error(Diag_Non_Con_Sel(Data = S20.Rainfall$Data$OsWL, x_lab = "O-sWL (ft NGVD 29)",
                              Omit = c("Gaus", "Gum", "Lapl", "Logis", "RGum"), Selected="Logis"),
                 "Cannot omit all distributions. At least one distribution must be tested.")
})


