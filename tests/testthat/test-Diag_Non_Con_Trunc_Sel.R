S20.OsWL<-Con_Sampling_2D(Data_Detrend=S20.Detrend.df[,-c(1,4)],   
                          Data_Declust=S20.Detrend.Declustered.df[,-c(1,4)],
                          Con_Variable="OsWL",Thres=0.97)

S20.OsWL$Data$Rainfall <- S20.OsWL$Data$Rainfall + runif(length(S20.OsWL$Data$Rainfall),0.001,0.01)       

# Test basic functionality for both distributions
# Group distribution tests
test_that("All distributions work correctly", {
  distributions <- c("BS", "Exp", "Gam(2)", "Gam(3)", "GamMix(2)", "GamMix(3)", "LNorm", "TNorm", "Twe", "Weib")
  for(dist in distributions) {
    expect_silent(Diag_Non_Con_Trunc_Sel(Data = S20.OsWL$Data$Rainfall, 
                                         x_lab = "Rainfall (Inches)", 
                                         y_lim_min = 0, y_lim_max = 1.5, 
                                         Selected = dist))
  }
})

# Test NA handling
test_that("Handles NA values", {
  
  data_with_na <- c(S20.OsWL$Data$Rainfall , NA, NA, NA)
  
  expect_warning(Diag_Non_Con_Trunc_Sel(Data = data_with_na,x_lab="Rainfall (Inches)",
                                        y_lim_min=0,y_lim_max=1.5, Selected="Twe"),
                 "Removed 3 NA values from Data")
  
  
  # Test all NA data - should error
  expect_error(Diag_Non_Con_Trunc_Sel(Data = rep(NA, 20), x_lab = "Rainfall (Inches)", Selected="Twe"),
               "Data must be numeric, got: logical")
})

# Test that invalid inputs gives errors
test_that("Invalid inputs produce errors", {
  
  expect_error(Diag_Non_Con_Trunc_Sel(Data=numeric(0),x_lab="Rainfall (Inches)",
                                  y_lim_min=0,y_lim_max=-1.5, Selected="Twe"),
               "Data is empty")
  
  expect_error(Diag_Non_Con_Trunc_Sel(Data=S20.OsWL$Data$Rainfall [1:5],x_lab="Rainfall (Inches)",
                                  y_lim_min=0,y_lim_max=1.5, Selected="Twe"),
               "Data must have at least 10 non-missing observations, got: 5")
  
  expect_error(Diag_Non_Con_Trunc_Sel(Data=S20.OsWL$Data$Rainfall ,x_lab="Rainfall (Inches)",
                                  y_lim_min=0,y_lim_max=-1.5, Selected="Twe"),
               "y_lim_min must be less than y_lim_max, got: y_lim_min = 0, y_lim_max = -1.5")
  
  expect_error(Diag_Non_Con_Trunc_Sel(Data=S20.OsWL$Data$Rainfall ,x_lab="Rainfall (Inches)",
                                  Omit= "Gaussian", y_lim_min=0,y_lim_max=1.5, Selected="Twe"),
               "Invalid distribution names in Omit: Gaussian")
  
  data_with_inf <- c(S20.OsWL$Data$Rainfall , Inf, -Inf)
  expect_warning(Diag_Non_Con_Trunc_Sel(Data = data_with_inf, x_lab = "Rainfall (Inches)", Selected="Twe"),
                 "Data contains 2 infinite values")
  
  expect_error(Diag_Non_Con_Trunc_Sel(Data = S20.OsWL$Data$Rainfall , x_lab = "Rainfall (Inches)", 
                                  Omit = c("BS","Exp","Gam(2)","Gam(3)","GamMix(2)","GamMix(3)","LNorm","TNorm","Twe","Weib"), Selected="Twe"),
               "Cannot omit all distributions")
  
  expect_error(Diag_Non_Con_Trunc_Sel(Data = S20.OsWL$Data$Rainfall, 
                                x_lab = "Rainfall (Inches)"),
               "argument \"Selected\" is missing")
})

