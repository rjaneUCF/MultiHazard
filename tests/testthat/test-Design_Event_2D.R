S22.Rainfall<-Con_Sampling_2D(Data_Detrend=S22.Detrend.df[,-c(1,4)],
                              Data_Declust=S22.Detrend.Declustered.df[,-c(1,4)],
                              Con_Variable="Rainfall",u=0.97)
S22.OsWL<-Con_Sampling_2D(Data_Detrend=S22.Detrend.df[,-c(1,4)],
                          Data_Declust=S22.Detrend.Declustered.df[,-c(1,4)],
                          Con_Variable="OsWL",u=0.97)
S22.Copula.Rainfall<-Copula_Threshold_2D(Data_Detrend=S22.Detrend.df[,-c(1,4)],
                                         Data_Declust=S22.Detrend.Declustered.df[,-c(1,4)],u1 =0.97,
                                         y_lim_min=-0.075,y_lim_max=0.25,
                                         Upper=c(2,9),Lower=c(2,10),GAP=0.15)$Copula_Family_Var1
S22.Copula.OsWL<-Copula_Threshold_2D(Data_Detrend=S22.Detrend.df[,-c(1,4)],
                                     Data_Declust=S22.Detrend.Declustered.df[,-c(1,4)],u2 =0.97,
                                     y_lim_min=-0.075, y_lim_max =0.25,
                                     Upper=c(2,9),Lower=c(2,10),GAP=0.15)$Copula_Family_Var2

test_that("Design_Event_2D works", {
  
  result <- Design_Event_2D(Data=S22.Detrend.df[,-c(1,4)],
                            Data_Con1=S22.Rainfall$Data, Data_Con2=S22.OsWL$Data,
                            u1=0.97, u2=0.97,
                            Copula_Family1=S22.Copula.Rainfall, Copula_Family2=S22.Copula.OsWL,
                            Marginal_Dist1="Logis", Marginal_Dist2="Twe",
                            RP=100,Interval=10000,N=10^4,N_Ensemble=10,
                            Plot_Quantile_Isoline=FALSE)
  
  # Checking type of output
  expect_type(result, 'list')
  expect_named(result, c("FullDependence", "MostLikelyEvent", "Ensemble", "Isoline","Contour",
                         "Quantile_Isoline_1","Quantile_Isoline_2","Threshold_1","Threshold_2"))
  
  
  #Checking length of outputs
  expect_true(nrow(result$FullDependence$`100`) == 1 & ncol(result$FullDependence$`100`)==2)
  expect_true(nrow(result$MostLikelyEvent$`100`) == 1 & ncol(result$MostLikelyEvent$`100`)==2)
  expect_equal(nrow(result$Ensemble$`100`),10)
  expect_equal(nrow(result$Isoline$`100`),10000)
  expect_equal(length(result$Contour$`100`),10000)
  expect_equal(nrow(result$Quantile_Isoline_1$`100`),1429)
  expect_equal(nrow(result$Quantile_Isoline_2$`100`),1212)
  expect_equal(length(result$Threshold_1),1)
  expect_equal(length(result$Threshold_2),1)
  
  #Checking column names of outputs
  expect_equal(colnames(S20.Detrend.df[,2:4]), colnames(result$FullDependence$`100`))
  expect_equal(colnames(S20.Detrend.df[,2:4]), colnames(result$MostLikelyEvent$`100`))
  expect_equal(colnames(S20.Detrend.df[,2:4]), colnames(result$Ensemble$`100`))
  expect_equal(colnames(S20.Detrend.df[,2:4]), colnames(result$Isoline$`100`))
  expect_equal(colnames(S20.Detrend.df[,2:4]), colnames(result$Quantile_Isoline_1$`100`))
  expect_equal(colnames(S20.Detrend.df[,2:4]), colnames(result$Quantile_Isoline_1$`100`))
})

test_that("Invalid inputs", {
  
  expect_error(
    Design_Event_2D(Data=S22.Detrend.df[,-c(1,4)],
                    Data_Con1=S22.Rainfall$Data, Data_Con2=S22.OsWL$Data,
                    u1=0.97, u2=0.97,
                    Copula_Family1=S22.Copula.Rainfall, Copula_Family2=S22.Copula.OsWL,
                    Marginal_Dist1="Logis", Marginal_Dist2="Twe",
                    RP=100,Interval=10000,N=10^4,N_Ensemble=10,
                    Isoline_Type = "Null",
                    Plot_Quantile_Isoline=FALSE),
    "Isoline_Type must be one of: Combined, Con1, Con2.")
  
  expect_error(
    Design_Event_2D(Data=S22.Detrend.df[,-c(1,4)],
                    Data_Con1=S22.Rainfall$Data, Data_Con2=S22.OsWL$Data,
                    u1=0.97, u2=0.97,
                    Copula_Family1=S22.Copula.Rainfall, Copula_Family2=25,
                    Marginal_Dist1="Logis", Marginal_Dist2="Twe",
                    RP=100,Interval=10000,N=10^4,N_Ensemble=10,
                    Plot_Quantile_Isoline=FALSE),
    "Invalid Copula_Family2.")
})
