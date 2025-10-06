#' Groundwater Levels at Well G-3355
#'
#' Time series of daily water elevation at Well G-3355.
#'
#' @format A data frame with 12,340 rows and 2 variables:
#' \describe{
#'   \item{Date}{Date of observation (YYYY-MM-DD)}
#'   \item{Value}{Groundwater elevation (ft NGVD29)}
#' }
#' @source South Florida Water Management District, via
#' [https://www.sfwmd.gov/science-data/dbhydro/](https://www.sfwmd.gov/science-data/dbhydro/)
#'
#' @usage data(G_3355)
"G_3355"

#' Groundwater Levels at Well G-3356
#'
#' Time series of daily water elevation at Well G-3356.
#'
#' @format A data frame with 12,272 rows and 2 variables:
#' \describe{
#'   \item{Date}{Date of observation (YYYY-MM-DD)}
#'   \item{Value}{Groundwater elevation (ft NGVD29)}
#' }
#' @source South Florida Water Management District, via
#' [https://www.sfwmd.gov/science-data/dbhydro/](https://www.sfwmd.gov/science-data/dbhydro/)
#'
#' @usage data(G_3356)
"G_3356"

#' Groundwater Levels at Well G-860
#'
#' Time series of daily water elevation at Well G-860.
#'
#' @format A data frame with 16,507 rows and 2 variables:
#' \describe{
#'   \item{Date}{Date of observation (YYYY-MM-DD)}
#'   \item{Value}{Groundwater elevation (ft NGVD29)}
#' }
#' @source South Florida Water Management District, via
#' [https://www.sfwmd.gov/science-data/dbhydro/](https://www.sfwmd.gov/science-data/dbhydro/)
#'
#' @usage data(G_860)
"G_860"

#' Groundwater Levels at Well G-580A
#'
#' Time series of daily water elevation at Well G-580A.
#'
#' @format A data frame with 12,301 rows and 2 variables:
#' \describe{
#'   \item{Date}{Date of observation (YYYY-MM-DD)}
#'   \item{Value}{Groundwater elevation (ft NGVD29)}
#' }
#' @source South Florida Water Management District, via
#' [https://www.sfwmd.gov/science-data/dbhydro/](https://www.sfwmd.gov/science-data/dbhydro/)
#'
#' @usage data(G_580A)
"G_580A"

#' Atlantic hurricane database (HURricane DATa 2nd generation) - HURDAT2)
#'
#' "Best-track" data at a six-hourly resolution for all known tropical cyclones and subtropical cyclones available from the b-decks in the Automated Tropical Cyclone Forecast (ATCF – Sampson and Schrader 2000) system database.
#'
#' @format A data frame with columns parsed from the fixed-width HURDAT2 format:
#' \describe{
#'    \item{V1}{See \url{https://www.nhc.noaa.gov/data/hurdat/hurdat2-format-atl-1851-2021.pdf} for description.}
#' }
#' @source National Hurricane Center, via
#' [https://www.nhc.noaa.gov/data/](https://www.nhc.noaa.gov/data/)
#'
#' @usage data(HURDAT2)
"HURDAT2"

#' Rainfall Totals at Miami International Airport
#'
#' Time series of daily rainfall totals from the gauge at Miami International Airport, FL (Network:ID	GHCND:USW00012839).
#'
#' @format A data frame with 25,959 rows and 2 variables:
#' \describe{
#'   \item{X}{Index}
#'   \item{Date}{Date of observation (YYYY-MM-DD)}
#'   \item{Value}{Rainfall totals (inches)}
#' }
#' @source NOAA National Centers for Environmental Information, via
#' [https://www.ncdc.noaa.gov/](https://www.ncdc.noaa.gov/)
#'
#' @usage data(Miami_Airport_df)
"Miami_Airport_df"

#' 2012 NOAA Sea Level Rise Projections for Miami Beach
#'
#' Sea level rise projections for Miami Beach under four emission scenarios from NOAA Technical Report OAR CPO-1 (2012).
#'
#' @format A data frame containing the projections with columns as follows:
#' \describe{
#'   \item{Year}{Year of Projection (YYYY)}
#'   \item{Low}{Projections for the low emission scenario (m)}
#'   \item{Int.Low}{Projections for the intermediate low emission scenario (m)}
#'   \item{Int.High}{Projections for the intermediate high emission scenario (m)}
#'   \item{High}{Projections for the high emission scenario (m)}
#' }
#' @source NOAA National Ocean Service. Report available at:
#' [https://repository.library.noaa.gov/view/noaa/11124](https://repository.library.noaa.gov/view/noaa/11124)
#'
#' @usage data(NOAAetal2012)
"NOAAetal2012"

#' 2017 NOAA Sea Level Rise Projections for Miami Beach
#'
#' Sea level rise projections for Miami Beach derived from NOAA Technical Report NOS CO-OPS 083 (2017). The report defines six global mean sea level rise scenarios, which are regionally downscaled on a 1-degree grid covering U.S. coastlines.
#'
#' @format A data frame containing the projections with columns as follows:
#' \describe{
#'   \item{Year}{Year of Projection (YYYY)}
#'   \item{VLM}{Virtical land motion (m)}
#'   \item{Low}{Projections for the low emission scenario (m)}
#'   \item{Int.Low}{Projections for the intermediate low emission scenario (m)}
#'   \item{Intermediate}{Projections for the intermediate emission scenario (m)}
#'   \item{Int.High}{Projections for the intermediate high emission scenario (m)}
#'   \item{High}{Projections for the high emission scenario (m)}
#'   \item{Extreme}{Projections for the extreme emission scenario (m)}
#' }
#' @source NOAA National Ocean Service. Report available at: [https://tidesandcurrents.noaa.gov/publications/techrpt83_Global_and_Regional_SLR_Scenarios_for_the_US_final.pdf](https://tidesandcurrents.noaa.gov/publications/techrpt83_Global_and_Regional_SLR_Scenarios_for_the_US_final.pdf)
#'
#' @usage data(NOAAetal2017)
"NOAAetal2017"

#' 2022 NOAA Global and Regional Sea Level Rise Scenarios for the United States
#'
#' Probabilistic sea level rise scenarios for Miami Beach from given in NOAA Technical Report NOS 01. In the report, five global mean sea level rise scenarios are used to derive probabilistic regional RSL responses on a 1-degree grid covering the coastlines of the U.S. mainland.
#'
#' @format A data frame containing the projections with columns as follows:
#' \describe{
#'   \item{psmsl_id}{PSML identifier}
#'   \item{process}{Process e.g total sea level}
#'   \item{Units}{Measurement units (mm)}
#'   \item{scenario}{Emission scenario (e.g., Low, Int_Medium, High)}
#'   \item{quantile}{Quantile of the projection distribution (e.g., 0.05, 0.50, 0.95)}
#'   \item{X2020}{Projected sea level rise in 2020 (mm)}
#'   \item{X2030}{Projected sea level rise in 2030 (mm)}
#'   \item{X2040}{...}
#'   \item{X2050}{...}
#'   \item{X2060}{...}
#'   \item{X2070}{...}
#'   \item{X2080}{...}
#'   \item{X2090}{...}
#'   \item{X2100}{...}
#'   \item{X2110}{...}
#'   \item{X2120}{...}
#'   \item{X2130}{...}
#'   \item{X2140}{...}
#'   \item{X2150}{...}
#' }
#' @source NOAA National Ocean Service, via [https://oceanservice.noaa.gov/hazards/sealevelrise](https://oceanservice.noaa.gov/hazards/sealevelrise)
#'
#' @usage data(NOAAetal2022)
"NOAAetal2022"

#' Automatic Peaks Over Threshold (POT) Threshold Selection
#'
#' P-values used in the Automatic Peaks Over Threshold (POT) threshold selection procedure in Solari et al. (2017) for a subset of shape and scale parameter values.
#'
#' @format A data frame or matrix of p-values, where rows correspond to XX and columns to XX.
#' @source Supplied by the authors of Solari et al. (2017): \doi{10.1002/2016WR019426}
#'
#' @references Solari, et al. (2017). [https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2016WR019426](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2016WR019426)
"PVAL_AU2_LMOM_1"

#' Automatic Peaks Over Threshold (POT) Threshold Selection
#'
#' P-values used in the Automatic Peaks Over Threshold (POT) threshold selection procedure in Solari et al. (2017) for a subset of shape and scale parameter values.
#'
#' @format A data frame or matrix of p-values, where rows correspond to XX and columns to XX.
#' @source Supplied by the authors of Solari et al. (2017): \doi{10.1002/2016WR019426}
#'
#' @references Solari, et al. (2017). [https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2016WR019426](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2016WR019426)
"PVAL_AU2_LMOM_2"

#' Automatic Peaks Over Threshold (POT) Threshold Selection
#'
#' P-values used in the Automatic Peaks Over Threshold (POT) threshold selection procedure in Solari et al. (2017) for a subset of shape and scale parameter values.
#'
#' @format A data frame or matrix of p-values, where rows correspond to XX and columns to XX.
#' @source Supplied by the authors of Solari et al. (2017): \doi{10.1002/2016WR019426}
#'
#' @references Solari, et al. (2017). [https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2016WR019426](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2016WR019426)
"PVAL_AU2_LMOM_3"

#' Automatic Peaks Over Threshold (POT) Threshold Selection
#'
#' P-values used in the Automatic Peaks Over Threshold (POT) threshold selection procedure in Solari et al. (2017) for a subset of shape and scale parameter values.
#'
#' @format A data frame or matrix of p-values, where rows correspond to XX and columns to XX.
#' @source Supplied by the authors of Solari et al. (2017): \doi{10.1002/2016WR019426}
#'
#' @references Solari, et al. (2017). [https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2016WR019426](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2016WR019426)
"PVAL_AU2_LMOM_4"

#' Rainfall Totals at Perrine, FL.
#'
#' Time series of daily rainfall totals from the gauge at Perrine 4W, FL (Network:ID GHCND:USC00087020).
#'
#' @format A data frame with 22,100 rows and 2 variables:
#' \describe{
#'   \item{Date}{Date of observation (YYYY-MM-DD)}
#'   \item{Value}{Rainfall totals (inches)}
#' }
#' @source NOAA National Centers for Environmental Information, via
#' [https://www.ncdc.noaa.gov/](https://www.ncdc.noaa.gov/)
#'
#' @usage data(Perrine_df)
"Perrine_df"

#' Measurement from Control Structure S-13
#'
#' Hourly time series of rainfall totals and ocean-side water levels (tailwater levels) measured at control structure S-13. The O-sWLs are detrended using a 3-month moving window.
#'
#' @format A data frame with 281,846 rows and 3 variables:
#' \describe{
#'   \item{Date_Time}{Time of observation (YYYY-MM-DD HH:MM:SS)}
#'   \item{Rainfall}{Rainfall totals (inches)}
#'   \item{OsWL}{Ocean-side water levels totals (ft NGVD29)}
#' }
#' @source South Florida Water Management District, via
#' [https://www.sfwmd.gov/science-data/dbhydro/](https://www.sfwmd.gov/science-data/dbhydro/)
#'
#' @usage data(S13.Detrend.df)
"S13.Detrend.df"


#' Rainfall Totals at S-13
#'
#' Hourly time series of rainfall totals from the gauge at control structure S-13.
#'
#' @format A data frame with 281,846 rows and 2 variables:
#' \describe{
#'   \item{Date_Time}{Time of observation (YYYY-MM-DD HH:MM:SS)}
#'   \item{Rainfall}{Rainfall totals (inches)}
#' }
#' @source South Florida Water Management District, via
#' [https://www.sfwmd.gov/science-data/dbhydro/](https://www.sfwmd.gov/science-data/dbhydro/)
#'
#' @usage data(S13_Rainfall)
"S13_Rainfall"

#' Declustered Time Series for Case Study Site S-20 in Jane et al. (2020)
#'
#' Time series of rainfall totals from the Perrine 4W gauge, ocean-side water levels at control structure S-20, and groundwater at Well G-3356 (with missing values interpolated using groundwater levels at Well G-3355).
#' All time series were declustered using a peaks-over-threshold approach (runs method with a 0.98 quantile threshold and a 3-day separation criterion).
#' Ocean-side water levels and groundwater levels were each detrended using a 3-month moving average window prior to declustering.
#'
#' @format A data frame with 22,100 rows and 4 variables:
#' \describe{
#'   \item{Date}{Date of observation (YYYY-MM-DD)}
#'   \item{Rainfall}{Rainfall totals (inches)}
#'   \item{OsWL}{Ocean-side water levels (ft NGVD29)}
#'   \item{Groundwater}{Groundwater levels (ft NGVD29)}
#' }
#'
#' @references Jane, R., Cadavid, L., Obeysekera, J., and Wahl, T.: Multivariate statistical modelling of the drivers of compound flood events in south Florida, Nat. Hazards Earth Syst. Sci., 20, 2681–2699. \url{https://doi.org/10.5194/nhess-20-2681-2020}
#'
#' @source Rainfall totals from the NOAA National Centers for Environmental Information: \url{https://www.ncdc.noaa.gov/}
#' Ocean-side water levels and groundwater levels from the South Florida Water Management District: \url{https://www.sfwmd.gov/science-data/dbhydro/}
#'
#' @usage data(S20.Detrend.Declustered.df)
"S20.Detrend.Declustered.df"

#' Time Series for Case Study Site S-20 in Jane et al. (2020)
#'
#' Time series of rainfall totals from the Perrine 4W gauge, ocean-side water levels at control structure S-20, and groundwater at Well G-3356 (with missing values interpolated using groundwater levels at Well G-3355).
#' Ocean-side water levels and groundwater levels were each detrended using a 3-month moving average window.
#'
#' @format A data frame with 22,100 rows and 4 variables:
#' \describe{
#'   \item{Date}{Date of observation (YYYY-MM-DD)}
#'   \item{Rainfall}{Rainfall totals (inches)}
#'   \item{OsWL}{Ocean-side water levels (ft NGVD29)}
#'   \item{Groundwater}{Groundwater levels (ft NGVD29)}
#' }
#'
#' @references Jane, R., Cadavid, L., Obeysekera, J., and Wahl, T.: Multivariate statistical modelling of the drivers of compound flood events in south Florida, Nat. Hazards Earth Syst. Sci., 20, 2681–2699. \url{https://doi.org/10.5194/nhess-20-2681-2020}
#'
#' @source Rainfall totals from the NOAA National Centers for Environmental Information: \url{https://www.ncdc.noaa.gov/}
#' Ocean-side water levels and groundwater levels from the South Florida Water Management District: \url{https://www.sfwmd.gov/science-data/dbhydro/}
#'
#' @usage data(S20.Detrend.df)
"S20.Detrend.df"

#' Declustered Time Series of the Ocean-side Water Level at control structure S-20
#'
#' Declustered Time Series of the Ocean-side Water Level at control structure S-20.
#' Declustering is implemented using a peaks-over-threshold approach (runs method with a 0.98 quantile threshold and a 3-day separation criterion).
#'
#' @format A data frame with 18,320 rows and 5 variables:
#' \describe{
#'   \item{X}{Index}
#'   \item{Date}{Date of observation (YYYY-MM-DD)}
#'   \item{Value}{Raw time series (ft NGVD29)}
#'   \item{ValueFilled}{Raw times series with missing values in-filled using record at structure S-21A_T (ft NGVD29)}
#'   \item{Detrend}{Detrended time series (ft NGVD29)}
#'   \item{Declustered}{Declustered time series (ft NGVD29)}
#' }
#'
#' @source South Florida Water Management District: \url{https://www.sfwmd.gov/science-data/dbhydro/}
#'
#' @usage data(S20_T_MAX_Daily_Completed_Detrend_Declustered)
"S20_T_MAX_Daily_Completed_Detrend_Declustered"

#' Time Series of the Ocean-side Water Level at control structure S-22
#'
#' Time Series of the Ocean-side Water Level at control structure S-22.
#'
#' @format A data frame with 12,137 rows and 4 variables:
#' \describe{
#'   \item{Date}{Date of observation (YYYY-MM-DD)}
#'   \item{ValueFilled}{Raw times series with missing values in-filled using record at structure S-21A_T (ft NGVD29)}
#'   \item{Detrend}{Detrended time series (ft NGVD29)}
#' }
#'
#' @source South Florida Water Management District: \url{https://www.sfwmd.gov/science-data/dbhydro/}
#'
#' @usage data(S22_T_MAX_Daily_Completed_Detrend)
"S22_T_MAX_Daily_Completed_Detrend"

#' Declustered Time Series for Case Study Site S-22 in Jane et al. (2020)
#'
#' Time series of rainfall totals from the Miami International Airport gauge, ocean-side water levels at control structure S-22, and groundwater at Well G-580A (with missing values interpolated using groundwater levels at Well G-860).
#' All time series were declustered using a peaks-over-threshold approach (runs method with a 0.98 quantile threshold and a 3-day separation criterion).
#' Ocean-side water levels and groundwater levels were each detrended using a 3-month moving average window prior to declustering.
#'
#' @format A data frame with 26,067 rows and 4 variables:
#' \describe{
#'   \item{Date}{Date of observation (YYYY-MM-DD)}
#'   \item{Rainfall}{Rainfall totals (inches)}
#'   \item{OsWL}{Ocean-side water levels (ft NGVD29)}
#'   \item{Groundwater}{Groundwater levels (ft NGVD29)}
#' }
#'
#' @references Jane, R., Cadavid, L., Obeysekera, J., and Wahl, T.: Multivariate statistical modelling of the drivers of compound flood events in south Florida, Nat. Hazards Earth Syst. Sci., 20, 2681–2699. \url{https://doi.org/10.5194/nhess-20-2681-2020}
#'
#' @source Rainfall totals from the NOAA National Centers for Environmental Information: \url{https://www.ncdc.noaa.gov/}
#' Ocean-side water levels and groundwater levels from the South Florida Water Management District: \url{https://www.sfwmd.gov/science-data/dbhydro/}
#'
#' @usage data(S22.Detrend.Declustered.df)
"S22.Detrend.Declustered.df"

#' Time Series for Case Study Site S-22 in Jane et al. (2020)
#'
#' Time series of rainfall totals from the Miami International Airport gauge, ocean-side water levels at control structure S-22, and groundwater at Well G-580A (with missing values interpolated using groundwater levels at Well G-860).
#' Ocean-side water levels and groundwater levels were each detrended using a 3-month moving average window.
#'
#' @format A data frame with 25,992 rows and 4 variables:
#' \describe{
#'   \item{Date}{Date of observation (YYYY-MM-DD)}
#'   \item{Rainfall}{Rainfall totals (inches)}
#'   \item{OsWL}{Ocean-side water levels (ft NGVD29)}
#'   \item{Groundwater}{Groundwater levels (ft NGVD29)}
#' }
#'
#' @references Jane, R., Cadavid, L., Obeysekera, J., and Wahl, T.: Multivariate statistical modelling of the drivers of compound flood events in south Florida, Nat. Hazards Earth Syst. Sci., 20, 2681–2699. \url{https://doi.org/10.5194/nhess-20-2681-2020}
#'
#' @source Rainfall totals from the NOAA National Centers for Environmental Information: \url{https://www.ncdc.noaa.gov/}
#' Ocean-side water levels and groundwater levels from the South Florida Water Management District: \url{https://www.sfwmd.gov/science-data/dbhydro/}
#'
#' @usage data(S22.Detrend.df)
"S22.Detrend.df"

#' Sea Level Projections for Fort Myers
#'
#' Probabilistic sea level rise projections for Fort Myers (PSMSL ID 1106) derived from the Interagency Sea Level Rise Scenario Tool.
#' Data are based on PSMSL records and include multiple emission scenarios (e.g., Low, Intermediate, High), each represented as quantile estimates.
#'
#' @format A data frame containing the projections with columns as follows:
#' \describe{
#'   \item{psmsl_id}{PSML identifier}
#'   \item{process}{Process e.g total sea level}
#'   \item{Units}{Measurement units (mm)}
#'   \item{scenario}{Emission scenario (e.g., Low, Int_Medium, High)}
#'   \item{quantile}{Quantile of the projection distribution (e.g., 0.05, 0.50, 0.95)}
#'   \item{X2020}{Projected sea level rise in 2020 (mm)}
#'   \item{X2030}{Projected sea level rise in 2030 (mm)}
#'   \item{X2040}{...}
#'   \item{X2050}{...}
#'   \item{X2060}{...}
#'   \item{X2070}{...}
#'   \item{X2080}{...}
#'   \item{X2090}{...}
#'   \item{X2100}{...}
#'   \item{X2110}{...}
#'   \item{X2120}{...}
#'   \item{X2130}{...}
#'   \item{X2140}{...}
#'   \item{X2150}{...}
#' }
#' @source Interagency Sea Level Rise Scenario Tool (PSMSL ID 1106 – Fort Myers), accessed via the National Sea Level Explorer.
#'
#' @usage data(sl_taskforce_scenarios_psmsl_id_1106_Fort_Myers)
"sl_taskforce_scenarios_psmsl_id_1106_Fort_Myers"

#' USACE Sea Level Projections for Virginia Key
#'
#' U.S. Army Corps. of Engineers (USACE) 2013 sea level rise projections for Virginia Key, as provided in Engineering Regulation ER 1100-2-8162.
#' These projections represent future sea level rise scenarios, resulting in global sea level increases of approximately 0.2 meters (Low scenario), 0.5 meters (Intermediate scenario), and 1.5 meters (High scenario) by the year 2100.
#'
#' @format A data frame containing the projections with columns as follows:
#' \describe{
#'   \item{Year}{Year of Projection (YYYY)}
#'   \item{Low}{Projections for the low emission scenario (m)}
#'   \item{Int}{Projections for the intermediate emission scenario (m)}
#'   \item{High}{Projections for the high emission scenario (m)}
#' }
#' @source U.S. Army Corps. of Engineers, via. USACE Sea Level Change Calculator: \url{http://www.corpsclimate.us/ccaceslcurves.cfm}
#'
#' @usage data(USACE2013)
"USACE2013"
