#' Implements a monthly bootstrap
#'
#' Months in which at least one variable exceeds the user-specified minimum proportion of non-missing values are sampled with replacement. February of leap years are treated as a 13th month.
#'
#' @param data Data frame of raw data detrended if necessary. First column should be of the \code{Date}.
#' @param boot_prop Numeric vector of length one specifying the minimum proportion of non-missing values of at least one of the variables for a month to be included in the bootstrap. Default is \code{0.8}.
#' @return Dataframe containing a bootstrap undertaken with replacement that accounts for monthly-scale seasonality.
#' @export
#' @examples
#' #Let's assess the sampling variability in kendall's tau
#' #correlation coefficient between rainfall and OsWL at S-22.
#'
#' #Data starts on first day of 1948
#' head(S22.Detrend.Declustered.df)
#'
#' #Dataframe ends on 1948-02-03
#' tail(S22.Detrend.Declustered.df)
#'
#' #Adding dates to complete final month of combined records
#' final.month = data.frame(seq(as.Date("2019-02-04"),as.Date("2019-02-28"),by="day"),NA,NA,NA)
#' colnames(final.month) = c("Date","Rainfall","OsWL","Groundwater")
#' S22.Detrend.Declustered.df = rbind(S22.Detrend.Declustered.df,final.month)
#'
#' #Generate 100 monthly bootstrap samples of rainfall and OsWL
#' cor = rep(NA,100)
#' for(i in 1:100){
#'  boot_df = bootstrap_month(S22.Detrend.df[,c(1:3)], boot_prop=0.8)
#'  boot_df = na.omit(boot_df)
#'  cor[i] = cor(boot_df$Rainfall, boot_df$OsWL, method="kendall")
#' }
#'
#' #Compare means of bootstrap samples with the mean of the observed data
#' hist(cor)
#' df = na.omit(S22.Detrend.df[,1:3])
#' abline(v=cor(df$Rainfall,df$OsWL, method="kendall"),col=2,lwd=2)
bootstrap_month <- function(data, boot_prop = 0.8) {

  # Pre-compute dates once
  data_months <- month(data$Date)
  data_years <- year(data$Date)

  # Create complete month-year grid
  year_range <- min(data_years):max(data_years)
  m.y.complete <- data.frame(
    month = rep(1:12, length(year_range)),
    year = rep(year_range, each = 12)
  )

  # Vectorized indicator calculation using data.table or dplyr (much faster)
  # Create a lookup of month-year combinations in data
  data$month_year <- paste(data_months, data_years, sep = "_")
  m.y.complete$month_year <- paste(m.y.complete$month, m.y.complete$year, sep = "_")

  # Calculate indicators vectorized
  indicator_calc <- function(my) {
    rows <- data$month_year == my
    if (sum(rows) == 0) return(NA)

    d <- data[rows, ]
    n <- nrow(d)

    x_prop <- sum(!is.na(d[, 2]) & is.na(d[, 3])) / n
    y_prop <- sum(is.na(d[, 2]) & !is.na(d[, 3])) / n
    b_prop <- sum(!is.na(d[, 2]) & !is.na(d[, 3])) / n

    if (b_prop > boot_prop) return("B")
    if (x_prop > boot_prop) return("x")
    if (y_prop > boot_prop) return("y")
    return(NA)
  }

  m.y.complete$indicator <- sapply(m.y.complete$month_year, indicator_calc)

  # Handle leap years - vectorized
  leap_years <- seq(1880, 2052, 4)
  leap_years <- leap_years[leap_years >= min(data_years) & leap_years <= max(data_years)]

  # Mark February in leap years as month 13
  leap_feb_idx <- which(m.y.complete$month == 2 & m.y.complete$year %in% leap_years)
  m.y.complete$month[leap_feb_idx] <- 13

  # Bootstrap only for rows with data
  valid_idx <- which(!is.na(m.y.complete$indicator))

  # Vectorized bootstrap sampling by group
  boot_idx <- integer(nrow(m.y.complete))

  for (i in valid_idx) {
    month_i <- m.y.complete$month[i]
    indicator_i <- m.y.complete$indicator[i]

    # Find all matching month-indicator combinations
    candidates <- which(m.y.complete$month == month_i &
                          m.y.complete$indicator == indicator_i)

    boot_idx[i] <- sample(candidates, 1)
  }

  # Get bootstrap sample
  m.y.boot <- m.y.complete[boot_idx, c("month", "year")]

  # Convert month 13 back to 2
  m.y.boot$month[m.y.boot$month == 13] <- 2
  m.y.complete$month[m.y.complete$month == 13] <- 2

  # Create bootstrap dataframe - VECTORIZED APPROACH
  boot_df <- data.frame(
    Date = data$Date,
    V2 = NA,
    V3 = NA
  )
  colnames(boot_df) <- colnames(data[,1:3])

  # Use match for faster lookup instead of loops
  data$orig_month_year <- paste(data_months, data_years, sep = "_")

  for (i in valid_idx) {
    orig_my <- paste(m.y.complete$month[i], m.y.complete$year[i], sep = "_")
    boot_my <- paste(m.y.boot$month[i], m.y.boot$year[i], sep = "_")

    orig_rows <- data$orig_month_year == orig_my
    boot_rows <- data$orig_month_year == boot_my

    boot_df[orig_rows, 2:3] <- data[boot_rows, 2:3]
  }

  # Clean up temporary columns
  data$month_year <- NULL
  data$orig_month_year <- NULL

  return(boot_df)
}



