#' Implements the version of the conditional multivariate approach of Heffernan and Tawn (2004) proposed in Keef et al. (2013) which incorporates lags between the variables.
#'
#' Implements the version of the conditional multivariate approach of Heffernan and Tawn (2004) proposed in Keef et al. (2013) which incorporates lags between the variables. Function utilizes the \code{mexDependence} and \code{predict.mex.conditioned} functions from the \code{texmex} package.
#'
#' @param data_Detrend_Dependence_df  A data frame with (n+1) columns, containing in column \itemize{
#' \item 1 - Continuous sequence of dates spanning the first to the final time of any of the variables are recorded.
#' \item 2:(n+1) - Values, detrended where necessary, of the variables to be modelled.
#' }
#' @param data_Detrend_Declustered_df A data frame with (n+1) columns, containing in column \itemize{
#' \item 1 - Continuous sequence of dates spanning the first to the final time of any of the variables are recorded.
#' \item 2:(n+1) - Declustered and if necessary detrended values of the variables to be modelled.
#' }
#' @param Lags Matrix specifying the lags. The no lag i.e. \code{0} lag cases need to be specified. Row n denotes the lags applied to the variable in the nth column of \code{data_Detrend_Dependence_df}. Column n corresponds to the nth largest lag applied to any variable. Default is \code{matrix(c(NA,0,1,NA,0,1,NA,0,NA),nrow=3,byrow = TRUE)}, which corresponds to a lag of 1 being applied to variables in the first and second columns of \code{data_Detrend_Dependence_df} and no lag being applied to the variable in the third column of \code{data_Detrend_Dependence_df}.
#' @param u_Dependence Dependence quantile. Specifies the (sub-sample of) data to which the dependence model is fitted, that for which the conditioning variable exceeds the threshold associated with the prescribed quantile. Default is \code{0.7}, thus the dependence parameters are estimated using the data with the highest \code{30\%} of values of the conditioning variables.
#' @param Migpd An \code{Migpd} object, containing the parameterized Pareto models fitted (independently) to each of the variables.
#' @param mu Numeric vector of length one specifying the (average) occurrence frequency of events in \code{data_Detrend_Dependence_df}. Default is \code{365.25}, daily data.
#' @param N Numeric vector of length one specifying the number of years worth of extremes to simulate. Default is \code{100} years.
#' @param Margins Character vector specifying the form of margins to which the data are transformed for carrying out dependence estimation. Default is \code{"gumbel"}, alternative is \code{"laplace"}. Under Gumbel margins, the estimated parameters \code{a} and \code{b} describe only positive dependence, while \code{c} and \code{d} describe negative dependence in this case. For Laplace margins, only parameters \code{a} and \code{b} are estimated as these capture both positive and negative dependence.
#' @param V See documentation for \code{mexDependence}.
#' @param Maxit See documentation for \code{mexDependence}.
#' @return List comprising the fitted HT04 models \code{Models}, proportion of the time each variable is most extreme, given at least one variable is extreme \code{Prop}, residuals \code{z}, as well as the simulated values on the transformed \code{u.sim} and original \code{x.sim} scales.
#' @seealso \code{\link{Dataframe_Combine}} \code{\link{Decluster}} \code{\link{GPD_Fit}} \code{\link{Migpd_Fit}}
#' @export
#' @examples
#' #' #Format date column
#' S20.Detrend.Declustered.df$Date =
#'       as.Date(as.character(S20.Detrend.Declustered.df$Date), format = "%m/%d/%Y")
#' #Fit GPD marginal distributions above the threshold
#' S20_GPD<-Migpd_Fit(Data=S20.Detrend.Declustered.df[,-1],
#'                      Data_Full=S20.Detrend.Declustered.df[,-1],
#'                      mqu =c(0.99,0.99,0.99))
#' #Fitting and simulating from the Heffernan and Tawn (2004) model
#' HT04_Lag(data_Detrend_Dependence_df = S20.Detrend.df,
#'      data_Detrend_Declustered_df = S20.Detrend.Declustered.df,
#'      Lags = matrix(c(NA,0,1,NA,0,1,NA,0,NA),nrow=3,byrow = TRUE),
#'      Migpd = S20_GPD, u_Dependence=0.7,Margins = "gumbel")
HT04_Lag<-function (data_Detrend_Dependence_df, data_Detrend_Declustered_df, Lags, u_Dependence, Migpd, mu = 365.25, N = 100, Margins = "gumbel",V = 10, Maxit = 10000){

  # Input Validation

  # Check if data frames exist and are valid
  if(missing(data_Detrend_Dependence_df) || is.null(data_Detrend_Dependence_df) || !is.data.frame(data_Detrend_Dependence_df)){
    stop("data_Detrend_Dependence_df must be a data frame.")
  }

  if(missing(data_Detrend_Declustered_df) || is.null(data_Detrend_Declustered_df) || !is.data.frame(data_Detrend_Declustered_df)){
    stop("data_Detrend_Declustered_df must be a data frame.")
  }

  # Check if data frames have matching dimensions after removing date/factor columns
  temp_dep <- data_Detrend_Dependence_df
  temp_declust <- data_Detrend_Declustered_df

  if(class(temp_dep[,1])[1] %in% c("Date", "factor", "POSIXct", "character")){
    temp_dep <- temp_dep[,-1]
  }
  if(class(temp_declust[,1])[1] %in% c("Date", "factor", "POSIXct", "character")){
    temp_declust <- temp_declust[,-1]
  }

  if(ncol(temp_dep) != ncol(temp_declust)){
    stop("Data frames must have the same number of numeric columns after removing date/factor columns.")
  }

  if(!all(colnames(temp_dep) == colnames(temp_declust))){
    stop("Column names must match between data frames after removing date/factor columns.")
  }

  # Check if Lag is a matrix
  if (!is.matrix(Lags)) {
    stop("Lag must be a matrix.")
  }

  # Check if Lag contains only numeric values and NA
  if (!all(is.numeric(Lags) | is.na(Lags))) {
    stop("Lag must contain only non-negative integer values or NA.")
  }

  # Check for negative values
  if (any(Lags < 0, na.rm = TRUE)) {
    stop("Lag values must be non-negative integers or NA. Negative lag values are not permitted.")
  }

  # Check for non-integer values
  if (any(Lags != round(Lags), na.rm = TRUE)) {
    stop("Lag values must be integers or NA. Non-integer lag values are not permitted.")
  }

  # Check dimension match with data
  if (nrow(Lags) != ncol(temp_dep)) {
    stop("Number of rows in Lags must match number of columns in data_Detrend_Dependence_df after removing any date/factor columns. Row n denotes lags applied to variable in nth column of data.")
  }

  # Validate u_Dependence
  if(missing(u_Dependence) || is.null(u_Dependence) || !is.numeric(u_Dependence) || length(u_Dependence) != 1 || u_Dependence <= 0 || u_Dependence >= 1){
    stop("u_Dependence must be a single numeric value between 0 and 1 (exclusive).")
  }

  # Validate Migpd
  if(missing(Migpd) || is.null(Migpd) || !is.list(Migpd)){
    stop("Migpd must be a list object (migpd class).")
  }

  if(is.null(Migpd$models)){
    stop("Migpd$models is NULL. Please ensure Migpd object contains fitted models.")
  }

  if(length(Migpd$models) != ncol(temp_dep)){
    stop("Number of models in Migpd$models (", length(Migpd$models),
         ") must match number of data columns (", ncol(temp_dep), ").")
  }

  # Validate numeric parameters
  if(!is.numeric(mu)){
    stop("mu must be numeric.")
  }

  if(length(mu) != 1 || mu <= 0){
    stop("mu must be a positive numeric value.")
  }

  if(!is.numeric(N)){
    stop("N must be a numeric.")
  }

  if(length(N) != 1 || N <= 0 || N != round(N)){
    stop("N must be a positive numeric value.")
  }

  if(!is.numeric(V)){
    stop("V must be a numeric.")
  }

  if(length(V) != 1 || V <= 0 || V != round(V)){
    stop("V must be a positive integer.")
  }

  if(!is.numeric(Maxit)){
    stop("Maxit must be a numeric value.")
  }

  if(length(Maxit) != 1 || Maxit <= 0 || Maxit != round(Maxit)){
    stop("Maxit must be a positive integer.")
  }

  # Validate Margins parameter
  valid_margins <- c("gumbel", "laplace", "exponential")
  if(!is.character(Margins) || length(Margins) != 1 || !Margins %in% valid_margins){
    stop("Margins must be one of: ", paste(valid_margins, collapse=", "), ".")
  }

  if(inherits(data_Detrend_Dependence_df[,1], c("Date", "factor", "POSIXct", "character"))){
    data_Detrend_Dependence_df <- data_Detrend_Dependence_df[,-1]
  }

  if(inherits(data_Detrend_Declustered_df[,1], c("Date", "factor", "POSIXct", "character"))){
    data_Detrend_Declustered_df <- data_Detrend_Declustered_df[,-1]
  }

  HT04_Model<-vector('list',ncol(data_Detrend_Declustered_df))
  u<-array(NA,dim=c(nrow(na.omit(data_Detrend_Dependence_df)),ncol(data_Detrend_Declustered_df)))
  Gumbel_df<-array(NA,dim=c(nrow(na.omit(data_Detrend_Dependence_df)),ncol(data_Detrend_Declustered_df)))
  Gumbel_df.Threshold<-rep(NA,ncol(data_Detrend_Declustered_df))
  u.extremes<-array(0,dim=c(nrow(na.omit(data_Detrend_Dependence_df)),ncol(data_Detrend_Declustered_df)))
  Prop<-rep(NA,ncol(data_Detrend_Declustered_df))
  HT04.Predict<-vector('list',ncol(data_Detrend_Declustered_df))
  HT04_z<-vector('list',ncol(data_Detrend_Declustered_df))

  Lag <- function(x, k) {
    if (k > 0) {
      return(c(rep(NA, k), x)[1:length(x)])
    }
    else {
      return(c(x[(-k + 1):length(x)], rep(NA, -k)))
    }
  }

  DF <- matrix( NA_real_,nrow = nrow(data_Detrend_Dependence_df),ncol = sum(!is.na(Lags)) )

  lag_map <- data.frame( df_col = integer(0),variable = character(0), lag = integer(0),name = character(0), stringsAsFactors = FALSE)

  L <- 0

  for (i in seq_len(ncol(data_Detrend_Dependence_df))) {

    variable <- colnames(data_Detrend_Dependence_df)[i]

    for (lag in sort(na.omit(Lags[i, ]))) {

      L <- L + 1

      DF[, L] <- if (lag == 0) {
        data_Detrend_Dependence_df[, i]
      } else {
        Lag(data_Detrend_Dependence_df[, i], lag)
      }

      name <- if (lag == 0) {
        variable
      } else {
        paste0(variable, "_", lag)
      }

      lag_map <- rbind(lag_map, data.frame(df_col = L, variable = variable, lag = lag, name = name, stringsAsFactors = FALSE))
    }
  }

  colnames(DF) <- lag_map$name

  vars <- colnames(data_Detrend_Dependence_df)

  for (i in seq_along(vars)) {

    con.var <- vars[i]

    # All lagged columns except those belonging to conditioning variable
    other_cols <- lag_map$df_col[lag_map$variable != con.var]

    HT04.Dataset <- matrix(NA_real_, nrow = nrow(data_Detrend_Declustered_df),ncol = length(other_cols) + 1)

    HT04.Dataset[, 1] <- data_Detrend_Declustered_df[, con.var]

    if (length(other_cols) > 0) {
      HT04.Dataset[, -1] <- DF[, other_cols]
    }

    colnames(HT04.Dataset) <- c( con.var,lag_map$name[other_cols])

    Migpd_1 <- Migpd

    Migpd_1$models <- Migpd$models[c(i,match(lag_map$variable[other_cols], vars))]

    Migpd_1$mth <- Migpd$mth[c(i,match(lag_map$variable[other_cols], vars))]

    Migpd_1$mqu <- Migpd$mqu[c(i,match(lag_map$variable[other_cols], vars))]

    names(Migpd_1$models) <- colnames(HT04.Dataset)
    names(Migpd_1$mth) <- colnames(HT04.Dataset)
    names(Migpd_1$mqu) <- colnames(HT04.Dataset)

    Migpd_1$data <- na.omit(as.data.frame(HT04.Dataset))

    HT04_Model[[i]] <- mexDependence(Migpd_1, which = con.var, dqu = u_Dependence, margins = Margins, constrain = FALSE, v = V, maxit = Maxit)

    Migpd$data <- na.omit(data_Detrend_Dependence_df)

    u[, i] <- as.numeric(transFun.HT04(x = Migpd$data[, i], mod = Migpd$models[[i]]))

    Gumbel_df[, i] <- -log(-log(u[, i]))

    Gumbel_df.Threshold[i] <-
      quantile(Gumbel_df[, i], u_Dependence)

    u.extremes[
      Gumbel_df[, i] > Gumbel_df.Threshold[i],
      i
    ] <- 1
  }

  Gumbel_df_Extremes<-Gumbel_df[which(apply(u.extremes, 1, sum)>0),]
  colnames(Gumbel_df_Extremes)<-colnames(data_Detrend_Dependence_df)

  for(i in 1:ncol(data_Detrend_Declustered_df)){
    Prop[i]<-length(which(apply(Gumbel_df_Extremes, 1, function(x) which.max(x))==i))/nrow(Gumbel_df_Extremes)
  }

  #Simulations
  n <- setNames(integer(length(vars)),vars)
  for (i in seq_along(vars)) {
    HT04.Predict[[i]]<-predict.mex.conditioned(HT04_Model[[i]], which=colnames(data_Detrend_Dependence_df)[i], pqu = u_Dependence, nsim = round((1-u_Dependence)*mu*N*Prop[i],0), trace=10)
    n[i]<-nrow(HT04.Predict[[i]]$data$transformed)
  }


  HT04_Predict_Transformed_list <- vector("list",length(vars))
  HT04_Predict_Simulated_list <- vector( "list",length(vars) )

  for (con.i in seq_along(vars)) {

    con.var <- vars[con.i]

    n.sim <- nrow(HT04.Predict[[con.i]]$data$transformed)

    HT04_Predict_Transformed_1 <- matrix(NA_real_, nrow = n.sim,  ncol = length(vars), dimnames = list(NULL, vars))
    HT04_Predict_Simulated_1 <- matrix(NA_real_, nrow = n.sim, ncol = length(vars),dimnames = list(NULL, vars))


    # Conditioning variable
    HT04_Predict_Transformed_1[, con.var] <-HT04.Predict[[con.i]]$data$transformed[, 1]
    HT04_Predict_Simulated_1[, con.var] <- HT04.Predict[[con.i]]$data$simulated[, 1]

    row_max <- function(x) {
      apply(x, 1, function(z) {
        if (all(is.na(z))) NA_real_ else max(z, na.rm = TRUE)
      })
    }


    # All non-conditioning variables
    for (other.var in vars[vars != con.var]) {

      # All lagged versions of this variable
      lag_names <- lag_map$name[lag_map$variable == other.var]

      HT04_Predict_Transformed_1[, other.var] <-row_max(HT04.Predict[[con.i]]$data$transformed[, paste0(lag_names, ".trans"), drop = FALSE])
      HT04_Predict_Simulated_1[, other.var] <-row_max(HT04.Predict[[con.i]]$data$simulated[,lag_names, drop = FALSE])
      }


    # Store simulations for this conditioning variable
    HT04_Predict_Transformed_list[[con.i]] <- HT04_Predict_Transformed_1
    HT04_Predict_Simulated_list[[con.i]] <-HT04_Predict_Simulated_1
  }

  # Combine simulations from all conditioning variables

  HT04_Predict_Transformed <- do.call(rbind,HT04_Predict_Transformed_list)
  HT04_Predict_Simulated <- do.call(rbind,HT04_Predict_Simulated_list)

  HT04_Predict_Transformed <- data.frame(HT04_Predict_Transformed)
  HT04_Predict_Simulated <- data.frame(HT04_Predict_Simulated)


  # Store residuals
  for (i in seq_along(data_Detrend_Dependence_df)) {
    HT04_z[[i]] <- HT04.Predict[[i]]$data$z
  }

  # Output
  res <- list(Model = HT04_Model, Prop = Prop, n =n, z = HT04_z,  u.sim = HT04_Predict_Transformed,x.sim = HT04_Predict_Simulated)

  return(res)
}
