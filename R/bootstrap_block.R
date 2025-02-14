#' Implements a block bootstrap 
#' 
#' This function performs block bootstrapping on a given dataset.
#'
#' @param data Data frame of raw data detrended if necessary. First column should be of the \code{Date}.
#' @param k Numeric vector of length one specifying the block length. Default is \code{14}.
#' @return Dataframe containing th eblock bootstap sample.
#' @export
#' @examples
#' #Boostrap the detrended data at site S-22
#' boot_df = bootstrap_block(S22.Detrend.df[,2:3])
#' #Calculate the mean of the drivers in bootstrapped sample
#' apply(boot_df,2,function(x) mean(x,na.rm=T))
bootstrap_block = function(data,k=14){ #function for performing block bootstrapping
  #data is bivariate dataset
  #k is block length
  n=length(as.matrix(data)[,1])
  data = as.matrix(data)
  no_blocks = ceiling(n/k)
  n_new = no_blocks*k
  new_data = matrix(NA,nrow=n_new,ncol=dim(data)[2])
  indices = 1:(n-k+1)
  start_points = sample(x=indices,size=no_blocks,replace=TRUE)
  for(i in 1:no_blocks){
    new_data[((i-1)*k+1):(i*k),] = data[(start_points[i]:(start_points[i]+k-1)),]
  }
  return(new_data)
}