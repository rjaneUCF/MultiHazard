#' @noRd
properties <- function(w, lambda){
  idx <- w <= 0.5
  idxhalf <- which.min(idx)
  lower <- pmax(w, 1 - w)
  viol <- which(lambda < lower)
  if(!identical(viol, integer(0))){
    idxl <- viol[1]
    idxu <- viol[length(viol)]
    if(idxl < idxhalf){
      ind_wl <- max(which((lambda < lower)[1:(idxhalf-1)]))
      lambda[1:ind_wl] <- lower[1:ind_wl]
    }
    if(idxu >= idxhalf){
      ind_wu <- min(which((lambda < lower)[idxhalf:length(w)])) + idxhalf - 1
      lambda[ind_wu:length(w)] <- lower[ind_wu:length(w)]
    }
  }
  for(i in length(w[w < 0.5]):1){
    if((w / lambda)[i + 1] < (w / lambda)[i]){
      lambda[i] <- w[i] * (lambda / w)[i + 1]
    }
    if(((1 - w) / lambda)[i + 1] > ((1 - w) / lambda)[i]){
      lambda[i] <- (1 - w)[i] * (lambda / (1 - w))[i + 1]
    }
  }
  for(i in length(w[w > 0.5]):length(w)){
    if((w / lambda)[i] < (w / lambda)[i - 1]){
      lambda[i] <- w[i] * (lambda / w)[i - 1]
    }
    if(((1 - w) / lambda)[i] > ((1 - w) / lambda)[i - 1]){
      lambda[i] <- (1 - w)[i] * (lambda / (1 - w))[i - 1]
    }
  }
  return(lambda)
}
