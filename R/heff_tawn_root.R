#' @noRd
heff_tawn_root = function(y,prob,Ypar,YZ,nsim){ #function for uniroot, aiming to find the point on the return curve where y=x
  q = 1 - Laplace_cdf(y)
  sample_y = rexp(nsim) + y
  sample_z = sample(YZ,nsim,replace = T)
  sample_x = (sample_y^Ypar[2])*sample_z + Ypar[1]*sample_y
  x = quantile(sample_x,(1-prob/q)) 
  return(x-y)
}