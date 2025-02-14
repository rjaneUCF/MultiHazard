#' @noRd
ConstraintsAreSatisfied <- function(a,b,z,zpos,zneg,v){
  C1e <- a <= min(1, 1 - b*min(z)*v^(b-1), 1 - v^(b-1)*min(z) + min(zpos)/v) &
    a <= min(1, 1 - b*max(z)*v^(b-1), 1 - v^(b-1)*max(z) + max(zpos)/v)
  
  C1o <- a <= 1 & 
    a > 1 - b * min(z) * v^(b-1) & 
    a > 1 - b * max(z) * v^(b-1) &
    (1 - 1/b)*(b*min(z))^(1/(1-b)) * (1-a)^(-b/(1 - b)) + min(zpos) > 0 &
    (1 - 1/b)*(b*max(z))^(1/(1-b)) * (1-a)^(-b/(1 - b)) + max(zpos) > 0
  
  C2e <- -a <= min(1, 1 + b*v^(b-1)*min(z), 1 + v^(b-1)*min(z) - min(zneg)/v) &
    -a <= min(1, 1 + b*v^(b-1)*max(z), 1 + v^(b-1)*max(z) - max(zneg)/v)
  
  C2o <- -a <= 1 & 
    -a > 1 + b*v^(b-1)*min(z) & 
    -a > 1 + b*v^(b-1)*max(z) &
    (1-1/b)*(-b*min(z))^(1/(1-b))*(1+a)^(-b/(1-b)) - min(zneg) > 0 &
    (1-1/b)*(-b*max(z))^(1/(1-b))*(1+a)^(-b/(1-b)) - max(zneg) > 0
  
  if (any(is.na(c(C1e, C1o, C2e, C2o)))) {
    warning("Strayed into impossible area of parameter space")
    C1e <- C1o <- C2e <- C2o <- FALSE
  }
  
  (C1e | C1o) && (C2e | C2o)
}
