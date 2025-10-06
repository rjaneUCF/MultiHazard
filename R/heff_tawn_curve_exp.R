#' @noRd
heff_tawn_curve_exp = function(data_exp_x,data_exp_y,prob,q,nsim){ #function for estimating return curve for data on standard exponential margins using HT model
  #prob is curve survival probability (small)
  #q is quantile cdf probability for fitting HT model
  #nsim is number of simulation for HT model
  nOptim = 25

  #transform data_exp to Laplace margins

  data_u_x = apply(data_exp_x,2,pexp)
  data_x = apply(data_u_x,2,Laplace_inverse)

  data_u_y = apply(data_exp_y,2,pexp)
  data_y = apply(data_u_y,2,Laplace_inverse)

  #fitting the HT model, conditioning on one variable

  ux = quantile(data_x[,1],q)
  uy = quantile(data_y[,2],q)
  dataX = data_x[data_x[,1]>ux,]
  dataY = data_y[data_y[,2]>uy,]

  Qpos <- function(param, yex, ydep, constrain, v, aLow) {

    a <- param[1]
    b <- param[2]

    res <- PosGumb.Laplace.negProfileLogLik(yex, ydep, a, b, constrain, v, aLow) # defined in file mexDependenceLowLevelFunctions
    res$profLik
  } # Close Qpos <- function

  #conditioning on large values of Y
  Yopt <- try(optim(par=c(.01, .01), fn=Qpos,
                    control=list(maxit=100000),
                    yex = dataY[,2], ydep = dataY[,1], constrain=TRUE, v=10, aLow=-1 + 10^(-10)), silent=TRUE)

  for( i in 2:nOptim ){
    Yopt <- try(optim(par=Yopt$par, fn=Qpos,
                      control=list(maxit=100000),
                      yex = dataY[,2], ydep = dataY[,1], constrain=TRUE, v=10, aLow=-1 + 10^(-10)), silent=TRUE)
  }


  Ypar = Yopt$par
  YZ = (dataY[,1] - Ypar[1]*dataY[,2])/(dataY[,2]^Ypar[2])

  #conditioning on large values of X
  Xopt <- try(optim(par=c(.01, .01), fn=Qpos,
                    control=list(maxit=100000),
                    yex = dataX[,1], ydep = dataX[,2], constrain=TRUE, v=10, aLow=-1 + 10^(-10)), silent=TRUE)

  for( i in 2:nOptim ){
    Xopt <- try(optim(par=Xopt$par, fn=Qpos,
                      control=list(maxit=100000),
                      yex = dataX[,1], ydep = dataX[,2], constrain=TRUE, v=10, aLow=-1 + 10^(-10)), silent=TRUE)
  }
  Xpar = Xopt$par
  XZ = (dataX[,2] - Xpar[1]*dataX[,1])/(dataX[,1]^Xpar[2])

  #Finding points on curve where y=x
  test = tryCatch(uniroot(heff_tawn_root,interval = c(Laplace_inverse(q),Laplace_inverse(1-prob-10^(-10))),prob=prob,Ypar=Ypar,YZ=YZ,nsim=10000),error = function(e){1})
  test2 = tryCatch(uniroot(heff_tawn_root,interval = c(Laplace_inverse(q),Laplace_inverse(1-prob-10^(-10))),prob=prob,Ypar=Xpar,YZ=XZ,nsim=10000),error = function(e){1})

  if(is.list(test) & is.list(test2)){ #if roots can be found for both conditioning variables

    y = test$root #extract values
    x = test2$root

    val = mean(c(x,y))#take the average of these values and use this as the crossover point from region 1 to region 2
    qvalsY = Laplace_cdf(val) #cdf probability value
    qvalsX = Laplace_cdf(val)

    qvalsY = exp(seq(log(1-prob-10^(-10)),log(qvalsY),length.out=100)) #create sequence of exceedance probabilities for Y variable in region 1
    y = Laplace_inverse(qvalsY) #find corresponding values of Y
    x = rep(NA,length.out=100)
    for(i in 1:length(y)){ #find corresponding values of X variable using fitted model
      sample_y = rexp(nsim) + y[i]
      sample_z = sample(YZ,nsim,replace = TRUE)
      sample_x = (sample_y^Ypar[2])*sample_z + Ypar[1]*sample_y
      x[i] = quantile(sample_x,(1-prob/(1-qvalsY[i])))
    }

  } else { #if rootfinder gives an error, we must try to find the exceedance probability manually or just consider y values above the q-th quantile

    qvalsY = exp(seq(log(1 - prob - 10^(-10)),log(min(1-sqrt(prob),q)),length.out=1000)) #exceedance probabilities
    y = Laplace_inverse(qvalsY) #correspoding values of y
    x = rep(NA,length.out = 1000)
    i = 1
    sample_y = rexp(nsim) + y[i]
    sample_z = sample(YZ,nsim,replace = TRUE)
    sample_x = (sample_y^Ypar[2])*sample_z + Ypar[1]*sample_y
    x[i] = quantile(sample_x,(1-prob/(1-qvalsY[i])))

    while(y[i]>x[i] & i<1000){ #finding the crossover point where we go from region 1 to region 2 (if it exists)
      #using model to obtain estimates of x values on curve
      i = i+1
      sample_y = rexp(nsim) + y[i]
      sample_z = sample(YZ,nsim,replace = TRUE)
      sample_x = (sample_y^Ypar[2])*sample_z + Ypar[1]*sample_y
      x[i] = quantile(sample_x,(1-prob/(1-qvalsY[i])))
    }

    qvalsX = Laplace_cdf(x[i])
    qvalsY = exp(seq(log(qvalsY[1]),log(qvalsY[i-1]),length.out=100)) #selecting probabilities up to the crossover point
    y = Laplace_inverse(qvalsY) #corresponding values of y
    x = rep(NA,length.out=100)
    for(i in 1:length(y)){#using model to obtain estimates of x values on curve
      sample_y = rexp(nsim) + y[i]
      sample_z = sample(YZ,nsim,replace = TRUE)
      sample_x = (sample_y^Ypar[2])*sample_z + Ypar[1]*sample_y
      x[i] = quantile(sample_x,(1-prob/(1-qvalsY[i])))
    }
  }

  qvalsX = exp(seq(log(max(qvalsX,q)),log(1 - prob - 10^(-10)),length.out=100)) #selecting exceedance probabilities of X variable such that curve point is in region 2
  x = c(x, Laplace_inverse(qvalsX))
  y = c(y,rep(NA,length.out = length(qvalsX)))

  for(i in (length(qvalsX)+1):(2*length(qvalsX))){#using model to obtain estimates of y values on curve
    sample_x = rexp(nsim) + x[i]
    sample_z = sample(XZ,nsim,replace = TRUE)
    sample_y = (sample_x^Xpar[2])*sample_z + Xpar[1]*sample_x
    y[i] = quantile(sample_y,(1-prob/(1-qvalsX[i-length(qvalsX)])),na.rm = TRUE)
  }

  #transforming the curve back to standard exponential margins
  curve = cbind(x,y)
  curve = apply(curve,2,Laplace_cdf)
  curve = apply(curve,2,qexp)

  #imposing property 3.1
  end = qexp(1-prob)
  curve = rbind(c(0,end),curve,c(end,0))

  #imposing property 3.4
  #change the way you do this - make it from centre
  for(i in 2:(length(curve[,1])-1)){
    if(curve[i+1,1]<curve[i,1]){
      curve[i+1,1] = curve[i,1]
    }
    if(curve[i+1,2]>curve[i,2]){
      curve[i+1,2] = curve[i,2]
    }
  }
  return(curve)#returns matrix containing curve estimate
}
