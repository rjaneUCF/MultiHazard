#' @noRd
predict.mex.conditioned <- function(object, which, pqu = .99, nsim = 1000, trace=10, smoothZdistribution=FALSE, ...){

    theCall <- match.call()



    # Class can be either mex or bootmex

    theClass <- class(object)[1]

    if (! theClass %in% c("mex", "bootmex")){

      stop("object must have class 'mex' or 'bootmex'")

    }



    if (theClass == "bootmex" ){

      which <- object$which

      migpd <- object$simpleMar

      margins <- object$margins

      constrain <- object$constrain

      dall <- mexDependence( migpd , which=which , dqu=object$dqu, margins = margins[[1]], constrain=constrain )

    } else {

      which <- object$dependence$which

      if(is.null(object$margins$referenceMargin)){

        migpd <- object$margins

      } else {

        migpd <- object$margins$referenceMargin

      }

      margins <- object$dependence$margins

      constrain <- object$dependence$constrain

      dall <- object

    }



    ################################################################

    MakeThrowData <- function(dco,z,coxi,coxmi,data){

      numb <- 1:nsim

      ui <- runif( nsim , min=pqu )

      z <- as.matrix(z[ sample( 1:( dim( z )[ 1 ] ), size=nsim, replace=TRUE ) ,])

      while(sum(numb)>0){

        nsim <- length(numb)

        ui[numb] <- runif( nsim , min=pqu )

        y <- margins$p2q(ui)

        distFun <- margins$q2p

        z[numb,] <- as.matrix(z[ sample( 1:( dim( z )[ 1 ] ), size=nsim, replace=TRUE ) ,])

        if(smoothZdistribution){

          z <- apply(z,2,function(x)x + rnorm(length(x),0,bw.nrd(x)))

        }

        ymi <- sapply( 1:( dim( z )[[ 2 ]] ) , makeYsubMinusI, z=z, v=dco , y=y )

        xmi <- apply( ymi, 2, distFun)

        xi <- u2gpd( ui, p = 1 - migpd$mqu[ which ], th=migpd$mth[ which ], sigma=coxi[ 1 ], xi = coxi[ 2 ] )

        for( i in 1:( dim( xmi )[[ 2 ]] ) ){

          xmi[, i ] <- revTransform( xmi[ ,i ], as.matrix(data[,-which])[, i ],

                                     th = migpd$mth[ -which ][ i ],

                                     qu = migpd$mqu[ -which ][ i ],

                                     sigma=coxmi[ 1,i ], xi=coxmi[ 2,i ] )

        }

        sim <- data.frame( xi , xmi , y, ymi, z)

        names( sim ) <- c( colnames( migpd$data )[ which ], colnames( migpd$data )[ -which ],

                           paste0(c(colnames( migpd$data )[ which ], colnames( migpd$data )[ -which ]),".trans"),

                           paste0(c(colnames( migpd$data )[ -which ]),".Z"))

        sim[,dim(sim)[2]+1] <- y > apply(ymi,1,max) # condlargest extra column

        ymi.max<-apply(ymi,1,max)

        numb <- which(y < ymi.max)

      }

      sim

    }



    ################################################################

    makeYsubMinusI <- function( i, z, v , y ){

      v <- v[ , i ]

      z <- z[ , i ]

      if ( !is.na( v[ 1 ] ) ){

        if( v[ 1 ] < 10^(-5) & v[ 2 ] < 0 ){

          if( v[ 4 ] < 10^(-5 ) ) d <- 0

          else d <- v[ 4 ]

          a <- v[ 3 ] - d * log( y )

        }

        else a <- v[ 1 ] * y

      } # close if( !is.na...

      else a <- NA

      a + ( y^v[ 2 ] ) * z

    }



    ###############################################################

    if (theClass == "bootmex"){

      # The function lfun does most of the work

      lfun <- function( i , bo, pqu, nsim , migpd, which ){

        if ( i %% trace == 0 ) cat( i, "sets done\n" )



        res <- MakeThrowData(dco=bo[[ i ]]$dependence,z=bo[[ i ]]$Z, coxi = bo[[i]]$GPD[,which],

                             coxmi = as.matrix(bo[[ i ]]$GPD[,-which]),

                             data = bo[[i]]$Y)

        res <- res[,1:((dim(res)[2]-1)/2)]

        res

      }



      bootRes <- lapply( 1:length( object$boot ) , lfun ,

                         migpd=migpd, pqu=pqu, bo = object$boot, nsim=nsim,

                         which = which )

      # bootRes contains the bootstrap simulated complete vectors X on the original

      # scale of the data, conditional on having the _which_ component above the pqu quantile.

    } else {

      bootRes <- NULL

    }



    ##########################################################################

    # Get a sample using the point estimates of the parameters

    # that are suggested by the data



    cox <- coef(migpd)[3:4, which]

    coxmi <- as.matrix(coef(migpd)[3:4, -which])



    sim <- MakeThrowData(dco=dall$dependence$coefficients,z=dall$dependence$Z,coxi=cox,coxmi=coxmi,data=migpd$data)

    CondLargest <- sim[,dim(sim)[2]]

    z<- sim[,(2*ncol(migpd$data)+1):(dim(sim)[2]-1)]

    transformed <- sim[,(ncol(migpd$data)+1):(2*ncol(migpd$data))]
    #sim <- sim[,1:((dim(sim)[2]-1)/2)]
    sim <- sim[,1:ncol(migpd$data)]


    m <- 1 / ( 1 - pqu ) # Need to estimate pqu quantile

    zeta <- 1 - migpd$mqu[ which ] # Coles, page 81

    pth <- migpd$mth[ which ] + cox[ 1 ] / cox[ 2 ] * ( ( m*zeta )^cox[ 2 ] - 1 )



    data <- list( real = data.frame( migpd$data[, which ], migpd$data[, -which] ), simulated = sim, z=z, pth=pth,CondLargest=CondLargest, transformed = transformed)

    names(data$real)[1] <- colnames(migpd$data)[which]



    res <- list( call = theCall , replicates = bootRes, data = data,

                 which = which, pqu = pqu,

                 mth=c( migpd$mth[ which ], migpd$mth[ -which ] ),

                 gpd.coef = coef(migpd)[,c(which,(1:dim(data$real)[2])[-which])])



    oldClass( res ) <- "predict.mex"

    res

  }
