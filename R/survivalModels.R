# This file contains the public functions associated with survivalMOdels object

##' Class storing the fitted survival models
##' @slot subject.data a data frame with X columns
##' see vignette and \code{SurvivalModel} for further iunformation
##' @export
setClass("SurvivalModel", 
         slots= list( models="list", 
                      surv.data = "SurvivalData", 
                      by="character", 
                      surv.formula = "formula", 
                      subgroups="character" ) )


# Fits the five standard models to the data
fitStandardModels <- function( f, my.data, by="", subgroups ){
  x <- my.data@subject.data
  
  # If subgroups!
  x.sub <- getSubgroupString( my.data, subgroups )
  if( x.sub != "" ){
    x <- x[ with( x, eval( parse( text=x.sub ) ) ), ]
  }  
  if( by=="arm" ){
    if( sum( grepl( "arm", as.character( f ) ) > 0 ) ) {
      warning( "The formulae should not contain arm if using option by=\"arm\"" )
    }
    x.ctrl   <- x[ x$arm == my.data@ctrl.arm, ]
    x.active <- x[ x$arm == my.data@active.arm, ]
    
    exponential <- list() 
    exponential[[ my.data@ctrl.arm ]]   <- survreg( f, data=x.ctrl, dist="exponential" )
    exponential[[ my.data@active.arm ]] <- survreg( f, data=x.active, dist="exponential" )
    
    weibull <- list() 
    weibull[[ my.data@ctrl.arm ]]   <- survreg( f, data=x.ctrl, dist="weibull" )
    weibull[[ my.data@active.arm ]] <- survreg( f, data=x.active, dist="weibull" )
    
    loglogistic <- list() 
    loglogistic[[ my.data@ctrl.arm ]]   <- survreg( f, data=x.ctrl, dist="loglogistic" )
    loglogistic[[ my.data@active.arm ]] <- survreg( f, data=x.active, dist="loglogistic" )
    
    lognormal <- list() 
    lognormal[[ my.data@ctrl.arm ]]   <- survreg( f, data=x.ctrl, dist="lognormal" )
    lognormal[[ my.data@active.arm ]] <- survreg( f, data=x.active, dist="lognormal" )
    
    gompertz <- list() 
    gompertz[[ my.data@ctrl.arm ]]   <- flexsurvreg( f, data=x.ctrl, dist="gompertz" )
    gompertz[[ my.data@active.arm ]] <- flexsurvreg( f, data=x.active, dist="gompertz" )

    new( "SurvivalModel", 
         models = list( 
           "exponential"=exponential,
          "weibull"=weibull,
          "loglogistic"=loglogistic,
          "lognormal"=lognormal,
          "gompertz" = gompertz ),
         surv.data = my.data, 
         by = by,
         surv.formula=f,
         subgroups = subgroups )
  }
  else{
    new( "SurvivalModel",
    models = list( 
        "exponential" = survreg( f, data=my.data@subject.data, dist="exponential" ),
        "weibull" = survreg( f, data=my.data@subject.data, dist="weibull" ),
        "loglogistic" = survreg( f, data=my.data@subject.data, dist="loglogistic" ),
        "lognormal" = survreg( f, data=my.data@subject.data, dist="lognormal" ),
        "gompertz" =  flexsurvreg( f, data=my.data@subject.data, dist="gompertz" )
        ),
    surv.data = my.data, 
    by = by,
    surv.formula=f,
    subgroups = subgroups )
  }
}




##' Method to create to check that covariates exist in data object
##' @rdname addModel
##' @param object survival data object
##' @param subgroups vector of subgroups
setGeneric( "addModel", function( object, ... ) standardGeneric( "addModel" ))


##' @rdname addModel
##' @export
setMethod( "addModel", 
           signature( object="SurvivalModel" ),
           function( object, name="gengamma" ){
             if( name != "gengamma" ){
               stop( "Can currently only add \"gengamma\"" )
             }
             models  <- object@models
             my.data <- object@surv.data
             if( object@by=="arm" ){
               if( sum( grepl( "arm", as.character( object@surv.formula ) ) > 0 ) ) {
                 warning( "The formulae should not contain arm if using option by=\"arm\"" )
               }
               x <- my.data@subject.data
               x.ctrl   <- x[ x$arm == my.data@ctrl.arm, ]
               x.active <- x[ x$arm == my.data@active.arm, ]
               gengamma <- list() 
               gengamma[[ my.data@ctrl.arm ]]   <- flexsurvreg( object@surv.formula, data=x.ctrl, dist="gengamma" )
               gengamma[[ my.data@active.arm ]] <- flexsurvreg( object@surv.formula, data=x.active, dist="gengamma" )
               models <- c( models, 
                            list( "gengamma" = gengamma ) )
             }
             else{
               models <- c( models, 
                            list( "gengamma" = flexsurvreg( object@surv.formula, data=my.data@subject.data, dist="gengamma" ) ) )
             }
             new( "SurvivalModel",
                  models = models,
                  surv.data = my.data, 
                  by = object@by,
                  surv.formula=object@surv.formula,
                  subgroups = object@subgroups )
})


##' Method to create to check that covariates exist in data object
##' @rdname printICtable
##' @param object survival data object
##' @param subgroups vector of subgroups
setGeneric( "printICtable", function( object,  ... ) standardGeneric( "printICtable" ))


##' @rdname printICtable
##' @export
setMethod( "printICtable", 
           signature( object="SurvivalModel"  ),
           function( object ){
              models <- object@models
              if( class( models[[1]] ) == "list" ){
                ic.per.model <- data.frame()
                for( i in seq_len( length( models ) ) ) {
                  xx <- sapply( models[[i]], AIC )
                  names(xx) <- paste0( "AIC:", names(xx) )
                  yy <- sapply( models[[i]], BIC )
                  names(yy) <- paste0( "BIC:", names(yy) )
                  row <- as.data.frame( c( xx, yy ) )
                  names(row) <- names(models)[i]
                  row <- t( row )
                  ic.per.model  <- rbind( ic.per.model, row ) 
                }
              }
              else{
                for( i in models ){
                  ic.per.model <- rbind( sapply( models, AIC ), sapply( models, BIC )  ) 
                  ic.per.model <- t( as.data.frame( ic.per.model ) )
                  colnames( ic.per.model ) = c( "AIC", "BIC" )
                }
        }
  print( ic.per.model ) 
})


survivalFormula <- function( my.data, armAsFactor=FALSE, covariates=NULL ){
  my.formula =  "Surv( time, has.event )"
  if( armAsFactor ){
    my.formula <- paste0( my.formula, " ~ as.factor( arm )" )
    
  } else{
    my.formula <- paste0( my.formula, " ~ 1" )
  }
  if( !is.null( covariates ) ){
    if( hasCovariates( my.data, covariates ) ){
      my.formula <- paste0( my.formula, " + ", paste0( covariates, collapse=" + " ) )
    }
    else{
      warning( paste0( "Covariates could not be added, please check names (", covariates, ")!" ) )
    }
  }
  as.formula( my.formula )
}

averageModel <- function( model, my.data, times, Nsim ){
  sim.params <- normboot.flexsurvreg( model, B=Nsim, newdata=my.data@subject.data )
  sim.params <- array( unlist(sim.params), 
                       dim=c( dim(sim.params[[1]] ), 
                              length(sim.params)), 
                       dimnames=dimnames(sim.params[[1]]))
  
  tmp.f <- model$dfns$p
  times <- list( times )
  cdf <- function( params ) { do.call( tmp.f, c( times, params )) }
  
  idx <- seq_len( nrow( my.data@subject.data ) )
  arm.levels <- my.data@subject.data$arm 
  
  
  curves <- by( idx, # get the indices of patients
                arm.levels,          # in each averaging group
                function( indices ) {
                  apply(sim.params[,, indices, drop=FALSE],
                        1, # compute the survival curves for each patient
                        function(x) {
                          ## "1 -" is here to turn into survival curves
                          ## the rowMeans averages over a Monte-Carlo
                          1 - rowMeans(apply(x, 2, cdf))
                        })
                })
  
  new( "AverageModelResult",
        subject.data=my.data@subject.data,
        curves=curves,
        arms=levels( arm.levels), 
        times=unlist( times ),
        nsim=Nsim )
}



##' @rdname plot-methods
##' @name plot
##' @aliases plot,SurvivalModel,missing-method
##' @export
setMethod( "plot",
         signature( x="SurvivalModel", y="missing" ),
         function( x,  ... ) {
           if( x@by == "arm" ){
             cat( "FLEXSURVMODELS NEED TO BE ADDED HERE!")
             models <- x@models
             data <- x@surv.data
             plot( data, loglogS=FALSE )
             xscale <- 1
             x.ctrl   <- predict( models$exponential[[data@ctrl.arm]], type="quantile", p=seq(.01,.99,by=.01))[1,]/xscale
             x.active <- predict( models$exponential[[data@active.arm]], type="quantile", p=seq(.01,.99,by=.01))[1,]/xscale
             y <- seq( .99,.01,by=-.01 )
             lines( x.active, y, col="gray", lwd=2 )
             lines( x.ctrl, y, col="gray", lwd=2 )
             
             x.ctrl   <- predict( models$weibull[[data@ctrl.arm]], type="quantile", p=seq(.01,.99,by=.01))[1,]/xscale
             x.active <- predict( models$weibull[[data@active.arm]], type="quantile", p=seq(.01,.99,by=.01))[1,]/xscale
             y <- seq( .99,.01,by=-.01 )
             lines( x.active, y, col="orange", lwd=2 )
             lines( x.ctrl, y, col="orange", lwd=2 )
             
             x.ctrl   <- predict( models$loglogistic[[data@ctrl.arm]], type="quantile", p=seq(.01,.99,by=.01))[1,]/xscale
             x.active <- predict( models$loglogistic[[data@active.arm]], type="quantile", p=seq(.01,.99,by=.01))[1,]/xscale
             y <- seq( .99,.01,by=-.01 )
             lines( x.active, y, col="brown", lwd=2 )
             lines( x.ctrl, y, col="brown", lwd=2 )
             
             x.ctrl   <- predict( models$lognormal[[data@ctrl.arm]], type="quantile", p=seq(.01,.99,by=.01))[1,]/xscale
             x.active <- predict( models$lognormal[[data@active.arm]], type="quantile", p=seq(.01,.99,by=.01))[1,]/xscale
             y <- seq( .99,.01,by=-.01 )
             lines( x.active, y, col="darkgreen", lwd=2 )
             lines( x.ctrl, y, col="darkgreen", lwd=2 )
             
             legend( "topright", c( "exponential", "weibull", "loglogistic", "lognormal" ), col=c( "gray", "orange", "brown", "darkgreen" ), lwd=rep(4,2) )
           }
           else {
             par( mfrow= c(1,2))
             cat( "CURVES NEED TO BE ADDED HERE!")
             for( i in 1:2 ) {
               models <- x@models[[i]]
               data <- x@surv.data
               plot( data, loglogS=FALSE )
             }
             par( mfrow= c(1,1))
           }
         })


##' Method to create to check that covariates exist in data object
##' @rdname printCoefficientEstimates
##' @param object survival data object
setGeneric( "printCoefficientEstimates", function( object,  ... ) standardGeneric( "printCoefficientEstimates" ))


##' @rdname printCoefficientEstimates
##' @export
setMethod( "printCoefficientEstimates", 
           signature( object="SurvivalModel"  ),
           function( object ){
  models <- object@models
  v_names <- names( models )
  cat( "\n---------------------------------------------------------------\n" )
  for( i in seq_along( v_names ) ){
    printCoef( models[[i]], v_names[i] )
  }
})

printCoef <- function( model, name ){
  cat( name )
  cat( "\n---------------------------------------------------------------\n" )
  if( length( model ) == 2 ) {  # Separate models per arm
      if( name == "exponential"){
        for( i in 1:2 ){
          if( i == 2 ) cat( "\n" )
          cat( toupper( names( model )[i] ) )
          cat( "\n")
          res <- summary( model[[i]] )[[9]]
          rownames(res)[1] = c( "Log(1/Rate)"  )
          print( res[,1:2] )
        }
    }
    else if( name == "weibull" || name == "loglogistic" || name == "lognormal" ){
      for( i in 1:2 ){
        if( i == 2 ) cat( "\n" )
        cat( toupper( names( model )[i] ) )
        cat( "\n" )
        res <- summary( model[[i]] )[[9]]
        rownames( res )[ 1 ] <- "Log.Scale"
        idx <- which( rownames( res ) == "Log(scale)" ) 
        rownames( res )[ idx ] = "1/Shape"
        print( res[,1:2] )
        if( length(res[,1]) > 2 ){
          rownames(res)[ -c(1,idx) ] <- paste0( "Log(", rownames(res)[ -c(1,idx) ], ")" )
        }
      }
    }
    else if( name=="gompertz" || name == "gengamma" ){
      for( i in 1:2 ){
        if( i == 2 ) cat( "\n" )
        cat( toupper( names( model )[i] ) )
        cat( "\n")
        res <- model[[i]]$res[,1:3]
        print( res )
      }
    }
  }
  else { # Arm is coefficient in model
    cat( "WARNING, OUTPUT NEEDS TO BE REFORMATED\nSEE PARAMETERIZATIONS FOR EACH DIST!\n")
    if( name == "exponential"){
      res <- summary( model )[[9]]
      rownames(res)[1:2] = c( "Log(1/Rate.1)", "Log(1/Rate.2 - 1/Rate.1)"  )
      print( res[,1:2] )
    }
    else if( name == "weibull" || name == "loglogistic" || name == "lognormal" ){
      res <- summary( model )[[9]]
      rownames( res )[ 1 ] <- "Log.Scale.1"
      rownames( res )[ 2 ] <- "Log.Scale.2 - Log.Scale.1"
      idx <- which( rownames( res ) == "Log(scale)" ) 
      rownames( res )[ idx ] = "1/Shape"
      print( res[,1:2] )
      if( length(res[,1]) > 3 ){
        rownames(res)[ -c(1,2,idx) ] <- paste0( "Log(", rownames(res)[ -c(1,2,idx) ], ")" )
      }
    }
    else if( name=="gompertz" || name == "gengamma" ){
      res <- model$res[,1:3]
      print( res )
    }
  }
  cat( "---------------------------------------------------------------\n" )
}
