

# Fits the five standard models to the data
fitStandardModels <- function( f, data ){
  list( 
    "exponential" = survreg( f, data=my.data@subject.data, dist="exponential" ),
    "weibull" = survreg( f, data=my.data@subject.data, dist="weibull" ),
    "loglogistic" = survreg( f, data=my.data@subject.data, dist="loglogistic" ),
    "lognormal" = survreg( f, data=my.data@subject.data, dist="lognormal" ),
    "gompertz" =  flexsurvreg( f, data=my.data@subject.data, dist="gompertz" )
  )
}

# Fits generalized gamma model
fitGeneralizedGammaModel <- function( f, data ){
  list( "gengamma" = flexsurvreg( f, data=my.data@subject.data, dist="gengamma" ) )
}

printICtable <- function( models ){ 
  ic.per.model <- rbind( sapply( models, AIC ), sapply( models, BIC )  ) 
  ic.per.model <- t( as.data.frame( ic.per.model ) )
  colnames( ic.per.model ) = c( "AIC", "BIC" )
  print( ic.per.model )
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


