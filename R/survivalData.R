# This file contains the public functions associated with survivalData object


##' @include common.R 

# Function that will check the validity of SurvivalData object
# when constructed
# @param object SurvivalData object
checkSurvivalData <- function( object ){
  errors <- ""
  
  required.colnames <- c( "subject", "arm", "censor", "time" )
  
  data <- object@subject.data
  
  if( nrow( data ) > 0 ){
  
     if( any( !data$censor %in% c(0,1) ) ) errors <- paste(errors,"censor must be 1 or 0" )
     if( any( !required.colnames %in% colnames( data ) ) ) errors <- paste( errors,"invalid column names")

     if( class(data$time)  != "numeric" ) {
         errors <- paste(errors,"times must be numeric")
       }
       else{
         if( all( data$time==0 ) ) errors <- paste( errors, "time values incorrect, all are 0!")  
         if( any( data$time < 0 | is.na( data$time ) ) ) {
           errors <- paste( errors,"subjects cannot have non-positive time on trial. Please check subjects ",
                            paste( data[data$time <= 0 | is.na( data$time),]$subject, collapse=", ") )
         }
       }
     
       if(any(duplicated(data$subject))) errors <- paste(errors,"subject ID must be unique")
     
       if( all( data$censor==1 ) ){
         warning("No events have occurred - a model cannot be fit!")
       }
       if( all( data$censor==0 ) ) errors <- paste(errors,"all events have occurred!")
  }
  if( errors == "" ) TRUE else errors
}
  
##' Class representing data to use for extrapolation
##' @slot subject.data a data frame with X columns
##' "subject",...
##' see vignette and \code{SurvivalData} for further iunformation
##' @export
setClass("SurvivalData", 
          slots= list( subject.data = "data.frame", ctrl.arm="character", active.arm="character", subgroups="character", covdef="list" ),
          validity = checkSurvivalData )


##' Constructor for event data object
##' 
##' All dates must be in one of the following formats:
##' YYYY-MM-DD, DD/MM/YY or DD Month YYYY
##'  
##' 
##' @param data A data frame 
##' @param subject string, the column name of subject identifiers
##' @return An \code{SurvivalData} object
##' @export
SurvivalData <- function( data, 
                          subject,
                          arm, 
                          covariates=NULL,
                          covdef=NULL,
                          subgroups=NULL,
                          ctrl.arm,
                          active.arm,
                          censor, 
                          time ){
  
  # Check columns are in data frame 
  for(x in c(subject,
             arm, 
             censor, 
             time, covariates, subgroups )){
    if( !is.null( x ) && !is.na( x ) && !x %in% colnames( data ) ){
      stop( paste( "Column name", x, "not found in data frame" ) )
    }
  }

  c.time  <- if( !is.na( time ) )  data[,time]  else NA
  c.censor  <- if( !is.na( censor ) )  data[,censor]  else NA
  c.has.event <-  if( !is.na( censor ) )  -1*( data[,censor]-1 ) else NA
  
  if( any( c.time == 0 ) ){
    warning( "Some times are zero, +1 will be added to all times! Please check data!" )
    c.time <- c.time + 1
  }
  
  
  subject.data <- data.frame(
    subject=data[,subject],
    arm=data[,arm],
    time =      c.time,
    censor =    c.censor,
    has.event = c.has.event
  )
  c_cov <- as.data.frame( data[, covariates ] )
  names( c_cov ) <- covariates
  c_sub <-  as.data.frame( data[, subgroups ] )
  names( c_sub ) <- subgroups
  
  # ADD SOME CHECKS FOR COVs and SUBs
  # 1. MAKE SURE THE COVARIATE IS DEFINED
  for( i_cov in covariates ) {
    if( !any( sapply( covdef, function(x) x@name==i_cov  ) ) ) {
      stop( paste0( "Covariate ", i_cov, " not defined!" ))
    }
  }
  # 2. IF CATEGORICAL VARIABLE, CHECK THAT VALID VALUES
  idx <- which( sapply( covdef, function(x){ x@type=="categorical" && x@name==i_cov } ) )
  if( length( idx > 0 ) ) {
    for( i in idx )
    name <- covdef[[ i ]]@name
    categories <- covdef[[ i ]]@categories
    if( !is.factor( data[, name] ) ) data[, name] = as.factor( data[, name] )
    if( !all( levels( data[, name ] ) %in% categories ) ) {
      stop( paste0( " Covariate ", name, " is categorical but values: ",  paste0( levels( data[, name ] ), collapse=", " ), 
                    " not all found in definitions: ", paste0( categories, collapse=", " ) ))
    }
  }
  subject.data <- cbind( subject.data, c_cov, c_sub )
  
  if( !any( data[,arm]==ctrl.arm ) )   { stop( "Control arm label is not found in the arm column!" )}
  if( !any( data[,arm]==active.arm ) ) { stop( "Active arm label is not found in the arm column!" )}
  
  # Validation occurs in the validity function of the class 
  return( new( "SurvivalData", subject.data = subject.data, ctrl.arm=ctrl.arm, active.arm=active.arm, subgroups=subgroups, covdef = covdef ) )
}




##' @name show
##' @rdname show-methods
##' @aliases show,SurvivalData-method
##' @export
setMethod("show",
          "SurvivalData",
    function(object) {
         cat("SurvivalData, use object@subject.data$param to access individual columns:\n")
         cat( str( object@subject.data ) )
     })


##' @name summary
##' @rdname summary-methods
##' @aliases summary,SurvivalData-method
##' @export
setMethod( "summary",
           "SurvivalData",
  function( object ){
    df <- object@subject.data 
    
    if( nrow( df ) > 0 ){
      cat( paste( "Number of subjects:", nrow( df ),"\n") )
      cat( paste( "Number of events:", sum( df$has.event ),"\n") )
    }
})

##' @rdname plot-methods
##' @name plot
##' @aliases plot,SurvivalData,missing-method
##' @export
setMethod( "plot",
           signature( x="SurvivalData", y="missing" ),
           function( x, xlab="log(t)", ylab="log(-log(S(t)))", main="", separate.plots=FALSE, loglogS=TRUE, ... ) {
             
             if( nrow( x@subject.data ) == 0 ) stop( "Empty data frame!" )
             if( sum( x@subject.data$has.event ) == 1 ){
               stop( "Cannot fit a model to a dataset with no events" )
             }
             
          if( loglogS ){
            loglogSplot( x@subject.data, 
                         arms=c(x@ctrl.arm, x@active.arm ),
                         xlab, 
                         ylab, 
                         main, 
                         separate.plots=separate.plots,
                         ... )
          }
          else {
            if( xlab=="log(t)" ){ xlab = "time" }
            if( ylab=="log(-log(S(t)))" ){ ylab = "S(t)" }
            kmPlot( x@subject.data, 
                    arms=c(x@ctrl.arm, x@active.arm ),
                    xlab, 
                    ylab, 
                    main, 
                    separate.plots=separate.plots,
                    ... )
          }
      }
)


##' @name loglogS plot for a given endpoint
##' @param time Time variable
##' @param has.event Has event variable
##' @param title text to be displayed above plot
##' @param xlab xlab
##' @param ylab ylab
##' @param ... additional plot params 
##' @export
loglogSplot = function( x, arms, xlab, ylab, main, separate.plots=FALSE, ... ){
  
  x1 <- x[ x$arm == arms[1], ]
  x2 <- x[ x$arm == arms[2], ]
  
  
  model.1 <- survfit( Surv( time, has.event ) ~ 1, data=x1 )
  model.2 <- survfit( Surv( time, has.event ) ~ 1, data=x2 )
   
  res.1 <- data.frame( t = model.1$time, s = model.1$surv )
  res.1 <- res.1[ res.1$t > 0 & res.1$s > 0 & res.1$s < 1, ]
  res.2 <- data.frame( t = model.2$time, s = model.2$surv )
  res.2 <- res.2[ res.2$t > 0 & res.2$s > 0 & res.2$s < 1, ]  

  df.1 <- data.frame( x = log( res.1$t ), y = log( -log( res.1$s ) ) )
  df.2 <- data.frame( x = log( res.2$t ), y = log( -log( res.2$s ) ) )
  
  if( separate.plots == TRUE ){
   p1 <- qplot(x, y,data=df.1) + 
          stat_smooth( method="lm", se=FALSE, color="red", size=1 ) +
          xlab( xlab ) + ylab( ylab ) +
          ggtitle( arms[1] ) + theme_bw() +
          theme( panel.grid.minor = element_line( colour="gray", size=0.5 ) ) 
   p2 <- qplot(x, y,data=df.2) + 
    stat_smooth( method="lm", se=FALSE, color="red", size=1 ) +
    xlab( xlab ) + ylab( ylab ) +
    ggtitle( arms[2] ) + theme_bw() +
    theme( panel.grid.minor = element_line( colour="gray", size=0.5 ) ) 
    multiplot( p1, p2, cols=1 )
  }
  else {
    df.1$arm=arms[1]
    df.2$arm=arms[2]
    x <- rbind( df.1, df.2 )
    ggplot( x, aes(x, y, color=arm) ) + 
      geom_point() + 
      stat_smooth( method="lm", se=FALSE, size=1, aes(color=arm) ) +
      xlab( xlab ) + ylab( ylab ) +
      ggtitle( "" ) + theme_bw() +
      theme( panel.grid.minor = element_line( colour="gray", size=0.5 ) ) 
  }
}


##' @name KM plot for a given endpoint
##' @param time Time variable
##' @param has.event Has event variable
##' @param title text to be displayed above plot
##' @param xlab xlab
##' @param ylab ylab
##' @param ... additional plot params 
##' @export
kmPlot = function( x, arms, xlab, ylab, main, separate.plots=FALSE, ... ){
  # Make base graphics plot for KM curve
  x1 <- x[ x$arm == arms[1], ]
  x2 <- x[ x$arm == arms[2], ]
  
  model.1 <- survfit( Surv( time, has.event ) ~ 1, data=x1 )
  model.2 <- survfit( Surv( time, has.event ) ~ 1, data=x2 )
  
  res.1 <- data.frame( t = model.1$time, s = model.1$surv )
  res.1 <- res.1[ res.1$t > 0 & res.1$s > 0 & res.1$s < 1, ]
  res.2 <- data.frame( t = model.2$time, s = model.2$surv )
  res.2 <- res.2[ res.2$t > 0 & res.2$s > 0 & res.2$s < 1, ]  
  
  df.1 <- data.frame( x = res.1$t, y = res.1$s ) 
  df.2 <- data.frame( x = res.2$t, y = res.2$s ) 
  
  if( separate.plots == TRUE ){
    par( mfrow=c(1,2 ))
    plot( df.1$x, df.1$y, col="red", "s", xlab=xlab,  ylab=ylab, main=main ) 
    plot( df.2$x, df.2$y, col="red", "s", xlab=xlab,  ylab=ylab, main=main )
    par( mfrow=c(1,1) )
  }
  else {
    plot( df.1$x, df.1$y, col="red", "s", xlab=xlab,  ylab=ylab, main=main )
    lines( df.2$x, df.2$y, col="blue", "s" )
  }
}



setMethod( "summary", 
           signature( object="SurvivalData" ),
           function( object, ... ) {
  
  ## MAY NEED MODIFICATION FOR A GENERAL SETTING
  output <- ddply( object@subject.data, 
                   .( eval(as.symbol( object@subgroups )), arm ), 
                   summarize, 
                   n = length( subject ),
                   events = sum( has.event ),
                   p.events = round( sum( has.event )/length( subject ), 2 ) )
  names( output )[ length(object@subgroups) ] <- object@subgroups
  output <- t( output ) 
  colnames( output ) = rep( "", ncol(output) )
  x <- table( object@subject.data[ object@subgroups ] )
  xx <- table( object@subject.data$arm )
  cat( "\n-------------------------------------------------------------\n")
  cat( paste( "arm", ":", names(xx), xx, collapse="\t" ) )
  cat("\n")
  cat( paste( object@subgroups, ":", names(x), x, collapse="\t" ) )
  cat( "\n-------------------------------------------------------------")
  print( output )
  cat( "\n-------------------------------------------------------------")
})


##' Method to create to check that covariates exist in data object
##' @rdname hasCovariates
##' @param object survival data object
##' @param covariates vector of covariates 
setGeneric( "hasCovariates", function( object, covariates, ... ) standardGeneric( "hasCovariates" ))


##' @rdname hasCovariates
##' @export
setMethod( "hasCovariates", 
  signature( object="SurvivalData", covariates="character" ),
  function( object, covariates ){
    xx <- sapply( covariates, function( i_cov ) {
      i_cov %in% names( object@subject.data ) &&
        i_cov %in% sapply( object@covdef, function( i ) i@name )
    })
    if( all( xx ) ) { return( TRUE ) }
    else { warning( paste0( covariates[ xx==FALSE ], " not found!" ) ); return( FALSE ) }
})


##' Method to create to check that covariates exist in data object
##' @rdname getSubgroupString
##' @param object survival data object
##' @param subgroups vector of subgroups
setGeneric( "getSubgroupString", function( object, subgroups, ... ) standardGeneric( "getSubgroupString" ))


##' @rdname getSubgroupString
##' @export
setMethod( "getSubgroupString", 
           signature( object="SurvivalData", subgroups="character" ),
           function( object, subgroups ){
             xx <- sapply( subgroups, function( i_sub ) {
               i_sub %in% object@subgroups
             })
             if( all( xx ) ) { 
               paste0( "(",
                 paste0( subgroups, "==1", collapse="&" ),
                 ")" )
             }
             else { warning( paste0( subgroups[ xx==FALSE ], " not found!" ) ); return( "" ) }
           })
