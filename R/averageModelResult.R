# This file contains the public functions associated with AverageModelResult object


setOldClass( "by" )

##' Class representing data to use for extrapolation
##' @slot subject.data a data frame with X columns
##' "subject",...
##' see vignette and \code{AverageModelResult} for further iunformation
##' @export
setClass( "AverageModelResult", 
          slots= list( subject.data = "data.frame", curves="by", times="numeric", nsim="numeric", arms="character"  ) )



##' @name summary
##' @rdname summary-methods
##' @aliases summary,AverageModelResult-method
##' @export
# setMethod( "summary",
#            "AverageModelResult",
#   function( object ){
#     df <- object@subject.data 
#     
#     if( nrow( df ) > 0 ){
#       cat( paste( "Number of subjects:", nrow( df ),"\n") )
#       cat( paste( "Number of events:", sum( df$has.event ),"\n") )
#     }
# })

##' @rdname plot-methods
##' @name plot
##' @aliases plot,AverageModelResult,missing-method
##' @export
setMethod( "plot",
           signature( x="AverageModelResult", y="missing" ),
           function( x, xlab="log(t)", ylab="log(-log(S(t)))", main="",  ... ) {
               colours <- c("red", "blue")
               km.curve <- survfit(Surv(time, has.event) ~ arm, data=x@subject.data )
               prob<- c( 0.025, 0.5, 0.975) 
               newDim <- c(length(prob), length(x@times)  )
               quantiles <- 
                 lapply( x@curves,
                         function( curveSet, newDim ) {
                           res <- apply(curveSet, 1, quantile,
                                        prob=prob, names=FALSE)
                           dim(res) <- newDim
                           res
                         },
                         newDim )
               
               myxscale <- 365.25 / 12
               plot( km.curve, xscale=myxscale, col=colours, xlab='time (months)',
                     ylab="Survival", main=main,... )
               lines( x@times/myxscale, quantiles$combination[2,], col=colours[1])
               lines( x@times/myxscale, quantiles$patchOnly[2,], col=colours[2] )
               lines( x@times/myxscale, quantiles$combination[1,], col=colours[1], lty="dashed" )
               lines( x@times/myxscale, quantiles$patchOnly[1,], col=colours[2], lty="dashed"  )
               lines( x@times/myxscale, quantiles$combination[3,], col=colours[1], lty="dashed"  )
               lines( x@times/myxscale, quantiles$patchOnly[3,], col=colours[2], lty="dashed"  )
               legend( "topright", legend=x@arms, col=colours, lty=c(1,1))
           }
)

