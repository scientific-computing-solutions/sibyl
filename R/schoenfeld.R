##' @include semiParametric.R
NULL


##' Method to calculate the Schoenfeld residuals for a SemiParametricModel
##' object.
##' @details Note there is no stratification used when calculating fitting 
##' the Cox model
##' @name schoenfeldResiduals
##' @rdname schoenfeldResiduals-methods
##' @param object (SemiParametricModel) The object for which to 
##' @param ... additional arguments passed to cox.zph  
##' @return An object of class \code{cox.zph}
##' @export
setGeneric("schoenfeldResiduals", function(object, ...) 
  standardGeneric("schoenfeldResiduals")
)


##' @name schoenfeldResiduals
##' @aliases schoenfeldResiduals,SemiParametricModel-method
##' @rdname schoenfeldResiduals-methods
setMethod("schoenfeldResiduals", "SemiParametricModel", 
  function(object, ...){
    cox.zph(object@cox, ...)  
  }
)
