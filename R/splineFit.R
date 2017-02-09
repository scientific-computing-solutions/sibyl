
##' Method to create SurvivalModel object with different spline models
##' @name fitSplines
##' @rdname fitSplines-methods
##' @param object \code{SurvivalModel} for which spline models are to be fitted
##' @param ... additional arguments for specific instances of this generic
##' @return A SurvivalModel object with different spline model fits
setGeneric( "fitSplines", function(object,...) standardGeneric("fitSplines"))

##' @rdname fitSplines-methods
##' @aliases fitSplines,SurvivalModel-methods
##' @param k (numeric vector) The number of knots to be used when fitting spline model
##' (if k=c(2,3) then models will be fit with both 2 and 3 knots)
##' @param scale (character) The scale arguement to be passed to flexsurvspline
##' @export
setMethod("fitSplines", signature(object="SurvivalModel"),
  function(object, k=0:5, scale=c("hazard", "odds", "normal")[1]){
    
    #validation
    if(length(scale)!=1 || ! scale %in% c("hazard","odds","normal")){
      stop("scale must be 'hazard', 'odds' or 'normal'")
    }
    if(length(k)==0 || any(!is.numeric(k) | k < 0)){
      stop("invalid k")
    }
    
    retVal <- NULL
    
    #for each k (for loop used as building up SurvivalModel object)
    for(thisK in k){
      modelOptions <- list(spline=list(k=thisK, scale=scale))
      
      #if no model successfully fit
      suppressWarnings(
      retVal <- if(is.null(retVal))
        fitModels(object@survData,
                  armAsFactor=object@armAsFactor, 
                  covariates=object@covariates, 
                  subgroup=object@subgroup, 
                  endPoint=object@endPoint,
                  model="spline", modelOptions=modelOptions) 
                else
        addModel(retVal, "spline",modelOptions=modelOptions)
      )
    }
    
    if(is.null(retVal)){
      warning("None of the spline models could be fit")  
    }
    
    retVal
  }
)

#return TRUE iff all model names are spline_k_scale for 
#the same scale AND there are at least 2 models
isSplineFit <- function(modelNames){
  splitNamesList <- strsplit(modelNames, split="_")
  scales <- unlist(lapply(splitNamesList, function(name){
    if(length(name)!=3) return(NA)
    name[3]
  }))
  #at least two, all are spline and all same scale
  !any(is.na(scales)) && length(scales)>1 && length(unique(scales))==1
}



#extract the knots/scale from the model names
#which are of the form spline_knots_scale (all scale are same)
extractKnots <- function(modelNames){
  extractForSpline(modelNames, 2)  
}

extractScale <- function(modelNames){
  extractForSpline(modelNames, 3)[1]
}


extractForSpline <- function(modelNames, index){
  vapply(modelNames, function(name){
    strsplit(name, "_")[[1]][index]
  }, FUN.VALUE=character(1))  
}