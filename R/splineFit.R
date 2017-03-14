
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




##' Given a SurvivalModel object, output a data frame
##' containing the knot locations (NOT on log-scale) of the given model name 
##' @param object (SurvivalModel) Survival model containing the
##' models fit
##' @param k (numeric) The number of knots of the spline model
##' whose knot locations are required
##' @param scale ("hazard", "odds", or "normal") The scale argument
##' of the spline model whose knots are required 
##' @param class ("data.frame" or "FlexTable") whether to output the table
##' as a data.frame or FlexTable
##' @param digits (numeric) The number of digits to round the locations
##' when class="FlexTable"  
##' @export
getSplineKnotLocations <- function(object, k, scale, class=c("data.frame", "FlexTable")[2],
                                   digits=5){
  
  if(class(object)!="SurvivalModel"){
    stop("Object must be of type SurvivalModel")
  }
  
  if(length(k)!=1 || length(scale) != 1){
    stop("k and scale must be arguments of length 1")
  }
  
  if(length(class)!=1 || !class %in% c("data.frame","FlexTable")){
    stop("class argument must be data.frame or FlexTable")
  }
  
  if(length(digits)!=1 || !is.numeric(digits) || !digits > 0 || is.infinite(digits) ||
     is.na(digits)){
    stop("Invalid digits argument")
  } 
  
  splineModelName <- paste("spline", k, scale, sep="_")
  if(!splineModelName %in% names(object@models)){
    stop("Spline model with k ", k, "and scale ", scale, " has not been fitted")
  }
  
  #list of the given spline model, one per arm
  splineModel <- object@models[[splineModelName]]
  
  #no knot locations
  if(k==0){
    stop("No knot locations as there are no knots!")
  }
  
  #for each model (one per arm) extract the knot locations
  knotLocations <- lapply(splineModel, function(oneArmModel){
    
    knots <- oneArmModel$knots
    #safe as knots includes boundary positions so is always a vector of length >2
    knots <- knots[2:(length(knots)-1)]
    
    data.frame(knots=exp(knots))
    
  })
  
  retVal <- do.call(cbind, knotLocations)
  rownames(retVal) <- NULL
  colnames(retVal) <- if(is.na(names(knotLocations)[1])) "Spline knot locations" else names(knotLocations)
  
  if(class=="data.frame"){
    return(retVal)
  }
  
  MyFTable <- FlexTable(round(retVal,digits=digits),
            body.par.props=parProperties(text.align="right"),
            header.text.props = textProperties(font.weight = "bold"),
            body.cell.props = cellProperties(padding.right=2))
  
  if(ncol(retVal) > 1){
    hR <- FlexRow("Spline knot locations",colspan = ncol(retVal),
                  par.properties=parProperties(text.align="center",padding=1),
                  text.properties = textProperties(font.weight = "bold"),
                  cell.properties = cellProperties(border.top.width=1, border.bottom.width=0,
                                                   border.left.width=1, border.right.width=1))
    
    MyFTable <- addHeaderRow(first=TRUE, MyFTable, hR)
  }
  
  MyFTable
}
