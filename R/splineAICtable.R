##' Method to create the AIC table for the different spline models
##' @name createSplineAICTable
##' @rdname createSplineAICTable-methods
##' @param object \code{SurvivalModel} contains data to be fitted
##' @param ... additional arguments for specific instances of this generic (
##' @return A data frame or FlexTable object
setGeneric( "createSplineAICTable", function(object,...) standardGeneric( "createSplineAICTable" ))


##' @rdname createSplineAICTable-methods
##' @aliases createSplineAICTable,SurvivalModel-methods
##' @param k (numeric vector) The number of knots to be used when fitting spline model
##' (if k=c(2,3) then models will be fit with both 2 and 3 knots)
##' @param scale (character vector) The scale arguement to be passed to flexsurvspline
##' (if scale=c("hazard", "odds") then models will be fitted with both hazard and odds)
##' @param class ('data.frame' or 'FlexTable') the output data type
##' @param digits (numeric) Number of digits to round the flexTable output to
##' @export
setMethod("createSplineAICTable", signature(object="SurvivalModel"),
  function(object,
           k=2:5, scale=c("hazard", "odds", "normal"),
           class=c("data.frame","FlexTable")[2], digits=2){
  
    internalCreateSplineAICTable(object@survData, k=k, scale=scale, class=class, digits=digits,
                         armAsFactor=object@armAsFactor, 
                         covariates=object@covariates, 
                         subgroup=object@subgroup, 
                         endPoint=object@endPoint)            
  }
)



internalCreateSplineAICTable <- function(object,
           k=2:5, scale=c("hazard", "odds", "normal"),
           class=c("data.frame","FlexTable")[2], digits=2, ...){
  
    if(length(class) != 1 || !class %in% c("data.frame","FlexTable")){
      stop("Invalid class argument, should be 'data.frame' or 'FlexTable")
    }
    
    if(length(digits)!=1 || !is.numeric(digits) || !digits > 0 || is.infinite(digits) ||
       is.na(digits)){
      stop("Invalid digits argument")
    } 
    if(length(scale)==0 || any(! scale %in% c("hazard","odds","normal"))){
      stop("scale must be a vector containing only hazard, odds and/or normal")
    }
    if(length(k)==0 || any(!is.numeric(k) | k < 0)){
      stop("invalid k")
    }
  
    retVal <- lapply(k, function(thisK){
      
      allScale <- lapply(scale, function(thisScale){
        tryCatch({
          suppressWarnings(mod <- fitModels(object, model="spline",  
             modelOptions=list(spline=list(k=thisK,scale=thisScale)), ... ))
        
          icTable <- createIcTable(mod, summaryFn="sum", class="data.frame")
        
        icTable$Model <- NULL
        icTable$BIC <- NULL
        icTable$Arm <- NULL
        icTable$k <- thisK
        icTable$scale <- thisScale
        icTable
        },error=function(cond){
          data.frame(AIC=NA,k=thisK,scale=thisScale)
        })
        
      })
      
      do.call("rbind",allScale)
      
    })

    retVal <- do.call("rbind", retVal)
    
    if(all(is.na(retVal$AIC))){
      warning("Could not fit any spline models!")
    }
    
    if(class=="data.frame") return(retVal)
    
    numRows <- length(k)
    numCols <- 1 + length(scale)
    
    #create table
    MyFTable <- FlexTable(numrow=numRows,numcol=numCols, 
                          body.par.props=parProperties(text.align="right"),
                          header.text.props = textProperties(font.weight = "bold"),
                          body.cell.props = cellProperties(padding.right=3))
    
    #set up borders
    MyFTable[1:numRows,1:numCols,side='bottom'] <- borderProperties(width=0)
    MyFTable[1:numRows,1:numCols,side='left'] <- borderProperties(width=0)
    MyFTable[1:numRows,1:numCols,side='top'] <- borderProperties(width=0)
    MyFTable[1:numRows,1:numCols,side='right'] <- borderProperties(width=0)
    
    MyFTable[numRows,1:numCols,side='bottom'] <- borderProperties(width=3)
    MyFTable[1,1:numCols,side='top'] <- borderProperties(width=3)
    
    MyFTable[1:numRows,1] <- textProperties(font.weight = "bold")
    
    MyFTable[1:numRows,1] <- k
    
    hR <- FlexRow(c("knots",paste0("AIC (", scale , ")")), 
                  par.properties=parProperties(text.align="left"),
                  text.properties = textProperties(font.weight = "bold"),
                  cell.properties = cellProperties(border.top.width=3, border.bottom.width=0,
                                                   border.left.width=0, border.right.width=0))
    
    MyFTable <- addHeaderRow(MyFTable,hR)
    MyFTable[1:numRows,2:numCols] <- round(retVal$AIC, digits=digits)
    MyFTable
}
