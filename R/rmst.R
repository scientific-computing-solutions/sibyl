##' @include semiParametric.R survivalModels.R
NULL


##' Method to calculate model restricted means
##' @name calcModelRmst
##' @rdname calcModelRmst-methods
##' @param object (SurvivalModel) A survival model - note there cannot be 
##' any covariates and armAsFactor must be FALSE
##' @param ... additional arguments for this generic  
##' @return (data.frame or FlexTable)
##' @export
setGeneric("calcModelRmst", function(object, ...){
  standardGeneric("calcModelRmst")
})


##' @name calcModelRmst
##' @aliases calcModelRmst,SurvivalModel-method
##' @rdname calcModelRmst-methods
##' @param model (character) The name of the model for which to calculate the restricted mean
##' @param times (nuermic vector) times to calculate the restricted mean
##' @param class ('data.frame' or "FlexTable' (default)) type of output required
##' @param digits (numeric default=3) if outputting a FlexTable then the number of digits
##' to round the entries to
##' @export
setMethod("calcModelRmst", "SurvivalModel", 
  function(object, model, times, class=c("data.frame","FlexTable")[2], digits=3, ...){
  
    #validation
    if(object@armAsFactor){
      stop("Cannot calculate restricted means if armAsFactor is TRUE")
    }
    if(length(object@covariates)!=0){
      stop("Cannot calculate restricted means if covariates fitted in model")
    }
    
    if(any(!is.numeric(times) | times < 0)){
      stop("Times must be numeric and non-negative")
    }
    
    if(length(class) != 1 || !class %in% c("data.frame","FlexTable")){
      stop("Invalid class argument, should be 'data.frame' or 'FlexTable")
    }
    
    if(length(digits)!=1 || !is.numeric(digits) || !digits > 0 || is.infinite(digits) ||
       is.na(digits)){
      stop("Invalid digits argument")
    } 
    
    if(length(model)!=1 || !model %in% names(object@models)){
      stop("Invalid model argument must be one of ",
           paste(names(object@models),collapse=", "))
    }
    
    
    #for each arm
    rmsts <- lapply(object@models[[model]],function(oneModel){
      
      #get the cdf function
      tempF <- oneModel$dfns$p
      
      #if spline need to add the knots argument for the dfns$p function to work
      args <- list()
      if(!is.null(oneModel$knots)){
        args$knots <- oneModel$knots  
      }  
      
      #survival function
      survFn <- function(x){
        1 - do.call("tempF", c(args, list(q=x), oneModel$res[,"est"]))
      }
      
      #calculate restricted means (an optimization would be to 
      #not handle times indpendently but calculate [0,t1], [t1, t2], ... and
      #sum them up as needed)
      vapply(times, function(time){
          tryCatch(
            integrate(survFn, 0, time)$value,
            error=function(cond) NA
          )  
        },
        FUN.VALUE = numeric(1))
      })
    
    rmsts <- as.data.frame(do.call("rbind",rmsts))
    colnames(rmsts) <- NULL
    
    #if two arms add a difference row
    if(nrow(rmsts)==2){
      rmsts <- rbind(rmsts,difference=as.numeric(rmsts[2,])-as.numeric(rmsts[1,]) )
    }
    
    #Add row of times
    rmsts <- rbind(time=times, rmsts)
    
    if(class=="data.frame") return(rmsts)
    
    #create FlexTable
    numRows <- nrow(rmsts)
    numCols <- 1+ncol(rmsts)
    
    MyFTable <- MyFTable <- FlexTable(numrow=numRows,numcol=numCols, 
                                      body.par.props=parProperties(text.align="right"),
                                      header.text.props = textProperties(font.weight = "bold"),
                                      body.cell.props = cellProperties(padding.right=1))
    
    #Set borders
    MyFTable[1:numRows,1:numCols,side='bottom'] <- borderProperties(width=0)
    MyFTable[1:numRows,1:numCols,side='left'] <- borderProperties(width=0)
    MyFTable[1:numRows,1:numCols,side='top'] <- borderProperties(width=0)
    MyFTable[1:numRows,1:numCols,side='right'] <- borderProperties(width=0)
    
    MyFTable[numRows,1:numCols,side='bottom'] <- borderProperties(width=3)
    MyFTable[1,1:numCols,side='top'] <- borderProperties(width=3)
    MyFTable[2,1:numCols,side='top'] <- borderProperties(width=3)
    
    #Add in data to table  
    MyFTable[2:numRows,2:numCols] <- round(rmsts[2:numRows,], digits=digits)
    MyFTable[1,2:numCols] <- times
    MyFTable[1:numRows, 1] <- rownames(rmsts)
    
    #Add header denoting which distribution
    hR <- FlexRow(paste(getDistributionDisplayNames(model),"\nrestricted means"),
                   colspan = numCols,
                   par.properties=parProperties(text.align="center",padding=1),
                   cell.properties = cellProperties(border.width = 0),
                   text.properties = textProperties(font.weight = "bold"))
    
    MyFTable <- addHeaderRow(MyFTable,hR)
    
    MyFTable
  }
)


##' Method to calculate RMST on subset of data contained in a
##' SemiParametricModel object
##' @name calcRmst
##' @rdname calcRmst-methods
##' @param object (SemiParametricModel) The object which was created when
##' fitting the Cox model
##' @param ... additional parameters needed for specific instances of this
##' generic
##' @return (rmst object) contains list of RMST values, differences and call or
##' a FlexTable for output into a word document (depending on the class variable) 
##' @export
setGeneric("calcRmst", function(object, ...){
  standardGeneric("calcRmst")
})

##' @name calcRmst
##' @aliases calcRmst,SemiParametricModel-method
##' @rdname calcRmst-methods
##' @param class ('rmst' or "FlexTable' (default)) type of output required
##' @param digits (numeric default=3) if outputting a FlexTable then the number of digits
##' to round the entries to
##' @export
setMethod("calcRmst", "SemiParametricModel", function(object,class=c("rmst","FlexTable")[2], digits=3, ...){
  
  if(length(class) != 1 || !class %in% c("rmst","FlexTable")){
    stop("Invalid class argument, should be 'rmst' or 'FlexTable")
  }
  
  if(isSingleArm(object)){
    stop("Cannot calculate rmst for a single arm trial")
  }
  
  # Create formula for Kaplan-Meier estimator
  formulaToFit <- survivalFormula(armAsFactor=!isSingleArm(object),
                                  covariates=character(0),
                                  timeCol = object@endPointDef[["timeCol"]],
                                  censorCol = object@endPointDef[["censorCol"]])
  
  
  # Call RMST 
  result <- rmst(formula = formulaToFit,
                 data = object@survData@subject.data, ...)
  
  if(class=="rmst"){
    return(result)
  }
  
  #create FlexTable (note rmst function only works with two arms so size of table is fixed)
  numRows <- 3
  numCols <- 6
  
  MyFTable <- FlexTable(numrow=numRows,numcol=numCols, 
                        body.par.props=parProperties(text.align="right"),
                        header.text.props = textProperties(font.weight = "bold"),
                        body.cell.props = cellProperties(padding.right = 1))
  
  #Add data
  MyFTable[3,2:6] <-   round(result$diff, digits)
  MyFTable[1:2,2:5] <- round(result$RMST[,2:5],digits)
  
  #Add 1st column (the arm names)
  MyFTable[1:3,1] <- c(as.character(getArmNames(object@survData)),"Difference")
  MyFTable[1:numRows,1] <- parProperties(text.align="left")
  MyFTable[1:numRows,1] <- textProperties(font.weight = "bold")
  
  #Add header
  hR <- FlexRow(c("Arm","RMST", "SE","Lower CI", "Upper CI", "p-value"),
                par.properties=parProperties(text.align="left"),
                cell.properties =cellProperties(padding.right = 1), 
                text.properties = textProperties(font.weight = "bold"))
  
  MyFTable <- addHeaderRow(MyFTable,hR)
  
  MyFTable
})