##' @include semiParametric.R
NULL

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
  
  # Create formula for Kaplan-Meier estimator
  formulaToFit <- survivalFormula(armAsFactor=TRUE,
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