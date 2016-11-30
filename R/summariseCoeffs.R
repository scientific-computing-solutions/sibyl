##' @include survivalModels.R
NULL

##' Method to create summary tables of coefficients for each model
##' @rdname summariseCoeffs
##' @param object (SurvivalModel object) contains fitted models to be summarised
##' @param ... additional arguments to summariseCoeffs
setGeneric( "summariseCoeffs", function(object, ...) standardGeneric( "summariseCoeffs" ))


##' @rdname summariseCoeffs
##' @param class ('matrix' or 'FlexTable' (default)) output format
##' for the summary tables
##' @param digits (numeric) The number of significant digits to round the FlexTable output to.
##' This option is ignored if class is not 'FlexTable'
##' @export
setMethod( "summariseCoeffs","SurvivalModel",
  function(object, class=c("matrix","FlexTable")[2], digits=3){
   
    #get covariate and arm column definitions for calculating the rownames of the output table
    columnDetails <- c(list(object@survData@armDef),object@survData@covDef)
    
    calcModelTables(object@models, class=class,summaryFn=getParamSummariser, digits=digits,
                    columnDetails)
  }
)


##' Method to return covariance matrix of estimated model parameters
##' @name vcov
##' @rdname vcov-methods
##' @param object (SurvivalModel object) contains fitted models with estimated
##'        parameters
##' @param ... additional arguments to vcov       
##' @export
setGeneric("vcov",function(object,...) standardGeneric("vcov"))

##' @rdname vcov-methods
##' @aliases vcov,SurvivalModel-methods
##' @param class ('matrix' or 'FlexTable' (default)) output format
##' for the summary tables
##' @param digits (numeric) The number of significant digits to round the FlexTable output to.
##' This option is ignored if class is not 'FlexTable'
##' @export
setMethod("vcov","SurvivalModel",
  function(object,  class=c("matrix","FlexTable")[2], digits=3){
    
    #get covariate and arm column definitions for calculating the rownames of the output table
    columnDetails <- c(list(object@survData@armDef),object@survData@covDef)
    
    calcModelTables(object@models, class=class, summaryFn=function(modelName, modelClass){vcov},
                    digits=digits, columnDetails, leftBorder = TRUE)
})

# For each model (weibull, loglogistic, etc.) and each data set (e.g. arm),
# return a table of model parameters obtained using a summary function
# returned when calling summaryFn using the model's name and class as parameters
# class 'data.frame' or 'FlexTable' - output format
# digits  The number of significant digits to round the FlexTable output to.
#leftBorder - should there be a border between the first and second columns
#of the FlexTable
#columnDetails is the list of treatment arm and covariate ColumnDef objects used to replace
#column names with the display names for the FlexTable outputs
calcModelTables <- function(models, class, summaryFn, digits, columnDetails,leftBorder=FALSE){

  #validate class
  if(length(class) != 1 || !class %in% c("matrix","FlexTable")){
    stop("Invalid class argument, should be 'matrix' or 'FlexTable")
  }
  
  # Initialise list of parameters for all results (all models, all data sets)
  dataTables <- vector("list", 0)

  allModelNames <- names(models)

  for (idxModel in seq_len(length(models))){

    # Extract model-fitting results for this model
    thisModelName <- tolower(allModelNames[[idxModel]])
    thisModelResults <- models[[idxModel]]

    # Look up summariser function for this distribution. Note that elements of
    # thisModelResults all come from the same model (applied to different data)
    # so the summariser can depend only on the first set of results.
    summariser <- summaryFn(thisModelName, class(thisModelResults[[1]]))

    # Look up parameter info for this model fitted to each data set
    thisModelData <- lapply(thisModelResults, summariser)

    # Name the data sets
    names(thisModelData) <- names(thisModelResults)

    #Convert to FlexTable if required
    if(class=="FlexTable"){
      
      #calculate the rownames for the output table
      #i.e. changing COV_racehispanic to race:hispanic
      #this needs to be done on a per arm basis as some arms 
      #may not have all factors
      updatedRowNames <- lapply(thisModelResults, getDisplayRowNames, columnDetails) 
      
      #convert matrices
      thisModelData <- mapply(convertSummariseToFlexTable, data=thisModelData, 
                              theRowNames=updatedRowNames, 
                              MoreArgs=list(modelName=thisModelName,
                                            digits=digits,
                                            leftBorder=leftBorder),SIMPLIFY=FALSE) 
      
    }
    
    # Put results for this model into list of all results
    dataTables[[1 + length(dataTables)]] <- thisModelData
      
  }

  # Name the models
  names(dataTables) <- allModelNames

  dataTables
}


#convert a summary table to a FlexTable
#data - the data frame to convert to FlexTable
#theRowNames the first column of the table should contain this vector
#modelName - the name of the model, to be used in header
#digits - number of significant digits to round the numbers
#leftBorder - should there be a border between the first and second columns
#of the FlexTable
convertSummariseToFlexTable <- function(data,theRowNames, modelName, digits, leftBorder){
  
  if(nrow(data)==ncol(data) && all(colnames(data)==rownames(data))){
    colnames(data) <- theRowNames
  }
  
  
  data <- formatC(data,digits=digits,format="g", preserve.width = "individual")
  #Add the row names 
  data <- cbind(data.frame(dispname=theRowNames), data)
  
  
  
  #Set the left column header to be the display name of the distribution
  colnames(data)[1] <- getDistributionDisplayNames(modelName)
  
  #create FlexTable
  MyFTable <- FlexTable(data,header.columns = TRUE,
                        body.par.props=parProperties(text.align="right"),
                        header.text.props = textProperties(font.weight = "bold"),
                        body.cell.props = cellProperties(padding.right=3,border.width = 0),
                        header.cell.props = cellProperties(border.left.width = 0,
                                                           border.top.width = 3,
                                                           border.right.width = 0,
                                                           border.bottom.width = 0,
                                                           padding.right=3))
  
  numRows <- MyFTable$numrow
  numCols <- MyFTable$numcol
  
  #set borders
  MyFTable[numRows,1:numCols,side='bottom'] <- borderProperties(width=3)
  MyFTable[1,1:numCols,side='top'] <- borderProperties(width=3) 
  
  if(leftBorder){ #border to right of first column
    MyFTable[1:numRows,1,side='right'] <- borderProperties(width=3)  
  }  
  
  MyFTable[1:numRows,1] <- parProperties(text.align="left")
  MyFTable[1:numRows,1] <- textProperties(font.weight = "bold")
  
  MyFTable
}


# Look up summary function from name and class of model
getParamSummariser <- function(modelName, modelClass){

  if (modelClass == "flexsurvreg"){
    return(flexsurvSummariser)
  }
  else if(modelClass == "survreg"){
    switch(tolower(modelName),
           "weibull"     = weibullSummariser,
           "loglogistic" = loglogisticSummariser,
           stop(paste0("Parameter summary for '", modelName, "' model is not supported for package survival")))
  }
}


# Summarise parameters of flexsurv models
flexsurvSummariser <- function(fittedModel){
  result <- fittedModel$res[, 1:3]
  #in the case with exactly one model parameter (which
  #will only occur if no covariates and using exponential distribution)
  #then R is useless and coerces the matrix to a numeric vector 
  if(class(result)=="numeric"){
    result <- t(data.frame(result))
    rownames(result) <- rownames(fittedModel$res)
  }
  
  result
}


# Summarise parameters of Weibull model using summary() method of Survival package
weibullSummariser <- function(fittedModel){
  paramTable <- summary(fittedModel)[[9]]

  rownames(paramTable)[1] <- "Log.Scale.1"
  rownames(paramTable)[2] <- "Log.Scale.2 - Log.Scale.1"

  idx <- which(rownames(paramTable) == "Log(scale)")
  rownames(paramTable)[idx] = "1/Shape"

  paramTable <- paramTable[, 1:2]

  if (length(paramTable[, 1]) > 3){
    rownames(paramTable)[-c(1,2,idx)] <- paste0("Log(", rownames(paramTable)[-c(1,2,idx)], ")")
  }

  return(paramTable)
}


# Summarise parameters of loglogistic model using summary() method of Survival package
loglogisticSummariser <- function(fittedModel){
  paramTable <- summary(fittedModel)[[9]]

  rownames(paramTable)[ 1 ] <- "Log.Scale"

  idx <- which(rownames(paramTable) == "Log(scale)")
  rownames(paramTable)[idx] = "1/Shape"

  paramTable <- paramTable[, 1:2]

  if (length(paramTable[,1]) > 2){
    rownames(paramTable)[-c(1,idx)] <- paste0("Log(", rownames(paramTable)[-c(1,idx) ], ")")
  }

  return(paramTable)
}


#' Convert the distribution names from FlexSurv (e.g. `lnorm')
#' to display names (e.g. Lognormal)
#' @details If name is not  lnorm, llogis, gengamma then the function returns
#' the name unchanged except the first character capitalized.
#' @param modelNames (character vector) A list of distribution names
#' @return A vector of display names
#' @export
getDistributionDisplayNames <- function(modelNames){
  
  convertFunc <- function(x){
    if(x=="lnorm") return("Lognormal")
    if(x=="llogis") return("Loglogistic")
    if(x=="gengamma") return("Generalized Gamma")
    if(nchar(x)==1) return(toupper(x))
    paste0(toupper(substr(x,1,1)),substr(x,2,nchar(x)))
  }
  
  vapply(modelNames,convertFunc,FUN.VALUE=character(1))
}


#take the rownames from model$res and output a vector
#of display names to be output in the FlexTable outputs
#for example   COV_racehispanic would be replaced with
#race:hispanic if the column name was COV_race and its display name
#was race
#model - an individual fitted model
#columnDetails is the list of the treatment arm and covariate ColumnDef objects used to replace
#column names with the display names for the FlexTable outputs
getDisplayRowNames <- function(model, columnDetails){
  
  originalNames <- rownames(model$res)
  
  #keep rate, scale, Q etc. unchanged
  newNames <- originalNames[model$basepars]
  
  #get columnNames
  columnNames <- listColumnDefSlot(columnDetails,"columnName")
  
  #loop over each covariate/coefficient
  for(covName in attributes(model$data$m)$covnames.orig){
    
    #get displayName
    displayName <- columnDetails[[which(columnNames==covName)]]@displayName
    
    #get type
    colType <- columnDetails[[which(columnNames==covName)]]@type
    
    numReplace <- 1
    #if factor replacing more 
    if(!is.null(attributes(model$data$m[[covName]]))){
      numReplace <- length(attributes(model$data$m[[covName]])$levels) - 1 
    }
    
    replacedNames <- originalNames[length(newNames)+(1:numReplace)]
    replacedNames <- substr(replacedNames,nchar(covName)+1,nchar(replacedNames))
    
    sepChar <- if(colType != "numeric") ":" else ""
    replacedNames <- paste(displayName,replacedNames,sep = sepChar)
    
    newNames <- c(newNames, replacedNames)
    
  }
  
  newNames
}
