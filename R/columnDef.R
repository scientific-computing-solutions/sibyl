# This file contains the public functions associated with Covariates object


##' Class representing a column of input data for the \code{SurvivalData} object
##'   
##' @slot columnName (character) The column name for which 
##' @slot displayName (character) The name to display when referring to this column in
##' output tables 
##' @slot type (One of "logical", "categorical", "numeric") The expected type of the column
##' @slot categories (NULL or factor) if \code{type=="categorical"} this contains the 
##' set of allowed categories, ordered as required by the user
##' @slot unit (character) The unit, if any, associated with the data in this column  
##' @seealso \code{\link{ColumnDef}}
##' @export
setClass("ColumnDef",
          slots = list(columnName = "character",
                       displayName = "character",
                       type = "character",
                       categories = "factor",
                       unit = "character"))


##' @name show
##' @rdname show-methods
##' @aliases show,ColumnDef-method
##' @export
setMethod("show","ColumnDef",
  function(object) {
    cat("Column Name:", object@columnName,"\n")
    cat("Display Name:", object@displayName,"\n")
    cat("Column Type:", object@type,"\n")
    
    if(object@type=="categorical"){
      cat("Categories:", levels(object@categories))
    }
    
    if(object@unit != ""){
      cat("Column unit:", object@unit,"\n")
    }
    
})



##' Constructor for \code{ColumnDef} object
##'
##' @param columnName (character) name of column as it appears in the raw data
##' @param type (character) how to interpret values in the column,
##'        e.g. "numeric" or "categorical" or "logical"
##' @param  categories (NULL or factor) expected levels in a categorical column. If this
##'        is entered as vector of character values, they will be coerced into factors     
##' @param  unit (character) unit of the column
##' @param displayName (character) name to use for column in outputs
##' @return An \code{ColumnDef} object
##' @seealso \code{\link{ColumnDef-class}} 
##' @examples
##' 
##' #Define columns corresponding to arm
##' #coercing categories into factors
##' armDef <- ColumnDef(columnName = "grp",
##'                    type = "categorical",
##'                    categories = c("patchOnly", "combination"))
##'
##' #Define columns corresponding to covariates
##' covariateDef <- list(
##'   ColumnDef(columnName = "age",
##'             type = "numeric",
##'             unit= "years"),
##'   
##'   ColumnDef(columnName = "race",
##'             type = "categorical",
##'             categories = factor(c("black", "hispanic", "other", "white"),
##'                                 levels=c("black", "white", "hispanic", "other"))))
##'
##' 
##' @export
ColumnDef <- function(columnName, type=NULL, categories=NULL, unit="", displayName = NULL){

  if (is.null(categories)){
    if (is.null(type) || type == "categorical" ){
      stop("If categories are not given, type must be defined and cannot be 'categorical'")
    }
  }
  else{
    # Categories are given, so must be categorical
    if (type!="categorical"){
      warning("Categories are given so type is assumed to be categorical")
    }
    type <- "categorical"
  }

  ensureIsCharacterLengthOne <- function(value, name){
    if(length(value)!=1){
      stop(paste0("'",name,"' should be a vector of length 1"))
    }
    if (class(value) != "character"){
      value <- as.character(value)
      warning(paste0("'", name, "' should be of type character. Coercing to character vector: ",
                     value))
    }
    return(value)
  }
  columnName <- ensureIsCharacterLengthOne(columnName, "columnName")
  unit <- ensureIsCharacterLengthOne(unit, "unit")
  

  # Type must be one of a defined set of possibilities
  supportedTypes <- c("numeric", "categorical", "logical")
  if (length(type)!= 1 || !type %in% supportedTypes){
    stop(paste0("Type must be one of: ",
                paste0(supportedTypes, collapse = ", ")))
  }

  # Check categories are correctly defined
  if (type == "categorical" && !is.factor(categories)){
      # Coerce categories to factor
      categories <- factor(categories,levels=categories)

      # Warn
      levelsAsChar <- paste0(levels(categories), collapse = ", ")
      message(paste0("If type is 'categorical', 'categories' should be a factor variable. Coercing to factor with levels: ",
                     levelsAsChar))
  }

  
  if (type != "categorical"){
    categories <- as.factor(NULL)
  }

  # If display not provided, default to column name
  if (is.null(displayName)){
    displayName <- columnName
  }
  else{
    displayName <- ensureIsCharacterLengthOne(displayName, "displayName")  
  }

  # Create object
  return(new("ColumnDef",
             columnName = columnName,
             displayName = displayName,
             type = type,
             categories = categories,
             unit = unit))
}


# Extract value of a slot from every element of list of columnDef objects
# Return as a vector.
listColumnDefSlot <- function(colDefs, thisSlot){
  if(length(colDefs)==0){
    character(0)
  }
  else{
    vapply(colDefs, function(x){slot(x, thisSlot)}, FUN.VALUE = slot(colDefs[[1]],thisSlot))
  }
}  
