# This file contains the public functions associated with survivalData object

##' @include columnDef.R
NULL

##' Class representing the data to be used by Sibyl for model fitting
##' and extrapolating.
##' @slot subject.data (data frame), one row per subject
##' with unique subject ID, their covariate values, treatment arm,
##' time to event and censor indicators for multiple endpoints and
##' indicators as to which subgroups the subject belongs.
##' Colums names are described by the data in other slots (although the treatment arm
##' column is named "arm").
##' @slot armDef (\code{ColumnDef}) object describing the treatment arm column and its categories
##' (the control group should be first) in the subject.data data frame
##' @slot endPoints (list of lists) each endpoint is the name of a two element list, "timeCol" - the time
##' to event for the given endpoint and "censorCol" whether the subject was censored (\code{TRUE}) or not (\code{FALSE}).
##' @slot endPointUnit (days, months, years) The units for the endpoints
##' @slot subgroupDef (list of \code{ColumnDef}) definitions of the columns in the subject.data data frame
##' describing subgroups
##' @slot covDef (list of \code{ColumnDef}) definitions of columns in the
##' subject.data data frame describing covariates
##' @seealso \code{\link{SurvivalData}}
##' @export
setClass("SurvivalData",
          slots = list(subject.data = "data.frame",
                       armDef = "ColumnDef",
                       endPoints = "list",
                       endPointUnit = "character",
                       subgroupDef = "list",
                       covDef = "list"))


##' Constructor for \code{SurvivalData} object
##'
##' @param data (data frame) raw data to be formatted into SurvivalData object
##' @param armDef (columnDef object) definition of the data column that
##'        specifies the trial arm each subject was on
##' @param subjectCol (string) name of the data column containing the unique subject
##'        identifier
##' @param covDef (list of columnDef objects, default=NULL) definition of each
##'        required covariate, including name, type, unit, etc.
##' @param subgroupDef (list of columnDef objects, default=NULL) definition of
##'        data columns that indicate membership of a subgroup.
##' @param endPointNames (vector of strings) names of the endpoints
##' @param censorCol (vector of strings) names of the data columns that specify
##'        which subjects are censored for a particular end point. The N-th
##'        element should correspond to the N-th end point.
##' @param timeCol (vector of strings) names of the data columns that specify
##'        the times at which subjects reached a particular end point. The N-th
##'        element should correspond to the N-th end point.
##' @param endPointUnit ("days", "months" or "years" - default "months") The unit of time
##' for the endPoint time columns
##' @return A \code{SurvivalData} object
##' @seealso \code{\link{SurvivalData-class}}
##' @export
SurvivalData <- function(data,
                         armDef,
                         subjectCol,
                         covDef=NULL,
                         subgroupDef=NULL,
                         endPointNames,
                         censorCol,
                         timeCol,
                         endPointUnit=c("days","months","years")[2]){

  if (class(data)!="data.frame"){
    stop("'data' must be a data frame")
  }

  # Where there is a column definition, ensure they are lists
  ensureIsList <- function(x){
    if (!is.null(x) && class(x) != "list"){
      x <- list(x)
    }
    return(x)
  }
  covDef <- ensureIsList(covDef)
  subgroupDef <- ensureIsList(subgroupDef)

  # Validation
  if(length(endPointUnit)> 1 || !endPointUnit %in% c("days","months","years")){
    stop("endPointUnit must be one of 'days', 'months' or 'years'")
  }
  
  validateArm(data, armDef)
  validateEndPoints(data, endPointNames, timeCol, censorCol)
  validateSubjects(data, subjectCol)
  validateCovariates(data, covDef)
  validateSubgroups(data, subgroupDef)

  # Arm is always unique, so shouldn't be a list
  if (is.list(armDef)){armDef <- armDef[[1]]}

  # Ensure values in subgroup columns are logical
  for (col in listColumnDefSlot(subgroupDef, "columnName")){
    data[, col] <- as.logical(data[, col])
  }

  # Extract relevant part of raw data
  subject.data <- data[, c(subjectCol,
                           timeCol,
                           censorCol,
                           armDef@columnName,
                           listColumnDefSlot(covDef, "columnName"),
                           listColumnDefSlot(subgroupDef, "columnName"))]

  # Rename relevant columns with standard names
  allNames = names(subject.data)
  names(subject.data)[allNames == subjectCol] <- "subject"
  names(subject.data)[allNames == armDef@columnName] <- "arm"
  #set armDef's column name to arm
  armDef@columnName <- "arm"


  # Convert arm column to factor variable, preserving order of categories
  subject.data$arm <- factor(subject.data$arm, levels = armDef@categories)

  # Zip up (time, censor) pairs for each type of end point
  endPoints <- lapply(seq_along(endPointNames),
                      function(i){list(timeCol = timeCol[i],
                                       censorCol = censorCol[i])})
  names(endPoints) <- endPointNames

  #set null to empty list to satisfy S4's demands
  if(is.null(subgroupDef)){
    subgroupDef <- list()
  }

  if(is.null(covDef)){
    covDef <- list()
  }

  # Create SurvivalData object
  return(new("SurvivalData",
             subject.data = subject.data,
             armDef = armDef,
             endPoints = endPoints,
             endPointUnit = endPointUnit,
             subgroupDef = subgroupDef,
             covDef = covDef))
}


checkColumnNameValid <- function(dataColNames, colNames){

  isInvalid <- !vapply(colNames,
                       function(x){!is.null(x) && !is.na(x) && (x %in% dataColNames)},
                       FUN.VALUE = FALSE)

  if (any(isInvalid)){
    stop(paste0("The following column names are invalid or not found in the raw data: '",
                paste0(colNames[isInvalid], collapse = ", ")))
  }
}


validateArm <- function(data, armDef){

  # Handle arm defined as list
  if (is.list(armDef)){
    if (length(armDef) == 1){
      armDef <- armDef[[1]]
    }
    else{
      stop(paste0("Multiple arm columns defined"))
    }
  }

  # Only 1 arm column defined
  armCol <- unique(armDef@columnName)
  if (length(armCol) > 1){
    stop(paste0("Multiple arm columns defined: ", paste0(armCol, collapse = ", "),
                ". Arms must be defined by categorical values within a single column."))
  }

  # Arm column name is valid
  checkColumnNameValid(colnames(data), armDef@columnName)

  # Arm definition has at least two levels (one per arm)
  if (length(armDef@categories) < 2){
    stop(paste0("Definition of arm column '", armDef@columnName,
                "' contains fewer than two levels. There must be one level per treatment group and at least two groups."))
  }

  # Values in arm column must match defined categories
  categoriesNotFound <- setdiff(armDef@categories, data[, armCol])
  if (length(categoriesNotFound) > 0){
    stop(paste0("The following arm categories are not found in the data: ",
                paste0(categoriesNotFound, collapse = ", ")))
  }
  unexpectedValues <- setdiff(data[, armCol], armDef@categories)
  if (length(unexpectedValues) > 0){
    stop(paste0("The following values in the arm column are not defined as arm categories: ",
                paste0(unexpectedValues, collapse = ", ")))
  }

}


validateEndPoints <- function(data, endPointNames, timeCol, censorCol){

  # Equal numbers of event time and censor columns
  numEndPoints <- length(endPointNames)
  if (length(censorCol) != numEndPoints){
    stop(paste("Different number of censor columns and end points"))
  }
  if (length(timeCol) != numEndPoints){
    stop(paste("Different number of time columns and end points"))
  }

  
  if(length(unique(endPointNames)) != numEndPoints){
    stop("Endpoint names must be unique")
  }
  
  # Unique event time and censor columns
  if (length(unique(censorCol)) != numEndPoints){
    stop(paste("Repeated censor columns. Each censor column must correspond to one type of end point"))
  }
  if (length(unique(timeCol)) != numEndPoints){
    stop(paste("Repeated time columns. Each time column must correspond to one type of end point"))
  }

  # Column names occur in the data
  checkColumnNameValid(colnames(data), timeCol)
  checkColumnNameValid(colnames(data), censorCol)

  # Values in censor and time columns
  for (idx in seq(numEndPoints)){

    # Censor
    thisCensorCol <- censorCol[idx]

    if (all(data[, thisCensorCol] %in% c(0,1))){
      # Ensure logical
      data[, thisCensorCol] <- data[, thisCensorCol] > 0.5
    }
    else{
      stop(paste0("Values in censor column '", thisCensorCol, "' must be either boolean or 0/1"))
    }

    # Time
    thisTimeCol <- timeCol[idx]
    if (!any(is.numeric(data[, thisTimeCol]))){
      stop("Event times must be numeric")
    }
    if (any( data[, thisTimeCol] <  0)){
      stop("Event times must be non-negative")
    }
  }
}


validateSubjects <- function(data, subjectCol){

  checkColumnNameValid(colnames(data), subjectCol)

  subjectIds <- data[, subjectCol]
  isDuplicate <- duplicated(subjectIds)

  if(any(isDuplicate)){
    stop(paste0("Subject ID must be unique. The following IDs are repeated: ",
                paste0(subjectIds[isDuplicate], collapse = ", ")))
  }
}


validateCovariates <- function(data, covDef){

  checkColumnNameValid(colnames(data), listColumnDefSlot(covDef, "columnName"))

  # Categorical covariates
  isCategorical <- listColumnDefSlot(covDef, "type") == "categorical"
  idxCategorical <- which(isCategorical)
  for (idx in idxCategorical){

    name <- covDef[[idx]]@columnName
    categories <- covDef[[idx]]@categories

    # Ensure factor
    if (!is.factor(data[, name])){
      data[, name] <- as.factor(data[, name])
    }

    # Values must match defined categories
    uniqueValues <- levels(data[, name])
    unmatchedValues <- setdiff(uniqueValues, categories)

    if (!all(uniqueValues %in% categories)){
      stop(paste0(" Covariate '", name, "' is categorical but values {",
                  paste0(unmatchedValues, collapse=", "),
                  "} do not match the defined categories {",
                  paste0(categories, collapse=", "),
                  "}"))
    }
  }
  
  #displayNames must be unique
  displayNames <- listColumnDefSlot(covDef,"displayName")
  if(length(displayNames) != length(unique(displayNames))){
    stop("Covariate display names must all be different")   
  }
}


validateSubgroups <- function(data, subgroupDef){

  # Column names occur in raw data
  subgroupCols <- listColumnDefSlot(subgroupDef, "columnName")
  checkColumnNameValid(colnames(data), subgroupCols)

  # Values in subgroup columns are logical or numeric 0/1
  containsNa <- vapply(subgroupCols,
                       function(col){any(is.na(data[, col]))},
                       FUN.VALUE = FALSE)
  isValid <- vapply(subgroupCols,
                    function(col){all(is.logical(data[, col])) || all(data[, col] %in% c(0,1))},
                    FUN.VALUE = FALSE)

  if (any(containsNa)){
    stop(paste0("NA values occur in the following subgroup columns: ",
                paste0(subgroupCols[containsNa], collapse = ", ")))
  }

  if (any(!isValid)){
    stop(paste0("Values in the following subgroup columns are neither logical nor numeric 0/1: ",
                paste0(subgroupCols[!isValid], collapse = ", ")))
  }
  
  #displayNames must be unique
  displayNames <- listColumnDefSlot(subgroupDef,"displayName")
  if(length(displayNames) != length(unique(displayNames))){
    stop("Subgroup display names must all be different")   
  }

}


##' Show methods for Sibyl objects
##' @name show
##' @rdname show-methods
##' @aliases show,SurvivalData-method
##' @param object (SurvivalData object)
##' @export
setMethod("show", "SurvivalData",
  function(object){
    cat(str(object@subject.data),"\n")

    if(length(object@covDef) > 0){
      cat("Covariate columns:\n")
      print(object@covDef)
      cat("\n")
    }

    if(length(object@armDef) > 0){
      cat("Arm column:\n")
      print(object@armDef)
      cat("\n\n")
    }

    if(length(object@subgroupDef) > 0){
      cat("Subgroup columns:\n")
      print(object@subgroupDef)
      cat("\n")
    }

    if(length(object@endPoints) > 0){
      cat("End-point columns:\n")
      print(object@endPoints)
      cat("\n")
    }
  }
)


#check that the covariates in the vector covariates exist as column names in the Survival Data object
#throw error if covariate not found
hasCovariates <-  function(object, covariates){
  xx <- vapply(covariates, function(i_cov){
          i_cov %in% colnames(object@subject.data)
          },FUN.VALUE = logical(1))
  if(!all(xx)){
    stop(paste0("Covariate ", paste(covariates[xx==FALSE], collapse=", "), " not found!"))
  }
  TRUE
}


#Get the armnames from a SurvivalData object
getArmNames <- function(survData){
  survData@armDef@categories
}



##' Calculate the maximum observed time for each arm
##' for a given endPoint and output the minimum
##' 
##' @param object (SurvivalData object)
##' @param endPointName (character) An endpoint associated with the
##' SurvivalData object 
##' @return (numeric) minimum of the maximum observed time
##' for each arm
##' @export
minOfMaxObserved <- function(object, endPointName){
  if(class(object)!="SurvivalData"){
    stop("object must be a SurvivalData object")  
  }
  
  if(length(endPointName)!=1 || !endPointName %in% names(object@endPoints)){
    stop(paste("Invalid endPointName argument should be one of",
               paste(names(object@endPoints),collapse=", ")))
  }
  
  dataSplitByArm <- split(object@subject.data, object@subject.data[,object@armDef@columnName])
  
  maxObserved <- vapply(dataSplitByArm, function(oneArmDf){
    max(oneArmDf[,object@endPoints[[endPointName]]$timeCol])
  },FUN.VALUE=numeric(1))
  
  min(maxObserved)
}
