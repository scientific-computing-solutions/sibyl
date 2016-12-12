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
##' @details See Vignette for further details
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
  validateSubjects(data, subjectCol)
  validateCovariates(data, covDef)
  validateSubgroups(data, subgroupDef, listColumnDefSlot(covDef, "columnName"))

  validateEndPointCols(data, endPointNames, timeCol, censorCol)
  validateEndPointVals(data, endPointNames, timeCol, censorCol)
  
  # Arm is always unique, so shouldn't be a list
  if (is.list(armDef)){armDef <- armDef[[1]]}

  # Ensure values in subgroup columns are logical
  for (col in listColumnDefSlot(subgroupDef, "columnName")){
    data[, col] <- as.logical(data[, col])
  }
  
  #Ensure timeCols are numeric 
  for(col in timeCol){
    data[,col] <- convertToNumeric(data, col)
  }
  
  #and censorCols are logical 
  for(col in censorCol){
    data[,col] <- convertToLogical(data, col)
  }
  
  #for each covariate convert to correct type.
  #Categorical
  for(cov in covDef){
    name <- cov@columnName
    data[,name] <- switch(cov@type,
                          "categorical"=factor(data[, name], levels=cov@categories),
                          "logical"=convertToLogical(data, name),
                          "numeric"=convertToNumeric(data, name))
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
  if(is.null(subgroupDef)) subgroupDef <- list()
  if(is.null(covDef)) covDef <- list()
  
  # Create SurvivalData object
  return(new("SurvivalData",
             subject.data = subject.data,
             armDef = armDef,
             endPoints = endPoints,
             endPointUnit = endPointUnit,
             subgroupDef = subgroupDef,
             covDef = covDef))
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
  
  dataSplitByArm <- split(object@subject.data, 
                          object@subject.data[,object@armDef@columnName])
  
  maxObserved <- vapply(dataSplitByArm, function(oneArmDf){
    max(oneArmDf[,object@endPoints[[endPointName]]$timeCol])
  },FUN.VALUE=numeric(1))
  
  min(maxObserved)
}

##' Calculate the maximum observed time 
##' for a given endPoint 
##' 
##' @param object (SurvivalData object)
##' @param endPointName (character) An endpoint associated with the
##' SurvivalData object 
##' @return (numeric) maximum observed time
##' @export
maxObservedTime <- function(object, endPointName){
  if(class(object)!="SurvivalData"){
    stop("object must be a SurvivalData object")  
  }
  
  if(length(endPointName)!=1 || !endPointName %in% names(object@endPoints)){
    stop(paste("Invalid endPointName argument should be one of",
               paste(names(object@endPoints),collapse=", ")))
  }
  max(object@subject.data[,object@endPoints[[endPointName]]$timeCol])
}