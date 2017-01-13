#validation routines for creation of SurvivalData object

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


validateEndPointCols <- function(data, endPointNames, timeCol, censorCol){
  
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
  
  
}


validateEndPointVals <- function(data, endPointNames, timeCol, censorCol){
  # Values in censor and time columns
  numEndPoints <- length(endPointNames)
  for (idx in seq(numEndPoints)){
    
    # Censor
    thisCensorCol <- censorCol[idx]
    
    #check can convert to logical
    thisCensorColData <- convertToLogical(data,thisCensorCol)
    
    
    # Time
    #check can convert to numeric
    thisTimeCol <- timeCol[idx]
    thisTimeColData <- convertToNumeric(data, thisTimeCol)
    
    if (any(thisTimeColData <  0 & !is.na(thisTimeColData))){
      stop("Event times must be non-negative")
    }
    
    #Check that either both or neither are NA
    if(any(xor(is.na(thisCensorColData), is.na(thisTimeColData)))){
      stop(paste("Subjects with an endpoint time must have a",
                 "censor indicator and vice versa. This is not the case",
                 "for endpoint",endPointNames[idx]))
    }
    
    #and that there is at least one non-NA
    if(all(is.na(thisCensorColData))){
      stop(paste("No data for endpoint:",endPointNames[idx]))
    }
    
  }
}

#convert the timeCol to numeric, it may be a character vector with empty strings
#throw error if cannot convert
convertToNumeric <- function(data, thisCol){
  thisColData <- data[, thisCol]
  
  if(class(thisColData)%in% c("factor","logical")){
    if(all(is.na(thisColData))) return(rep(as.numeric(NA,length(thisColData))))   
    
    stop(paste0("Values in '", thisCol,
                "' must be either numeric or empty and cannot be",
                " logical or a factor"))   
  }
  
  if(class(thisColData)=="character"){
    suppressWarnings(
      if(!all(nchar(trimws(thisColData))==0 | is.na(thisColData) |
              !is.na(as.numeric(thisColData)))){
        stop(paste0("Values in column '", thisCol,
                    "' must be either numeric or empty"))   
      }  
    )
    thisColData <- as.numeric(thisColData)
  }
  thisColData
}

#used to convert censorCol to logical, it may be numeric (0/1) or 
#a character vector with empty strings for missing data
#throw error if cannot convert
convertToLogical <- function(data,thisCol){ 
  thisColData <- data[, thisCol]
  
  if(class(thisColData)== "factor"){
    stop(paste0("Values in column '", thisCol,
                "' must be either boolean or 0/1 or empty not a factor"))   
  }
  
  if(class(thisColData)=="character"){
    if(!all(thisColData %in% c("0","1","TRUE","FALSE","") | is.na(thisColData))){
      stop(paste0("Values in column '", thisCol,
                  "' must be either boolean or 0/1 or empty"))  
    }
    
    thisColData <- ifelse(thisColData %in% c("0","FALSE"),FALSE,
                          ifelse(thisColData %in% c("1","TRUE"),TRUE,NA))
    
  }
  
  if(class(thisColData)%in% c("integer","numeric")){
    if(!all(thisColData %in% c(0,1) | is.na(thisColData))){
      stop(paste0("Values in column '", thisCol,
                  "' must be either boolean or 0/1 or empty"))  
    }
    thisColData <- as.logical(thisColData)
  }
  
  thisColData
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
  
  for(cov in covDef){
    if(cov@type == "categorical"){
      
      name <- cov@columnName
      categories <- cov@categories
      
      # convert to factor
      convertedValues <- factor(data[, name], levels=categories)
      
      #make sure all <NA>s are empty strings or NA 
      if(any(is.na(convertedValues) & !(is.na(data[,name]) | nchar(trimws(data[,name]))==0))){
        stop(paste0("Covariate '",name,"' is categorical but some values do not match the",
                    " defined categories {", paste0(categories, collapse=", "),"}"))
      }
    }
    
    if(cov@type == "logical") convertToLogical(data, cov@columnName)
    if(cov@type == "numeric") convertToNumeric(data, cov@columnName)
  }
  
  #displayNames must be unique
  displayNames <- listColumnDefSlot(covDef,"displayName")
  if(length(displayNames) != length(unique(displayNames))){
    stop("Covariate display names must all be different")   
  }
}


validateSubgroups <- function(data, subgroupDef, covarColNames){
  
  # Column names occur in raw data
  subgroupCols <- listColumnDefSlot(subgroupDef, "columnName")
  checkColumnNameValid(colnames(data), subgroupCols)
  
  #cannot have subgroup columns also as covariates
  if(any(subgroupCols %in% covarColNames)){
    stop("Cannot have a column as both a subgroup and a covariate column")
  }
  
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