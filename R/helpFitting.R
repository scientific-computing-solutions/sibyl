# Select rows of data frame with a given value in the subgroup column
extractSubgroup <- function(data, subgroup){

  # Pull out data set for this subgroup
  if (!is.na(subgroup)){
    # Select subset corresponding to the required subgroup
    isThisSubgroup <- data[, subgroup]
    data <- data[isThisSubgroup, ]
  }
  # else{
  # Don't need to do anything: use the whole data set
  # }

  return(data)
}

# Create the formula for fitting survival models to
# @param armAsFactor (logical) TRUE if the formula is to regress on arms,
#        FALSE otherwise
# @param covariates (vector of strings) names of covariates on which the
#        formula is to regress
# @param timeCol, censorCol - the time and censor columns associated with the 
# endpoint used to fit the model
# @param strata (vector of strings) names of covariates for stratification
# @return A formula object
survivalFormula <- function(armAsFactor,
                            covariates=character(0),
                            timeCol,
                            censorCol,
                            strata=character(0)){
  
  # Create default formula: regress (time, !censor) on a constant only
  my.formula <- formula(paste0("Surv(", timeCol, ", ", "!", censorCol, ") ~ 1"))
  
  # Add arms into regression model if specified
  if (armAsFactor){
    my.formula <- update.formula(my.formula, .~arm)
  }
  
  #append covariates to formula
  my.formula <- update.formula(my.formula,paste0(c(".~.",covariates),collapse="+ "))
  
  #append strata to formula
  if(length(strata)>0){
    strata_string <- paste("strata(",strata,")")
    my.formula <- update.formula(my.formula,paste0(c(".~.",strata_string),collapse="+ "))
  }
  
  return(my.formula)
}


removeNonPositiveTimes <- function(data, timeCol){
  timeIsPositive <- data[, timeCol] > 0
  if(!any(timeIsPositive)){
    stop("All subjects have time to event=0")  
  }
  
  if(!all(timeIsPositive)){
    warning("When fitting models subjects with time = 0 are not included")
    data <- data[timeIsPositive, ]
  }
  
  data
}

#check each arm contains an event
checkEventsPerArm <- function(data, censorCol){
  # Check that each arm sees at least one event 
  dataSplitByArm <- split(data, data[, "arm"])
  
  vapply(dataSplitByArm,function(x){
    sum(x[,censorCol])==nrow(x)
  },FUN.VALUE=logical(1))
}

# Check each arm contains valid data and group as needed for model fitting
checkValidDataPerArm <- function(data){
  
  # Check that each arm has subject with positive time
  dataSplitByArm <- split(data, data[, "arm"])
  
  vapply(dataSplitByArm,function(x){
    nrow(x)==0
  },FUN.VALUE=logical(1))
  
}
