#' @include survivalData.R
NULL

##' Class storing the fitted survival models for a specific subgroup, end point
##' and armAsFactor
##' @slot models (list) names of model distributions fitted
##' @slot covariates (character vector) names of covariates used
##' @slot survData (SurvivalData object) contains data for the selected subgroup
##'       only
##' @slot armAsFactor (logical) TRUE if arm is included in the model, FALSE if
##'       separate models are fit for each arm
##' @slot survFormula (formula) formula used to fit the model
##' @slot subgroup (character string or NA) name of subgroup fitted
##' @slot endPointDef (list) defines time and censor column for the fitted
##'       endpoint
##' @slot endPoint (character) name of the end-point to be fitted
##' see vignette and \code{SurvivalModel} for further iunformation
##' @export
setClass("SurvivalModel",
         slots = list(models = "list",
                      covariates = "character",
                      survData = "SurvivalData",
                      armAsFactor = "logical",
                      survFormula = "formula",
                      subgroup = "character",
                      endPointDef = "list",
                      endPoint="character"))

##' Method to create the \code{SurvivalModel} object
##' @name fitModels
##' @rdname fitModels-methods
##' @param object \code{SurvivalModel} contains data to be fitted
##' @param ... additional arguments for specific instances of this generic
##' @return A \code{SurvivalModel} object
setGeneric( "fitModels", function(object,...) standardGeneric( "fitModels" ))

##' @rdname fitModels-methods
##' @aliases fitModels,SurvivalModel-methods
##' @param armAsFactor (logical) TRUE if arm is to be included in the model,
##'        FALSE if separate models are to be fitted for each arm
##' @param covariates (character vector) names of covariates to use
##' @param endPoint (character) name of the end-point to be fitted
##' @param subgroup (character string or NA) name of subgroup to fit
##' @param model (character vector) names of model distributions to fit
##' @param modelOptions (named list) extra options to be passed to the fitting functions.
##' For example modelOptions=list(spline=list(k=4)) will pass the argument k=4 to flexsurvSpline
##' @param preferredPackage (character) name of the preferred package in which
##'        to look for the specified model
##' @export
setMethod("fitModels", signature(object="SurvivalData"),
function(object,
         armAsFactor = TRUE,
         covariates = character(0),
         endPoint,
         subgroup = as.character(NA),
         model = c("exponential", "weibull", "llogis", "lnorm", "gompertz"),
         modelOptions = NULL,
         preferredPackage = getDefaultPackage()){

  # Check the required covariates are present in the data
  hasCovariates(object, covariates)

  if (!is.logical(armAsFactor)){
    stop("armAsFactor must be logical")
  }

  if (!all(vapply(covariates, is.character, FUN.VALUE = TRUE))){
    stop("covariates must be a vector of character strings")
  }

  if (!is.character(endPoint)){
    stop("endPoint must be a character string")
  }

  if (!is.na(subgroup) && !is.character(subgroup)){
    stop("subgroup must be either NA or a character string")
  }

  #  subgroup matches one of the subgroups defined in the data
  allowedSubgroups <- listColumnDefSlot(object@subgroupDef,"columnName")
  if(!is.na(subgroup) && !subgroup %in% allowedSubgroups){
    stop(paste("subgroup argument must be one of: ",paste(allowedSubgroups, collapse = ", ")))
  }

  if (any(vapply(list(armAsFactor, endPoint, subgroup), function(x){length(x) != 1}, FUN.VALUE = TRUE))){
    stop("Can only fit for single values of subgroup, endPoint or armAsFactor")
  }

  
  #warn if modelOptions given for a model not being fitted
  if(any(!names(modelOptions) %in% model)){
    warning("modelOptions given for a model which is not being fitted")
  }
  

  #endPoint matches one of the known endpoints (object@endPoints)
  # Look up relevant column names for the specified endpoint
  if(!endPoint %in% names(object@endPoints)){
    stop("Invalid endpoint argument, must be one of: ", paste(names(object@endPoints), collapse=", "))
  }
  endPointDef <- object@endPoints[[endPoint]]

  # Create new SurvivalData object with data for this subgroup only
  survData <- object
  survData@subject.data <- extractSubgroup(object@subject.data, subgroup)

  # Create formula to fit model
  formulaToFit <- survivalFormula(armAsFactor,
                                  covariates,
                                  timeCol = endPointDef[["timeCol"]],
                                  censorCol = endPointDef[["censorCol"]])

  # Fit the specified models
  fittedModels <- internalFitModels(model,
                                    modelOptions,
                                    preferredPackage,
                                    formulaToFit,
                                    survData,
                                    armAsFactor,
                                    covariates,
                                    endPointDef)
  if(length(fittedModels)==0){
    stop("No models could be fitted!")
  }
  

  # Create output object
  new("SurvivalModel",
      models = fittedModels,
      covariates = covariates,
      survData = survData,
      armAsFactor = armAsFactor,
      survFormula = formulaToFit,
      subgroup = subgroup,
      endPointDef = endPointDef,
      endPoint=endPoint)
})


# Helper function to do model fitting: see fitModels and addModel for examples of use
internalFitModels <- function(model, modelOptions, preferredPackage, formulaToFit, survData, armAsFactor, covariates, endPointDef){

  # Create one data set per grouping required for model fitting (e.g. per arm,
  # or all data together)
  dataSets <- groupDataForFitting(survData@subject.data, endPointDef[["timeCol"]], 
                                  endPointDef[["censorCol"]], armAsFactor, getArmNames(survData))

  # Create model-fitting functions for supported distributions
  modelFitters <- lapply(model, function(x){getModelFitter(x, modelOptions[[x]], preferredPackage, formulaToFit)})

  # Fit each model to each data set
  fittedModels <- lapply(modelFitters, function(f){
    
    #set warnings as errors to refuse to fit models if there are warnings
    options("warn"=2)
    on.exit(options("warn"=0))
    
    tryCatch(fitModelToData(f, dataSets),
               error=function(msg){msg$message})
  
  })

  # Name the models with the matching distribution
  names(fittedModels) <- model

  #Next remove models which produced errors
  #successful models
  successfulModelNames <- NULL
  successfulModels <- list()
  
  #for loop not lapply as need to build up successfulModelNames
  for(modelName in model){
    if(class(fittedModels[[modelName]]) !="character"){
      successfulModelNames <- c(successfulModelNames, modelName)
      successfulModels[[modelName]] <- fittedModels[[modelName]]
    } 
    else warning(paste("Sibyl could not fit the model", modelName, "because:", fittedModels[[modelName]]))
  }
  
  names(successfulModels) <- successfulModelNames
  return(successfulModels)
}



##' Method to remove fitted models from a SurvivalData object
##' @name removeModel
##' @rdname removeModel-methods
##' @param object (SurvivalModel object) contains list of already fitted models
##' @param ... additional arguments for specific instances of this generic
##' @return a \code{SurvivalModel} object
setGeneric("removeModel", function(object, ...) standardGeneric("removeModel"))

##' @rdname removeModel-methods
##' @aliases removeModel,SurvivalModel-methods
##' @param modelName (list of character strings) names of models to be removed
##' @export
setMethod("removeModel", signature(object = "SurvivalModel"),
function(object, modelName){

  # Remove each specified model in turn from list of models
  models <- object@models
  for (thisModelName in modelName){
    models <- removeOneModel(models, thisModelName)
  }
  object@models <- models

  return(object)
})

#remove one model (with name modelName) from list of models, if suppressError is
#false then do not complain if last model is removed
removeOneModel <- function(models, modelName, suppressError=FALSE){
  # Match model names case-insensitively
  modelName <- tolower(modelName)
  
  if (modelName %in% names(models)){
    
    if (length(models) == 1 && !suppressError){
      stop(paste0("'", modelName, "' is the only model in this SurvivalModel object. Removing it would leave an empty object"))
    }
    else{
      # Remove the model from the object
      models[[modelName]] <- NULL
    }
  }
  else{
    warning(paste0("Trying to remove non-existent model '", modelName,
                   "'from SurvivalModel object has no effect."))
  }
  return(models)
}


##' Method to add extra models to be fitted
##' @name addModel
##' @rdname addModel-methods
##' @param object (SurvivalModel object) contains list of already fitted models
##' @param ... additional arguments for specific instances of this generic
##' @return a \code{SurvivalModel} object
setGeneric("addModel", function(object, ...) standardGeneric("addModel"))

##' @rdname addModel-methods
##' @aliases addModel,SurvivalModel-methods
##' @param modelName (list of character strings) names of models to be added
##' @param modelOptions (named list) extra options to be passed to the fitting functions.
##' For example modelOptions=list(spline=list(k=4)) will pass the argument k=4 to flexsurvSpline
##' @param preferredPackage (character string) name of the preferred package
##'        in which to look for the specified model
##' @param suppressOverwriteWarning (logical) should be FALSE so that a warning is displayed
##' when overwriting existing model  
##' @export
setMethod("addModel", signature = c("SurvivalModel"),
function(object, modelName, modelOptions=NULL, preferredPackage = getDefaultPackage(),suppressOverwriteWarning=FALSE){

  # Match model names case-insensitively
  modelName <- lapply(modelName, tolower)

  # Warn if over-write existing models
  isDuplicate <- vapply(modelName, function(n){n %in% names(object@models)}, FUN.VALUE = FALSE)
  if (any(isDuplicate) && !suppressOverwriteWarning){
    duplicateNames <- modelName[isDuplicate]
    warning(paste0("The following models have already been fitted and (if new models are successfully 
                   fitted then) they will be over-written: ",
                   paste0(duplicateNames, collapse = ", ")))
  }

  # Fit the specified models
  fittedModels <- internalFitModels(modelName,
                                    modelOptions,
                                    preferredPackage,
                                    object@survFormula,
                                    object@survData,
                                    object@armAsFactor,
                                    object@covariates,
                                    object@endPointDef)

  # Add existing models to fittedModels removing those to be overwritten
  for(x in names(fittedModels)){
    if(x %in% names(object@models)){
      object@models <- removeOneModel(object@models, x, suppressError=TRUE)
    }
    object@models[[x]] <- fittedModels[[x]]
  }
  
  return(object)
})


getDefaultPackage <- function(){
  "flexsurv"
}


# Return model-fitting function from the specified package
getModelFitter <- function(model, modelOptions, preferredPackage, formulaToFit){

  # Define which distributions are supported by which fitting package
  supportedFitters <- list(survival = c("weibull",
                                        "exponential",
                                        "gaussian",
                                        "logistic",
                                        "lognormal",
                                        "loglogistic"),
                           flexsurv = c("exponential",
                                        "gengamma",
                                        "gengamma.orig",
                                        "genf",
                                        "genf.orig",
                                        "weibull",
                                        "gamma",
                                        "exp",
                                        "llogis",
                                        "lnorm",
                                        "gompertz",
                                        "spline"))

  # Standardise name of preferred package
  preferredPackage <- tolower(preferredPackage)
  nameLength <- nchar(preferredPackage)
  if (preferredPackage == substr("flexsurv", 1, nameLength)){
    preferredPackage <- "flexsurv"
  }

  else if(preferredPackage == substr("survival", 1, nameLength)){
    preferredPackage <- "survival"
  }

  else{
    stop(paste0("'", preferredPackage, "' is not a recognised model-fitting package"))
  }

  # Check if model is supported by this package
  if (model %in% supportedFitters[[preferredPackage]]){
    # Model is supported by preferred package
    packageToUse <- preferredPackage
  }

  else{
    # Model isn't supported by the specified package. Try the default package.
    defaultPackage <- getDefaultPackage()

    if (model %in% supportedFitters[[defaultPackage]]){

      # Model is available in default package, but not in specified one
      packageToUse <- defaultPackage
      warning(paste0("Model '", model, "' is not part of package '", preferredPackage,
                     "'; using default package '", defaultPackage, "' instead"))

    }
    else{
      # Model isn't available in either package
      stop(paste0("Model '", model, "' is not part of package '", preferredPackage,
                  "' or default package '", defaultPackage, "'"))
    }
  }

  # Get fitting function from package
  if (packageToUse == "flexsurv"){
    if(model == "spline"){
      packageFunc <- flexsurvspline
    }
    else{
      packageFunc <- flexsurvreg
    }
  }

  else if (packageToUse == "survival"){
    packageFunc <- survreg
  }

  else{
    stop(paste0("'", packageToUse, "' is not a recognised model-fitting package"))
  }

  # Construct model-fitting function
  args <- list(formula=formulaToFit)
  
  if(model != "spline"){
    args$dist <- model
  }
  
  args <- c(args, modelOptions)
  
  fitter <- function(d){
    args$data <- d
    do.call(packageFunc, args)
  }
  
  return(fitter)
}


# Apply 1 model fitter function to all data sets, and name the list of results
fitModelToData <- function(modelFitter, dataSets){

    # Fit the model to each data set in turn
    fittedModel <- lapply(dataSets, modelFitter)

    # Name the data sets
    names(fittedModel) <- names(dataSets)

  return(fittedModel)
}


# Group data either by arm or all together
groupDataForFitting <- function(data, timeCol, censCol, armAsFactor, armNames){

  # Remove data with zero end point times
  data <- removeNonPositiveTimes(data,timeCol)

  # Check that all arms contain  data
  armHasNoData <- checkValidDataPerArm(data)
  
  if (any(armHasNoData)){
    stop(paste0("The following arms contain no data: ",
                paste(armNames[armHasNoData], collapse = ", ")))
  }
  
  #Check that all arms contain events
  armHasNoEvents <- checkEventsPerArm(data, censCol) 
  
  if (any(armHasNoEvents)){
    stop(paste0("The following arms have no events: ",
                paste(armNames[armHasNoEvents], collapse = ", ")))
  }
  
  
  if (armAsFactor){
    # Use whole data set
    dataSets <- list(data)

    # There is no arm name
    names(dataSets) <- c(NA)
  }
  else{
    # Separate the data set by arm
    dataSets <- split(data, data[, "arm"])
  }

  dataSets
}

