#' @include survivalModels.R
NULL

##' Method to create life table
##' @name createLifeTable
##' @rdname createLifeTable-methods
##' @param object (SurvivalModel object) contains fitted models
##' @param ... additional arguments for specific instances of this generic
##' @return The lifetable
setGeneric("createLifeTable", function(object,  ...) standardGeneric("createLifeTable"))

##' @rdname createLifeTable-methods
##' @aliases createLifeTable,SurvivalModel-methods
##' @param times (numeric vector) times at which to evaluate life table
##' @param modelToUse (character) name of model to base life table upon. If
##'        NULL, the model with lowest AIC is chosen
##' @param Nsim (numeric) number of simulations used to estimate the average
##'        curve - only used when models includes covariates
##' @param class (data.frame or default flexTable) the class of the created table
##' @param digits (numeric default 3) the number of decimal places to round the
##' FlexTable output to 
##' @param seed (numeric, default NULL) if not NULL then set random seed (although it will only be
##' used if models include covariates)
##' @param B (integer) Only used when no covariates in model. See summary.flexsurvreg 
##' Number of simulations from the normal asymptotic distribution of the estimates used 
##' to calculate confidence intervals. Decrease for greater speed at the expense of accuracy, 
##' or set B=0 to turn off calculation of CIs.
##' @details If the models include covariates a simulation procedure is required to generate averaged survival
##' curves. If the models do not include covariates then \code{summary.flexsurvreg} is used 
##' @export
setMethod("createLifeTable", signature(object="SurvivalModel"),
  function(object, times, modelToUse = NULL, Nsim=500, class=c("data.frame","FlexTable")[2], digits=3,
           seed=NULL, B=1000){

    if(!is.null(seed)){
      set.seed(seed)
    }
    
    if(is.null(modelToUse)){
      modelToUse <- getBestModel(object)
    }
    
    tableHeading <- modelToUse
    

    # Check required model has been fitted
    if (!any(vapply(names(object@models),
       function(x){x == modelToUse}, FUN.VALUE = TRUE))){
         stop(paste0("Model '", modelToUse, "' has not been fitted"))
       }

    # Validate times
    if (any(!is.numeric(times) | times < 0)){
      stop("Times must be numeric and non-negative")
    }

    # Extract fits of this model to each arm
    modelFits <- object@models[[modelToUse]]

    # Construct life table for this model
    lifeTable <- calcLifeTable(modelFits, object@survData@subject.data, object@endPointDef, times, modelToUse, Nsim,
                               object@armAsFactor, length(object@covariates)!=0, getArmNames(object@survData), B)

    if(class=="data.frame"){
      return(lifeTable)
    }
    
    #create FlexTable
    arms <- as.character(getArmNames(object@survData))
    numRows <- length(times)
    numCols <- 1+2*length(arms)
    
    MyFTable <- MyFTable <- FlexTable(numrow=numRows,numcol=numCols, 
                                      body.par.props=parProperties(text.align="right"),
                                      header.text.props = textProperties(font.weight = "bold"),
                                      body.cell.props = cellProperties(padding.right=1))
    
    #Add data to table
    MyFTable[1:numRows,1] <- times
    
    oneTime <- split(lifeTable, lifeTable$t)
    lapply(seq_along(oneTime), function(i){
      df <- oneTime[[i]]
      ans <- round(c(df$KM,df$S), digits=digits)
      MyFTable[i,2:numCols] <- ifelse(is.na(ans),"", ans)
    })
    
    #Add headers
    hR <- FlexRow(c("time", rep(arms,2)),
                  par.properties=parProperties(text.align="center",padding=1),
                  text.properties = textProperties(font.weight = "bold"))
    
    hR2 <- FlexRow(c("", "Kaplan Meier",getDistributionDisplayNames(tableHeading)),
                   colspan = c(1,length(arms), length(arms)),
                   par.properties=parProperties(text.align="center",padding=1),
                   text.properties = textProperties(font.weight = "bold"))
    
    MyFTable <- addHeaderRow(MyFTable,hR2)
    MyFTable <- addHeaderRow(MyFTable,hR)

    MyFTable
})



##' Output the model distribution with the lowest AIC
##' 
##' @details If armAsFactor is FALSE and separate models are
##' fitted then for each distribution calculate the sum of the
##' AIC and returns the distribution with the lowest sum of AIC 
##' @param object (SurvivalModel)
##' @return The name of the model with the lowest AIC
##' @export 
getBestModel <- function(object){
  
  if(class(object)!="SurvivalModel"){
    stop("Argument must be a Survival Model object")
  }
  
  icTable <- internalGetICTable(object)
  #generate AIC sums
  AICSums <- by(icTable,icTable$Model,function(x){sum(x$AIC)})
  #find minimum
  AICSumMin <- min(AICSums)
  #find name of model which model has the minimum value
  names(which(as.list(AICSums)== AICSumMin)[1])
}

# modelFits: list of model fit to each arm
# times: vector, times at which to evaluate the life table
# modelToUse: string; name of the parametric model used
#Nsim number of simulations (if simulating)
#armAsFactor: logical is the model included in the arm?
#useCovariates: logical does the model include covariates
#B: see summary.flexsurvreg 
calcLifeTable <- function(modelFits, subjectData, endPointDef, times,
                          modelToUse, Nsim, armAsFactor, useCovariates, armNames,
                          B){

  # Ensure times in ascending order
  times <- sort(times)

  # Calculate KM-lifetable with times given by survfit
  kmValuesByArm <- calcKMLifeTable(subjectData,endPointDef)

  #Next convert the standard KM table to one
  #whose time points match the times input by the user
  KMLifeTables <- calcInterpolatedKMLifeTable(kmValuesByArm, times)

  #Now calculate life table for given modelfits:

  #Apply calcParametricLifeTable to each {modelFits, dataByArm} pair,
  #if there is only one modelFit (as it is was fitted with arm as factor)
  #it is reused for each arm
  if(useCovariates){
    #split data by arm
    dataByArm <- split(subjectData, subjectData[, "arm"])
    parametricLifeTables <- mapply(calcParametricLifeTable, mod=modelFits, oneArmData=dataByArm,
                              MoreArgs=list(times=times,Nsim=Nsim), SIMPLIFY = FALSE)
  }
  else{
    parametricLifeTables <- mapply(calcLifeTableNoCovariates, mod=modelFits, armName=armNames,
                                   MoreArgs=list(times=times,armAsFactor=armAsFactor,B=B), SIMPLIFY = FALSE)  
  }
    
  names(parametricLifeTables) <- names(KMLifeTables)

  #merge the parametric and KM lifetables together
  mergedTables <- mapply(merge, KMLifeTables, parametricLifeTables,
                          MoreArgs=list(by="t"), SIMPLIFY=FALSE)

  #bind the different arms into the same data frame
  results <- data.frame(do.call("rbind",mergedTables))
  rownames(results) <- NULL

  results
}


# Calculate KM-lifetables with times given by survfit
calcKMLifeTable <- function(subjectData,endPointDef, outputCI=FALSE){
  # Create formula for fitting KM curv. Note:
  #   - This is *not* the formula used to originally fit the model, because
  #     that (may) depend on covariates, whereas for the KM curve we don't
  #     want these.
  #   - Arm is always factor because then survfit will fit each arm
  #     separately, without having to explicitly loop over arms. We can then
  #     separate out results for different arms.
  formulaToFit <- survivalFormula(armAsFactor = TRUE,
                                  timeCol = endPointDef[["timeCol"]],
                                  censorCol = endPointDef[["censorCol"]])

  # Estimate KM curve (for all arms at once)
  kmFit <- survfit(formulaToFit, data = subjectData)

  # Use extractCumHazData to split up data points by arm
  kmValuesByArm <- extractCumHazData(kmFit,levels(subjectData$arm),outputCI)

  # Name values according to arm
  names(kmValuesByArm) <- levels(subjectData$arm)

  kmValuesByArm
}


#given a list of KM life tables (one per arm) with times from survfit
#output interpolated KM table with times given by the input parameter times
calcInterpolatedKMLifeTable <- function(kmValuesByArm, times){

  KMLifeTables <- lapply(kmValuesByArm,function(oneArmDataFrame){

  #interpolate for each time point t, the value of S:
  #P(have not had event by time t)
  S <- vapply(times, function(t){
    getLastPoint(t, oneArmDataFrame)
  }, FUN.VALUE=numeric(1))

  #convert to data frame
  data.frame(t=times,
             KM=S,
             Arm=rep(oneArmDataFrame$Arm[1],length(S)))

  })
}


# Interpolate KM curve to a given time, taking the last point
getLastPoint <- function(t, km){

  # Find the last data point before x
  if(t > max(km$t)){
    return(as.numeric(NA))
  }
  
  idxBefore <- which(km$t <= t)
  if (length(idxBefore) > 0){
    idx <- max(idxBefore)
  }
  else{
    # If there is no previous data point, return the first one
    idx <- 1
  }

  return(km$S[idx])
}


#Calculate the parametric life table when the model has no covariates
#mod: The (flexSurv) model used to fit this arm, from which to generate the lifetable
#armName: The name of a single arm
#times: The times
#armAsFactor: logical whether arm is a covariate in the model
#conf.int: whether to output the confidence interval?
#level: The percentiles to calculate for the confidence interval of S
#returns a data frame with columns:
# t:times (with 0 added if not part of input parameter)
# S: median of the Nsim average cumulative surival curves at the given times
# if outputCI=True then output lower and upper confidence intervals at 100*(1-conf.int) and
# 100*conf.int  percentiles of the Nsim S values
calcLifeTableNoCovariates <- function(mod, armName, times, armAsFactor, outputCI=FALSE, conf.int=0.95, B){
  if(!class(mod)=="flexsurvreg"){
    stop("Cannot create average lifetable unless using a flexsurv model")
  } 
  
  #Do not need B if not outputting confidence intervals
  if(!outputCI) B <- 0
  
  if(armAsFactor){
    result <- summary.flexsurvreg(mod,t=times, cl=conf.int, newdata=data.frame(arm=armName), ci=outputCI, B=B)
  }
  else{
    result <- summary.flexsurvreg(mod,t=times, cl=conf.int, ci=outputCI, B=B)
  }
  
  result <- result[[1]]
  
  cols <- c("t","S")
  if(outputCI) cols <- c(cols,"lower", "upper")
  colnames(result) <- cols
  rownames(result) <- NULL
  
  
  if(!0 %in% times){
    result <-rbind(c(0,rep(1,1+2*outputCI)),result)
  }
  result
}

#Calculate the parametric life table when the model has covariates
#mod: The (flexSurv) model used to fit this arm, from which to generate the lifetable
#oneArmData: The data associated with a single arm
#times: The times
#Nsim: The number of simulations from which to generate the average
#conf.int: whether to output the confidence interval?
#level: The percentiles to calculate for the confidence interval of S
#returns a data frame with columns:
# t:times (with 0 added if not part of input parameter)
# S: median of the Nsim average cumulative surival curves at the given times
# if outputCI=True then output lower and upper confidence intervals at 100*(1-conf.int) and
# 100*conf.int  percentiles of the Nsim S values
calcParametricLifeTable <- function(mod, oneArmData, times, Nsim, outputCI=FALSE, conf.int=0.95){

  #check we are using flexsurv
  if(!class(mod)=="flexsurvreg"){
    stop("Cannot create average lifetable unless using a flexsurv model")
  }

  #Create a 3D array of model parameters for subjects
  #which takes into account their covariate values and uncertainties
  #sim.params has 3 dimensions, the first is which simulation
  #number (1 to Nsim) the second is which model parameters (e.g. for
  #weibull it is shape and scale) and the third is which subject
  #so sim.params[6, ,125] is the model parameters for subject 125 for simulation 6

  #create parameters
  sim.params <- normboot.flexsurvreg(mod, B=Nsim, newdata=oneArmData)
  if(!is.list(sim.params)){
    #This should never actually be called (unless Nsim=1) as
    #in the case of no covariates the calcLifeTableNoCovariates
    #should be used instead
    sim.params <- rep(list(sim.params),nrow(oneArmData))
  }
  
  #coerce into form described above
  sim.params <- array(unlist(sim.params),
        dim=c(dim(sim.params[[1]]),length(sim.params)),
        dimnames=dimnames(sim.params[[1]]))
  
  
  #Set up cdf function:  As an input
  #this function will take a (named) vector of arguments specifying
  #the parameters of the model and will output the values
  #of the cdf of the model with given parameters at the user inputted times

  ##get the unparameterized-model cdf
  tmp.f <- mod$dfns$p
  ##create our cdf function - splines need knots argument
  cdf <- function(params){do.call(tmp.f, c(list(times), params))}
  if(!is.null(mod$knots)){
    cdf <- function(params){do.call(tmp.f, c(list(times, knots=mod$knots), params))} 
  }
  
  
  
  #for each simulation:
  simulationResults <- apply(sim.params, 1, function(oneSimulation){

    #for each person get prob had event by given input times
    oneSimulation <- apply(oneSimulation, 2, cdf)

    #oneSimulation is a matrix with each column a subject
    #and each time a column
    1 - rowMeans(oneSimulation, na.rm=TRUE)

  })
  #simulation results is a matrix with columns for each simulation (the indivdual
  #values of S) and one row for each time requested

  #calculate median at each time point and wrap into data frame
  result <- data.frame(t=times,
                       S=apply(simulationResults,1,median))

  #if calculating confidence intervals calculate them and
  # add to result data frame
  if(outputCI){
    result <- cbind(result,
                    lower=apply(simulationResults,1,quantile,probs=1-conf.int),
                    upper=apply(simulationResults,1,quantile,probs=conf.int))
  }

  #if time = 0 not included at it t=0, S=1, confidence intervals (if used) 1
  if(!0 %in% times){
    result <-rbind(c(0,rep(1,1+2*outputCI)),result)
  }

  result
}
