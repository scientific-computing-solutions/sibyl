##' @include survivalModels.R
NULL

##' Class storing the dataframe used to plot the
##' parametric averaged survival curves
##' @slot plotdata (data frame) to be used to plot
##' the survival data object it has the following columns
##' \itemize{
##'   \item{t} survival time 
##'   \item{S} probability of survival 
##'   \item{Arm} which arm the data in this row is associated with
##'   \item{lower} lower (95%) confidence interval for S
##'   \item{upper} upper (95%) confidence interval for S
##'   \item{model} which model was used to generate the data (or "KM" for the Kaplan-Meier estimator)
##' } 
##' @export
setClass("AvCurvePlotData",
         slots = list(plotdata="data.frame"))


##' Method to create the \code{AvCurvePlotData} object
##' @name createAvCurvePlotData
##' @rdname createAvCurvePlotData-methods
##' @param object \code{SurvivalModel} object
##' @param ... Additional arguments for specific instances of this generic
##' @return A \code{AvCurvePlotData} object
setGeneric( "createAvCurvePlotData",
            function(object,...) standardGeneric("createAvCurvePlotData"))

##' @rdname createAvCurvePlotData-methods
##' @aliases createAvCurvePlotData,SurvivalModel-methods
##' @param maxTime (numeric) the time to extrapolate the curves to (default=NULL implies no extrapolation)
##' @param Npoints (numeric) the number of time points for which the survival curves are to be evaluated at
##' @param Nsim (numeric) the number of simulations used to generate the averaged curves
##' @param models (character vector) which models from \code{names(object@@models)} are to be used when 
##' calculating averaged survival curves - default NULL implies use all
##' @param seed (numeric, default NULL) if not NULL then set random seed (although it will only be
##' used if models include covariates) - random numbers are used when models includes covariates 
##' or when confidence intervals are required  
##' @param B (integer) Only used when no covariates in model. See summary.flexsurvreg 
##' Number of simulations from the normal asymptotic distribution of the estimates used 
##' to calculate confidence intervals. Decrease for greater speed at the expense of accuracy, 
##' or set B=0 to turn off calculation of CIs.
##' @param conf.type ("none", "plain", "log" [default], "log-log") argument passed to survfit
##' @details If the models include covariates a simulation procedure is required to generate averaged survival
##' curves. If the models do not include covariates then \code{summary.flexsurvreg} is used
##' @export
setMethod("createAvCurvePlotData", signature(object="SurvivalModel"),
  function(object, maxTime=NULL, Npoints=201, Nsim=500, models=NULL, seed=NULL, B=1000,
           conf.type=c("none", "plain", "log", "log-log")[3]){

    if(length(conf.type)!=1 || !conf.type %in% c("none", "plain", "log", "log-log")){
      stop("Invalid conf.type argument")
    }
    
    validateCreateAvCurvePLotDataArgs(maxTime, Npoints, Nsim, seed)

    if(!is.null(seed)){
      set.seed(seed)
    }
    
    #If models is NULL then use all models
    if(is.null(models)){
      models <- names(object@models)
    }
    else{ #Otherwise check 
      models <- tolower(models)
      if(!all(models %in% names(object@models))){
        stop(paste("Invalid models argument - the vector can only contain elements from:",
                   paste(names(object@models),collapse=", ")))
      }
    }
    

    endPointDef <- object@endPointDef
    subjectData <- object@survData@subject.data

    # Calculate KM-lifetable with times given by survfit
    KMLifeTables <- calcKMLifeTable(subjectData,endPointDef,outputCI=TRUE, conf.type)

    KMLifeTables <- do.call("rbind",KMLifeTables)
    KMLifeTables$model <- "KM"

    #if maxTime is not set then take last timepoint in km lifetable
    if(is.null(maxTime)){
      maxTime <- max(KMLifeTables$t)
    }

    # Ensure parametric plots extend at least as far as the KM plot
    maxTime <- max(maxTime, max(KMLifeTables$t))

    #calculate times for lifetables
    stepSize <- maxTime/(Npoints-1)
    times <- seq(0,maxTime,stepSize)

    #Next create lifetables for each parametric model at the given times

    #split data by arm
    dataByArm <- split(subjectData, subjectData[, "arm"])

    #generate a list, one element per model each of which
    #contains a data frame containing the lifetables for
    #each arm, concatenated
    parametricLifeTables <- lapply(models, function(modelName){

      if(length(object@covariates)!=0){
      
      #calculate the lifetable for a single model, for each arm
      lifeTableOneModel <-
        mapply(calcParametricLifeTable, mod=object@models[[modelName]],
               oneArmData=dataByArm,
               MoreArgs=list(times=times,Nsim=Nsim, outputCI=TRUE), SIMPLIFY = FALSE)
      }
      else{
        lifeTableOneModel <- mapply(calcLifeTableNoCovariates, mod=object@models[[modelName]], 
                                    armName=getArmNames(object@survData),
                                    MoreArgs=list(times=times,outputCI=TRUE,B=B,
                                                  armAsFactor=object@armAsFactor),
                                    SIMPLIFY = FALSE)  
      }
    
      lifeTableOneModel <- do.call("rbind",lifeTableOneModel)
      lifeTableOneModel$Arm <- rep(names(dataByArm), each=length(times))
      lifeTableOneModel$model <- modelName

      lifeTableOneModel
    })

    parametricLifeTables <- do.call("rbind",parametricLifeTables)

    result <- rbind(KMLifeTables, parametricLifeTables)
    rownames(result) <- NULL

    new("AvCurvePlotData",
        plotdata=result)
})


validateCreateAvCurvePLotDataArgs <- function(maxTime, Npoints, Nsim, seed){
  
  if(!is.null(seed)){
    if(length(seed)>1 || !is.numeric(seed)){
      stop("Invalid seed")
    }
  }
  
  isWholeNumber <- function(x, tol = .Machine$double.eps^0.5){ 
    abs(x - round(x)) < tol
  }

  if(!is.null(maxTime)){
    if(length(maxTime)!= 1 || !is.numeric(maxTime) || !maxTime > 0 || is.infinite(maxTime)){
      stop("Invalid maxtime argument, must be positive")
    }
  }
  
  if(length(Npoints) != 1 || !is.numeric(Npoints) || !Npoints > 1 || !isWholeNumber(Npoints) || is.infinite(Npoints)){
    stop("Invalid Npoints argument, must be positive finite integer greater than 1")
  } 
  
  if(length(Nsim) != 1 || !is.numeric(Nsim) || !Nsim > 0 || !isWholeNumber(Nsim) || is.infinite(Nsim)){
    stop("Invalid Nsim argument, must be positive finite integer")
  } 

}  


##' @rdname plot-methods
##' @name plot
##' @aliases plot,AvCurvePlotData,missing-method
##' @param use.facet (logical) should the treatment arms be split into different facets
##' @param outputModel (character) which model's survival curves should be plotted (default NULL
##' implies use all) if =character(0) then output only the KM curve
##' @param xMax (numeric or default NULL) the x-axis limit for the graph. If not included then
##' all data is displayed
##' @param useCI (logical) TRUE if including CI on graph, FALSE otherwise 
##' @export
setMethod("plot", signature(x="AvCurvePlotData", y="missing"),
  function(x, use.facet=TRUE, outputModel=NULL, xMax=NULL, useCI=FALSE){

    #R-cmd-check thinks t, s, model, ... are global 
    #variables inside the ggplot commands so complains about them
    #they are not global variables but for them to pass R-cmd-check we 
    #need to create dummy variables
    t <- NULL
    S <- NULL
    model <- NULL
    Arm <- NULL
    upper <- NULL
    lower <- NULL
    
    
    #setup legend for distributions
    modelNames <- unique(x@plotdata$model)
    modelNames <- modelNames[modelNames!="KM"]
    
    isModelSplineFit <- isSplineFit(modelNames)
    if(isModelSplineFit){
      legendLabel <- extractKnots(modelNames)
      legendLabel <- c("KM", legendLabel)
      names(legendLabel) <- c("KM", modelNames)
      cols <- getSplineDistColours(c("KM",modelNames))
      legendTitle <- paste0("Spline\n(scale=", extractScale(modelNames[1]), ")\nknots")
    }
    else{
      cols <- getDistColours(unique(x@plotdata$model))
      legendLabel <- getDistributionDisplayNames(unique(x@plotdata$model))
      legendTitle <- "Distribution"
    }  
   
    # Pull out data for KM curve - this is plotted differently from the others
    kmData <- x@plotdata[x@plotdata$model == "KM", ]
    
    #if model given then only plot KM and the given models:
    if(!is.null(outputModel)){
     
      outputModel <- tolower(outputModel)
      
      #check model exists
      if(any(!outputModel %in% x@plotdata$model)){
        warning(paste0("There is no data for some models so they cannot be plotted!"))
      }
      
      modelData <- x@plotdata[x@plotdata$model %in% outputModel, ]
      
    }
    else{
      modelData <- x@plotdata[x@plotdata$model != "KM", ]
    }
    

    # Create plot, with KM line shown step-wise and set colours
    #and lgend names correctly
    p <- ggplot(modelData, aes(x = t, y = S, colour = model)) +
      scale_colour_manual(name=legendTitle, values=cols, labels=legendLabel)

    # Add labels
    p <- p + xlab("Time")
    p <- p + ylab("P(survival)")

    # Show arms on separate plots
    if (use.facet){

      # Separate plots for each arm
      p <- p + facet_grid(Arm ~ .)

      # Use the same line type on each plot
      p <- p +
           geom_line() +
           geom_step(data = kmData, aes(x = t, y = S))
      # Add lines for confidence intervals
      if(useCI){
        p <- p +
          geom_line(aes(x = t, y = lower, colour = model), linetype = "dashed") +
          geom_step(data = kmData, aes(x = t, y = lower), linetype = "dashed") +
          geom_line(aes(x = t, y = upper, colour = model), linetype = "dashed") +
          geom_step(data = kmData, aes(x = t, y = upper), linetype = "dashed")
      }
      
      
    }
    else{
      # Use different line types for each arm, on a single plot
      p <- p +
           geom_line(aes(linetype = Arm)) +
           geom_step(data = kmData, aes(linetype = Arm))

      if(useCI){
        p <- p +
          geom_line(aes(x = t, y = lower, colour = model, linetype = Arm)) +
          geom_step(data = kmData, aes(x = t, y = lower, linetype = Arm)) +
          geom_line(aes(x = t, y = upper, colour = model,linetype = Arm)) +
          geom_step(data = kmData, aes(x = t, y = upper,linetype = Arm))
      }
      
    }
    
    #if xmax is set then set xlim
    if(!is.null(xMax)){
      p <- p + coord_cartesian(xlim=c(0, xMax), ylim=c(0,1.05),expand=FALSE)
    }
    
    # Format background and borders
    p <- p + theme(panel.background = element_blank(),
                   panel.border = element_rect(colour = "black", fill = NA),
                   panel.grid.major = element_line(color="black", linetype="dotted"),
                   panel.grid.minor = element_blank())
    
    # Display
    p
  }
)


#modelNames should be a vector of unique model
#names in the avCurvePlotData object
getDistColours <- function(modelNames){
  
  retVal <- vapply(modelNames, function(x){
    if(nchar(x) > 6 && substr(x,1,6)=="spline"){
      x <- "spline"
    }
    
    switch(x,
           KM="black",
           exp="red",
           exponential="red",
           weibull="darkblue",
           lnorm="darkgreen",
           lognormal="darkgreen",
           llogis="orange",
           loglogistic="orange",
           gompertz="purple",
           gengamma="cyan",
           gengamma.orig="cyan",
           spline="brown",
           logistic="blue",
           gaussian="yellow",
           genf="pink",
           genf.orig="pink",
           "red" #default if other
    )
  },FUN.VALUE = character(1))
  names(retVal) <- modelNames
  retVal
}

#different colours are used when all models are spline
#modelNames = c("KM", "0", "3", "4") if splines with 0, 3, 4
#knots have been fitted
getSplineDistColours <- function(modelNames){
  retVal <- vapply(modelNames, function(x){
    if(nchar(x) > 6 && substr(x,1,6)=="spline"){
      x <- strsplit(x,"_")[[1]][2]
    }
    switch(x,
         KM="black",
         "0"="red",
         "1"="darkblue",
         "2"="purple",
         "3"="orange",
         "4"="cyan",
         "5"="brown",
         "6"="blue",
         "7"="yellow",
         "8"="pink",
         "red")
  },FUN.VALUE = character(1))
  names(retVal) <- modelNames
  retVal
}