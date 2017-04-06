##' @include survivalModels.R
NULL

##' Method to create to create the information criteria table 
##' @name createIcTable
##' @rdname createIcTable-methods
##' @param object SurvivalModel object
##' @param ... additional arguments for the generic 
##' @return A table contain the AIC and BIC for each model
setGeneric("createIcTable", function( object,  ... ) standardGeneric( "createIcTable" ))


##' @rdname createIcTable-methods
##' @aliases createIcTable,SurvivalModel-methods
##' @param summaryFn ("sum" or "identity") if "sum" the AIC and BIC are summed
##' over each arm producing a single value per model, it "identity" then the AIC 
##' and BIC are output for each arm separately
##' @param class ('data.frame' or (default) 'FlexTable') whether to output a data frame or FlexTable 
##' @param digits (numeric default 2) The number of decimal places to round the AIC/BIC to
##' @export
setMethod("createIcTable","SurvivalModel",
  function(object, summaryFn=c("sum","identity")[1], class=c("data.frame","FlexTable")[2], digits=2){
    if(!summaryFn %in% c("sum","identity")){
      stop("Invalid summaryFn argument should be either 'sum' or 'identity'")
    }
      
    if(length(class) != 1 || !class %in% c("data.frame","FlexTable")){
      stop("Invalid class argument, should be 'data.frame' or 'FlexTable")
    }
    if(length(digits)!=1 || !is.numeric(digits) || !digits > 0 || is.infinite(digits) ||
       is.na(digits)){
      stop("Invalid digits argument")
    }  
    
    #create ICTable data frame
    icValues <- internalGetICTable(object)
    
    if(is.na(icValues$Arm[1])) icValues$Arm  <- rep("ALL",nrow(icValues))
    
    #combine ICtable as a sum 
    if(summaryFn == "sum"){
     
      aics <- as.numeric(by(icValues,icValues$Model, function(y){sum(y$AIC)}))
      bics <- as.numeric(by(icValues,icValues$Model, function(y){sum(y$BIC)}))
      icValues <- data.frame(Model=names(by(icValues,icValues$Model, function(y){sum(y$AIC)})),
                             Arm=rep("All",length(aics)),
                             AIC=aics, BIC=bics)
      icValues <- icValues[order(icValues$AIC),]
    }  
    
    if(class=="data.frame"){ 
      return(icValues)
    }
    
    isModelSplineFit <- isSplineFit(names(object@models))
    
    if(isModelSplineFit){
      scale <- extractScale(names(object@models)[1])
      rowText <- paste0("Spline\n(scale=", scale, ")\nknots")
    }
    else{
      rowText <-  "Model"
    }
      
    #if splitting by arm
    if(!object@armAsFactor && summaryFn=="identity"){
      arms <- as.character(getArmNames(object@survData))
      icValues$Arm <- factor(icValues$Arm,levels = rev(arms)) 
      numArms <- length(arms) 
      
      numCols <- 2*numArms + 1
    }
    else{
      numCols <- 3
      arms <- "ALL"
    }
    
    numRows <- length(unique(icValues$Model))
    
    #create table
    MyFTable <- FlexTable(numrow=numRows,numcol=numCols, 
                          body.par.props=parProperties(text.align="right"),
                          header.text.props = textProperties(font.weight = "bold"),
                          body.cell.props = cellProperties(padding.right=3))
    
    #set up borders
    MyFTable[1:numRows,1:numCols,side='bottom'] <- borderProperties(width=0)
    MyFTable[1:numRows,1:numCols,side='left'] <- borderProperties(width=0)
    MyFTable[1:numRows,1:numCols,side='top'] <- borderProperties(width=0)
    MyFTable[1:numRows,1:numCols,side='right'] <- borderProperties(width=0)
    
    MyFTable[numRows,1:numCols,side='bottom'] <- borderProperties(width=3)
    MyFTable[1,1:numCols,side='top'] <- borderProperties(width=3)
    
    MyFTable[1:numRows,1] <- parProperties(text.align="left")
    MyFTable[1:numRows,1] <- textProperties(font.weight = "bold")
    
    #for each arm 
    oneArmData <- split(icValues,icValues$Arm)
    col <- 2
   
    
    #order the data 
    df <- oneArmData[[length(oneArmData)]]
    #order by knots if spline fit
    if(isModelSplineFit){
      theOrderNames <- df$Model[order(as.numeric(extractKnots(as.character(df$Model))))]
    }
    else{#sort by control arm (or all if not split by arm)
      theOrderNames <- names(sort(by(df,df$Model, function(y){sum(y$AIC)})) )
    }
    
    #add model name column to flextable
    MyFTable[1:numRows,1] <- if(isModelSplineFit) extractKnots(as.character(theOrderNames))
      else getDistributionDisplayNames(theOrderNames)
    
    theOrder <- vapply(theOrderNames, function(y){which(df$Model==y)},FUN.VALUE = numeric(1))
    
    for(df in oneArmData){
      #Add data to FlexTable
      MyFTable[1:numRows,col] <- round(df$AIC[theOrder],digits)
      MyFTable[1:numRows,col+1] <- round(df$BIC[theOrder],digits)
      
      #if spline model fits then bold the minimum of each column 
      if(isModelSplineFit){
        MyFTable[which(df$AIC[theOrder]==min(df$AIC)),col] <- textBold()
        MyFTable[which(df$BIC[theOrder]==min(df$BIC)),col+1] <- textBold()
      }
      col <- col+2
    }
    
    #Add header rows
    topBorder <- 3
    if(length(oneArmData) > 1){
      hR2 <- FlexRow(c("",rev(arms)),colspan = c(1,rep(2,length(arms))), 
                     par.properties=parProperties(text.align="left"),
                     text.properties = textProperties(font.weight = "bold"),
                     cell.properties = cellProperties(border.top.width=topBorder, border.bottom.width=0,
                                                      border.left.width=0, border.right.width=0))
      
      MyFTable <- addHeaderRow(MyFTable,hR2)
      topBorder <- 0
    }
    
    hR <- FlexRow(c(rowText,rep(c("AIC","BIC"),length(oneArmData))), 
                  par.properties=parProperties(text.align="left"),
                  text.properties = textProperties(font.weight = "bold"),
                  cell.properties = cellProperties(border.top.width=topBorder, border.bottom.width=0,
                                                   border.left.width=0, border.right.width=0))
    
    MyFTable <- addHeaderRow(MyFTable,hR)
    MyFTable
})


internalGetICTable <- function(object){
  # object@models is a list of lists with structure:
  #
  #     x[[distributionName]][[dataSet]]
  #
  # where a data set might be either the entire data set or subset by arm
  models <- object@models
  
  # Create table of results as a list of lists, then convert to data frame at the end
  icValues <- vector("list", 0)
  
  allDistNames <- names(models)
  
  for (idxDist in seq_len(length(models))){
    
    thisDistResults <- models[[idxDist]]
    thisDistName <- allDistNames[[idxDist]]
    
    dataSetNames <- names(thisDistResults)
    
    # IC values for all data sets, with this distribution
    icValuesThisDist <- lapply(seq_len(length(thisDistResults)),
                               function(i){list(thisDistName,
                                                dataSetNames[[i]],
                                                akaike_ic(thisDistResults[[i]]),
                                                bayes_ic(thisDistResults[[i]]))})
    
    # Append to list
    icValues <- append(icValues, icValuesThisDist)
  }
  
  # Convert list of lists to data frame
  icValues <- data.frame(do.call("rbind", icValues))
  
  data.frame(Model=unlist(icValues[,1]),
             Arm=unlist(icValues[,2]),
             AIC=unlist(icValues[,3]),
             BIC=unlist(icValues[,4]))  
}

# Support function in case AIC unavailable 
akaike_ic <- function(modelFit){
  
  if(class(modelFit)=="survreg"){
    return(extractAIC(modelFit)[2])
  }
  
  tryCatch(expr = {AIC(modelFit)},
           error = function(x){2*(modelFit$npar - modelFit$loglik)})
}

# Support function in case BIC unavailable (works for flexsurv models only)
bayes_ic <- function(modelFit){
  
  if(class(modelFit)=="survreg"){
    return(NA)
  }
  
  tryCatch(expr = {BIC(modelFit)},
           error = function(x){
             n <- nrow(x$data$m)
             k <- x$npar
             k*log(n) - 2*x$loglik})
}