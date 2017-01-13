#' @include survivalData.R
NULL

# Create the patient characteristics table
# This outputs two tables: one for categorical and one for numeric data
covariatesSummary <- function(object, htmlEncoding, meanOrMedian){
  
  #Need to replace <NA> in categorical and logical with a "missing data" level
  #for output in table and converts logical covariates to factors for output
  object <- convertMissingFactorsToOwnLevel(object)
  
  numericCovariatesTable <- createCovariateSummarySubTable(object,
                                                           "numeric", digits=1,
                                                           htmlEncoding, meanOrMedian)
  
  categoricalCovariatesTable <- createCovariateSummarySubTable(object,
                                                               "categorical",
                                                               digits=2,
                                                               htmlEncoding, meanOrMedian)
  
  return(list(numeric=numericCovariatesTable,
              categorical=categoricalCovariatesTable))
}


#Create the covariatesMaturity table
covariatesMaturitySummary <- function(object){
  
  #Need to replace <NA> in categorical and logical with a "missing data" level
  #for output in table and converts logical covariates to factors for output
  object <- convertMissingFactorsToOwnLevel(object)
  
  #if no categorical covariates then return NULL
  if(length(object@covDef)==0 || !any(listColumnDefSlot(object@covDef,"type") =="categorical")){
    return(NULL)
  }
  
  
  retVal <- lapply(names(object@endPoints),function(endPointName){
    createCovariateSummarySubTable(object,
                                   "categorical",
                                   digits=1, 
                                   htmlEncoding=NULL,
                                   meanOrMedian=NULL,   
                                   endPoint=object@endPoints[[endPointName]]) 
  })
  names(retVal) <- names(object@endPoints)
  retVal
}


#create either the categorical or the numeric summary table
#if endPoint is not NULL then we are calculating a covariate-Maturity summary table for the
#given endpoint, otherwise we are calculating the standard covariates table
createCovariateSummarySubTable <- function(object, requiredTypes, digits, htmlEncoding, meanOrMedian, endPoint=NULL){
  
  numRowsFromCov <- vapply(object@covDef, numberOfRowsNeeded, requiredTypes, FUN.VALUE = numeric(1))
  
  #if no covariates for table return NULL
  if(length(numRowsFromCov)==0|| sum(numRowsFromCov)==0) return(NULL)
  
  #get objects which depend on whether this is the numeric or categorical
  #summary table
  typeSpecificValues <- getTypeSpecificValues("numeric" %in% requiredTypes, digits, htmlEncoding,
                                              meanOrMedian, endPoint)
  
  
  #create FlexTable
  numSubgroups <- length(object@subgroupDef)
  numArms <- length(object@armDef@categories)
  
  numRows <- sum(numRowsFromCov)
  numCols <- 2 + (1+numSubgroups)*numArms # for each (subgroup + all data) x for each arm
  
  MyFTable <- FlexTable(numrow=numRows,numcol=numCols, 
                        body.par.props=parProperties(text.align="right"),
                        header.text.props = textProperties(font.weight = "bold"),
                        body.cell.props = cellProperties(padding.right=1))
  
  #calculate table values
  ans <- extractCovariateOutput(object,typeSpecificValues$summaryFunc,requiredTypes, endPoint)
  MyFTable[1:numRows,3:numCols] <- ans
  
  #firstColumn
  firstCol <- unlist(lapply(seq_along(numRowsFromCov),function(i){
    rep(object@covDef[[i]]@displayName,numRowsFromCov[i])
  }))
  MyFTable[1:numRows,1] <- firstCol
  #merge cells in first column with same value
  MyFTable <- spanFlexTableRows(MyFTable, j = 1, runs = as.character(firstCol))
  
  #second Column
  MyFTable[1:numRows,2] <- unlist(lapply(object@covDef,typeSpecificValues$secondColFunction))
  
  MyFTable[1:numRows,1:2] <-  parProperties(text.align="left")
  MyFTable[1:numRows,1] <- textProperties(font.weight = "bold")
  
  #align text to top of cells
  MyFTable[1:numRows,3:numCols] <- cellProperties(vertical.align="top")
  
  #sort out headers  
  subgroupDetails <- extractSubgroupTable(object)
  headers <- getHeaders(subgroupDetails, leftCol1="",
                        leftCol2=c("Covariate", typeSpecificValues$leftCol2Header))
  
  MyFTable <- addHeaderRow(MyFTable,headers[[1]])
  MyFTable <- addHeaderRow(MyFTable,headers[[2]])
  
  
  #Add footer
  if(typeSpecificValues$outputFooterString != ""){
    fR3 <- FlexRow(c("",typeSpecificValues$outputFooterString),
                 colspan=c(2,numCols-2),
                 par.properties=parProperties(text.align="center",padding=1))
  
    MyFTable <- addFooterRow(MyFTable,fR3)
  }
  
  MyFTable
}

#How many rows of the tables are needed for covariate defined by 
#ColumnDef covDef - and if its type is not in requiredTypes
#vector then no rows needed
numberOfRowsNeeded <- function(covDef, requiredTypes){
  type <- covDef@type
  if(!type %in% requiredTypes) return(0)
  if(type=="numeric") return(1)
  length(covDef@categories) #categorical type (logical has been converted to categorical)
}  


#Given a survivalData object output a dataframe, rows for covariates
#if numeric or rows for single category of categorical covariate 
#columns for treatment arms, first all data then subgroup 1 then subgroup 2
#with values calculated using the summary function func
#which takes in the values of the covariate (and the censoring indicator for the
#given endPoint if endPoint is not NULL) 
#for subjects in the appropriate arm and subgroup
extractCovariateOutput <- function(object, func, requiredTypes, endPoint){
  
  #For each covariate
  retVal <- lapply(object@covDef, function(cov){
    
    if(!cov@type %in% requiredTypes) return(NULL)
    
    covVals <- object@subject.data[,cov@columnName]
    theArms <- object@subject.data$arm
    cens <- if(!is.null(endPoint)) object@subject.data[,endPoint$censorCol] else NULL
  
    
    #for each subgroup and "ALL" for everyone
    ans <- lapply(c("ALL",listColumnDefSlot(object@subgroupDef,"columnName")),function(subgroup){
      
      if(subgroup != "ALL"){
        covVals <- covVals[object@subject.data[,subgroup]]
        theArms <- theArms[object@subject.data[,subgroup]]
        if(!is.null(cens)){
          cens <- cens[object@subject.data[,subgroup]]
        }
      }
      
      #for each arm
      oneArmRes <- vapply(as.character(getArmNames(object)),function(arm){
        covVals <- covVals[theArms==arm]
        cens <- if(!is.null(cens)) cens[theArms==arm] else NULL
        func(covVals, cens)  
      }, FUN.VALUE = character(numberOfRowsNeeded(cov,requiredTypes)))
      
      if(class(oneArmRes)=="character"){
        oneArmRes <- matrix(oneArmRes, nrow=1)
      }
      
      oneArmRes  
        
    })
    
    ans <- do.call("cbind",ans)
    
    #reorder data frame as vapply messes it up
    if(cov@type == "categorical"){
      ans <- ans[order(as.numeric(cov@categories)),]
    }
    ans
  })
  
  do.call("rbind", retVal) 
}

#function which takes a SurvivalData object and replaces missing factor variables
#with their own level "(no data)" for output in covariate table and converts
#logical variables to factors TRUE, FALSE
convertMissingFactorsToOwnLevel <- function(object){
  
  for(idx in seq_along(object@covDef)){
    
    cov <- object@covDef[[idx]]
    name <- cov@columnName
    
    if(cov@type == "logical"){
      object@covDef[[idx]]@type <- "categorical"
      object@covDef[[idx]]@categories <- factor(c("TRUE","FALSE"),levels = c("TRUE","FALSE"))
      object@subject.data[,name] <- factor(object@subject.data[,name], levels=c(TRUE,FALSE)) 
      
    }
    
    
    if(cov@type %in% c("logical", "categorical")){
      #if missing data
      if(any(is.na(object@subject.data[,name]))){
        
        newCategories <- as.character(object@covDef[[idx]]@categories)
        
        #change column definition
        newCategories <- factor(c(newCategories, "(no data)"),
                                levels= c(newCategories, "(no data)"))
        
        object@covDef[[idx]]@categories <- newCategories
      
        #and data
        object@subject.data[,name] <-as.character(object@subject.data[,name])
        object@subject.data[,name] <- ifelse(is.na(object@subject.data[,name]), "(no data)",
                                             object@subject.data[,name] )
        
        object@subject.data[,name] <- factor(object@subject.data[,name],levels= newCategories)
      }
      
    }
  }
  
  object
}
