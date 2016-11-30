#' @include survivalData.R
NULL

# Create the patient characteristics table
# This outputs two tables: one for categorical and one for numeric data
covariatesSummary <- function(object, htmlEncoding){
  
  numericCovariatesTable <- createCovariateSummarySubTable(object,
                                                           "numeric", digits=1,
                                                           htmlEncoding)
  
  categoricalCovariatesTable <- createCovariateSummarySubTable(object,
                                                               c("categorical", "logical"),
                                                               digits=2,
                                                               htmlEncoding)
  
  return(list(numeric=numericCovariatesTable,
              categorical=categoricalCovariatesTable))
}


#create either the categorical or the numeric summary table
createCovariateSummarySubTable <- function(object, requiredTypes, digits, htmlEncoding){
  
  numRowsFromCov <- vapply(object@covDef, numberOfRowsNeeded, requiredTypes, FUN.VALUE = numeric(1))
  
  #if no covariates for table return NULL
  if(length(numRowsFromCov)==0|| sum(numRowsFromCov)==0) return(NULL)
  
  #get objects which depend on whether this is the numeric or categorical
  #summary table
  typeSpecificValues <- getTypeSpecificValues("numeric" %in% requiredTypes, digits, requiredTypes, htmlEncoding)
  
  
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
  ans <- extractCovariateOutput(object,typeSpecificValues$summaryFunc,requiredTypes)
  MyFTable[1:numRows,3:numCols] <- ans
  
  #firstColumn
  firstCol <- unlist(mapply(rep,listColumnDefSlot(object@covDef,"displayName"), each=numRowsFromCov))
  MyFTable[1:numRows,1] <- firstCol
  #merge cells in first column with same value
  MyFTable <- spanFlexTableRows(MyFTable, j = 1, runs = as.character(firstCol))
  
  #second Column
  MyFTable[1:numRows,2] <- unlist(lapply(object@covDef,typeSpecificValues$secondColFunction ))
  
  MyFTable[1:numRows,1:2] <-  parProperties(text.align="left")
  MyFTable[1:numRows,1] <- textProperties(font.weight = "bold")
  
  
  #sort out headers  
  subgroupDetails <- extractSubgroupTable(object)
  headers <- getHeaders(subgroupDetails, leftCol1="",
                        leftCol2=c("Covariate", typeSpecificValues$leftCol2Header))
  
  MyFTable <- addHeaderRow(MyFTable,headers[[1]])
  MyFTable <- addHeaderRow(MyFTable,headers[[2]])
  
  
  #Add footer
  fR3 <- FlexRow(c("",typeSpecificValues$outputFooterString),
                 colspan=c(2,numCols-2),
                 par.properties=parProperties(text.align="center",padding=1))
  
  MyFTable <- addFooterRow(MyFTable,fR3)
  
  MyFTable
}

#How many rows of the tables are needed for covariate defined by 
#ColumnDef covDef - and if its type is not in requiredTypes
#vector then no rows needed
numberOfRowsNeeded <- function(covDef, requiredTypes){
  type <- covDef@type
  if(!type %in% requiredTypes) return(0)
  if(type=="numeric") return(1)
  if(type=="logical") return(2)
  length(covDef@categories) #categorical type
}  


#set appropriate summary function,
#footer row string (to be displayed in last row of table)
#secondColFunction - the function used to generate the second column
#of the table (the units or category values)
getTypeSpecificValues <- function(isNumericTable, digits, requiredTypes, htmlEncoding){
  if(isNumericTable){
    
    pm <- if(htmlEncoding) "&plusmn;" else "\U00B1" 
    
    summaryFunc <- function(covVals){
      if(length(covVals)==0){
        return("NA")
      }
      mu <- round(mean(covVals),digits)
      ste <- round(se(covVals), digits)
      mi <- round(min(covVals), digits)
      ma <- round(max(covVals), digits)
      media <- round(median(covVals), digits)
      paste(mu," (", pm, ste,")\n [", mi, ",", ma,"] ", media, sep="")
    }
    
    outputFooterString <- paste("Output: mean (", pm, "se)\n [min, max] median",sep="")
   
    secondColFunction <- function(covDef){
      if(covDef@type %in% requiredTypes ) return(covDef@unit)
    }
    
    leftCol2Header <-"Unit"
  }
  else{ #if categorical table
    summaryFunc <- function(covVals){
      #coerce logical to factor for the analysis
      if(is.logical(covVals)){
        covVals <- factor(covVals, levels=c(TRUE,FALSE)) 
      }
      
      #named vector of results one per category
      vapply(levels(covVals),function(x){
        if(length(covVals)==0){
          return("NA (0)")
        }
        perCent <- round(100*sum(covVals==x)/length(covVals),digits)
        paste(perCent," (",sum(covVals==x),")",sep="")
      },character(1))
    }
    
    outputFooterString <- "Output: % (n)" 
    
    secondColFunction <- function(covDef){
      if(covDef@type=="logical") return(c("TRUE","FALSE"))
      if(covDef@type=="categorical"){
        return(levels(covDef@categories))
      } 
    }
    
    leftCol2Header <-""
  }
  
  return(
    list(leftCol2Header=leftCol2Header,
         secondColFunction=secondColFunction,
         outputFooterString=outputFooterString,
         summaryFunc=summaryFunc)
  )
  
}



#Given a survivalData object output a dataframe, rows for covariates
#if numeric or rows for single category of categorical covariate 
#columns for treatment arms, first all data then subgroup 1 then subgroup 2
#with values calculated using the summary function func
#which takes in the values of the covariate 
#for subjects in the appropriate arm and subgroup
extractCovariateOutput <- function(object, func, requiredTypes){
  
  #For each covariate
  retVal <- lapply(object@covDef, function(cov){
    
    if(!cov@type %in% requiredTypes) return(NULL)
    
    covVals <- object@subject.data[,cov@columnName]
    theArms <- object@subject.data$arm
    
    #for each subgroup and "ALL" for everyone
    ans <- lapply(c("ALL",listColumnDefSlot(object@subgroupDef,"columnName")),function(subgroup){
      
      if(subgroup != "ALL"){
        covVals <- covVals[object@subject.data[,subgroup]]
        theArms <- theArms[object@subject.data[,subgroup]]
      }
      
      #for each arm
      vapply(as.character(getArmNames(object)),function(arm){
        covVals <- covVals[theArms==arm]
        func(covVals)  
      }, FUN.VALUE = character(numberOfRowsNeeded(cov,requiredTypes)))
      
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
