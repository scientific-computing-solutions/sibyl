#set appropriate summary function for the covariate/covariate maturity summary tables
#return a list with the following:
#footer row string (to be displayed in last row of table)
#secondColFunction - the function used to generate the second column
#of the table (the units or category values)
#summaryFunc - a function which takes in covariate values (and censor indicators if endPoint not NULL)
#and outputs a vector of character strings one per row of the table
#leftCol2Header the column name for the second column in the table
getTypeSpecificValues <- function(isNumericTable, digits, htmlEncoding, meanOrMedian, endPoint){
  if(isNumericTable){
    pm <- if(htmlEncoding) "&plusmn;" else "\U00B1" 
    return(numericTypeSpecific(digits, pm, meanOrMedian))
  }  
  
  if(!is.null(endPoint)){ #covariate-maturity table
    return(categoricalMaturityTypeSpecific(digits))  
  }
  
  #if categorical table (logical has been set to categorical)
  categoricalTypeSpecific(digits)

}

#typeSpecific values for the covariate maturity table
categoricalMaturityTypeSpecific <- function(digits){
  summaryFunc <- function(covVals, cens){
    
    #named vector of results one per category
    vapply(levels(covVals),function(x){
      
      #if no subjects in the given {subgroup,arm}
      if(length(covVals)==0){
        return("NA")
      }
      
      #censor indicators of subjects with given covariate value
      censThisCovValue <- cens[covVals==x]
      
      #if no subjects with given {endpoint, covariate value}
      if(length(censThisCovValue)==0 || all(is.na(censThisCovValue))){
        return("-/0\n-")
      }
      
      #Number censored
      numCens <- sum(censThisCovValue,na.rm=TRUE)
      
      #Number with events
      numEvent <- sum(!censThisCovValue,na.rm=TRUE)  
      
      #maturity
      maturity <- 100*numEvent/(numEvent+numCens)
      
      paste(numEvent,"/",numEvent+numCens,"\n",round(maturity, digits) ,"%",sep="")
    
    },character(1))
  }
  
  outputFooterString <- "#Events / Total" 
  
  secondColFunction <- function(covDef){
    levels(covDef@categories)
  }
  
  leftCol2Header <-""
  
  typeSpecificConstructor(leftCol2Header, secondColFunction,
                          outputFooterString, summaryFunc)
  
}


#type specific values for the categorical covariates summary table
categoricalTypeSpecific <- function(digits){
  summaryFunc <- function(covVals, cens){
    
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
    levels(covDef@categories)
  }
  
  leftCol2Header <-""
  
  typeSpecificConstructor(leftCol2Header, secondColFunction,
                          outputFooterString, summaryFunc)
  
}

#function to create the type specific list needed for the numeric covariates
#summary table
#digits - for rounding
#pm - the unicode/html character for displaying the plus/minus symbol
#meanOrMedian - should the mean or median of the covariates be calculated?
numericTypeSpecific <- function(digits, pm, meanOrMedian){
  summaryFunc <- function(covVals, cens){
    if(length(covVals)==0 || sum(!is.na(covVals))==0){
      return("NA")
    }
    
    numberMissing <- sum(is.na(covVals))
    if(numberMissing > 0){
      missingString <- paste0("\n[", numberMissing,"]")
    }
    else{
      missingString <- ""
    }
    
    if(meanOrMedian=="mean"){
      mu <- round(mean(covVals,na.rm = TRUE),digits)
      ste <- round(se(covVals, na.rm = TRUE), digits)  
      return(paste(mu," (", pm, ste,")", missingString, sep=""))
    }
    
    media <- round(median(covVals,na.rm = TRUE), digits)
    q1 <- round(quantile(covVals,probs = 0.25, na.rm = TRUE), digits)
    q3 <- round(quantile(covVals,probs = 0.75, na.rm = TRUE), digits)
    paste(media, " [", q1, ",", q3,"]",missingString, sep="")
  }
  
  if(meanOrMedian=="mean"){
    outputFooterString <- paste("Output: mean (", pm, "se)\n [#missing - if any]",sep="")
  }
  else{
    outputFooterString <- "Output: median [Q1, Q3]\n [#missing - if any]"                      
  }
  
  secondColFunction <- function(covDef){
    if(covDef@type == "numeric") return(covDef@unit)
  }
  
  leftCol2Header <-"Unit"
  
  typeSpecificConstructor(leftCol2Header, secondColFunction,
                          outputFooterString, summaryFunc)
}


#factory to create the type specific list
typeSpecificConstructor <- function(leftCol2Header, secondColFunction,
                                    outputFooterString, summaryFunc){
  list(leftCol2Header=leftCol2Header,
       secondColFunction=secondColFunction,
       outputFooterString=outputFooterString,
       summaryFunc=summaryFunc)

}
