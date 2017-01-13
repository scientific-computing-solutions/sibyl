#This file contains a function used to check SurvivalData objects against a set of rules
#for specific named endpoints (i.e. OS, PFS) - if there are additional rules
#they should be added into this function and the roxygen comments updated, 
#some tests added to the testthat file test-specificEndpointValidation.R and
#numberRules variable updated below

##' Produce warning messages if SurvivalData object endpoints
##' do not satisfy a set of rules
##' @details Rules:
##' 1) If endpoint PFS2 and PFS exists then for all subjects time to PFS must be < time
##' to PFS2 
##' 2) If endpoint OS exists then for all subjects who have an event for OS, their event
##' time for OS must be >= their time for all other endpoints
##'
##' @param object (SurvivalData object)
##' @param rulesToCheck (numeric vector) - which rules should be checked
##' NULL if check all rules 
##' @return A string of warning messages or NULL if no warnings
##' @export
specificEndpointRules <- function(object, rulesToCheck=NULL){
  
  numberRules <- 2 #if number of rules changes update this number
  if(is.null(rulesToCheck)) rulesToCheck <- 1:numberRules 
  
  if(class(object)!="SurvivalData"){
    stop("object must be a SurvivalData object")  
  }
  
  if(any(!rulesToCheck %in% 1:numberRules)){
    stop("Invalid rulesToCheck argument")
  }
  
  
  data <- object@subject.data
  
  #Any rules which are applicable to this Survival object are
  #added into this list
  rules <- list()
  
  #helper function to get times and censor columns
  getTime <- function(endPoint, oneSubjectData){
    oneSubjectData[,object@endPoints[[endPoint]]$timeCol]  
  }
  getCensor <- function(endPoint, oneSubjectData){
    oneSubjectData[,object@endPoints[[endPoint]]$censorCol]  
  }
  
  
  #Rule 1, no PFS2 time can be less than PFS time
  #function takes the single row of the SurvivalData object's data frame
  #and returns TRUE if a warning should be output and FALSE otherwise
  pfs2Func <- function(oneSubjectData){
    #get times
    pfs <- getTime("PFS", oneSubjectData) 
    pfs2 <- getTime("PFS2", oneSubjectData)
    #if both not NA and pfs2 < pfs then return TRUE as need a warning output
    #otherwise return FALSE as no warning needed
    !is.na(pfs) && !is.na(pfs2) && pfs2 < pfs  
  }
  
  #this rule is applicable if PFS and PFS2 are endpoints in the data file
  if(1 %in% rulesToCheck){
    if("PFS" %in% names(object@endPoints) && "PFS2" %in% names(object@endPoints)){  
      rules[[length(rules)+1]] <- list(func=pfs2Func,warningMsg="has PFS2 time < PFS time")
    }
  }
  
  #Rule 2, no times can occur after the event for endpoint "OS" 
  #function takes the single row of the SurvivalData object's data frame
  #and returns TRUE if a warning should be output and FALSE otherwise
  osFunc <- function(oneSubjectData){
    #extract OS time/indicator
    osCens <- getCensor("OS", oneSubjectData)
    os <- getTime("OS", oneSubjectData) 
   
    #if no OS data or OS was censored do not need a warning
    if(is.na(osCens) || osCens) return(FALSE)
    
    #for each endpoint if its time > OS then return TRUE as a warning is needed 
    for(endPoint in names(object@endPoints)){
      if(endPoint != "OS"){
        time <- getTime(endPoint, oneSubjectData)
        if(!is.na(time) && time > os) return(TRUE)
      }
      
    }
    return(FALSE)
  } 
  
  #rule is applicable if OS and at least one other endpoint are included in the data file
  if(2 %in% rulesToCheck){
    if("OS" %in% names(object@endPoints) && length(object@endPoints)>1){
      rules[[length(rules)+1]] <- list(func=osFunc,
                                     warningMsg="has, for an endpoint other than OS, a time > OS event time")
    }  
  }
 
  if(length(rules)==0) return(NULL)
  
  evaluateRules(data, rules)
}


evaluateRules <- function(data, rules){
  
  #for each subject (use lapply not apply to preserve data types)
  retVal <- unlist(lapply(seq_len(nrow(data)), function(i){
    
    oneSubjectData <- data[i,]  
    
    #apply each rule to the subject
    oneSubjectWarning <- lapply(rules, function(rule){
      if(rule$func(oneSubjectData)){
        paste("WARNING: Subject", oneSubjectData$subject, rule$warningMsg)  
      }        
    })
    
    paste(unlist(oneSubjectWarning), collapse="\n")
    
  }))
  
  if(all(retVal=="")) return(NULL)
  return(paste(retVal[retVal!=""], collapse="\n"))
}
