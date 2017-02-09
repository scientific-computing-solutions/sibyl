##' @include survivalData.R 
NULL

#create the Numer of events per subgroup dataframe
#for SurvivalData object
endPointsSummary <- function(object, digits){
  
  #First create values for table 
  maturityFunc <- function(time, cens){
    
    #If no subjects with data then return NA
    if(length(time)==0 || sum(!is.na(cens))==0) return("NA")
    
    #function to extract maturity
    mat <- function(censorIndicators){
      maturity <- 100*(length(censorIndicators)-sum(censorIndicators))/length(censorIndicators)
      as.character(round(maturity, digits=digits))
    }
    
    #if no missing data
    if(!any(is.na(cens))) return(mat(cens))
    
    #if have missing data 
    notMissing <- cens[!is.na(cens)]
    fraction <- paste0("(", sum(!notMissing), "/",
                            length(notMissing), ")")
      
    return(paste(mat(notMissing),fraction,sep="\n"))
  }
  
  KMFunc <- function(time, cens){
    if(length(time)==0 || sum(!is.na(cens))==0) return("NA")
    s <- survfit(Surv(time,!cens)~1)
    as.character(round(100*tail(s$surv, 1),digits=digits))
  }
  
  maturityVals <- extractEndPointOutput(object, func= maturityFunc) 
  KMs <-  extractEndPointOutput(object, func=KMFunc)  
  
  ######################################
  
  #calculate size of table
  numSubgroups <- length(object@subgroupDef)
  numArms <- length(object@armDef@categories)
  
  numRows <- 2*length(object@endPoints) #maturity + KM for each endpoint
  numCols <- 2 + (1+numSubgroups)*numArms # for each (subgroup + all data) x for each arm
  
  #Then create table
  MyFTable <- FlexTable(numrow=numRows,numcol=numCols, 
                        body.par.props=parProperties(text.align="right"),
                        header.text.props = textProperties(font.weight = "bold"),
                        body.cell.props = cellProperties(padding.right=1))
  
  #Add values into table
  MyFTable[seq(1,numRows,2),3:numCols] <- maturityVals
  MyFTable[seq(2,numRows,2),3:numCols] <- KMs
  
  #Add first two columns (endpoints and outcome columns)
  MyFTable[1:numRows,2] <- rep( c("Maturity", "Kaplan Meier\n% end trial"), numRows/2    )
  MyFTable[1:numRows,2] <- parProperties(text.align="left")
  
  MyFTable[seq(1,numRows,2),1] <- names(object@endPoints) 
  MyFTable[1:numRows,1] <- parProperties(text.align="left")
  MyFTable[1:numRows,1] <- textProperties(font.weight = "bold")
  
  #merge cells in first column
  for(i in seq(1,numRows,2)){
    MyFTable <- spanFlexTableRows(MyFTable, j = 1, from=i, to=i+1)
  }
  
  #Add headers
  subgroupDetails <- extractSubgroupTable(object)
  headers <- getHeaders(subgroupDetails, leftCol1=c("Events \n(%)",""),
                        leftCol2=c("Endpoint", "Outcome"))
  
  MyFTable <- addHeaderRow(MyFTable,headers[[1]])
  MyFTable <- addHeaderRow(MyFTable,headers[[2]])
  
  MyFTable
}



#Given a survivalData object output a dataframe, rows for endpoints,
#columns for treatment arms, first all data then subgroup 1 then subgroup 2
#with values calculated using the summary function func
#which takes in the survival times and censor for the given endpoint
#for subjects in the appropriate arm and subgroup
extractEndPointOutput <- function(object, func){
  
  #For each endpoint
  retVal <- lapply(object@endPoints, function(eP){
    
    time <- object@subject.data[,eP$timeCol]
    cens <- object@subject.data[,eP$censorCol]
    theArms <- object@subject.data$arm
    
    #for each subgroup and "ALL" for everyone
    unlist(lapply(c("ALL",listColumnDefSlot(object@subgroupDef,"columnName")),function(subgroup){
      
      if(subgroup != "ALL"){
        time <- time[object@subject.data[,subgroup]]
        cens <- cens[object@subject.data[,subgroup]]
        theArms <- theArms[object@subject.data[,subgroup]]
      }
      
      #for each arm
      vapply(rev(as.character(getArmNames(object))),function(arm){
        time <- time[theArms==arm]
        cens <- cens[theArms==arm]
        func(time, cens)
      }, FUN.VALUE = character(1))
      
    }))
    
  })
  do.call("rbind", retVal) 
}