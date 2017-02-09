##' Summary methods for Sibyl objects
##' @name summary
##' @rdname summary-methods
##' @aliases summary,SurvivalData-method
##' @param object (SurvivalData or SemiParametricModel object) data to be summarised
##' @param type (one of "subgroups", "endpoints", "covariates", "covarMaturity") which output table
##' should be generated for the \code{SurvivalData} object OR (one of "medianTTE" or "KM") if 
##' object is of type SemiParametricModel
##' @param digits (numeric, default 1) number of decimal places for the output
##' @param htmlEncoding (logical) if TRUE use &plusmn instead of unicode plus minus sign 
##' @param meanOrMedian ("median" or "mean") output summary function for the numeric covariates
##' summary table
##' @export
setMethod("summary", signature(object="SurvivalData"),
  function(object, type=c("subgroups","endPoints","covariates", "covarMaturity")[1], 
           digits=1, htmlEncoding=FALSE, meanOrMedian=c("median","mean")[1]){

  if(length(meanOrMedian) != 1 || !meanOrMedian %in% c("median","mean")){
    stop("meanOrMedian argument must be 'mean' or 'median")
  }  
    
  if(length(digits)!=1 || !is.numeric(digits) || !digits > 0 || is.infinite(digits) ||
     is.na(digits)){
    stop("Invalid digits argument")
  }  
    
  switch(tolower(type),
                "subgroups" = subgroupsSummary(object) ,
                "endpoints" = endPointsSummary(object, digits=digits),
                "covariates" = covariatesSummary(object, htmlEncoding, meanOrMedian),
                "covarmaturity" = covariatesMaturitySummary(object, digits=digits),
                 stop("Invalid type argument must be one of 'subgroups', 'endpoints'
                      or 'covariates' or 'covarMaturity'"))
})


#create the raw data for the subgroups summary object,
#input SurvivalData object, output a data frame with columns
#for each arm and subgroups for each  row cells contain #subjects in given
#arm/subgroup pair 
extractSubgroupTable <- function(object){
  
  dataSplitByArm <- split(object@subject.data, object@subject.data[, "arm"])
  subgroupCols <- listColumnDefSlot(object@subgroupDef,"columnName")
  
  #for each arm
  retVal <- lapply(dataSplitByArm, function(df){
  
    if(length(subgroupCols)==0){
      return(cbind(n=nrow(df)))
    }
    
    #for each subgroup
    ans <- lapply(subgroupCols,function(x){
      sum(df[,x])
    })
    names(ans) <- subgroupCols
    cbind(n=nrow(df),data.frame(ans))
  })

  results <- do.call("rbind",retVal)
  rownames(results) <- as.character(object@armDef@categories)
  
  colnames(results) <- c("Total",listColumnDefSlot(object@subgroupDef,"displayName"))
  t(results)
}



#create the Number of patients per subgroup dataframe
#for SurvivalData object
subgroupsSummary <- function(object, colWidth=1.5){
  
  #create the summary data frame
  df <- extractSubgroupTable(object)
  
  #Then convert to FlextTable output
  numRows <- nrow(df)
  numCols <- length(object@armDef@categories) + 1
  
  MyFTable <- FlexTable(numrow=numRows,numcol=numCols, 
                        body.par.props=parProperties(text.align="right"),
                        header.text.props = textProperties(font.weight = "bold"))
  
  #Add in data
  for(thisRow in 1:numRows){
    MyFTable[thisRow, 2:numCols] <- df[thisRow,(numCols-1):1]
    MyFTable[thisRow,1] <- rownames(df)[thisRow]
  } 
  
  #Set properties
  MyFTable[1:numRows,1] <- parProperties(text.align="left")
  MyFTable[1:numRows,1] <- textProperties(font.weight = "bold")
  
  MyFTable[1:numRows,1:numCols,side='bottom'] <- borderProperties(width=0)
  MyFTable[1:numRows,1:numCols,side='left'] <- borderProperties(width=0)
  MyFTable[1:numRows,1:numCols,side='top'] <- borderProperties(width=0)
  MyFTable[1:numRows,1:numCols,side='right'] <- borderProperties(width=0)
  
  MyFTable[numRows,1:numCols,side='bottom'] <- borderProperties(width=3)
  MyFTable[1,1:numCols,side='top'] <- borderProperties(width=3)
  
  #Add headers
  hR <- FlexRow(c("","arm"),colspan=c(1,numCols-1), par.properties=parProperties(text.align="center"),
                text.properties = textProperties(font.weight = "bold"),
                cell.properties = cellProperties(border.top.width=3, border.bottom.width=0,
                                                 border.left.width=0, border.right.width=0)
  )
  hR2 <- FlexRow(c("",rev(as.character(object@armDef@categories))),
                 par.properties=parProperties(text.align="center",padding=1),
                 text.properties = textProperties(font.weight = "bold"),
                 cell.properties = cellProperties(border.width=0)
  )
  
  MyFTable <- addHeaderRow(MyFTable,hR)
  MyFTable <- addHeaderRow(MyFTable,hR2)
  
  setFlexTableWidths(MyFTable,rep(colWidth,numCols))
  MyFTable
  
}


#Function to output the header rows for the covariate and endpoint summary tables
#inputs subgroupDetails - the output from extractSubgroupTable function above
#leftCol1 a single character or vector of two characters for the two left-most column
#headings on top header row
#leftCol2 a length 2 character vector for the two left-most column
#headings on second header row
getHeaders <- function(subgroupDetails, leftCol1, leftCol2){
  
  numArms <- ncol(subgroupDetails)
  numSubgroups <- nrow(subgroupDetails)-1
  
  leftSpan <- if(length(leftCol1)==1) 2 else c(1,1)
  
  rownames(subgroupDetails)[1] <- "All" 
  
  #First row
  headerNames <- paste("\n(total=",rowSums(subgroupDetails),")",sep="")
  headerNames <- paste(rownames(subgroupDetails), headerNames)
  
  hR <- FlexRow(c(leftCol1,headerNames),colspan=c(leftSpan,rep(numArms,1+numSubgroups)), 
                par.properties=parProperties(text.align="center"),
                text.properties = textProperties(font.weight = "bold"))
  
  #reverse order so that control arm is last (done this way so works if 1 row) 
  subgroupDetails[,1:ncol(subgroupDetails)] <- subgroupDetails[,ncol(subgroupDetails):1]
  colnames(subgroupDetails) <- rev(colnames(subgroupDetails))
  
  #then arm header (second header row)
  armHeaders <- vapply(seq_len(nrow(subgroupDetails)),
                       function(x){
                         paste(colnames(subgroupDetails),"\n(total=",subgroupDetails[x,],")",sep="")  
                       },FUN.VALUE=character(numArms))
  
  hR2 <- FlexRow(c(leftCol2,armHeaders),
                 par.properties=parProperties(text.align="center",padding=1),
                 text.properties = textProperties(font.weight = "bold"))
  
  return(list(hR,hR2))
}

