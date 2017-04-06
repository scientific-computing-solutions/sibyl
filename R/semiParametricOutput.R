##' @include semiParametric.R
NULL 

##' @name summary
##' @rdname summary-methods
##' @aliases summary,SemiParametricModel-method
##' @param class ('data.frame' or "FlexTable' (default)) type of output for summary table
##' @export
setMethod("summary", signature(object="SemiParametricModel"),
  function(object, type=c("medianTTE", "KM")[1], class=c("data.frame","FlexTable")[2], 
           digits=if(type=="medianTTE") 2 else 5){
    
    if(length(class) != 1 || !class %in% c("data.frame","FlexTable")){
      stop("Invalid class argument, should be 'data.frame' or 'FlexTable")
    }
    
    if(length(type)!=1 || !type %in% c("medianTTE", "KM") ){
      stop("Invalid type argument must be one of 'medianTTE' or 'KM'")
    }
    
    if(length(digits)!=1 || !is.numeric(digits) || !digits > 0 || is.infinite(digits) ||
       is.na(digits)){
      stop("Invalid digits argument")
    } 
    
    if(type=="KM"){
      return(kmsummary(object, class, digits))
    }
    
    #calculate the median time to event with confidence interval
    results <- as.data.frame(quantile(object@km,prob=0.5,conf.int = TRUE))
            
    #bind on the number of events
    summaryKM <- summary(object@km)
    #notice the comma in the line below to handle different cases
    eventDetails <- if(isSingleArm(object)) summaryKM$table["events"] else summaryKM$table[,"events"] 
    results <- cbind(eventDetails, results)
            
    #set row/column names
    rownames(results) <- getArmNames(object@survData)
    colnames(results) <- c("Total Events","Median time to event","Lower.CI","Upper.CI")
            
    #calculate ratios
    if(nrow(results) > 2){
      warning("Not calculating ratios/differences when there is more than one arm")
    }
    else if(nrow(results) == 2){
      ratio <- apply(results,2,function(x){x[2]/x[1]})
      difference <- apply(results, 2, function(x){x[2]-x[1]})
      results <- rbind(results,ratio=ratio)
      results <- rbind(results,difference=difference)
    }
            
    #transpose to get arms as columns
    results <- t(results)
    
    #reorder
    armNames <- as.character(getArmNames(object@survData))
    numArms <- length(armNames)
    
    results[,1:numArms] <- results[,numArms:1]
    colnames(results)[1:numArms] <- colnames(results)[numArms:1]
            
    if(class=="data.frame"){
      return(results)
    }
    
    #Now create Flex Table        
    numRows <- 4
    numCols <- 1 + ncol(results)
   
    MyFTable <- FlexTable(numrow=numRows,numcol=numCols, 
                          body.par.props=parProperties(text.align="right"),
                          header.text.props = textProperties(font.weight = "bold"),
                          body.cell.props = cellProperties(padding.right=1))
    
    #Add number of events (want to be integers so have to do this separately)        
    MyFTable[1,1+(1:numArms)] <- results[1,1:numArms]
            
    #Add ratio of #events if 2 arms
    if(numArms==2){
      MyFTable[1,(numCols-1)] <- round(results[1,numArms+1],digits=digits)
      MyFTable[1,numCols] <- round(results[1,numArms+2],digits=digits)
    }
    
    #Add in rest of data        
    MyFTable[2:numRows,2:numCols] <- round(results[2:numRows,],digits = digits)
    
    #set left column to include appropriately formatted text        
    MyFTable[1:numRows,1] <- c("Total number of events", "Median time to event",
                               "95% lower CI", "95% upper CI")
            
    MyFTable[1:numRows,1] <- parProperties(text.align="left")
    MyFTable[1:numRows,1] <- textProperties(font.weight = "bold")
            
          
    #Set borders
    MyFTable[1:numRows,1:numCols,side='bottom'] <- borderProperties(width=0)
    MyFTable[1:numRows,1:numCols,side='left'] <- borderProperties(width=0)
    MyFTable[1:numRows,1:numCols,side='top'] <- borderProperties(width=0)
    MyFTable[1:numRows,1:numCols,side='right'] <- borderProperties(width=0)
            
    MyFTable[numRows,1:numCols,side='bottom'] <- borderProperties(width=3)
    MyFTable[1,1:numCols,side='top'] <- borderProperties(width=1)
            
    
    #create text for header row        
    armCounts <- vapply(armNames,function(x){
      sum(as.character(object@survData@subject.data$arm)==x)
    },FUN.VALUE=numeric(1))
            
    armHeaders <- paste(armNames,"\n(total=",armCounts,")",sep="")
            
    if(length(armNames)==2){
      armHeaders <- c(rev(armHeaders),"Ratio", "Difference")
    }
            
    hR <- FlexRow(c("",armHeaders),
                  par.properties=parProperties(text.align="center",padding=1),
                  text.properties = textProperties(font.weight = "bold"),
                  cell.properties = cellProperties(border.top.width=3, border.bottom.width=0,
                                                   border.left.width=0, border.right.width=0))
            
    MyFTable <- addHeaderRow(MyFTable,hR)
            
    MyFTable
})


#function to output the KM summary for the SemiParametricModel
kmsummary <- function(object, class, digits){
  
  #get summary
  s <- summary(object@km, censored=TRUE)
  #coerce to data frame
  
  arm <- if(isSingleArm(object)) getArmNames(object@survData) else s$strata
  
  summaryData <- data.frame(arm=arm,time=s$time,survival=s$surv,
                            n.risk=s$n.risk,n.event=s$n.event, std.err=s$std.err,
                            x=s$lower, y=s$upper)
  colnames(summaryData)[7:8] <- paste(c("lower","upper"), "95% CI")
  
  #split by arm
  listOfData <- split(summaryData, summaryData$arm)
  
  #set names of list elements as arms
  names(listOfData) <-  as.character(getArmNames(object@survData))
  
  #remove arm column from data frames
  listOfData <- lapply(names(listOfData), function(arm){
     listOfData[[arm]]$arm <- NULL
     listOfData[[arm]]
  })
  
  #set names of list elements as arms 
  names(listOfData) <-  as.character(getArmNames(object@survData))
  
  if(class=="data.frame") return(listOfData)
  
  #create flextable
  retVal <- lapply(listOfData, function(oneArmData){
    FlexTable(round(oneArmData,digits=digits),
              body.par.props=parProperties(text.align="right"),
              header.text.props = textProperties(font.weight = "bold"),
              body.cell.props = cellProperties(padding.right=1))
  })
  
  names(retVal) <-  as.character(getArmNames(object@survData))
  retVal 
}

##' Plot methods for Sibyl package
##' @rdname plot-methods
##' @aliases plot,SemiParametricModel,missing-method
##' @param x (SemiParametricModel object) contains data to be plotted
##' @param type (character) the type of plot to be created; one of "KM",
##'        "CumHaz", "LoglogS", "LogoddS" or "InvNormS"
##' @param logTime (logical) determines if x axis of plot is time or log(time)
##' @param armColours (vector of colours) the colours for the treatment arms for all
##' graphs except the KM curve (use the col argument to set the colour for the KM graph)
##' @param ... arguments to be passed to azplot.km when type = "KM"
##' @export
setMethod("plot", signature(x="SemiParametricModel", y="missing"),
  function(x, type=c("KM","CumHaz","LoglogS","LogoddS","InvNormS")[1], 
           use.facet=TRUE, logTime=NULL, 
           armColours=c("black", "red", "blue", "green", "yellow", "orange"), ...){
            
    #set defualt logTime
    if(is.null(logTime)){
      logTime <- tolower(type) !="cumhaz"
    }
            
    switch(tolower(type),
           "km"=kmPlot(x, ...),
           "cumhaz"=diagnosticPlot(x,logTime, yval="-log(S)", use.facet, cumHaz=TRUE, armColours=armColours),
           "loglogs"=diagnosticPlot(x, logTime, yval="log(-log(S))", use.facet, armColours=armColours),
           "logodds"=diagnosticPlot(x, logTime, yval="log(S/(1-S))",  use.facet, armColours=armColours),
           "invnorms"=diagnosticPlot(x, logTime, yval="qnorm(1-S)", use.facet, armColours=armColours),
           stop("type must equal one of KM, CumHaz, LoglogS or Logodds or InvNorms")
    )
  }
)


#Output KM plot
kmPlot <- function(x, ...){
  azplot.km(x@km, ...)
}


#output a ggplot object for diagnostic plots
#x is SemiParametricModel
#logTime logical, true if using log(time) on x-axis, false otherwise
#yval is a string representing the function to be plotted on the y-axis
#as a function of S, the cumulative hazard, for example "log(S/(1-S))"
#use.facet is TRUE if splitting the data into separate graphs per arm
#cumHaz is true if we are outputting cumulative hazard plot
diagnosticPlot <- function(x, logTime, yval, use.facet, cumHaz=FALSE, armColours){
  
  #cannot use facet if this is a single arm trial
  if(use.facet && isSingleArm(x)){
    stop("use.facet must be FALSE when trial is a single arm trial")
  }
  
  #R-cmd-check thinks t, s, model, ... are global
  #variables inside the ggplot commands so complains about them
  #they are not global variables but for them to pass R-cmd-check we
  #need to create dummy variables
  Arm <- NULL
  
  #x-axis value for the graphs
  xval <- if(logTime) "log(t)" else "t"
  
  #xlab
  xlabel <- if(logTime) "log(time)" else "time"
  
  #extract the cumulative hazard data from the survfit object
  armNames <- getArmNames(x@survData)
  data <- extractCumHazData(x@km, armNames, isSingleArm=isSingleArm(x))
  
  # Concatenate list of data frames row-wise into one big data frame
  data <- do.call("rbind", data)
 
  #remove the edge cases
  if(logTime){
    data <- data[data$t > 0,]
  }
  
  if(cumHaz){
    data <- data[data$S > 0,]
  }
  else{
    data <- data[data$S > 0 & data$S < 1,]
  }
  
  # Create the plot
  p <- ggplot(data, aes_string(color="Arm",x=xval,y=yval))
  
  # Set colour for each arm
  if (length(armNames) <= length(armColours)){
    p <- p + scale_colour_manual(values = armColours)
  }
  else{
    stop("The current value of armColours parameter supports a maximum of ",length(armColours)," arms. ",
         "Use armColours argument to choose colours for each arm")
  }
  
  #Add facet
  if(use.facet){
    p <- p + facet_grid(. ~ Arm)
  }
  
  if(cumHaz){
    #plot the data
    p <- p + geom_step()
  }
  else{
    #plot the data and best fit line
    p <- p + geom_point()
    p <- p + stat_smooth(method="lm", se=FALSE, size=1, aes(color=Arm))
  }
  
  
  #Add xlabel
  p <- p + xlab(xlabel)
  
  # Format background and borders
  p <- p + theme(panel.background = element_blank(),
                 panel.border = element_rect(colour = "black", fill = NA),
                 panel.grid.major = element_line(colour="grey", linetype="dotted"),
                 panel.grid.minor = element_blank())
  
  p
}