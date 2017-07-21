#Code given to Sibyl from an internal package (hence the different convention for naming arguments etc.)

##' Create survival plot with numbers at risk. 
##' @param x An objects of class 'survfit', created by the 'survfit' function.
##' @param times The timepoints at which the number of patients at risk should be
##'        displayed.
##' @param fun An arbitrary function defining a transformation of the survival curve.
##'        See the plot.survfit help for details.
##' @param main The main title for the plot. ("")
##' @param xlab The x-axis label. ("Time")
##' @param ylab The y-axis label. ("Proportion with event")
##' @param lty The types of lines for the Kaplan-Meier curves. Valid values are 
##'        1, 2, 3, ... (1:n)
##' @param pch Plot symbols used to label the Kaplan-Meier curves. Set to NA if 
##'        plot symbols are not to be included. Valid values are 1, 2, 3, ...(1:n)
##' @param col Color codes for the Kaplan-Meier curves. Valid values are 1, 2, 3, 
##'        ..., with a character string giving the color name (e.g., "red").
##' @param mark.time Default TRUE; Controls the labeling of the curves. If set to 
##'        FALSE, no labeling is done. If TRUE, then curves are marked at each 
##'        censoring time as defined by \code{x$time[x$n.censor > 0]}. If mark.time 
##'        is a numeric vector, then curves are marked at the specified time points.
##' @param cex Character expansion for plot symbols. The bigger the value, the 
##'        larger the characters on the plot. (1)
##' @param cex.legend Character expansion for the legend information. (0.9)
##' @param cex.axis Character expansion for axis labels. (1)
##' @param cex.nrisk Character expansion for number at risk table. (0.8)
##' @param cex.lty The types of lines for the Kaplan-Meier curves (used for legend). 
##'        Valid values are 1, 2, 3, ... (1:n).
##' @param labels Labels for the number at risk information. Defaults to the same 
##'        labels retrieved from the survfit object.
##' @param legend The location of the legend. Should be one of "bottomright", "bottom", 
##'        "bottomleft", "left", "topleft", "top", "topright", "right", "center". 
##'        Set to FALSE if legend is not to be included. ("topright")
##' @param n.at.risk Logical value specifying whether or not the Number at risk 
##'        should be included below the plot. (FALSE)
##' @param count.overlap Logical value specifying whether or not to display the 
##'        number of overlapping censored patients in the plot. (FALSE)
##' @param ... Arguments to be passed to the plot.survfit method.
##' @examples
##' lung.km <- survfit( Surv(time, status) ~ sex, data=lung)
##' kmPlotWrapper( lung.km, main="Kaplan-Meier plot", xlab="Time (days)", pch=NA, 
##' legend="topright", times=seq(0,1000,by=200), n.at.risk=TRUE )
##' @export
kmPlotWrapper <- function(x,
                      main="",
                      xlab="Time",
                      ylab="Survival probability",
                      lty=NULL,
                      col=NULL,
                      pch=NULL,
                      labels=NULL,
                      legend="topright",
                      mark.time=TRUE,
                      cex = 1,
                      cex.legend = 0.7,
                      cex.axis=1,
                      cex.nrisk=0.8,
                      cex.lty = NULL,
                      times=NULL,
                      n.at.risk=FALSE,
                      count.overlap=FALSE,
                      fun=function(x) x,
                      ...
){
  
  if(class(x) != "survfit") stop(paste("object x not of class survfit"))
  if(is.null(times)) times = round(quantile(c(0, max(x$time))))
  grps <- x$strata
  n.grps <- max( 1, length(grps) )
  
  # Line types, symbols, colors and labels for treatment groups
  if(is.null(lty)) lty <- 1:n.grps
  if(is.null(cex.lty)) cex.lty <- lty
  if(is.null(pch)) pch <- 1:n.grps
  if(is.null(labels)) {
    labels <- if( n.grps == 1 ) "Single arm" else names(grps) 
  }
  if(is.null(col)) {
    col <- c("#F8766D", "#C49A00", "#53B400", "#02C094", "#00B6EB", "#A58AFF", "#FB61D7")
    if( n.grps == 1 ) col[3] <- col[2]
  }
  # This is to enforce that all censor events are showed in the KM plot.
  #if(mark.time==TRUE){mark.time=x$time[x$n.censor > 0]}
  
  xaxt <- NULL  # Flag for keeping x-axis
  if(!is.null(times)) xaxt <- "n"  # Flag for removing x-axis in the plot if user specified time
  
  # Reserve space for number at risk in plot
  if(n.at.risk==TRUE){
    plt <- par("plt")
    plt.old <- par("plt")
    mar.old <- par("mar")
    alpha <- 0.2 + 0.06 * n.grps
    beta <- 1
    if(alpha > 0.6) {
      beta <- max((0.6/alpha), 0.6)
      alpha <- 0.6
    }
    plt[4] <- 0.95
    plt[3] <- alpha * 0.95
    par(list(plt = plt, mar=c(8+1.6*n.grps,4,4,1)))
    on.exit(par(list(plt = plt.old, mar=mar.old)))
  }
  
  plot(x, main=main, xlab=xlab, ylab=ylab, lty=lty, col=col, mark=pch, cex=cex, 
       xaxt=xaxt, yaxt="n", yaxs="i", fun=fun, mark.time=mark.time, ...)
  
  # Display time points specified by user (or at quartiles)
  axis(1, at=times, cex.axis=cex.axis)
  
  # Make y-labels horizontal
  axis(2, las=1, cex.axis=cex.axis)
  
  if(legend != FALSE){
    legend(legend, legend=rev(labels), 
           col=rev(col), lty=cex.lty, pt.cex=cex, cex=cex.legend, bty="n", 
           text.col=rev(col), 
           pch=rev(pch))
  }
  
  if(n.at.risk == TRUE){
    # Group names
    group.name.pos <- 0
    
    # Number at risk
    kms <- summary(x, time=times, extend=TRUE)
    if( !is.null( kms$strata )) {
      d1 <- data.frame(time = kms$time, n.risk = kms$n.risk, strata = c(kms$strata))
      d2 <- split(d1, f=d1$strata)
    } else {
      d1 <- data.frame( time = kms$time, n.risk = kms$n.risk )
      d2 <- NULL
    }
    
    # Right justify numbers
    if( is.null(d2) ){
      n.digits <- nchar( d1$n.risk )
      max.len <- length(n.digits)
      L <- data.frame( n.digits )
      nd <- n.digits
    } else{
      ndigits <- lapply(d2, function(x) nchar(x[,2]))
      max.len <- max(sapply(ndigits, length))
      L <- do.call('rbind', lapply(ndigits, function(x){length(x) <- max.len; x} ))
      nd <- apply(L, 2, max, na.rm=T)
    }
    
    mtext(side=1, text="Number at risk", at=group.name.pos, line=3, cex=cex.nrisk, adj=0, col=1, las=1)
    
    if( n.grps == 1 ) d2[[1]] <- d1    
    for(i in 1:n.grps){
      this <- d2[[i]]
      w.adj <- strwidth('0', font=par('font'))/5 * nd[1:nrow(this)]
      mtext(side=1, at=group.name.pos, text=labels[i], line=4+2*(n.grps-i), cex=cex.nrisk, adj=0, col=col[i], las=1)
      mtext(side=1, at=this$time+w.adj, text=this$n.risk, line=5+2*(n.grps-i), cex=cex.nrisk, adj=1, col=col[i], las=1)
    }
  }
  
  if(count.overlap == TRUE){
    text(x$time[x$n.censor>1], sapply(x$surv[x$n.censor>1], fun), labels=x$n.censor[x$n.censor>1], pos=3, cex=0.7)
  }
}
