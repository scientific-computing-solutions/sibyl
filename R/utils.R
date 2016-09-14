##' Some miscellaneous utility code for wrangling flexsurv models
##' which should probably belong in the flexsurv package.
##'
##' This code is densely written; there's not a lot of waste lines.
NULL

##' Bootstrap flexsurvreg model and get new parameters for each
##' subject.
##'
##' The return type of normboot.flexsurvreg is unhelpful for future
##' use; this function adapts it into a slightly more useful form.
##' 
##' @param model a flexsurvreg model
##' @param newdata a data frame defining the subjects for whom we want
##' parameters
##' @param nsim the number of simulations to do
##' @return an appropriately sized array
flexsurvParms <- function(model, newdata, nsim) {
  ## first get the parameters to use from the flexsurvreg bootstrap
  sim.params <- normboot.flexsurvreg(model, B=nsim, newdata=newdata)
  ## now roll them all up into an array to simplify future manipulations
  array(unlist(sim.params),
        dim=c(dim(sim.params[[1]]), length(sim.params)),
        dimnames=dimnames(sim.params[[1]]))
}

##' Build a callable CDF from a flexsurv model
##' @param model the flexsurv model
##' @param times the times at which we want output
##' @return a callable function
buildCdf <- function(model, times) {
  cdf <- model$dfns$p
  times <- list(times)
  ## and here's the return value; a function
  function(params) {
    do.call(cdf, c(times, params))
  }
}


##' Average a flexsurv model over provided data
##' @param model a flexsurvreg object
##' @param newdata a data frame compatible with the model
##' @param averageLevels a factor defining the rows of \code{newdata}
##' to average over
##' @param times the times at which predictions are required
##' @param nsim the number of simulations to do
##' @return a \code{flexsurvPred} object
averageModel <- function( model, newdata, averageLevels, times, nsim=1000 ) {
  the.call <- match.call()
  ## let's do a little checking
  if (!is.factor(averageLevels)) {
    stop("the levels to average over isn't a factor")
  }
  
  if (length(averageLevels) != nrow(newdata)) {
    stop("Levels to average over don't match data")
  }
  
  ## get our hands on the CDF
  cdf <- buildCdf( model, times )
  
  ## let's get our sample parameters so we can then predict the
  ## curves
  sim.params <- flexsurvParms( model, newdata, nsim )
  
  ## and now produce survival curves averaged over the appropriate
  ## groups; this code is dense
  
  curves <-
    by(seq_len(nrow(newdata)), # get the indices of patients
       averageLevels,          # in each averaging group
       function(indices) {
         ## for one averaging group, extract the parameters
         apply(sim.params[,, indices, drop=FALSE],
               1, # compute the survival curves for each patient
               function(x) {
                 ## "1 -" is here to turn into survival curves
                 ## the rowMeans averages over a Monte-Carlo
                 1 - rowMeans(apply(x, 2, cdf))
               })
       })
  
  ## and package everything up for return
  result <-
    list(call=the.call,
         curves=curves,
         times=times,
         nsim=nsim)
  ## that was easy. :-)
  class(result) <- "flexsurvPred"
  result
}

print.flexsurvPred <- function(x) {
  cat("\nCall:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  
  cat("Averaging over:",
      paste0("    ", names(x$curves), collapse="\n"),
      sep="\n")
  
  invisible(x)
}

summary.flexsurvPred <- function(object, prob=c(0.025, 0.5, 0.975)) {
  newDim <- c(length(prob), length(object$times))
  quantiles <- 
    lapply(object$curves,
           function(curveSet, newDim) {
             res <- apply(curveSet, 1, quantile,
                          prob=prob, names=FALSE)
             dim(res) <- newDim
             res
           },
           newDim)
  res <-
    list(call=object$call,
         quantiles=quantiles,
         prob=prob,
         times=times)
  class(res) <- "summary.flexsurvPred"
  res
}

print.summary.flexsurvPred <- function(x) {
  cat("\nCall:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  
  cat("Averaging over:",
      paste0("    ", names(x$quantiles), collapse="\n"),
      sep="\n")
  
  cat("\nQuantiles: ",
      paste0(x$prob, collapse=", "), "\n",
      sep="")
  
  invisible(x)
}

## and now let's draw some pretty pictures

lines.summary.flexsurvPred <- function(x, col,
                                       lty=c('dashed', 'solid', 'dashed'),
                                       xscale=1, ...) {
  ## warn if we've got silly inputs
  if (length(lty) != length(x$prob)) {
    warning("lty does not match quantiles given")
  }
  
  if (length(col) != length(x$quantiles)) {
    warning("col does not match average groups")
    col <- rep_len(col, length(x$quantiles))
  }
  
  ## these are the times at which to plot
  times <- x$times / xscale
  quants <- x$quantiles
  
  for (ind in seq_along(quants)) {
    matlines(times, t(quants[[ind]]),
             lty=lty,
             col=col[[ind]])
  }
  
  invisible(NULL)
}

lines.flexsurvPred <- function(x,
                               prob=c(0.025, 0.5, 0.975),
                               col,
                               lty=c('dashed', 'solid', 'dashed'),
                               xscale=1, ...) {
  lines(summary(x, prob=prob),
        col=col, lty=lty, xscale=xscale, ...)
}