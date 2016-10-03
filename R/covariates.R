# This file contains the public functions associated with Covariates object


##' @include common.R 

# # Function that will check the validity of Covariates object
# # when constructed
# # @param object Covariates object
 checkCovariates <- function( object ){
   errors <- ""
   
   def.types <- c( "numeric", "categorical" )
   if( !object@type %in% def.types  ){
     errors <- paste0( errors, "type must be: ", def.types  )
   }
   if( object@type=="categorical" && ( object@categories == "" || is.null( object@categories ) ) ){
     errors <- paste0( errors, "categories need to be defined for categorical variables: ", object@categories  )
   }
   if( errors == "" ) TRUE else errors 
}
   
##' Class representing data to use for extrapolation
##' @slot subject.data a data frame with X columns
##' "subject",...
##' see vignette and \code{Covariates} for further iunformation
##' @export
setClass("Covariate", 
          slots= list( name = "character", 
                       type = "character", 
                       categories = "character",
                       unit = "character"
                       ),
          validity = checkCovariates )


##' Constructor for event data object
##' 
##' All dates must be in one of the following formats:
##' YYYY-MM-DD, DD/MM/YY or DD Month YYYY
##'  
##' 
##' @param data A data frame 
##' @param subject string, the column name of subject identifiers
##' @return An \code{Covariates} object
##' @export
Covariate <- function( name, type, categories="", unit="" ){
  
  return( new( "Covariate", name=name, type=type, categories=categories, unit=unit ) )
}




