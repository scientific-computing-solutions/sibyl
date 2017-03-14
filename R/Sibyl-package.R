#' @importFrom graphics legend lines par
#' @importFrom methods new show slot
#' @importFrom stats AIC BIC as.formula predict quantile formula median sd update.formula extractAIC resid integrate
#' @importFrom utils str tail
#' @import flexsurv
#' @import ggplot2
#' @import survival
#' @importFrom azGraphics azplot.km
#' @import azRMST
#' @import ReporteRs
NULL



##' Sibyl: The package fits parametric models to time-to-event data and
##' extrapolates these curves to larger surivival times.
##' @name Sibyl
##' @docType package
NULL


#' Example time to event data set for subjects trying to stop smoking
#' @details This data set is loosely based on an example
#' data frame from the asaur package
#' 
#' @format A data frame with the following columns:
#' \describe{
#'   \item{ID}{the unique subject id}
#'   \item{age}{age of the subject, in years}
#'   \item{race}{'white', 'black', 'hispanic' or 'other'}
#'   \item{grp}{'combination' or 'patchOnly'; the name of the treatment arm to which the subject was assigned}
#'   \item{sub.isMale}{ndicator for membership of the subgroup 'isMale'}
#'   \item{sub.isHeavySmoker}{indicator for membership of the subgroup 'isHeavySmoker'}
#'   \item{ttr}{time to relapse, in days}
#'   \item{ttr.cens}{indicator for relapse}
#'   \item{end.2}{additional endpoint}
#'   \item{cens.2}{indicator for additional endpoint}
#' } 
"sibylData"


##' Method to extract the endpoint units for a given object
##' @rdname getEndpointUnits-methods
##' @name getEndpointUnits
##' @param object The object whose endpoint units is required
##' @return The units of the endpoints for the given object
##' @export
setGeneric( "getEndpointUnits", function(object)
  standardGeneric("getEndpointUnits"))




