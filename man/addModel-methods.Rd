% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/survivalModels.R
\docType{methods}
\name{addModel}
\alias{addModel}
\alias{addModel,SurvivalModel-method}
\alias{addModel,SurvivalModel-methods}
\title{Method to add extra models to be fitted}
\usage{
addModel(object, ...)

\S4method{addModel}{SurvivalModel}(object, modelName, modelOptions = NULL,
  preferredPackage = getDefaultPackage(), suppressOverwriteWarning = FALSE)
}
\arguments{
\item{object}{(SurvivalModel object) contains list of already fitted models}

\item{...}{additional arguments for specific instances of this generic}

\item{modelName}{(list of character strings) names of models to be added}

\item{modelOptions}{(named list) extra options to be passed to the fitting functions.
For example modelOptions=list(spline=list(k=4)) will pass the argument k=4 to flexsurvSpline}

\item{preferredPackage}{(character string) name of the preferred package
in which to look for the specified model}

\item{suppressOverwriteWarning}{(logical) should be FALSE so that a warning is displayed
when overwriting existing model}
}
\value{
a \code{SurvivalModel} object
}
\description{
Method to add extra models to be fitted
}
