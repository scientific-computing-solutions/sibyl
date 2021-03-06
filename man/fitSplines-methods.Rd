% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/splineFit.R
\docType{methods}
\name{fitSplines}
\alias{fitSplines}
\alias{fitSplines,SurvivalModel-method}
\alias{fitSplines,SurvivalModel-methods}
\title{Method to create SurvivalModel object with different spline models}
\usage{
fitSplines(object, ...)

\S4method{fitSplines}{SurvivalModel}(object, k = 0:5, scale = c("hazard",
  "odds", "normal")[1])
}
\arguments{
\item{object}{\code{SurvivalModel} for which spline models are to be fitted}

\item{...}{additional arguments for specific instances of this generic}

\item{k}{(numeric vector) The number of knots to be used when fitting spline model
(if k=c(2,3) then models will be fit with both 2 and 3 knots)}

\item{scale}{(character) The scale arguement to be passed to flexsurvspline}
}
\value{
A SurvivalModel object with different spline model fits
}
\description{
Method to create SurvivalModel object with different spline models
}
