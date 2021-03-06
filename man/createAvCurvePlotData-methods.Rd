% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/averageCurvePlotData.R
\docType{methods}
\name{createAvCurvePlotData}
\alias{createAvCurvePlotData}
\alias{createAvCurvePlotData,SurvivalModel-method}
\alias{createAvCurvePlotData,SurvivalModel-methods}
\title{Method to create the \code{AvCurvePlotData} object}
\usage{
createAvCurvePlotData(object, ...)

\S4method{createAvCurvePlotData}{SurvivalModel}(object, maxTime = NULL,
  Npoints = 201, Nsim = 500, models = NULL, seed = NULL, B = 1000,
  conf.type = c("none", "plain", "log", "log-log")[3])
}
\arguments{
\item{object}{\code{SurvivalModel} object}

\item{...}{Additional arguments for specific instances of this generic}

\item{maxTime}{(numeric) the time to extrapolate the curves to (default=NULL implies no extrapolation)}

\item{Npoints}{(numeric) the number of time points for which the survival curves are to be evaluated at}

\item{Nsim}{(numeric) the number of simulations used to generate the averaged curves}

\item{models}{(character vector) which models from \code{names(object@models)} are to be used when 
calculating averaged survival curves - default NULL implies use all}

\item{seed}{(numeric, default NULL) if not NULL then set random seed (although it will only be
used if models include covariates) - random numbers are used when models includes covariates 
or when confidence intervals are required}

\item{B}{(integer) Only used when no covariates in model. See summary.flexsurvreg 
Number of simulations from the normal asymptotic distribution of the estimates used 
to calculate confidence intervals. Decrease for greater speed at the expense of accuracy, 
or set B=0 to turn off calculation of CIs.}

\item{conf.type}{("none", "plain", "log" [default], "log-log") argument passed to survfit}
}
\value{
A \code{AvCurvePlotData} object
}
\description{
Method to create the \code{AvCurvePlotData} object
}
\details{
If the models include covariates a simulation procedure is required to generate averaged survival
curves. If the models do not include covariates then \code{summary.flexsurvreg} is used
}
