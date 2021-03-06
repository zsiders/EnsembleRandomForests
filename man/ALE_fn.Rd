% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ALE_fn.R
\name{ALE_fn}
\alias{ALE_fn}
\title{Accumulated Local Effects workhorse}
\usage{
ALE_fn(X, X.model, pred.fun, J, K = 40)
}
\arguments{
\item{X}{the data.frame to get the covariate from}

\item{X.model}{the name of the response variable}

\item{pred.fun}{a function to calculate new predictions from the model}

\item{J}{the column index of the covariate of interest}

\item{K}{an integer value that determines the number of "windows" or breaks to calculate the model predictions over. More increase computational time but serves smooths the ALE predictions.}
}
\value{
A list that contains:
\itemize{
\item \strong{K}: the number of realized breaks
\item \strong{x.values}: the break values trialed
\item \strong{class}: the class of the covariate
\item \strong{quantile}: the quantile of the breaks
\item \strong{fJ}: the ALEs evaluated at a given x value
}
}
\description{
Calculates the Accumulated Local Effects (ALE) from a given data.frame, model, predictions, and covariate. This riffs on the ALEPlot function available in the ALEPlot package.
}
