% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/erf_impt.R
\name{erf_impt}
\alias{erf_impt}
\title{Visualizes variable importance from ERF}
\usage{
erf_impt(fit, var.names, pal)
}
\arguments{
\item{fit}{The fitted object returned from calling ens_random_forests()}
}
\value{
A plot of the variable importance distributions ordered by the mean variable importance.
}
\description{
Plots the variable importance from the Ensemble Random Forests model
}
\examples{
#run an ERF with 10 RFs and 
ens_rf_ex <- ens_random_forests(df=simData$samples, var="obs", covariates=grep("cov", colnames(simData$samples),value=T), save=FALSE, cores=1)

erf_impt(ens_rf_ex)

}
