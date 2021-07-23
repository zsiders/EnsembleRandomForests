#' @title Example dataset for use with ERF
#'
#' @description A data set one dependent variable and five covariates
#'
#' @format A list with 3 objects:
#' \describe{
#' 	\item{samples}{A data.frame with the sampled locations to be used in the ERF model}
#' 	\item{betas}{A vector containing the beta coefficients used to generate the log odds}
#' 	\item{grid}{A raster brick containing the full gridded covariates, log odds, and the probability}
#' }
#' @source Generated using the `data_sim()` function
"simData"