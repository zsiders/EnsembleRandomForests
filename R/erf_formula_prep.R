#' Make a formula for ERF use
#' 
#' @description Takes a column names and makes the RF model formula
#' 
#' @param var A character string indicating the column name of the data frame that contains the number of interactions for the ERF to model; column should be a numeric column
#' @param covariates A character vector indicating the column name(s) of the data frame that contain the covariates
#' 
#' @return A formula class object
#' 
#' @export
#' 
erf_formula_prep <- function(var, covariates){
	form <- formula(paste0(var," ~ ",paste(covariates, collapse=" + ")))
	return(form)
}