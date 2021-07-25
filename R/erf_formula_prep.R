#' Make a formula for ERF use
#' 
#' @description Takes a column names and makes the RF model formula
#' 
#' @param var A character string indicating the column name of the data frame that contains the number of interactions for the ERF to model; column should be a numeric column
#' @param covariates A character vector indicating the column name(s) of the data frame that contain the covariates
#' @param add.random A logical flag for adding a random variable named 'random' to the formula; needed for the ERF models
#' 
#' @return A formula class object
#' 
#' @export
#' 
#' @examples
#' erf_formula_prep('death', c('sarlaac','lightsaber','wookie','stormtrooper'))
#' 
erf_formula_prep <- function(var, covariates, add.random=TRUE){
	if(add.random) covariates <- c(covariates, 'random')
	form <- formula(paste0(var," ~ ",paste(covariates, collapse=" + ")))
	return(form)
}