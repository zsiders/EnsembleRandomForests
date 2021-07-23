#' Determines the minimum sample size in the first RF bag
#' 
#' @description Takes a data.frame object determines the minimum split
#' 
#' @param v A data frame object returned from erf_data_prep()
#' @param p A numeric value between (0,1), default is just below 90%
#' @param nmax A numeric value > 0 specifying the maximum number of observations per bag, default is 1e4
#' 
#' @return A numeric value specifying the minimum split
#' @export
#' 
#' @examples
#' v <- data.frame(var=rbinom(1000,1,0.2), covar=rnorm(1000))
#' min_split(v)
#' 
min_splitter <- function(v, p=0.89, nmax=1e4){
	t <- table(v[,1]) #assumes variable of interest is in the first column
	#p is raised to the power of 2 for the two bagging events
	min_split <- pmin(nmax,as.numeric(floor(min(t*p^2))))
	return(min_split)
}