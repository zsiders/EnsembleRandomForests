#' Prepare a data.frame for ERF
#' 
#' @description Prepares a data.frame object for use with the ERF (not needed if using erf())
#' 
#' @param df A data.frame object
#' @param var A character string indicating the column name of the data frame that contains the number of interactions for the ERF to model; column should be a numeric column
#' @param covariates A character vector indicating the column name(s) of the data frame that contain the covariates
#' @param header A character vector indicating the column name(s) of the data frame that contain the additional columns you wish appended to the output
#' @param duplicate A logical flag that indicates whether to duplicate observations with more than one interaction. Default is TRUE to duplicate all records that interacted with more than one individual (i.e. a fishing set that caught two of the same species)
#' 
#' @description A data.frame with a new first column of var as a binary factor (duplicated if duplicate=TRUE), the header and covariates columns, and a random variable column
#' 
#' @export
#' 
erf_data_prep <- function(df, var, covariates, header, duplicate=TRUE){
	v <- cbind(df[,var], df[,c(header,covariates)])
	keep <- apply(v,2,function(x)all(!is.na(x)))
	v <- v[keep]
	colnames(v)[1] <- var #change the name of the first column to the var name
	# Handle whether to duplicate or not
		if(duplicate){
			dup.id <- which(v[,var] > 1)
			dup.rows <- v[rep.int(dup.id,v[dup.id,var]-1),]
			v <- rbind(v, dup.rows)
		}
	# Convert interactions to factor
		v[,var] <- factor(as.integer(v[,var]>0), levels=c(0,1))
	# Add a random variable
		v$random <- rnorm(nrow(v), 0, 1)
	return(v)
}