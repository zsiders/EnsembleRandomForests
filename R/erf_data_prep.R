#' Prepare a data.frame for ERF
#' 
#' @description Prepares a data.frame object for use with the ERF (not needed if using erf())
#' 
#' @param df A data.frame object
#' @param var A character string indicating the column name of the data frame that contains the number of interactions for the ERF to model; column should be a numeric column
#' @param covariates A character vector indicating the column name(s) of the data frame that contain the covariates
#' @param header A character vector indicating the column name(s) of the data frame that contain the additional columns you wish appended to the output
#' @param duplicate A logical flag that indicates whether to duplicate observations with more than one interaction. Default is TRUE to duplicate all records that interacted with more than one individual (i.e. a fishing set that caught two of the same species)
#' @param weights a vector equal in length to nrow(df) of weights, NULL by default
#' 
#' @description A data.frame with a new first column of var as a binary factor (duplicated if duplicate=TRUE), the header and covariates columns, and a random variable column
#' 
#' @export
#' 
#' @examples
#' data <- erf_data_prep(df = simData$samples, var = 'obs', covariates = grep('cov', colnames(simData$samples), value=TRUE), header = c('prob.raw','prob'))
#' head(data)
#' 
erf_data_prep <- function(df=NULL, var=NULL, covariates=NULL, header=NULL, weights = NULL, duplicate=TRUE, mode='bin'){
	if(is.null(df)) stop("Need to supply data.frame")
	if(is.null(var)) stop("Need to supply variable column")
	if(is.null(covariates)) stop("Need to supply covariate columns")
	if(is.null(header)){
		if(is.null(weights)){
			v <- cbind(df[,var], df[,c(covariates)])
		}else{
			v <- cbind(df[,var], df[,c(covariates)], weights = weights)
		}		
	}else{
		if(is.null(weights)){
			v <- cbind(df[,var], df[,c(header,covariates)])
		}else{
			v <- cbind(df[,var], df[,c(header,covariates)], weights = weights)
		}
	}
	
	keep <- apply(v,2,function(x)all(is.na(x)))
	v <- v[,!keep]
	keep <- apply(v,1,function(x)all(!is.na(x)))
	v <- v[keep,]
	colnames(v)[1] <- var #change the name of the first column to the var name
	# Handle whether to duplicate or not
		if(duplicate){
			if(any(v[,var]>1)){
				dup.id <- which(v[,var] > 1)
				dup.rows <- v[rep.int(dup.id,v[dup.id,var]-1),]
				v <- rbind(v, dup.rows)
			}
		}
	# Convert interactions to factor
		if(mode=='bin'){
			v[,var] <- factor(as.integer(v[,var]>0), levels=c(0,1))
		}else if(mode=='multi'){
			if(!is(v[,var],'factor')) v[,var] <- factor(v[,var])
		}
	# Add a random variable
		v$random <- rnorm(nrow(v), 0, 1)
	return(v)
}