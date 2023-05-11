#' Accumulated Local Effects for ERF
#' 
#' @description Calculates the Accumulated Local Effects (ALE) from an ERF object
#'
#' @param fit The fitted object returned from calling ens_random_forests()
#' @param var The name of the response variable
#' @param save A logical flag to save the output as an RData object, default is TRUE.
#' @param out.folder A path to the folder to write out too. If NULL then a folder is generated in the working directory
#' @param cores An integer value that either indicates the number of cores to use for parallel processing or a negative value to indicate the number of cores to leave free. Default is to leave two cores free.
#' @param type is either 'response' or 'prob' from predict.randomForest; if 'prob' then n sets of predictions are returned for the n levels in var; if "response" then the factorized predicted response values are returned
#' 
#' @return A list that contains a data.frame for each variable, ordered by the mean variable importance, and a vector of the covariate values (used for rug plot in plot_ALE). The columns in each data.frame are as follows:
#' \itemize{
#'      \item \strong{x}: the covariate values that the ALE was calculated for
#'      \item \strong{class}: the class of the covariate; used by subsequent plot_ALE function
#'      \item \strong{q}: the quantile of the x value of the covariate
#'      \item \strong{f.X}: the ALEs evaluated at a given x value
#' }
#' 
#' @export
#' 
#' @examples
#' #run an ERF with 10 RFs and 
#' ens_rf_ex <- ens_random_forests(df=simData$samples, var="obs", covariates=grep("cov", colnames(simData$samples),value=T), save=FALSE, cores=1)
#' 
#' ALEdf <- calc_ALE(ens_rf_ex, save=FALSE)
#' head(ALEdf[[1]]$df)
#' 
calc_ALE <- function(fit, var, save=TRUE, out.folder=NULL, cores=parallel::detectCores()-4, type='response'){
	if(missing(fit)) stop("Supply fit object")
	if(missing(var)){
        message("No name of response variable, making one")
        var <- "var"
    }
	if(!is.integer(cores)){
        message("rounding n.forests to the nearest one")
        cores <- floor(cores)
    } 

	cv.roc <- sd(sapply(fit$model,function(x)x$roc_test$auc))/mean(sapply(fit$model,function(x)x$roc_test$auc))
	#vi
	if(is.null(fit$var.imp)){
		message("No variable importance calculated in fit object >> running out of order")
		vi <- attr(fit$model[[1]]$mod$terms,'term.labels')
	}else{
		vi  <- rownames(fit$var.imp)[fit$var.imp$ord]
	}
	vi.ind <- match(vi,colnames(fit$data))


	#ALE wrapper

	ob <- as.numeric(object.size(fit)/1e9)

	if(cv.roc < 0.01){
		samp <- sort(sample(seq(1,length(fit$model)), 5))
	}else{
		nsamp <- pmin(length(fit$model),ceiling(exp(-2.041 + -1.068*log(ob))/0.01*5))
		samp <- sort(sample(seq(1,length(fit$model)), nsamp))
	}

	model <- fit$model[samp]
	data.df <- fit$data[c(1,vi.ind)]

	ob <- as.numeric(object.size(model)/1e9)

	UseCores <- floor(cores*as.numeric(as.character(cut(ob,breaks=c(Inf,20,15,10,6,0),labels=rev(c(0,0.125,0.25,0.5,0.8))))))

	cl <- makeCluster(UseCores)

	registerDoParallel(cl)

	ALEdf <- foreach(i = 2:ncol(data.df), .packages=c('randomForest'), .export=c("ALE_fn","yhat")) %dopar% {
		ex <- lapply(model, function(x) {ALE_fn(data.df, x$mod, yhat, J = i, K=50, type=type)})

		df <- data.frame(x = ex[[1]]$x.values, 
                         class = ex[[1]]$class, 
                         q = as.numeric(ex[[1]]$quantile), 
                         sapply(ex,function(x) {as.numeric(x$f.values)}))
		names(df)[4:ncol(df)] <- paste0("f.",1:length(ex))
		return(list(df=df, X=data.df[,i]))
	}
	stopCluster(cl)

	names(ALEdf) <- vi
	if(save){
		if(is.null(out.folder)){
			dir.create('Output/')
			out.folder <- "Output/"
			warning("No output folder provided; creating in working directory")
		}
		save(ALEdf, file=paste0(out.folder,"/ALE_",var,".Rdata"))
	}
	return(ALEdf)
}
yhat <- function(X.model, newdata, type='response'){
	if(type=='response'){
		as.numeric(predict(X.model, newdata, type=type))
	}else{
		predict(X.model, newdata, type=type)
	}
}