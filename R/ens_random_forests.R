#' Ensemble Random Forests
#' 
#' @description Execute Ensemble Random Forests model on a given dataset
#' 
#' @inheritParams erf_data_prep
#' @param out.folder A path to the folder to write out too. If NULL then a folder is generated in the working directory
#' @param duplicate A logical flag that indicates whether to duplicate observations with more than one interaction. Default is TRUE to duplicate all records that interacted with more than one individual (i.e. a fishing set that caught two of the same species)
#' @param n.forests An integer value indicating how many Random Forests to generate in the ensemble, default is 100
#' @param importance A logical flag for the randomForest model to calculate the variable importance
#' @param cores A integer value that either indicates the number of cores to use for parallel processing or a negative value to indicate the number of cores to leave free. Default is to leave two cores free.
#' @param save A logical flag to save the output as an RData object, default is TRUE.
#' @param ntree The number of decision trees to use in each RF, default is 1000
#' @param mtry The number of covariates to try at each node split, default is 5
#' @param var.q The quantiles for the distribution of the variable importance; only exectuted if importance=TRUE
#' 
#' @return A list containing the fitted ERF model and associated output. 
#' \itemize{
#' 	\item \strong{data}: the exact dataset used to fit each Random Forests within the ensemble. At the minimum, it will contain \code{var}, \code{covariates}, and \code{header} if provided. If \code{duplicate}==\code{TRUE}, then it will also contain the duplicated presence records. 
#' \item \strong{model}: the returned fitted individual Random Forests in the ensemble
#' \item \strong{ens.pred}: the ensemble model predictions. This is generally the prediction set to use.
#' \item \strong{ens.perf}: the ensemble model performance metrics including ROC curve.
#' \item \strong{mu.tr.perf}: the mean training set performance across all Random Forests for the AUC, TSS, and RMSE metrics. These are generally useless.
#' \item \strong{mu.te.perf}: the mean test set performance across all Random Forests for the AUC, TSS, and RMSE metrics. These can be informative.
#' \item \strong{roc_train}: the training performance metrics for each Random Forests in the ensemble
#' \item \strong{roc_test}: the test performance metrics for each Random Forests in the ensemble
#' \item \strong{pred}: a list with two objects: 
#' \itemize{
#' 		\item \strong{p}: predictions from each Random Forests in the ensemble to the dataset
#' 		\item \strong{resid}: residuals between the observed presence/absences and the predictions from each Random Forests in the ensemble to the dataset
#' 	}
#' }
#' 
#' @export
#' 
#' @examples
#' #run an ERF with 10 RFs and 
#' ens_rf_ex <- ens_random_forests(df=simData$samples, var="obs",covariates=grep("cov", colnames(simData$samples),value=T), save=FALSE, cores=1)
#' 
#' # view the dataset used in the model
#' head(ens_rf_ex$data) 
#' 
#' #view the model predictions
#' head(ens_rf_ex$ens.pred) 
#' 
#' #view the mean test threshold-free performance metrics
#' ens_rf_ex$mu.te.perf 
#' 
#' #view the threshold-free ensemble performance metrics
#' unlist(ens_rf_ex$ens.perf[c('auc','rmse','tss')]) 
#' 
ens_random_forests <- function(df, var, covariates, header=NULL, out.folder=NULL, duplicate=TRUE, n.forests=10L, importance=TRUE, cores=parallel::detectCores()-2, save=TRUE, ntree=1000, mtry=5, var.q = c(0.1,0.5,0.9)){
	#Prep
		if(missing(df)) stop("Supply a data.frame")
		if(missing(var)) stop("Supply a variable to model")
		if(missing(covariates)) stop("Supply covariates to use")
		if(!is.integer(n.forests)){
			n.forests <- floor(n.forests)
			message("rounding n.forests to the nearest one")
		}
		if(!is.integer(cores)){
			message("rounding n.forests to the nearest one")
			cores <- floor(cores)
		} 
		
		form <- erf_formula_prep(var, covariates) #Prepare the model formula
		if(!is.null(header)){
			v <- erf_data_prep(df, var, covariates, header, duplicate=duplicate)
		}else{
			v <- erf_data_prep(df, var, covariates, duplicate=duplicate)
		}
		
		max_split <- max_splitter(v)
	
	# Ensemble Random Forest
		if(cores!=1){
			UseCores <- ifelse(n.forests < cores, n.forests, cores) #use less cores if n.forests < cores
			UseCores <- ifelse(max_split < 7e3, pmin(UseCores,floor(cores * (1-max_split/7e3))), 2) #use less cores if # of max_split exceeds a certain value (need some RAM available for calcs)

			cl <- makeCluster(UseCores) #make clusters
			registerDoParallel(cl) #designate cores
			rf.ens <- foreach(i=1:n.forests, .packages=c('randomForest','ROCR')) %dopar%{
				rf_ens_fn(v, form, max_split, 
				          ntree=ntree, mtry=mtry, importance=TRUE)
			}
			stopCluster(cl)
		}else{
			rf.ens <- lapply(1:n.forests, function(i) rf_ens_fn(v, form, max_split, ntree=ntree, mtry=mtry, importance=TRUE))
		}
		
	# Get all the output from the ensemble
		pred_ens_p <- sapply(rf.ens, function(x) {x$preds$P.1}) #NULL line 
		pred_ens_z <- sapply(rf.ens, function(x) {x$preds$P.0})
		pred_ens_trAUC <- sapply(rf.ens, function(x) {x$roc_train$auc})
		pred_ens_teAUC <- sapply(rf.ens, function(x) {x$roc_test$auc})
		pred_ens_trRMSE <- sapply(rf.ens, function(x) {x$roc_train$rmse})
		pred_ens_teRMSE <- sapply(rf.ens, function(x) {x$roc_test$rmse})
		pred_ens_trTSS <- sapply(rf.ens, function(x) {x$roc_train$tss})
		pred_ens_teTSS <- sapply(rf.ens, function(x) {x$roc_test$tss})
		pred_ens_resid <- sapply(rf.ens, function(x) {as.numeric(as.character(x$preds$PRES)) - x$preds$P.1})
	# Get variable importance
		if(importance){
			var_ens <- sapply(rf.ens, function(x) {abs(x$mod$importance[,3]/x$mod$importanceSD[,3])})
			var_ens_df <- as.data.frame(t(apply(var_ens, 1, quantile, prob=var.q)))
			colnames(var_ens_df) <- c(paste0('q_',var.q*100))
			var_ens_df$mu <- rowMeans(var_ens)
			var_ens_df$sd <- apply(var_ens, 1, sd)
			var_ens_df$ord <- order(var_ens_df[,2], decreasing=TRUE)
		}
		

	# Generate ensemble predictions
		pred_ens <- data.frame("P.0"=rowMeans(pred_ens_z), 
		                       "P.1"=rowMeans(pred_ens_p))
		pred_ens$PRES <- v[,var]
		pred_ens$resid <- as.integer(as.character(v[,1])) - pred_ens[,2]
		rownames(pred_ens) <- rownames(v)
	#ROC on ensemble
		roc_ens <- rocr_ens(pred_ens$P.1, pred_ens$PRES)

	# Package up the output
		pack <- list(data = v, 
		             model = rf.ens, 
		             ens.pred = pred_ens,
		             ens.perf = roc_ens,
		             mu.tr.perf = c(trAUC=mean(pred_ens_trAUC),
		                             trRMSE=mean(pred_ens_trRMSE),
		                             trTSS=mean(pred_ens_trTSS)),
		             mu.te.perf = c(teAUC=mean(pred_ens_teAUC),
		                             teRMSE=mean(pred_ens_teRMSE),
		                             teTSS=mean(pred_ens_teTSS)),
		             roc_train = lapply(rf.ens, function(x)x$roc_train),
		             roc_test = lapply(rf.ens, function(x)x$roc_test),
		             pred = list(p = pred_ens_p,
		                         resid = pred_ens_resid))
		if(importance){
			pack$var.imp <- var_ens_df
			pack$var.imp.raw <- var_ens
		}
		if(save){
			if(is.null(out.folder)){
				dir.create('Output/')
				out.folder <- "Output/"
				warning("No output folder provided; creating in working directory")
			}
			save(pack, file=paste0(out.folder,"/ERF_",var,".Rdata"))
		}

	return(pack)
}
	