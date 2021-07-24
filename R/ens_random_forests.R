#' Ensemble Random Forests
#' 
#' @description Execute Ensemble Random Forests model on a given dataset
#' 
#' @inheritParams erf_data_prep
#' @param out.folder A path to the folder to write out too. If NULL then a folder is generated in the working directory
#' @param duplicate A logical flag that indicates whether to duplicate observations with more than one interaction. Default is TRUE to duplicate all records that interacted with more than one individual (i.e. a fishing set that caught two of the same species)
#' @param n.forests A numeric value indicating how many Random Forests to generate in the ensemble, default is 100
#' @param cores A numeric value that either indicates the number of cores to use for parallel processing or a negative value to indicate the number of cores to leave free. Default is to leave two cores free.
#' @param save A logical flag to save the output as an RData object, default is TRUE.
#' 
#' @return A fitted ERF model as a list object. Contains the model data, the fitted RF models, the ensemble predictions, the ensemble training/test performance, the ROC training/test performance for each RF, and the individual RF predictions.
#' 
#' @export
#' 
ens_random_forests <- function(df, var, covariates, out.folder=NULL, duplicate=TRUE, 
	                n.forests=100, cores=parallel::detectCores()-2, save=TRUE){
	#Prep
		if(is.null(out.folder)){
			dir.create('Output/')
			out.folder <- "Output/"
		}
		form <- erf_formula_prep(var, covariates) #Prepare the model formula
		v <- erf_data_prep(df, var, covariates, header, duplicate=duplicate)
		max_split <- max_splitter(v)
	
	# Ensemble Random Forest
		UseCores <- ifelse(n.forests < cores, n.forests, cores) #use less cores if n.forests < cores
		UseCores <- ifelse(max_split < 7e3, pmin(UseCores,floor(cores * (1-max_split/7e3))), 2) #use less cores if # of max_split exceeds a certain value (need some RAM available for calcs)

		cl <- makeCluster(UseCores) #make clusters
		registerDoParallel(cl) #designate cores
		rf.ens <- foreach(i=1:n.forests, .packages=c('randomForest','ROCR')) %dopar%{
			rf.ens.fn(v, var, form, max_split)
		}
		stopCluster(cl)
		
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

	# Generate ensemble predictions
		pred_ens <- data.frame("P.0"=rowMeans(pred_ens_z), 
		                       "P.1"=rowMeans(pred_ens_p))
		pred_ens$PRES <- v[,var]
		pred_ens$resid <- as.integer(as.character(v[,1])) - pred_ens[,2]
		rownames(pred_ens) <- rownames(v)

	# Package up the output
		pack <- list(data = v, 
		             model = rf.ens, 
		             ens.pred = pred_ens,
		             ens.tr.perf = c(trAUC=mean(pred_ens_trAUC),
		                             trRMSE=mean(pred_ens_trRMSE),
		                             trTSS=mean(pred_ens_trTSS)),
		             ens.te.perf = c(teAUC=mean(pred_ens_teAUC),
		                             teRMSE=mean(pred_ens_teRMSE),
		                             teTSS=mean(pred_ens_teTSS)),
		             roc_train = lapply(rf.ens, function(x)x$roc_train),
		             roc_test = lapply(rf.ens, function(x)x$roc_test),
		             pred = list(p = pred_ens_p,
		                         resid = pred_ens_resid))
		if(save) save(pack, file=paste0(out.folder,"/ERF_",var,".Rdata"))

		return(pack)
}
	