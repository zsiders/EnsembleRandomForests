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
#' @param ntree The number of decision trees to use in each RF, default is 1000
#' @param mtry The number of covariates to try at each node split, default is 5
#' 
#' @return A fitted ERF model as a list object. Contains the model data, the fitted RF models, the ensemble predictions, the ensemble model performance, the mean training/test performance, the ROC training/test performance for each RF, and the individual RF predictions.
#' 
#' @export
#' 
#' @examples
#' #run an ERF with 10 RFs and 
#' ens_rf_ex <- ens_random_forests(df=simData$samples, var="obs", covariates=grep("cov",colnames(simData$samples),value=T), save=FALSE, cores=1)
#' head(ens_rf_ex$data) # view the dataset used in the model
#' head(ens_rf_ex$ens.pred) #view the model predictions
#' ens_rf_ex$mu.te.perf #view the mean test threshold-free performance metrics
#' unlist(trial$ens.perf[c('auc','tss','rmse')]) #view the threshold-free ensemble performance metrics
#' 
ens_random_forests <- function(df, var, covariates, header=NULL, out.folder=NULL, duplicate=TRUE, n.forests=10, cores=parallel::detectCores()-2, save=TRUE, ntree=1000, mtry=5){
	#Prep
		if(missing(df)) stop("Supply a data.frame")
		if(missing(var)) stop("Supply a variable to model")
		if(missing(covariates)) stop("Supply covariates to use")
		if(is.null(out.folder) & save==TRUE){
			dir.create('Output/')
			out.folder <- "Output/"
			warning("No output folder provided; creating in working directory")
		}
		form <- erf_formula_prep(var, covariates) #Prepare the model formula
		if(is.null(header)){
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
				rf_ens_fn(v, form, max_split, ntree=ntree, mtry=mtry)
			}
			stopCluster(cl)
		}else{
			rf.ens <- lapply(1:n.forests, function(i) rf_ens_fn(v, form, max_split, ntree=ntree, mtry=mtry))
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
		if(save) save(pack, file=paste0(out.folder,"/ERF_",var,".Rdata"))

	return(pack)
}
	