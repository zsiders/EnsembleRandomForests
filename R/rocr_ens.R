#' Receiver Operator Characteristic Curve calculations
#' 
#' @description A function to calculate various ROC curve performance metrics given model predictions and the oberved samples
#' 
#' @param predictions Random Forest predictions (or any probability)
#' @param samples The original observations (must be a factor)
#' 
#' @return A list containing all the various ROC curve performance metrics
#' @export
#' 
#' @examples
#' predictions <- rbeta(100, 4, 4)
#' samples <- rbinom(100, 1, 0.5)
#' rocr <- rocr_ens(predictions, samples)
#' str(rocr)
#' rocr$auc
#' 
rocr_ens <- function(predictions, samples){
	#This function generates Receiving Operator Characteristic Cruves and model metrics
	#make predictions from ROCR
	pred <- ROCR::prediction(predictions, samples)
	# TPR/FPR
	perf <- ROCR::performance(pred, "tpr", "fpr")
	# TNR/FNR
	nperf <- ROCR::performance(pred, "tnr", "fnr")
	# SENS/SPEC
	ss <- ROCR::performance(pred, "sens", "spec")
	#PPV/NPV
	ppv <- ROCR::performance(pred, "ppv", "npv")
	#1-SENSITIVITY
	ss@x.values[[1]] <- 1-ss@x.values[[1]]
	#PHI
	phi <- ROCR::performance(pred, "phi")
	#ACCURACY
	acc <- ROCR::performance(pred, "acc")
	#ERROR
	err <- ROCR::performance(pred, "err")
	#FPR
	fpr <- ROCR::performance(pred, "fpr")
	#SENSITIVITY
	sens <- ROCR::performance(pred, "sens")
	#SPECIFICITY
	spec <- ROCR::performance(pred, "spec")
	#AUC
	auc <- ROCR::performance(pred, measure='auc')@y.values[[1]]
	
	return(list(ss=ss, auc=auc, phi=phi, 
	            acc=acc, err=err, fpr=fpr, 
	            sens=sens, spec=spec))
}