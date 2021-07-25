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
#' rocr$auc
#' str(rocr)
#' 
#' 
rocr_ens <- function(predictions, samples){
	#This function generates Receiving Operator Characteristic Cruves and model metrics
	#make predictions from ROCR
	pred <- ROCR::prediction(predictions, samples)
	# TPR/FPR
	tpr <- ROCR::performance(pred, "tpr")
	# TNR/FNR
	tnr <- ROCR::performance(pred, "tnr")
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
	tss <- max(tpr@y.values[[1]] + tnr@y.values[[1]]-1)
	rmse <- sqrt(sum((predictions-as.integer(as.character(samples)))^2)/length(predictions))
	
	return(list(tpr=tpr, tnr=tnr, auc=auc, phi=phi, 
	            acc=acc, err=err, fpr=fpr, 
	            sens=sens, spec=spec, tss=tss, rmse=rmse))
}