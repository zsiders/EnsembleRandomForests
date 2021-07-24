#' Random Forest model with additional bagging
#' 
#' @description Runs a single Random Forest model with an additional bagging layer and calculates performance metrics
#' 
#' @param v A data frame object created by `erf_data_prep()` or internally in `erf()`
#' @param var A character string specifying the column name to use as the dependent variable
#' @param form A formula class object specifying the RF model formulation (created by `erf_formula_prep()` or internal in `erf()`)
#' @param min_split The minimum number of samples in the RF bagging procedure (created internally by `erf()`)
#' @param ntree The number of decision trees to use in each RF, default is 1000
#' @param mtry The number of covariates to try at each node split, default is 5
#' @param importance A logical flag for the randomForest model to calculate the variable importance
#' 
#' @return A list containing mod (the randomForest model), preds (the predictions), roc_train (the Receiver Operator Characteristic Curve performance metrics calculated by rocr_ens() on the training set), roc_test (the Receiver Operator Characteristic Curve performance metrics calculated by rocr_ens() on the test set) 
#' @export
#' 
#' @examples
#' form <- erf_formula_prep('obs', grep('cov',colnames(simData$samples),value=TRUE))
#' data <- erf_data_prep(simData$samples, 'obs', grep('cov', colnames(simData$samples), value=TRUE))
#' max_split <- max_splitter(data)
#' 
#' #fit a single RandomForest
#' rf_ex <- rf_ens_fn(data, 'obs', form, max_split, ntree=200)
#' 
#' #see the training/test auc value
#' rf_ex$roc_train$auc
#' rf_ex$roc_test$auc
#' 
#' #see the distribution of predictions
#' plot(density(rf_ex$preds[,2],from=0,to=1,adj=2))
#' 
rf_ens_fn <- function(v, var, form, max_split, ntree=1000, mtry=5, importance=TRUE){
	#first bagging
	zero_sub_ens <- sample(1:nrow(v[v[,var]=="0",]),
	                       floor(0.9*nrow(v[v[,var]=="0",])))
	one_sub_ens <- sample(1:nrow(v[v[,var]!="0",]),
	                      floor(0.9*nrow(v[v[,var]!="0",])))

	#build a train and test set
	train_ens <- rbind(v[v[,var]=="0",][zero_sub_ens,], v[v[,var]=="1",][one_sub_ens,])
	test_ens <- rbind(v[v[,var]=="0",][-zero_sub_ens,], v[v[,var]=="1",][-one_sub_ens,])

	#run the RandomForests
	#second bagging internal in RF
	mod <- randomForest::randomForest(form, 
	                    data=train_ens, 
	                    ntree=1000, 
	                    mtry=5, 
	                    importance=TRUE, 
	                    sampsize=c(max_split,max_split))

	#predictions
	preds <- as.data.frame(randomForest::predict(mod, 
	                               newdata=v, 
	                               type='prob'))
	colnames(preds) <- paste0('P.',colnames(preds))
	preds$PRES <- v[,dv]
	preds$type <- 'train'
	preds$SetID <- v$SetID
	preds$type[preds$SetID %in%test_ens$SetID] <- 'test'

	roc_train <- rocr_ens(preds[preds$type=='train',2], preds$PRES[preds$type=='train'])
	roc_test <- rocr_ens(preds[preds$type=='test',2], preds$PRES[preds$type=='test'])

	pack <- list(mod = mod, 
	             preds = preds, 
	             roc_train = roc_train, 
	             roc_test = roc_test)
	return(pack)
}