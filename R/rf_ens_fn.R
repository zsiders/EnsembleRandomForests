#' Random Forest model with additional bagging
#' 
#' @description Runs a single Random Forest model with an additional bagging layer and calculates performance metrics
#' 
#' @param v A data frame object created by `erf_data_prep()` or internally in `ens_random_forests()`
#' @param form A formula class object specifying the RF model formulation (created by `erf_formula_prep()` or internal in `ens_random_forests()`)
#' @param max_split The maximum number of samples in the RF bagging procedure (created internally by `ens_random_forests()`)
#' @param ntree The number of decision trees to use in each RF, default is 100
#' @param mtry The number of covariates to try at each node split, default is 5
#' @param importance A logical flag for the randomForest model to calculate the variable importance
#' 
#' @return A list containing mod (the randomForest model), preds (the predictions), roc_train (the Receiver Operator Characteristic Curve performance metrics calculated by rocr_ens() on the training set), roc_test (the Receiver Operator Characteristic Curve performance metrics calculated by rocr_ens() on the test set) 
#' @export
#' 
#' @examples
#' form <- erf_formula_prep(var='obs', covariates=grep('cov',colnames(simData$samples),value=TRUE))
#' data <- erf_data_prep(df=simData$samples, var='obs', covariate=grep('cov', colnames(simData$samples), value=TRUE))
#' max_split <- max_splitter(data)
#' 
#' #fit a single RandomForest
#' rf_ex <- rf_ens_fn(v=data, form=form, max_split=max_split, ntree=50)
#' 
#' #see the training/test auc value
#' rf_ex$roc_train$auc
#' rf_ex$roc_test$auc
#' 
#' #see the distribution of predictions
#' par(mar=c(4,4,1,1))
#' plot(density(rf_ex$preds[,2],from=0,to=1,adj=2), main="", las=1)
#' 
rf_ens_fn <- function(v, form, max_split, ntree=100, mtry=5, importance=TRUE){
	if(missing(v)) stop("Supply a ERF data.frame")
	if(missing(form)) stop("Supply a model formula")
	var <- as.character(form)[2]
	if(missing(max_split)){
		max_split <- max_splitter(v)
		warning("max_split not provided, assuming default")
	}
	v$rf.ID <- 1:nrow(v)
	ncovar <- length(strsplit(gsub("[[:punct:]]+\\s","",as.character(form)[3]),"\\s")[[1]])
	if(mtry >= ncovar) mtry <- ncovar-1 #permutation
	#first bagging
	# if(mode=='bin'){
	# 	zero_sub_ens <- sample(1:nrow(v[v[,var]=="0",]),
	#                        floor(0.9*nrow(v[v[,var]=="0",])))
	# 	one_sub_ens <- sample(1:nrow(v[v[,var]!="0",]),
	# 	                      floor(0.9*nrow(v[v[,var]!="0",])))

	# 	#build a train and test set
	# 	train_ens <- rbind(v[v[,var]=="0",][zero_sub_ens,], v[v[,var]=="1",][one_sub_ens,])
	# 	test_ens <- rbind(v[v[,var]=="0",][-zero_sub_ens,], v[v[,var]=="1",][-one_sub_ens,])
	# }else{

		sub_ens <- do.call('c',sapply(levels(v[,var]),
		                  function(x){
		                  	n <- nrow(v[v[,var]==x,]);
		                  sample(v$rf.ID[v[,var]==x],
	                       floor(0.9*n))
		                  }))
		train_ens <- v[v$rf.ID %in% sub_ens,]
		test_ens <- v[!v$rf.ID %in% sub_ens,]
	# }
	

	#run the RandomForests
	#second bagging internal in RF
	mod <- randomForest(form, 
	                    data=train_ens, 
	                    ntree=ntree, 
	                    mtry=mtry, 
	                    importance=TRUE, 
	                    sampsize=rep(max_split,
	                                 nlevels(v[,var])))

	#predictions
	preds <- as.data.frame(predict(mod, 
	                               newdata=v, 
	                               type='prob'))
	colnames(preds) <- paste0('P.',colnames(preds))
	preds$PRES <- v[,var]
	preds$type <- 'train'
	preds$rf.ID <- 1:nrow(v)
	preds$type[preds$rf.ID %in% test_ens$rf.ID] <- 'test'
	preds <- preds[,-grep('rf.ID',colnames(preds))]
	if(nlevels(v[,var])==2){
		roc_train <- rocr_ens(preds[preds$type=='train',2], preds$PRES[preds$type=='train'])
		roc_test <- rocr_ens(preds[preds$type=='test',2], preds$PRES[preds$type=='test'])
	}else{
		roc_train <- lapply(1:nlevels(v[,var]),function(x) rocr_ens(preds[preds$type=='train',x], as.integer(preds$PRES[preds$type=='train']==levels(v[,var])[x])))
		roc_test <- lapply(1:nlevels(v[,var]),function(x) rocr_ens(preds[preds$type=='test',x], as.integer(preds$PRES[preds$type=='test']==levels(v[,var])[x])))
	}
	

	pack <- list(mod = mod, 
	             preds = preds, 
	             roc_train = roc_train, 
	             roc_test = roc_test)
	return(pack)
}