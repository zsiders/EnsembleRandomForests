# library(EnsembleRandomForests)

# ens_rf_ex <- ens_random_forests(df=simData$samples, var="obs",
#                                 covariates=grep("cov",colnames(simData$samples),value=T),
#                                 header = NULL,
#                                 save=FALSE,
#                                 out.folder=NULL,
#                                 duplicate = TRUE,
#                                 n.forests = 10L,
#                                 importance = TRUE,
#                                 ntree = 1000,
#                                 mtry = 5,
#                                 var.q = c(0.1,0.5,0.9),
#                                 cores = parallel::detectCores()-2)
# drawFun <- function(data, holdout=0.1, return.ho=TRUE){
#     #data is data.frame
#     #holdout is proportion to holdout
#     #return.ho is logical for returning the holdout or the rest
#     nho <- floor(nrow(data)*holdout)
#     samps <- sample.int(nrow(data),nho)
#     if(return.ho){
#         return(data[samps,])
#     }else{
#         return(data[-samps,])
#     }
# }

# #first value is number of holdout datasets to make
# newdat <- replicate(5, erf_data_prep(df = drawFun(data=simData$samples), var = "obs",covariates = grep('cov',colnames(simData$samples), value=T),header = c('prob.raw','prob'), duplicate = TRUE),simplify=FALSE)

# ho.pred <- lapply(newdat, function(dat)sapply(ens_rf_ex$model,function(mod.l) predict(mod.l$mod,dat,type='prob')[,2]))

# ho.ens.pred <- sapply(ho.pred, rowMeans)