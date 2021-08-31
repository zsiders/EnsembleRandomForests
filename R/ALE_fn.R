#' Accumulated Local Effects for ERF
#' 
#' @description Calculates the Accumulated Local Effects (ALE) from an ERF object
#'
#' @param fit The fitted object returned from calling ens_random_forests()
#' @param var The name of the response variable
#' @param save A logical flag to save the output as an RData object, default is TRUE.
#' @param out.folder A path to the folder to write out too. If NULL then a folder is generated in the working directory
#' @param cores An integer value that either indicates the number of cores to use for parallel processing or a negative value to indicate the number of cores to leave free. Default is to leave two cores free.
#' 
#' @return A list that contains the ALEs for each variable, ordered by the mean variable importance
#' 
#' @export
#' 
#' @examples
#' #run an ERF with 10 RFs and 
#' ens_rf_ex <- ens_random_forests(df=simData$samples, var="obs", covariates=grep("cov", colnames(simData$samples),value=T), save=FALSE, cores=1)
#' 
#' ALEdf <- ALE_fn(ens_rf_ex, save=FALSE)
#' head(ALEdf[[1]])
#' 
ALE_fn <- function(fit, var, save=TRUE, out.folder=NULL, cores=parallel::detectCores()-4){
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
		vi  <- rownames(fit$var.imp[fit$var.imp$ord])
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

	ALEdf <- foreach(i = 2:ncol(data.df), .packages=c('randomForest'), .export=c("ALEfun.J1","yhat")) %dopar% {
		ex <- lapply(model, function(x) {ALEfun.J1(data.df, x$mod, yhat, J = i, K=50)})

		df <- as.data.frame(cbind(ex[[1]]$x.values, sapply(ex,function(x) {x$f.values})))
		names(df)[1] <- "x"
		names(df)[2:ncol(df)] <- paste0("f.",1:length(ex))
		return(df)
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
yhat <- function(X.model, newdata) as.numeric(predict(X.model, newdata))
ALEfun.J1 <- function (X, X.model, pred.fun, J, K = 40){
    N = dim(X)[1]
    d = dim(X)[2]
    if (class(X[, J]) == "factor") {
        X[, J] <- droplevels(X[, J])
        x.count <- as.numeric(table(X[, J]))
        x.prob <- x.count/sum(x.count)
        K <- nlevels(X[, J])
        D.cum <- matrix(0, K, K)
        D <- matrix(0, K, K)
        for (j in setdiff(1:d, J)) {
            if (class(X[, j]) == "factor") {
              A = table(X[, J], X[, j])
              A = A/x.count
              for (i in 1:(K - 1)) {
                for (k in (i + 1):K) {
                  D[i, k] = sum(abs(A[i, ] - A[k, ]))/2
                  D[k, i] = D[i, k]
                }
              }
              D.cum <- D.cum + D
            }
            else {
              q.x.all <- quantile(X[, j], probs = seq(0, 
                1, length.out = 100), na.rm = TRUE, names = FALSE)
              x.ecdf = tapply(X[, j], X[, J], ecdf)
              for (i in 1:(K - 1)) {
                for (k in (i + 1):K) {
                  D[i, k] = max(abs(x.ecdf[[i]](q.x.all) - 
                    x.ecdf[[k]](q.x.all)))
                  D[k, i] = D[i, k]
                }
              }
              D.cum <- D.cum + D
            }
        }
        D1D <- cmdscale(D.cum, k = 1)
        ind.ord <- sort(D1D, index.return = T)$ix
        ord.ind <- sort(ind.ord, index.return = T)$ix
        levs.orig <- levels(X[, J])
        levs.ord <- levs.orig[ind.ord]
        x.ord <- ord.ind[as.numeric(X[, J])]
        row.ind.plus <- (1:N)[x.ord < K]
        row.ind.neg <- (1:N)[x.ord > 1]
        X.plus <- X
        X.neg <- X
        X.plus[row.ind.plus, J] <- levs.ord[x.ord[row.ind.plus] + 
            1]
        X.neg[row.ind.neg, J] <- levs.ord[x.ord[row.ind.neg] - 
            1]
        y.hat <- pred.fun(X.model = X.model, newdata = X)
        y.hat.plus <- pred.fun(X.model = X.model, newdata = X.plus[row.ind.plus, 
            ])
        y.hat.neg <- pred.fun(X.model = X.model, newdata = X.neg[row.ind.neg, 
            ])
        Delta.plus <- y.hat.plus - y.hat[row.ind.plus]
        Delta.neg <- y.hat[row.ind.neg] - y.hat.neg
        Delta <- as.numeric(tapply(c(Delta.plus, Delta.neg), 
            c(x.ord[row.ind.plus], x.ord[row.ind.neg] - 1), 
            mean))
        fJ <- c(0, cumsum(Delta))
        fJ = fJ - sum(fJ * x.prob[ind.ord])
        x <- levs.ord
    }
    else if (class(X[, J]) == "numeric" | class(X[, J]) == 
        "integer") {
        z = c(min(X[, J]), as.numeric(quantile(X[, J], seq(1/K, 
            1, length.out = K), type = 1)))
        z = unique(z)
        K = length(z) - 1
        fJ = numeric(K)
        a1 = as.numeric(cut(X[, J], breaks = z, include.lowest = TRUE))
        X1 = X
        X2 = X
        X1[, J] = z[a1]
        X2[, J] = z[a1 + 1]
        y.hat1 = pred.fun(X.model = X.model, newdata = X1)
        y.hat2 = pred.fun(X.model = X.model, newdata = X2)
        Delta = y.hat2 - y.hat1
        Delta = as.numeric(tapply(Delta, a1, mean))
        fJ = c(0, cumsum(Delta))
        b1 <- as.numeric(table(a1))
        fJ = fJ - sum((fJ[1:K] + fJ[2:(K + 1)])/2 * b1)/sum(b1)
        x <- z
    }
    list(K = K, x.values = x, f.values = fJ)
}