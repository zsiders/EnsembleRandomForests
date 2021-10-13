#' Plots Accumulated Local Effects for ERF
#' 
#' @description Plots the Accumulated Local Effects (ALE) from an ERF object
#'
#' @param ALE an ALE object for a given variable
#' @param name the name of the variable
#' @param cex.axis the cex of the axis
#' @param cex.lab the cex of the labels
#' 
#' @return A plot of the ALEs
#' 
#' @export
#' 
#' @examples
#' #run an ERF with 10 RFs and 
#' ens_rf_ex <- ens_random_forests(df=simData$samples, var="obs", covariates=grep("cov", colnames(simData$samples),value=T), save=FALSE, cores=1)
#' 
#' ALEdf <- ALE_fn(ens_rf_ex, save=FALSE)
#' 
#' plot_ALE(ALEdf[1])
#' 
plot_ALE <- function(ALE, name, cex.axis=1, cex.lab=1){
	if(missing(name) & length(ALE)==1) name <- names(ALE)
	if(length(ALE)==1) ALE <- ALE[[1]]
	if(class(ALE[,1])=='character'){
		ALE[,2:ncol(ALE)] <- sapply(ALE[,2:ncol(ALE)],as.numeric)
		y.range <- range(unlist(ALE[,2:ncol(ALE)]))
		x.seq <- c(0,nrow(ALE)+1)
		x.range <- range(x.seq)
		quant <- t(apply(ALE[,2:ncol(ALE)], 1, quantile, probs=c(0.1, 0.5, 0.9)))
		plot(1:nrow(ALE), quant[,2], 
	     type='n', 
	     xlim= range(pretty(x.range)), 
	     ylim = range(pretty(range(y.range))), 
	     las=1, xaxt='n',
	     xlab=parse(text=name), 
	     ylab="", 
	     cex.axis = cex.axis, 
	     cex.lab = cex.lab)
		axis(1, at=1:nrow(ALE), labels=ALE[,1], cex.axis=cex.axis)
		segments(x0 = 1:nrow(ALE),
		         x1 = 1:nrow(ALE), 
		        y0 = quant[,1],
		        y1=quant[,3], 
		        col = "gray30",
		        lwd=3)
		points(1:nrow(ALE),quant[,2], pch=16, cex=2)
		abline(h=0, lty=3)
	}else{
		x.range <- range(ALE[,1])
		y.range <- range(unlist(ALE[,2:ncol(ALE)]))
		quant <- t(apply(ALE[,2:ncol(ALE)], 1, quantile, probs=c(0.1, 0.5, 0.9)))
		plot(ALE[,1], quant[,2], 
	     type='n', 
	     xlim= range(pretty(x.range)), 
	     ylim = range(pretty(range(y.range))), 
	     las=1, 
	     xlab=parse(text=name), 
	     ylab="", 
	     cex.axis = cex.axis, 
	     cex.lab = cex.lab)
		polygon(x = c(ALE[,1], rev(ALE[,1])), 
		        y = c(quant[,1], rev(quant[,3])), 
		        border = NA, 
		        col = "gray80")
		lines(ALE[,1], quant[,2], lwd=3)
		abline(h=0, lty=3)
	}
}