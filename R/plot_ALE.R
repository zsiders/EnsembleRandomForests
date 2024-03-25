#' Plots Accumulated Local Effects for ERF
#' 
#' @description Plots the Accumulated Local Effects (ALE) from an ERF object
#'
#' @param ALE an ALE object for a given variable; indexed as a list for full functionality
#' @param xquantiles lower and upper quantile bounds to limit x values to
#' @param yquantiles lower, middle, and upper quantiles to plot the confidence interval around the ALE predictions
#' @param name the x-axis label of the variable; self-supplied if indexed as a list
#' @param cex.axis the cex of the axis
#' @param cex.lab the cex of the labels
#' @param rug logical for turning the rug plot on/off (defaults to table values for class=='factor')
#' @param rug.col color of rug ticks
#' @param rug.tick length of rug ticks
#' @param rug.lwd width of rug ticks
#' @param rug.alpha transparency of rug ticks
#' @param rug.max maximum number of ticks to plot (useful for big data)
#' 
#' @return A plot of the ALEs
#' 
#' @export
#' 
#' @examples
#' #run an ERF with 10 RFs and 
#' ens_rf_ex <- ens_random_forests(df=simData$samples, var="obs", covariates=grep("cov", colnames(simData$samples),value=T), save=FALSE, cores=1)
#' 
#' ALE_df <- ALE_fn(ens_rf_ex, save=FALSE)
#' 
#' plot_ALE(ALE_df[1])
#' 
plot_ALE <- function(ALE, xquantiles=c(0.025,0.975), yquantiles = c(0.1, 0.5, 0.9), name, cex.axis=1, cex.lab=1, rug=TRUE, rug.col='gray50', rug.tick = 0.02, rug.lwd=0.5, rug.alpha=0.2, rug.max=1000){
	if(missing(name) & length(ALE)==1) name <- names(ALE)
	if(length(ALE)==1){
		ALEdf <- ALE[[1]]$df
		if(is.null(ALE[[1]]$X)){
			rug = FALSE
		}else { X <- ALE[[1]]$X }
	}else{
		stop('ALE must have structure list(df=..., X=...)')
	}
	if(ALEdf$class[1]=='factor'){
		ALEdf <- ALEdf[match(levels(X),ALEdf$x),]
		ALEdf[,4:ncol(ALEdf)] <- sapply(ALEdf[,4:ncol(ALEdf)],as.numeric)
		y.range <- range(unlist(ALEdf[,4:ncol(ALEdf)]))
		x.seq <- c(0,nrow(ALEdf)+1)
		x.range <- range(x.seq)
		quant <- t(apply(ALEdf[,4:ncol(ALEdf)], 1, quantile, probs=yquantiles))
		plot(1:nrow(ALEdf), quant[,2], 
	     type='n', 
	     xlim= range(pretty(x.range)), 
	     ylim = range(pretty(range(y.range))), 
	     las=1, xaxt='n',
	     xlab=parse(text=name), 
	     ylab="", 
	     cex.axis = cex.axis, 
	     cex.lab = cex.lab)
		axis(1, at=1:nrow(ALEdf), labels=ALEdf[,1], cex.axis=cex.axis)
		segments(x0 = 1:nrow(ALEdf),
		         x1 = 1:nrow(ALEdf), 
		        y0 = quant[,1],
		        y1=quant[,3], 
		        col = "gray30",
		        lwd=3)
		points(1:nrow(ALEdf),quant[,2], pch=16, cex=2)
		abline(h=0, lty=3)
		if(rug){
			tab <- table(X)
			text(1:nrow(ALEdf),par('usr')[3],tab,
			     font=3,adj=c(0.5,-0.1),col=rug.col)
		}
		
	}else{
		ALEdf$x <- as.numeric(ALEdf$x)
		if(length(unique(ALEdf$q))<5){
			keep <- rep(TRUE,nrow(ALEdf))
		}else{
			keep <- ALEdf$q>xquantiles[1] & ALEdf$q<xquantiles[2]
		}
		
		x.range <- ax.range(ALEdf$x[keep],ntry=2,rug=FALSE)
		y.range <- ax.range(unlist(ALEdf[keep,4:ncol(ALEdf)]),
		                    ntry=2,rug=rug)
		quant <- t(apply(ALEdf[keep,4:ncol(ALEdf)], 1, quantile, probs=yquantiles))

		plot(ALEdf[keep,1], quant[,2], 
	     type ='n', 
	     xlim = x.range$range, 
	     ylim = y.range$range, 
	     las = 1, 
	     xlab = parse(text=name), 
	     ylab = "", 
	     cex.axis = cex.axis, 
	     cex.lab = cex.lab,
	     xaxs = 'i', yaxs = 'i',
	     xaxt = 'n', yaxt = 'n')
		axis(1, at=x.range$pretty)
		axis(2, at=y.range$pretty, las=1)
		polygon(x = c(ALEdf[keep,1], rev(ALEdf[keep,1])), 
		        y = c(quant[,1], rev(quant[,3])), 
		        border = NA, 
		        col = "gray80")
		lines(ALEdf[keep,1], quant[,2], lwd=3)
		abline(h=0, lty=3)
		if(rug){
			Xrug <- X[sample.int(length(X),rug.max)]
			Axis(side = 1, at = Xrug, labels = FALSE, 
			     lwd = 0, lwd.ticks = rug.lwd, 
			     col.ticks = col2rgbA(rug.col,rug.alpha), 
			     tck = rug.tick)
		}
	}
}
ax.range <- function(val,tol=0.05,ntry=1,rug=TRUE){
	vr <- range(val)
	vr.tol <- c(vr[1]-abs(diff(vr))*tol,
	            vr[2]+abs(diff(vr))*tol)
	for(i in 1:ntry){
		if(i==1){
			pr <- pretty(vr)
		}else{
			pr <- pretty(pr)
		}
		pr <- pr[pr > vr.tol[1] & pr < vr.tol[2]]
	}
	rpr <- range(pr)
	if(rug){
		vr[1] <- pmin(vr.tol[1],rpr[1])
	}else{
		vr[1] <- pmin(vr[1],rpr[1])
	}
	vr[2] <- pmax(vr[2],rpr[2])
	return(list(range = vr,
	            pretty = pr))
}
col2rgbA<-function(color,transparency){
  rgb(t(col2rgb(color))/255,alpha=transparency)
}