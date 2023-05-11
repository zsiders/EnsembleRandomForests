#' Plots Accumulated Local Effects for ERF for multivariate factors
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
#' @param poly.colp is a vector of colors equal to n levels in response metric
#' @param poly.alpha is a transparency value between 0 and 1 for the confidence interval polygon from the ERF predictions
#' 
#' @return A plot of the ALEs
#' 
#' @export
#' 
#' @examples
#' #run an ERF with 10 RFs and 
#' logit <- function(x){log(x/(1-x))}
#' inv_logit <- function(x){exp(x)/(exp(x)+1)}
#' 
#' x_mat <- as.data.frame(replicate(4, rnorm(1e4)))
#' x_mat_bin <- t(rmultinom(1e4,1,prob=c(0.33,0.33,0.33)))
#' x_mat$V5 <- factor(apply(x_mat_bin,1,function(x)which(x==1)))
#' 
#' y_gin <- function(){
#' 	eff <- rnorm(ncol(x_mat)-1)
#'	eff_bin <- rnorm(ncol(x_mat_bin))
#'	y_logit <- t(eff %*% t(x_mat[,-5])) + t(eff_bin %*% t(x_mat_bin))
#'	y <- inv_logit(y_logit)
#'	return(y)
#' }
#' 
#' y_poss <- data.frame(y1 = y_gin(),
#'                     y2 = y_gin(),
#'                     y3 = y_gin())
#' 
#' y_bin <- apply(y_poss,1,function(x) rmultinom(1,1,x))
#' y <- apply(y_bin,2,function(x)which(x==1))
#' df <- data.frame(y = factor(y),x_mat)
#' 
#' ens_rf_ex <- ens_random_forests(df=df, var="obs", covariates=colnames(df[,-1]), save=FALSE, cores=1)
#' 
#' 
#' ALE_df <- ALE_fn(ens_rf_ex, save=FALSE, multi=TRUE, type='prob')
#' 
#' plot_ALE_multi(ALE_df[1])
#' 
plot_ALE_multi <- function(ALE, xquantiles=c(0.025,0.975), yquantiles = c(0.1, 0.5, 0.9), name, cex.axis=1, cex.lab=1, rug=TRUE, rug.col='gray50', rug.tick = 0.02, rug.lwd=0.5, rug.alpha=0.2, rug.max=1000, poly.colp = c("#fcde9c","#e34f6f","#7c1d6f"),poly.alpha=0.5, level.names){
	if(missing(name) & length(ALE)==1) name <- names(ALE)
	if(length(ALE)==1){
		ALEdf <- ALE[[1]]$df
		if(is.null(ALE[[1]]$X)){
			rug = FALSE
		}else { X <- ALE[[1]]$X }
	}else{
		stop('ALE must have structure list(df=..., X=...)')
	}
	if(missing(level.names)) level.names <- as.character(seq(1:length(ALEdf)))
	if(ALEdf[[1]]$class[1]=='factor'){
		ALEdf <- lapply(ALEdf, function(X){X[,4:ncol(X)] <- sapply(X[,4:ncol(X)],as.numeric); return(X)})
		y.range <- range(unlist(lapply(ALEdf,function(X)X[,4:ncol(X)])))
		x.tick <- 1:nrow(ALEdf[[1]])
		x.seq <- c(0,nrow(ALEdf[[1]])+1)
		x.range <- range(x.seq)
		quant <- lapply(ALEdf, function(X) t(apply(X[,4:ncol(X)], 1, quantile, probs=yquantiles)))
		plot(x.tick, quant[[1]][,2], 
	     type='n', 
	     xlim= range(pretty(x.range)), 
	     ylim = range(pretty(range(y.range))), 
	     las=1, xaxt='n',
	     xlab=parse(text=name), 
	     ylab="", 
	     cex.axis = cex.axis, 
	     cex.lab = cex.lab)
		axis(1, at=x.tick, 
		     labels=ALEdf[[1]][,1], cex.axis=cex.axis)
		x.adj <- (x.tick - median(x.tick)) * 0.1
		for(i in 1:length(quant)){
			segments(x0 = x.tick+x.adj[i],
		         x1 = x.tick+x.adj[i], 
		        y0 = quant[[i]][,1],
		        y1=quant[[i]][,3], 
		        col = poly.colp[i],
		        lwd=3)
			points(x.tick+x.adj[i],quant[[i]][,2], 
			       pch=16, cex=2, col=poly.colp[i])
		}
		
		
		abline(h=0, lty=3)
		if(rug){
			tab <- table(X)
			text(x.tick,par('usr')[3],tab,
			     font=3,adj=c(0.5,-0.1),col=rug.col)
		}
		
	}else{
		ALEdfx <- as.numeric(ALEdf[[1]]$x)
		if(length(unique(ALEdf[[1]]$q))<5){
			keep <- rep(TRUE,nrow(ALEdf[[1]]))
		}else{
			keep <- ALEdf[[1]]$q>xquantiles[1] & ALEdf[[1]]$q<xquantiles[2]
		}
		
		x.range <- ax.range(ALEdfx[keep],ntry=2,rug=FALSE)
		y.range <- ax.range(unlist(lapply(ALEdf,function(X)X[keep,4:ncol(X)])),
		                    ntry=2,rug=rug)
		quant <- lapply(ALEdf, function(X) t(apply(X[keep,4:ncol(X)], 1, quantile, probs=yquantiles)))

		plot(ALEdfx[keep], quant[[1]][,2], 
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
		for(i in 1:length(quant)){
			polygon(x = c(ALEdfx[keep], rev(ALEdfx[keep])), 
		        y = c(quant[[i]][,1], rev(quant[[i]][,3])), 
		        border = poly.colp[i], 
		        col = col2rgbA(poly.colp[i],poly.alpha))
			lines(ALEdfx[keep], quant[[i]][,2], lwd=3,
			      col = poly.colp[i])
		}
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