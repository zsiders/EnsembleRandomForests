#' Visualizes variable importance from ERF
#' 
#' @description Plots the variable importance from the Ensemble Random Forests model
#'
#' @param fit The fitted object returned from calling ens_random_forests()
#' 
#' @return A plot of the variable importance distributions ordered by the mean variable importance.
#' 
#' @export
#' 
#' @examples
#' #run an ERF with 10 RFs and 
#' ens_rf_ex <- ens_random_forests(df=simData$samples, var="obs", covariates=grep("cov", colnames(simData$samples),value=T), save=FALSE, cores=1)
#' 
#' erf_impt(ens_rf_ex)
#' 
erf_impt <- function(fit, var.names, pal){
	if(missing(fit)) stop('Supply fit object')
	if(is.null(fit$var.imp)) stop('Fit object is missing variable importance')
	var.imp <- fit$var.imp
	if(missing(var.names)) var.names <- rownames(var.imp)
	nvars <- nrow(var.imp)
	max.v <- max(var.imp[,1:3])
	var.imp.o <- var.imp[var.imp$ord,][nvars:1,]
	rand.i <- which(var.names[var.imp$ord][nvars:1]=='random')

	par(mar=c(4,max(nchar(var.names))/1.5,1,1))
	plot(seq(0, max.v, length.out=nvars), 
	     seq(1,nvars), 
	     type='n', 
	     xlab="Mean Decrease in Accuracy", 
	     ylab="", bty='l', 
	     xaxt='n', yaxt='n', 
	     ylim=c(0,nvars+1))
	axis(2, seq(1,nvars), parse(text=var.names[var.imp$ord][nvars:1]), las=1)
	
	axis(1, pretty(c(0,max.v)))
	
	for(i in 1:nvars){
		d <- den.sc(rnorm(10000, 
		                   var.imp.o$mu[i], 
		                   var.imp.o$sd[i]),
					to.y = c(0,0.8)+i,
					probs=c(0.5,0.1,0.9,0.25,0.75),adj=3,
					from=0)
		
		med_seq <- seq(min(var.imp.o[,2]),
		               max(var.imp.o[,2]), 
		               length.out=100)
		med_id <- which.min(abs(med_seq - var.imp.o[i,2]))

		greys <- colorRampPalette(rev(c('white','grey50','black')))
		if(missing(pal)){
			pal <- colorRampPalette(c('#fcde9c', '#faa476', '#f0746e', '#e34f6f', '#dc3977', '#b9257a', '#7c1d6f'))
		}
		if(i==rand.i){
			cols <-'gray50'
		}else if(i < rand.i){
			cols <- "gray95"
		}else{
			cols <- pal(100)[med_id]
		}
		polygon(x=c(d$d$x.q, rev(d$d$x.q)), 
		        y=c(rep(i,length(d$d$y.q)), rev(d$d$y.q)), 
		        col=cols, border=FALSE)
		
		segments(x0=d$d$q[2,1], 
		         y0=i, 
		         x1=d$d$q[3,1], 
		         y1=i, lty=3)

		segments(x0=d$d$q[1,1], 
		         y0=i, 
		         x1=d$d$q[1,1], 
		         y1=d$d$q[1,2],
		         lty=1, col=greys(100)[med_id])

		lines(d$d$x, d$d$y, lwd=2, col='black')
	}
	legend("bottomright", 
	       legend=c("80% CI","50% CI","Median"), 
	       pch=c(NA,15,124), 
	       col=c('black',pal(5)[3],'black'), 
	       lty=c(3, NA, NA), 
	       lwd=c(3, NA, NA), 
	       pt.cex=c(NA,2,1),
	       inset=c(0,0), 
	       bty='n', seg.len = 1.1)
}
den.sc <- function(x,to.x,to.y,from,ret.x=FALSE,probs,...){
	x <- na.omit(unlist(x))

	if(!missing(from)){
		if(length(from)>1){
			d <- density(x, from=from[1],to=from[2],...)
		}else{
			d <- density(x, from=from[1],...)
		}
		
	}else{
		d <- density(x, ...)
	}
	dx <- d$x
	d$y <- d$y/max(d$y)
	px <- pretty(d$x)
	if(!missing(to.y)){
		d$y <- scales::rescale(d$y,to=to.y)
	}
	if(!missing(to.x)){
		x.ax <- scales::rescale(px,from=range(dx),to=to.x)
		d$x <- scales::rescale(d$x,to=to.x)
		if(!missing(probs)){
			q <- quantile(x,probs)
			q.id <- sapply(q, function(x)which.min(abs(dx-x)))
			d$q <- cbind(d$x[q.id],d$y[q.id])
		}
		if(ret.x){
			return(list(d=d,ax=data.frame(v=px,
	                              p=x.ax),x=dx))
		}else{
			return(list(d=d,ax=data.frame(v=px,
	                              p=x.ax)))
		}
	}else{
		if(!missing(probs)){
			q <- quantile(x,probs)
			q.id <- sapply(q, function(x)which.min(abs(dx-x)))
			d$q <- cbind(d$x[q.id],d$y[q.id])
			d$x.q <- d$x[q.id[4]:q.id[5]]
			d$y.q <- d$y[q.id[4]:q.id[5]]
		}
		if(ret.x){
			return(list(d=d,x=dx))
		}else{
			return(list(d=d))
		}
		
	}
}