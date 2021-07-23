#' data simulation for ERF example Data
#' 
#' @description A function to calculate various ROC curve performance metrics given model predictions and the oberved samples
#' 
#' @param ncell number of spatial grid cells
#' @param ncov number of covariates to simulate
#' @param prob_missing probability of nondetection (1-detection)
#' @param unif.bnds bounds of the uniform distribution draws of the scale parameter for the Matérn autocorrelation
#' @param mat.var variance of the Matérn autocorrelation
#' @param nsamp number of random samples to draw
#' 
#' @return A list containing the samples, the beta coefficients, and the raster brick
#' @export
#' 
#' @examples
#' df <- data_sim()
#' head(df$samples)
#' table(df$samples$obs)
#' 
#' colr <- colorRampPalette(c('dodgerblue4',dodgerblue2','ivory','firebrick2','firebrick4'))
#' plot(df$grid$prob, col=colr(100), xaxt='n', yaxt='n', asp=NA)
#' with(df$samples[df$samples$obs==1,], points(x,y,pch=16))
#' 
data_sim <- function(ncell=100, ncov=5, prob_missing=0.95, unif.bnd=c(10,30), mat.var=0.05, nsamp=1e4){
	#grid
		Dim <- c("n_x"=ncell, "n_y"=ncell) #grid dimensions
		loc_xy <- expand.grid("x"=1:Dim['n_x'], "y"=1:Dim['n_y']) #coordinates

	# Autocorrelation parameters
		Scale <- runif(ncov, unif.bnd[1],unif.bnd[2])*(ncell/100) #scale parameters
	
	# Matern spatial autocorrelation models
		model <- sapply(Scale, function(x) {RMmatern(nu=1, var=mat.var, scale=x)})
		loc_xy[,paste0("cov",1:ncov)] <- sapply(model, function(d) array(RFsimulate(model=d, x=loc_xy[,'x'], y=loc_xy[,'y'])@data[,1], dim=Dim))
		cov.id <- grep('cov',colnames(loc_xy))
	#detrend
		loc_xy[,cov.id] <- apply(loc_xy[,cov.id],2,function(x) x-mean(x))
	# Define covariance matrix
	# 	intercept <- 0
	# 	corr.mat <- matrix(rbeta(ncov^2,2,2)*2-1, ncov, ncov) #simulate correlations between parameters; 0.25 is mean correlation, 1 is sd, generation occurs in logit space and then is back-transformed
	# 	corr.mat[lower.tri(corr.mat)] <- t(corr.mat)[lower.tri(corr.mat)] #make symmetric
	# 	diag(corr.mat) <- 1 #diagonal equals 1
		
	# 	# mu <- rnorm(ncov) #means of the multivariate normal
	# 	sigma <- apply(loc_xy[,cov.id],2,sd) #sd of the multivariate normal between 0 and 1

	# 	cov.mat <- cor2cov(corr.mat, sigma) #convert correlation matrix to covariance matrix using sd
	# 	cov.mat[lower.tri(cov.mat)] <- t(cov.mat)[lower.tri(cov.mat)] #make covariance matrix symmetric
	# #Draw Beta coefficients
	# 	betas <- MASS::mvrnorm(1, rep(-1,ncov), Matrix::nearPD(cov.mat)$mat)
		betas <- rnorm(5,-1,5)
	#Calculate probabilities
		loc_xy$prob.raw <- rowSums(sapply(1:ncov, function(x) {loc_xy[,x+2]*betas[x]}))
		loc_xy$prob <- inv_logit(loc_xy$prob.raw)
		prob <- matrix(loc_xy$prob, Dim)
	#Sample	
		samp <- sample(1:nrow(loc_xy), nsamp, replace=TRUE)
		#extract data at random locations
		p_xy_df <- loc_xy[samp,]
		#extract probability of presence at random locations
		p_xy_df$prob <- prob[cbind(p_xy_df$x, p_xy_df$y)]
		#draw process component
		p_xy_df$pred <- rbinom(nsamp, 1, p_xy_df$prob)
		#draw observation component
		p_xy_df$obs <- rbinom(nsamp, 1, (1-prob_missing)*p_xy_df$prob) * p_xy_df$pred
	#clean up
		rownames(p_xy_df) <- 1:nrow(p_xy_df)
		loc_arr <- raster::brick(array(as.matrix(loc_xy[,-c(1:2)]), dim=c(Dim, ncol(loc_xy)-2)), xmx=ncell, ymx=ncell)
		loc_arr <- t(loc_arr)
		loc_arr <- flip(loc_arr,'y')
		names(loc_arr) <- colnames(loc_xy)[-c(1:2)]
	return(list(samples = p_xy_df,
	            betas = betas,
	            grid = loc_arr))
}
cor2cov <- function(V, sigma) {
  V * tcrossprod(sigma)
}
logit <- function(x){
	log(x/(1-x))
}
inv_logit <- function(x){
	exp(x)/(exp(x)+1)
}