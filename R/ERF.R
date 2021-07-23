#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
#			ERF Executable
#		
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Zach Siders

###-----------------------------------------------------
#		Initialization
###-----------------------------------------------------
	library(viridisLite)
	library(sf)
	library(raster)
	library(gridExtra)
	library(grid)
	library(gridBase)
	library(gtable)
	library(RANN)
	library(latticeExtra)
	library(doParallel)
	library(randomForest)
	library(rgeos)
	#also need MASS but lazyload

###-----------------------------------------------------
#		Paths
###-----------------------------------------------------
	
	setwd(path)
	
	dir.create(paste0(path,'Output/'))
	sub.dir <- paste0(path,"Output/ERF-",substr(Sys.time(),1,10))
	dir.create(file.path(sub.dir))

###-----------------------------------------------------
#		Functions
###-----------------------------------------------------
	source(paste0(path,"HelperFn.ERF.r"))
	source(paste0(path,"HelperFn.Plot.ERF.r"))
	## A helper function that tests whether an object is either NULL _or_ a list of NULLs
	is.NullOb <- function(x) is.null(x) | all(sapply(x, is.null)) && !is.data.frame(x)

	## Recursively step down into list, removing all such objects 
	rmNullObs <- function(x) {
	   x <- Filter(Negate(is.NullOb), x)
	   lapply(x, function(x) if (is.list(x)) rmNullObs(x) else x)
	}
	bbox.adj <- function(bbox,fac=0.1,hard.bnd=4){
		adj <- bbox
		xrange <- bbox$xmax - bbox$xmin
		yrange <- bbox$ymax - bbox$ymin
		xrange <- ifelse(xrange < 1, 1, xrange)
		yrange <- ifelse(yrange < 1, 1, yrange)
		xfac <- ifelse(xrange > yrange, fac, (yrange/xrange) * fac)
		yfac <- ifelse(yrange > xrange, fac, (xrange/yrange) * fac)
		xmulti <- ifelse(xfac*xrange < hard.bnd, hard.bnd, xfac*xrange)
		ymulti <- ifelse(yfac*yrange < hard.bnd, hard.bnd, yfac*yrange)
		adj[c(2,4)] <- bbox[c(2,4)] + c(-1,1)*ymulti
		adj[c(1,3)] <- bbox[c(1,3)] + c(-1,1)*xmulti
		adj[2] <- pmax(-90,adj[2])
		adj[4] <- pmin(90,adj[4])
		adj[1] <- pmax(0,adj[1])
		adj[3] <- pmin(360,adj[3])
		return(adj)
	}
###-----------------------------------------------------
#		Data Read In
###-----------------------------------------------------
	#CHECK THAT THE RData OBJECT NAME IS THE RIGHT ONE, AS WELL AS THE OBJECT IT LOADS.

	#-------
	# load set dataset
		load(paste0(data.path,"pt.covars.df.RData"))
		set.covars.df <- pt.covars.df
		pt.sf <- st_as_sf(pt.covars.df, coords=c('lon','lat'), 
		                  crs="+proj=longlat +datum=WGS84 +no_defs")
		bbox <- bbox.adj(st_bbox(pt.sf),fac=0.05)
		bbox.kde <- bbox.adj(st_bbox(pt.sf),fac=0.1)
		rm(pt.covars.df, pt.sf); gc()
	
	#-------
	# changes the names for plotting
		translator <- read.csv(paste0(path,'translator.csv'), as.is=TRUE, header=TRUE)

###-----------------------------------------------------
#		Add RV
###-----------------------------------------------------
	set.covars.df$random <- rnorm(nrow(set.covars.df), 0, 1)
###-----------------------------------------------------
#		covariates
###-----------------------------------------------------
	header <- c("SetID","TRIP_NUM","PERMIT_NUM",
	            "cent.lon", "cent.lat", 
				"area_km2", "soak_time_mins",
				"EFFORT.hr", "NUM_HKS_SET",
				'time','day','julian','month')
	env_simp <- c('bathymetry',
	              'lunar_rad', 'sst', 'chla',
	              'current.zonal', 'current.meridonal',
	              'current.speed', 'OkuboWeiss', 'EKE',
	              'wind.midpoint.zonal','wind.midpoint.meridonal',
	              'wind.midpoint.speed','wind.midpoint.vorticity',
	              'wind.midpoint.divergence',"seamt_ses", "random")
	env_moderate <- paste(env_simp,
	                      c('sst.front','dist.sst.front.ses',
	                        'chla.front','dist.current.front.ses',
	                        'current.front','SLA',
	                        'wind.front','dist.wind.front.ses'))
###-----------------------------------------------------
#		Execute
###-----------------------------------------------------
	
	duplicate = TRUE #repeat obs with more than 1 
#' Ensemble Random Forests
#' 
#' @description Execute Ensemble Random Forests model on a given dataset
#' 
#' @inheritParams erf_data_prep
#' @param out.folder A path to the folder to write out too. If NULL then a folder is generated in the working directory
#' @param duplicate A logical flag that indicates whether to duplicate observations with more than one interaction. Default is TRUE to duplicate all records that interacted with more than one individual (i.e. a fishing set that caught two of the same species)
#' @param n.forests A numeric value indicating how many Random Forests to generate in the ensemble, default is 100
#' @param cores A numeric value that either indicates the number of cores to use for parallel processing or a negative value to indicate the number of cores to leave free. Default is to leave two cores free.
#' 
#' @return A fitted ERF model 
#' 
#' @export
#' 
	erf <- function(df, var, covariates, out.folder=NULL, duplicate=TRUE, 
	                n.forests=100, cores=parallel::detectCores()-2)

	
		form <- erf_formula_prep(var, covariates) #Prepare the model formula
		v <- erf_data_prep(df, var, covariates, header, duplicate=duplicate)
		min_split <- min_splitter(v)
	
	# Ensemble Random Forest
		UseCores <- ifelse(n.forests < cores, n.forests, cores) #use less cores if n.forests < cores
		UseCores <- ifelse(min_split < 7e3, pmin(UseCores,floor(cores * (1-min_split/7e3))), 2) #use less cores if # of min_split exceeds a certain value (need some RAM available for calcs)

		cl <- makeCluster(UseCores) #make clusters
		registerDoParallel(cl) #designate cores
		rf.ens <- foreach(i=1:n.forests, .packages=c('randomForest','ROCR')) %dopar%{
			rf.ens.fn()
		}
		stopCluster(cl)
		
		
		pred_ens_p <- sapply(rf.ens, function(x) {x$preds$P.1}) #NULL line 
		pred_ens_z <- sapply(rf.ens, function(x) {x$preds$P.0})
		pred_ens_trAUC <- sapply(rf.ens, function(x) {x$roc_train$auc})
		pred_ens_teAUC <- sapply(rf.ens, function(x) {x$roc_test$auc})
		pred_ens_resid <- sapply(rf.ens, function(x) {as.numeric(as.character(x$preds$PRES)) - x$preds$P.1})

		pred_ens_l <- list(pred=pred_ens_p, resid = pred_ens_resid)

		pred_ens_pfin <- rowMeans(pred_ens_p)
		pred_ens_zfin <- rowMeans(pred_ens_z)

		pred_ens <- data.frame("P.0"=pred_ens_zfin, 
		                       "P.1"=pred_ens_pfin)
		pred_ens$PRES <- v[,var]
		pred_ens$resid <- as.integer(as.character(v[,1])) - pred_ens[,2]
		rownames(pred_ens) <- rownames(v)

		pack <- list(data = v, 
		             model = rf.ens, 
		             pred = pred_ens)
		save(pack, file=paste0(data.folder,"/randfor_ens_simple_",var,".Rdata"))

		rm(pred_ens_teAUC)
	 	rm(pred_ens_trAUC)
	 	rm(pred_ens_p)
	 	rm(pred_ens_z)
		rm(pred_ens_resid)
		rm(pred_ens_zfin)
	 	rm(pred_ens_pfin)
		gc()
	
	},silent=FALSE)
	