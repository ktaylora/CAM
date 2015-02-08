#
# author: Kyle Taylor (kyle.a.taylor@gmail.com)
#

require(raster)
require(rgdal)

##
## Various helper functions 
##

# 
# SWPtoVWC() 
# convert SWP to VWC, e.g., to calculate field capacity and wilting point
# Author: D. Schlaepfer (drschlaepfer@uwyo.edu) [2014]
#

SWPtoVWC <- function(swp, sand, clay) {
#Cosby, B. J., G. M. Hornberger, R. B. Clapp, and T. R. Ginn. 1984. A statistical exploration of the relationships of soil moisture characteristics to the physical properties of soils. Water Resources Research 20:682-690.
	
	#1. SWP in MPa [single value] + sand and clay in fraction [single values] --> VWC in fraction [single value]
	#2. SWP in MPa [single value] + sand and clay in fraction [vectors of length d] --> VWC in fraction [vector of length d]
	#3. SWP in MPa [vector of length l] + sand and clay in fraction [single values] --> VWC in fraction [vector of length l]
	#4. SWP in MPa [vector of length l] + sand and clay in fraction [vectors of length d] --> VWC in fraction [matrix with nrow=l and ncol=d, SWP vector repeated for each column]: probably not used
	#5. SWP in MPa [matrix with nrow=l and ncol=d] + sand and clay in fraction [single values] --> VWC in fraction [matrix with nrow=l and ncol=d]
	#6. SWP in MPa [matrix with nrow=l and ncol=d] + sand and clay in fraction [vectors of length d] --> VWC in fraction [matrix with nrow=l and ncol=d, sand/clay vector repeated for each row]
	
	stopifnot(length(sand) == length(clay))
	na.act <- na.action(na.exclude(apply(data.frame(sand, clay), MARGIN=1, FUN=sum)))
	
	if(length(sand) > length(na.act)){
		na.index <- as.vector(na.act)
		
		if(length(na.index) > 0){
			sand <- sand[-na.index]
			clay <- clay[-na.index]
		}
		
		thetas <- -14.2 * sand - 3.7 * clay + 50.5
		psis <- 10 ^ (-1.58 * sand - 0.63 * clay + 2.17)
		b <- -0.3 * sand + 15.7 * clay + 3.10
		if(any(b <= 0)) stop("b <= 0")
		
		bar_conversion <- 1024
		MPa_toBar <- -10
		
		get_vector <- function(swp, sand, clay, thetas=thetas, psis=psis, b=b, do.na=TRUE){#either swp or sand/clay needs be a single value
			vwc <- ifelse(!is.na(swp) & swp <= 0 & sand <= 1 & sand >= 0 & clay <= 1 & clay >= 0, thetas * (psis / (swp * MPa_toBar * bar_conversion))^(1/b) / 100, NA)
			if(do.na & length(na.index) > 0){
				vwc <- napredict(na.act, vwc)
			}
			return(vwc)
		}
		
		if(is.null(dim(swp))){
			if(length(swp) == 1 & length(sand) >= 1 | length(swp) >= 1 & length(sand) == 1){ #cases 1-3		
				vwc <- get_vector(swp, sand, clay, thetas=thetas, psis=psis, b=b)
			} else if(length(swp) > 1 & length(sand) > 1){ #case 4
				vwc <- t(sapply(1:length(swp), FUN=function(d) get_vector(swp[d], sand, clay, thetas=thetas, psis=psis, b=b)))
			}
		} else {
			if(length(sand) == 1){ #case 5
				vwc <- sapply(1:ncol(swp), FUN=function(d) get_vector(swp[, d], sand, clay, thetas=thetas, psis=psis, b=b))
			} else { #case 6
				sand <- napredict(na.act, sand)
				clay <- napredict(na.act, clay)
				stopifnot(ncol(swp) == length(sand))
				psis <- napredict(na.act, psis)
				thetas <- napredict(na.act, thetas)
				b <- napredict(na.act, b)
				vwc <- sapply(1:ncol(swp), FUN=function(d) get_vector(swp[, d], sand[d], clay[d], thetas=thetas[d], psis=psis[d], b=b[d], do.na=FALSE))
			}
		}
	} else {
		vwc <- swp
		vwc[!is.na(vwc)] <- NA
	}
	return(vwc) #fraction m3/m3 [0, 1]
}

#
# HWSDGridPts_parse()
# Accepts the HWSD raster and a SpatialPointsDataFrame containing arbitrary site_id's and populates the DataFrame with MU_GLOBAL values associated with each point that
# can be used by HWSDToGridPts_populate() to fill in the soil characterists for each corresponding site.
#
# Author: Kyle Taylor (kyle.a.taylor@gmail.com)
#

HWSDToGridPts_extract <- function(r=NULL,p=NULL){
  require(raster)
  require(rgdal)
  # sanity-checks
  if(is.null(r)) r <- raster(paste(Sys.getenv("HOME"), "Products/soils/world/hwsd.bil", sep="/"))
  # extract
  p<-spTransform(p, CRS(projection(r)))
    p$site_id <- 1:nrow(p)
    p$HWSD <- raster::extract(r,p)
      return(p[,c("site_id","HWSD")])
}

#
# HWSDToGridPts_populate()
#
# Author: Kyle Taylor (kyle.a.taylor@gmail.com)
#

HWSDToGridPts_populate <- function(gridPts=NULL, t=NULL){
  #
  # sanity-check our input
  #
  if(is.null(gridPts)) gridPts <- try(readOGR(paste(Sys.getenv("HOME"),"Products/uw/soilwat_wna_runs_10_km_grid/", sep="/"),"grid.pts")) 
    if(class(gridPts) == "Try-Error") stop("please pass a SpatialPointsDataFrame to gridPts= that we can use for our HWSD extraction.\n")
  # did the user pass an HWSD translation table?  
  if(is.null(t)) t <- read.csv(paste(Sys.getenv("HOME"), "Products/soils/world/HWSD_DATA.csv", sep="/"))
    if(class(t) == "Try-Error") stop("please pass a data.frame containing HWSD MU_GLOBAL codes and corresponding soil data to t= that we can use for calculating SOILWAT site parameters.\n") 
  # ensure we have an HWSD column in our input shapefile
  if(sum(names(gridPts) == "HWSD") < 1) gridPts<-HWSDToGridPts_extract(p=gridPts) 

  results <- data.frame()

  cat("\n\n(PROCEDURAL) SOIL TEXTURE PARSING INTERFACE FOR THE HWSD\n")
  cat("\n  -- processing: ")

  for(i in 1:nrow(gridPts)){
    cat(".");
    # calculate THIS result
    result <- data.frame(site_id=i, depth=NA, topsoil_impermeabilityFraction=NA, topsoil_sandFraction=NA,
					     topsoil_clayFraction=NA, subsoil_impermeabilityFraction=NA, subsoil_sandFraction=NA, subsoil_clayFraction=NA, 
					     topsoil_bulkDensity=NA, subsoil_bulkDensity=NA);

    # prime-the-pump and make sure we are working with data that at least has a reference depth ~ covering the depth of our soilwat run.  
    # otherwise throw it out

	focal <- t[t$MU_GLOBAL == as.numeric(gridPts[i,2]@data),]
	  focal <- focal[which(!is.na(focal$REF_DEPTH)),]

    if(nrow(focal)>0) { # sanity check the record length
		if(sum(((focal$SHARE)/100)*focal$REF_DEPTH) > 90){ # second sanity check : make sure that the reference depth is > 90 cm
		  # calculate depth 
		  if(sum(!is.na(floor((focal$SHARE)/100*focal$ROOTS))) > 0){ # any non-NA values in the ROOT OBS. columns?
		  	rootClasses <- c(0,80,70,50,30,40,10) # units: cm (yes, they are ordered irregularly)
		  	result$depth <- sum(as.numeric(na.omit(floor((focal$SHARE)/100*focal$ROOTS)))); 
		  	  if(result$depth != 7) { result$depth <- result$depth+1; }
		  	    result$depth <- rootClasses[result$depth]
		  } else if(sum(!is.na(floor((focal$SHARE)/100*focal$IL))) > 0){  # any non-NA values in the IL columns?
		  	ilClasses <- c(0,150,115,60,40) # units: cm (yes, they are ordered irregularly)
		    result$depth <- sum(as.numeric(na.omit(floor((focal$SHARE)/100*focal$IL))))
		      if(result$depth != 5) { result$depth <- result$depth+1; }
		        result$depth <- ilClasses[result$depth]
		  } else { # default to using the unit reference depth
		  	result$depth <- sum(((focal$SHARE)/100)*focal$REF_DEPTH) # units: cm
		  }
		  # calculate weighted sand fraction
		  result$topsoil_sandFraction <- sum(((focal$SHARE)/100)*focal$T_SAND)/100
		  result$subsoil_sandFraction <- sum(((focal$SHARE)/100)*focal$S_SAND)/100
		  # calculate weighted clay fraction
		  result$topsoil_clayFraction <- sum(((focal$SHARE)/100)*focal$T_CLAY)/100
		  result$subsoil_clayFraction <- sum(((focal$SHARE)/100)*focal$S_CLAY)/100
		  # calculate field capacity
		  result$topsoil_fieldCapacity <- SWPtoVWC(-0.033, result$topsoil_sandFraction, result$subsoil_sandFraction)
		  result$subsoil_fieldCapacity <- SWPtoVWC(-0.033, result$topsoil_sandFraction, result$subsoil_sandFraction)
		  # calculate wilting point
		  result$topsoil_wiltingPoint <- SWPtoVWC(-1.5, result$topsoil_sandFraction, result$subsoil_sandFraction)
		  result$subsoil_wiltingPoint <- SWPtoVWC(-1.5, result$topsoil_sandFraction, result$subsoil_sandFraction)		  
		  # calculate impermability fraction
		  if(result$depth < 100){
		  	if(result$depth <= 30){
		      result$topsoil_impermeabilityFraction <- 1
		  	} else {
		  	  result$topsoil_impermeabilityFraction <- 0; 
		  	}
		  	if(result$depth <= 90){
		  	  result$subsoil_impermeabilityFraction <- 1
		  	} else {
		  	  result$subsoil_impermeabilityFraction <- 0;
		  	}
		  } else {
		  	result$topsoil_impermeabilityFraction <- 0; 
		  	result$subsoil_impermeabilityFraction <- 0;
		  }
		  # calculate bulk density
		  result$topsoil_bulkDensity <- sum((focal$SHARE)/100*focal$T_REF_BULK_DENSITY)
		  result$subsoil_bulkDensity <- sum((focal$SHARE)/100*focal$S_REF_BULK_DENSITY)
		  # assign our calculated result to the table
		  result <- round(result,3)
		  if(nrow(results) > 0){ 
		  	results <- rbind(results,result); 
		  } else { 
		  	results <- result; 
		  }
	    } 
    };   
  }; cat("\n");

  gridPts <- gridPts[which(gridPts@data$site_id %in% results$site_id),];
    gridPts@data <- results
      return(gridPts);
}

#
# getSwInputDataBySiteID()
#
# Helper function to parse a processed HWSD grid points SpatialPointsDataFrame into a sw_inputData object that Rsoilwat
# can work with. Accepts s=SpatialPointsDataFrame containing processed site information (as provided by HWSDToGridPts_populate()), which
# can be processed according to the site_id=site specified by the user.
#
# Kyle Taylor
# Author: Kyle Taylor (kyle.a.taylor@gmail.com)
#

getSwInputDataBySiteID <- function(s=NULL, site_id=NULL, db=NULL){
  t<-s@data
    t<-t[which(t$site_id == site_id),]

  template <- db
  if(is.null(template)){ 
  	template <- Rsoilwat::sw_inputData()
  } else if(class(template) != "swInputData"){
    stop(" -- unknown template object passed as input.")
  }

  ## parse out the top and bottom soil layers associated with the HWSD to our SOILWAT soils layers
  
  # assign corresponding layers
  topsoils <- template@soils@Layers[template@soils@Layers[,1] < 30,]
  subsoils <- template@soils@Layers[template@soils@Layers[,1] > 30,]
  
  # populate
  topsoils[,2] <- t$tpsl_bD  # bulk denisty
  topsoils[,3] <- t$tpsl_fC  # field capacity
  topsoils[,4] <- t$tpsl_wP  # wilting point

  topsoils[,9]  <- t$tpsl_sF  # sand fraction
  topsoils[,10] <- t$tpsl_cF # clay fraction
  topsoils[,11] <- t$tpsl_mF # impermeability fraction
  #topsoils[,12] <- t$tpsl_sT # soil temperature

  subsoils[,2] <- t$sbsl_bD  # bulk denisty
  subsoils[,3] <- t$sbsl_fC  # field capacity
  subsoils[,4] <- t$sbsl_wP  # wilting point
  
  subsoils[,9]  <- t$sbsl_sF  # sand fraction
  subsoils[,10] <- t$sbsl_cF # clay fraction
  subsoils[,11] <- t$sbsl_mF # impermeability fraction
  #subsoils[,12] <- t$sbsl_sT # soil temperature

  # merge and assign table to @soils@Layers
  template@soils@Layers <- rbind(topsoils,subsoils)

  return(template)
}
