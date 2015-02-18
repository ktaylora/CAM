#
# parse output data generated from each model run into a composit table and dump this table into a SpatialPointsDataFrame
#

# default includes
require(rgdal)

argv <- commandArgs(trailingOnly=T)

files <- list.files(pattern="csv")                           # CSV files in CWD
  files <- files[(!grepl(files,pattern="processed.pop"))]    # Sanity check.  Don't use processed.pop.csv if it exists in CWD from previous run

out <- data.frame()                                          # Output data.frame (holds processed CSV data)

cat(" -- processing .")

# globals
if(length(list.files(pattern=argv[1]))<1) {
  SP_PTS_GRID = "gridPts"
} else {
  SP_PTS_GRID = argv[1]  
}

# iterate over csv files in CWD, calculating statistics and aggregating into our output data.frame as we go
for(f in files){
  t<-read.csv(f,stringsAsFactors=F)
    t<-na.omit(t) # strip NA values
      if(sum(is.na(t[nrow(t),])) > 1) { t<-t[1:(nrow(t)-1),]; } # bug fix: sometimes the last row in the CSV contains a NULL value.
        for(n in names(t)){ t[,n]<-as.numeric(t[,n]) } # force numeric for everything in the table
  
  site <- as.numeric(gsub("([0-9]+).*$", "\\1", f)) # GREP out the site_id from the CSV filename
  
  # calculate the site mean over time period and normalized variance for population size and seedbank size [~dispersion index] (mean/sd) 
  abg_mn <- mean(t[,'p.size']*t[,'ag.mdBiomass'], na.rm=T) # p.size * meanBiomass
    abg_mn <- round(abg_mn, 2)
  abg_id <- abg_mn/sd(t[,5]*t[,6], na.rm=T)
    abg_id <- round(abg_id,2)
  sb.size_mn <- mean(t[,14], na.rm=T) 
    sb.size_mn <- round(sb.size_mn)
  sb.size_id <- sb.size_mn/sd(t[,14], na.rm=T) # seedbank size
    sb.size_id <- round(sb.size_id)
  persistence <- nrow(t) # number of days across simulation the population survived (indicates population crashes)

  focal<-data.frame(site_id=site,abg_mn=abg_mn,abg_id=abg_id,sb.size_mn=sb.size_mn,sb.size_id=sb.size_id,pers=persistence);
  
  if(nrow(out)>0){
    out <- rbind(out,focal)
  } else {
    out <- focal
  }
  cat(".");
}; cat("\n");

cat(" -- writing to processed.pop.metrics.csv.\n")
write.csv(out, "processed.pop.metrics.csv", row.names=F)
cat(" -- done.\n")
cat(" -- writing to SpatialPointsDataFrame.\n")
s<-readOGR(".",SP_PTS_GRID, verbose=F)
  s.subs<-s[s@data$site_id %in% out$site,]

# order our tables by site_id 
s.subs <- s.subs[order(s.subs@data$site_id,decreasing=F),] 
out <- out[order(out$site_id,decreasing=F),]

if(sum(s.subs@data$site_id != out$site_id) == 0){ # sanity check: do all of our sites match between our tables?
  s.subs@data <- out
  writeOGR(s.subs,".","modelRunResults", driver="ESRI Shapefile", overwrite=T)
  cat(" -- done.\n")
} else {
  stop(" -- error: site_id's in gridPts and processed table don't match.  quitting.")
}
