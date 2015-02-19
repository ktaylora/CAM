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

# helper functions
processChunks <- function(x,var=NULL,size=30) { as.vector(split(as.numeric(x[,var]), ceiling(seq_along(x[,var])/size))) }

# iterate over csv files in CWD, calculating statistics and aggregating into our output data.frame as we go
for(f in files){
  t<-read.csv(f,stringsAsFactors=F)
    t<-na.omit(t) # strip NA values
      if(sum(is.na(t[nrow(t),])) > 1) { t<-t[1:(nrow(t)-1),]; } # bug fix: sometimes the last row in the CSV contains a NULL value.
        for(n in names(t)){ t[,n]<-as.numeric(t[,n]) } # force numeric for everything in the table
  
  site <- as.numeric(gsub("([0-9]+).*$", "\\1", f)) # GREP out the site_id from the CSV filename
  
  # calculate the site mean over time period and normalized variance for population size and seedbank size [~dispersion index] (mean/sd) 
  t_years <- split(t, f=t$year)
  t_years_ordered <- order(as.numeric(names(t_years)))
  
  # months of the year definitions
  winter <- c(12,13,1,2)
  spring <- c(3,4,5)
  summer <- c(6,7,8)
    fall <- c(9,10,11) 


  monthly_aggregated_winterBiomass <- rep(NA,length(t_years_ordered))
  monthly_aggregated_springBiomass <- rep(NA,length(t_years_ordered))
  monthly_aggregated_summerBiomass <- rep(NA,length(t_years_ordered))
    monthly_aggregated_fallBiomass <- rep(NA,length(t_years_ordered))
  monthly_aggregated_winterSeedbank <- rep(NA,length(t_years_ordered))
  monthly_aggregated_springSeedbank <- rep(NA,length(t_years_ordered))
  monthly_aggregated_summerSeedbank <- rep(NA,length(t_years_ordered))
    monthly_aggregated_fallSeedbank <- rep(NA,length(t_years_ordered))
  monthly_aggregated_winterRtLen <- rep(NA,length(t_years_ordered))
  monthly_aggregated_springRtLen <- rep(NA,length(t_years_ordered))
  monthly_aggregated_summerRtLen <- rep(NA,length(t_years_ordered))
    monthly_aggregated_fallRtLen <- rep(NA,length(t_years_ordered))
  
  # chunk each year into monthly intervals
  mon_rootLen <- lapply(t_years, FUN=processChunks, var='p.mnRtLen', size=30)
     mon_rootLen <- mon_rootLen[(names(mon_rootLen) != "NULL")]
  mon_age <- lapply(t_years, FUN=processChunks, var='p.meanAge', size=30)
     mon_age <- mon_rootLen[(names(mon_age) != "NULL")]
  mon_pop_size <- lapply(t_years, FUN=processChunks, var='p.size', size=30)
     mon_pop_size <- mon_pop_size[(names(mon_pop_size) != "NULL")]
  mon_abg_bio  <- lapply(t_years, FUN=processChunks, var='ag.mdBiomass', size=30) 
    mon_abg_bio <- mon_abg_bio[(names(mon_abg_bio) != "NULL")]
  mon_sb_size  <- lapply(t_years, FUN=processChunks, var='sb.size', size=30) 
    mon_sb_size <- mon_sb_size[(names(mon_sb_size) != "NULL")]

  # iterate over each year in the aggregation
  for(i in 1:length(mon_pop_size)){
    # take the monthly mean across all months for year i
    popSize <- round(unlist(lapply(mon_pop_size[[i]], FUN=mean, na.rm=T)))
    abgBiomass <- unlist(lapply(mon_abg_bio[[i]], FUN=mean, na.rm=T))
    seedbankSize <- round(unlist(lapply(mon_sb_size[[i]], FUN=mean, na.rm=T)))
<<<<<<< HEAD
    rootLength <- round(unlist(lapply(mon_rootLen[[i]], FUN=mean, na.rm=T)),2)
=======
    rootLength <- round(unlist(lapply(mon_rootLen[i]], FUN=mean, na.rm=T)),2)
>>>>>>> aae792cd5e63004c8bb755b670e08fd377071545

    # aggregate our monthly means by season
    monthly_aggregated_winterBiomass[i] <- round(mean(popSize[winter],na.rm=T))*mean(abgBiomass[winter],na.rm=T)
    monthly_aggregated_springBiomass[i] <- round(mean(popSize[spring], na.rm=T))*mean(abgBiomass[spring],na.rm=T)
    monthly_aggregated_summerBiomass[i] <- round(mean(popSize[summer], na.rm=T))*mean(abgBiomass[summer],na.rm=T)
    monthly_aggregated_fallBiomass[i]   <- round(mean(popSize[fall], na.rm=T))*mean(abgBiomass[fall],na.rm=T)
    monthly_aggregated_winterSeedbank[i] <- round(mean(seedbankSize[winter],na.rm=T))
    monthly_aggregated_springSeedbank[i] <- round(mean(seedbankSize[spring],na.rm=T))
    monthly_aggregated_summerSeedbank[i] <- round(mean(seedbankSize[summer],na.rm=T))
    monthly_aggregated_fallSeedbank[i]   <- round(mean(seedbankSize[fall],na.rm=T))
    monthly_aggregated_winterRtLen[i] <- round(mean(rootLength[winter],na.rm=T))
    monthly_aggregated_springRtLen[i] <- round(mean(rootLength[spring],na.rm=T))
    monthly_aggregated_summerRtLen[i] <- round(mean(rootLength[summer],na.rm=T))
    monthly_aggregated_fallRtLen[i]   <- round(mean(rootLength[fall],na.rm=T))
  }

  # statistics to record to attribute table
  abg_mn_mon <- round(mean(unlist(lapply(ls(pattern=glob2rx("monthly_*Biomass")), FUN=get)),na.rm=T),2)
  abg_se_mon <- round(sd(unlist(lapply(ls(pattern=glob2rx("monthly_*Biomass")), FUN=get)),na.rm=T),2)
  abg_mn_spr <- round(mean(monthly_aggregated_springBiomass, na.rm=T),2)
  abg_mn_sum <- round(mean(monthly_aggregated_summerBiomass, na.rm=T),2)
  abg_mn_win <- round(mean(monthly_aggregated_winterBiomass, na.rm=T),2)
  abg_mn_fal <- round(mean(monthly_aggregated_fallBiomass, na.rm=T),2)

  sb_mn_mon <-  round(mean(unlist(lapply(ls(pattern=glob2rx("monthly_*Seedbank")), FUN=get)),na.rm=T),2)
  sb_se_mon <-  round(sd(unlist(lapply(ls(pattern=glob2rx("monthly_*Seedbank")), FUN=get)),na.rm=T),2)
  sb_mn_spr <-  round(mean(monthly_aggregated_springBiomass, na.rm=T),2)
  sb_mn_sum <-  round(mean(monthly_aggregated_summerBiomass, na.rm=T),2)
  sb_mn_win <-  round(mean(monthly_aggregated_winterBiomass, na.rm=T),2)
  sb_mn_fal <-  round(mean(monthly_aggregated_fallBiomass, na.rm=T),2)

  rl_mn_mon <-  round(mean(unlist(lapply(ls(pattern=glob2rx("monthly_*RtLen")), FUN=get)),na.rm=T),2)
  rl_se_mon <-  round(sd(unlist(lapply(ls(pattern=glob2rx("monthly_*RtLen")), FUN=get)),na.rm=T),2)
  rl_mn_spr <-  round(mean(monthly_aggregated_springRtLen, na.rm=T),2)
  rl_mn_sum <-  round(mean(monthly_aggregated_summerRtLen, na.rm=T),2)
  rl_mn_win <-  round(mean(monthly_aggregated_winterRtLen, na.rm=T),2)
  rl_mn_fal <-  round(mean(monthly_aggregated_fallRtLen, na.rm=T),2)

  persistence <- nrow(t) # number of days across simulation the population survived (indicates population crashes)

  focal<-data.frame(site_id=site,abg_mn_mon=abg_mn_mon,abg_se_mon=abg_se_mon,abg_mn_spr=abg_mn_spr,abg_mn_sum=abg_mn_sum,abg_mn_win=abg_mn_win,abg_mn_fal=abg_mn_fal,
                    sb_mn_mon=sb_mn_mon,sb_se_mon=sb_se_mon,sb_mn_spr=sb_mn_spr,sb_mn_sum=sb_mn_sum,sb_mn_win=sb_mn_win,sb_mn_fal=sb_mn_fal,
                    rl_mn_mon= rl_mn_mon,rl_se_mon=rl_se_mon,rl_mn_spr=rl_mn_spr,rl_mn_sum=rl_mn_sum,rl_mn_win=rl_mn_win,rl_mn_fal=rl_mn_fal,
                    pers=persistence)
  
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
