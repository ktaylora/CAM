#!/bin/bash

#
# CAM Bash call script.  
# Usage: Accepts 1 argument at run-time containing <site_id> (e.g., # bash cam-soilwat-unix-wrapper.sh 25).
# cam-soilwat-unix-wrapper.R will automatically parse the CWD for a soildatabase shapefile (grid.pts.processed.shp) containing 
# the requisite soil information for SOILWAT and a field (site_id) specifying all of the sites across the extent of the run.  
# <site_id> is parsed from this list.
#
# The script is intended to be generalized to multiple sites / cores by the user in-shell, for example:
# for file in `ls -f -1 sites/split/*`; do #
#   for site in `cat $file`; do 
#     bash cam-soilwat-unix-wrapper.sh $site & 
#   done; 
#   while [[ `ps aux | grep "R --" | grep -v "grep"` > 3 ]]; do sleep 10; done; 
#  done
# 
# 

site=$1; rm -rf $site.log; 
oTime=`{ time nohup R --no-save --vanilla --slave --args $site < cam-soilwat-unix-wrapper.R >> $site.log & } 2>&1 | grep real` 
echo $oTime | awk '{ print "runtime: " $2 }' >> $site.log
grep -v " -" $site.log | grep -v "runtime" >> $site.csv
