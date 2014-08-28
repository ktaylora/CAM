#!/bin/bash
site=$1; rm -rf $site.log; 
oTime=`{ time nohup R --no-save --vanilla --slave --args $site < cam-soilwat-unix-wrapper.R >> $site.log & } 2>&1 | grep real` 
echo $oTime | awk '{ print "runtime: " $2 }' >> $site.log
