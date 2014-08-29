# callable via something like: 
for s in `cat drySites.lst`; 
  do bash cam-soilwat-unix-wrapper.sh $s & 
done
