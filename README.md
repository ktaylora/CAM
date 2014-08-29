# A Quick Walkthrough to the CAM
CAM is an interface for interpretting output from Rsoilwat, so it requires a local copy of a SQLite database file
containing requisite climate data in the working directory.
<pre>
# callable via something like: 
for s in `cat drySites.lst`; 
  do bash cam-soilwat-unix-wrapper.sh $s & 
done
</pre>
