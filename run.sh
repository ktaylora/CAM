for s in `cat us_*`; do
  while [[ `ps aux | grep "R --no-save" | grep -v "grep" | wc -l` -ge 10 ]]; do
    sleep 5;
  done
  sleep 1;
  if [[ `ps aux | grep "R --no-save" | grep -v "grep" | wc -l` -le 9 ]]; then
    nohup bash ./cam-soilwat-unix-wrapper.sh $s usClimateDb_1971_1998.sqlite >> /dev/null &
  fi
  sleep 1;
done
