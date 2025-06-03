#!/bin/bash

cd /exp/dune/app/users/wwang/

tar cvf to_grid.tar to_grid

source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup python v3_9_15
setup justin
justin get-token
# htgettoken -a htvaultprod.fnal.gov -i dune #$HTGETTOKENOPTS

USERF=$USER
FNALURL='https://fndcadoor.fnal.gov:2880/dune/scratch/users'

echo "to_grid.tar"
INPUT_TAR_DIR_LOCAL=`justin-cvmfs-upload to_grid.tar`
if [ $? -ne 0 ]; then
  echo "Could not upload. Exiting"
  exit 1
fi

sleep 90

ls -l $INPUT_TAR_DIR_LOCAL

#take a look if the tar file is in cvms for grid submission
#this function will take as first argument the name of your file
check_cvmfs() {
  a=2
  i=0
  while [ $i -le ${maxiter:-10} ]; do
    stat $1 > /dev/null 2>&1 
    if [ $? == 0 ]; then
      echo "Found in RCDS"
      break
    fi

    echo "iter $i"
    echo "File not found in RCDS"
    echo "sleeping $a"
    sleep $a
    a=$((2*a))
    i=$((i+1))
  done
}

#scope=${1}
#extra=${2:-""}
#echo "scope $scope"
#echo "extra $extra"

hour=$(date -u +"%H")
#convert to central (by hand because I couldn't get it to work)
hour=$(( $hour - 5 ))
echo "Hour: $hour"

today="$(date -u +"%Y%m%d")"

if [ "$hour" -ge 12 ]; then
  roundhour=12
  iter=1
else
  roundhour=00
  iter=0
fi

endhr=$(( $roundhour + 12 ))
begin="$(date -d '-60 day' -u +"%Y-%m-%d")T${roundhour}00"
end="$(date -d '-60 day' -u +"%Y-%m-%d")T${endhr}00"
endfind="$(date -d '-62 day' -u +"%Y%m%d")"

echo now: $(date)
echo begin: $begin
echo end: $end
echo endfind: $endfind

# query="files where created_timestamp >= ${begin} and created_timestamp < ${end} and core.file_type=detector and core.run_type=hd-protodune and core.data_tier=raw skip 80 ordered"

desc="1x2x6 semi-ana light sim Edinburgh (${begin} -- ${end})"
scope="usertests"
maxwall="$((3600*5))"
recopattern="*_reco_*.root:$FNALURL/$USERF"
rssmb=2000

echo "Checking cvmfs"
maxiter=10
check_cvmfs $INPUT_TAR_DIR_LOCAL
if [ $? -ne 0 ]; then
  echo "Could not find in cvmfs. Exiting"
  exit 1
fi

echo "Found in cvmfs"

#justin-test-jobscript --mql "$query" --jobscript $jobscript --env INPUT_TAR_DIR_LOCAL="$INPUT_TAR_DIR_LOCAL"

echo
echo "Will execute"
echo justin simple-workflow --monte-carlo 1 \
--jobscript $jobscript --rss-mb $rssmb --max-distance 30 --scope "usertests" \
--description "$desc" \
--refind-end-date $endfind \
--output-pattern $recopattern \
--env INPUT_TAR_DIR_LOCAL="$INPUT_TAR_DIR_LOCAL" \
--env NUM_EVENTS=300 \
--refind-interval-hours 1 --wall-seconds $maxwall \
--lifetime-days 90


justin simple-workflow --mql "$query" \
--jobscript $jobscript --rss-mb $rssmb --max-distance 30 --scope $scope \
--description "Light Simulation Edinburgh" \
--refind-end-date $endfind \
--output-pattern $recopattern \
--env INPUT_TAR_DIR_LOCAL="$INPUT_TAR_DIR_LOCAL" \
--env NUM_EVENTS=300 \
--refind-interval-hours 1 --wall-seconds $maxwall \
--lifetime-days 90

echo $?
