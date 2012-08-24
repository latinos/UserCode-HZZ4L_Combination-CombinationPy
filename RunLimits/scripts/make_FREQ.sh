#!/bin/bash
FORCE=0; if [[ "$1" == "-f" ]]; then FORCE=1; shift; fi;
STRICT=0; if [[ "$1" == "-s" ]]; then STRICT=1; shift; fi;
if test -d $1; then MASS=$1; else echo "Usage: $0 mass [what ]"; exit 1; fi; 
cd $MASS
MATCH=$2;

OPTIONS="-M HybridNew --freq --optimizeSim=1 --rAbsAcc=0 --rRelAcc=0"
if [[ "$FORCE" == "1" ]]; then OPTIONS="$OPTIONS --fullGrid"; fi;

function run {
    WHAT=$1; shift
    NAM=$(echo $1 | sed -e s/comb_// -e s/.root//   | tr '[a-z]' '[A-Z]')
    BSP=${1/.root};
    if [[ "$MATCH" == "" || "$MATCH" == "$1" || "$MATCH" == "${NAM}" ]]; then
        if test -f $1; then
             if ls crab_0_${BSP}_[lweS]*/res/*root > /dev/null 2>&1; then
                 echo "Runnining frequentist limits for $NAM at $MASS. ";
                 ../hadd2 -f grid-$BSP.root crab_0_${BSP}_[lweS]*/res/*root > /dev/null 2>&1;
                 combine $* -n ${NAM}_${WHAT} -m $MASS --grid=grid-$BSP.root $OPTIONS > ${1/.root/.log.$WHAT} 2>&1 
                 if [[ "$FORCE" == 0 ]]; then
                     for E in 50 16 025 84 975; do
                        combine $* -n ${NAM}_${WHAT} -m $MASS --grid=grid-$BSP.root --expectedFromGrid 0.$E $OPTIONS  > ${1/.root/.log.$WHAT}_E0$E 2>&1
                     done;
                 fi;
                 grep '^Limit:.*CL' ${1/.root/.log.$WHAT}* | sed 's/:/\t/';
             else
                 echo "No grid ready for $NAM at $MASS";
             fi;
        fi;
    fi;
}
if [[ "$MATCH" == "" ]]; then
    run FREQ comb_hgg.root
    run FREQ comb_hww.root
    run FREQ comb_htt.root
    run FREQ comb_hzz4l.root
    run FREQ comb_hzz2l2nu.root
    run FREQ comb_hzz2l2q.root
    run FREQ comb.root
else
    run FREQ $MATCH
fi;
