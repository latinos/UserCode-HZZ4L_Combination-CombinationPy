#!/bin/bash
RULE="CLs"; if [[ "$1" == "--clsb" ]]; then RULE="CLsplusb"; shift; fi;
if test -d $1; then MASS=$1; else echo "Usage: $0 mass [what ]"; exit 1; fi; 
cd $MASS
MATCH=$2;

OPTIONS="-M HybridNew --freq --optimizeSim=1 --singlePoint 1 --rule=$RULE"

function run {
    WHAT=$1; shift
    NAM=$(echo $1 | sed -e s/comb_// -e s/.root//   | tr [a-z] [A-Z])
    BSP=${1/.root};
    if [[ "$MATCH" == "" || "$MATCH" == "$1" || "$MATCH" == "${NAM}" ]]; then
        if test -f $1; then
             if ls crab_0_${BSP}_SM*/res/*root > /dev/null 2>&1; then
                 echo "Merging CLs grids for $NAM at $MASS. ";
                 ../hadd2 -f grid-$BSP.root crab_0_${BSP}_SM*/res/*root > /dev/null 2>&1;
             fi;
             if test -f grid-$BSP.root; then 
                 echo "Runnining CLs for SM Higgs for $NAM at $MASS. ";
                 GRID="--readHybridResult --toysFile=grid-$BSP.root"
                 combine $* -n ${NAM}_${WHAT} -m $MASS $GRID $OPTIONS > ${1/.root/.log.$WHAT} 2>&1 
                 for E in 50 16 025 84 975; do
                    combine $* -n ${NAM}_${WHAT} -m $MASS $GRID --expectedFromGrid 0.$E $OPTIONS  > ${1/.root/.log.$WHAT}_E0$E 2>&1
                 done;
                 grep "^$RULE =" ${1/.root/.log.$WHAT}* | sed 's/:/\t/';
             else
                 echo "No grid ready for $NAM at $MASS";
             fi;
        fi;
    fi;
}
WHAT="SMCLS";
if [[ "$RULE" == "CLsplusb" ]]; then WHAT="SMCLSB"; fi;
if [[ "$MATCH" == "" ]]; then
    run $WHAT comb_hgg.root
    run $WHAT comb_hww.root
    run $WHAT comb_htt.root
    run $WHAT comb_hzz4l.root
    run $WHAT comb_hzz2l2nu.root
    run $WHAT comb_hzz2l2q.root
    run $WHAT comb.root
else
    run $WHAT $MATCH
fi;
