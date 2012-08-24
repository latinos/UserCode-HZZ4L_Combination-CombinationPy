#!/bin/bash
TEST="LHC"; 
if [[ "$1" == "--lep" ]]; then TEST="LEP"; shift; fi;
if [[ "$1" == "--tev" ]]; then TEST="TEV"; shift; fi;
if test -d $1; then MASS=$1; else echo "Usage: $0 mass [what ]"; exit 1; fi; 
cd $MASS
MATCH=$2;

OPTIONS="-M HybridNew --optimizeSim=1 --singlePoint 1 --onlyTestStat --testStat=$TEST"

function run {
    WHAT=$1; shift
    NAM=$(echo $1 | sed -e s/comb_// -e s/.root//   | tr [a-z] [A-Z])
    BSP=${1/.root};
    if [[ "$MATCH" == "" || "$MATCH" == "$1" || "$MATCH" == "${NAM}" ]]; then
        if test -f $1; then
             combine $* -n ${NAM}_${WHAT} -m $MASS $OPTIONS 2>&1 | tee ${1/.root/.log.$WHAT} | tail -n 2
        fi;
    fi;
}
WHAT="Q$TEST";
if [[ "$RULE" == "CLsplusb" ]]; then WHAT="SMCLSB"; fi;
if [[ "$MATCH" == "" ]]; then
    run $WHAT comb_hgg.root
    run $WHAT comb_hww.root
    run $WHAT comb_htt.root
    run $WHAT comb_hzz4l.root
    run $WHAT comb_hzz2l2nu.root
    run $WHAT comb_hzz2l2q.root
    run $WHAT comb.root
    run $WHAT combs.root
else
    run $WHAT $MATCH
fi;
