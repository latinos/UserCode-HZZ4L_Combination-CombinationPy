#!/bin/bash
STRICT=0; if [[ "$1" == "-s" ]]; then STRICT=1; shift; fi;
MINIM=0; if [[ "$1" == "-1" ]]; then MINIM="--minimizerAlgo=Minuit"; shift; fi;
if test -d $1; then MASS=$1; else echo "Usage: $0 mass [what ]"; exit 1; fi; 
cd $MASS
WHAT="ML"
MATCH=$2;

OPTIONS=" --optimizeSim=1 --minimizerStrategy=2"
#OPTIONS=" --optimizeSim=1 $MINIM "
if [[ "$STRICT" == 1 ]]; then
   OPTIONS=" --optimizeSim=1 --minimizerStrategy=2 --minimizerTolerance=0.00001" 
fi
function run {
    WHAT=$1; shift
    NAM=$(echo $1 | sed -e s/comb_// -e s/.root//   | tr [a-z] [A-Z])
    if [[ "$MATCH" == "" || "$MATCH" == "$1" ]]; then
        if test -f $1; then
             combine -M MaxLikelihoodFit $* -n ${NAM}_${WHAT} -m $MASS $OPTIONS --out . 2>&1 | tee ${1/.root/.log.$WHAT} 
        fi;
    fi;
}


if [[ "$MATCH" == "" ]]; then
#    run $WHAT comb_hgg.root
#    run $WHAT comb_hww.root
#    run $WHAT comb_htt.root
    run $WHAT comb_hzz4l.root
#    run $WHAT comb_hzz2l2nu.root
#    run $WHAT comb_hzz2l2q.root
#    run $WHAT comb.root
else
    run $WHAT $MATCH
fi

