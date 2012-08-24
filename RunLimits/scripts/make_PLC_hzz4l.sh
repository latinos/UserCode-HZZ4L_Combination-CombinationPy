#!/bin/bash
STRICT=0; if [[ "$1" == "-s" ]]; then STRICT=1; shift; fi;
WHAT="PLC"
if [[ "$1" == "-S" ]]; then WHAT="PLS"; shift; fi;
if [[ "$1" == "-P" ]]; then WHAT="PLP"; shift; fi;
if [[ "$1" == "-E" ]]; then WHAT="PLE"; shift; fi;
if [[ "$1" == "--PE" ]]; then WHAT="PLPE"; shift; fi;
if test -d $1; then MASS=$1; else echo "Usage: $0 mass [what ]"; exit 1; fi; 
cd $MASS
MATCH=$2;

OPTIONS=" --optimizeSim=1 --minimizerTolerance=0.0001"
#OPTIONS="--useMinos --optimizeSim=1"
if [[ "$STRICT" == 1 ]]; then
    #OPTIONS="$OPTIONS --maxTries=100 --tries=20 --maxRelDeviation 0.1 --maxOutlierFraction=0.3 --maxOutliers=7 --minimizerTolerance=0.0001 --minimizerAlgo=Minuit"
    #OPTIONS="--optimizeSim=1 --maxTries=100 --tries=20 --maxRelDeviation 0.1 --maxOutlierFraction=0.3 --maxOutliers=7 --minimizerTolerance=0.0001"
    OPTIONS="--optimizeSim=1 --maxTries=100 --tries=7 --maxRelDeviation 0.1 --maxOutlierFraction=0.3 --maxOutliers=3 --minimizerTolerance=0.0001"
fi;
if [[ "$WHAT" == "PLS" ]]; then
    OPTIONS="$OPTIONS --signif"
elif [[ "$WHAT" == "PLP" ]]; then
    OPTIONS="$OPTIONS --signif --pvalue"
elif [[ "$WHAT" == "PLPE" ]]; then
    OPTIONS="$OPTIONS --signif --pvalue --expectSignal=1 -t 10 -s -1"
elif [[ "$WHAT" == "PLE" ]]; then
    OPTIONS="$OPTIONS -t 10 -s -1"
fi;

function run {
    WHAT=$1; shift
    NAM=$(echo $1 | sed -e s/comb_// -e s/.root//   | tr [a-z] [A-Z])
    if [[ "$MATCH" == "" || "$MATCH" == "$1" ]]; then
        if test -f $1; then
             combine $* -n ${NAM}_${WHAT} -m $MASS $OPTIONS 2>&1 | tee ${1/.root/.log.$WHAT} 
        fi;
    fi;
    echo ""
    echo "================================================================================"
    echo "Done $WHAT for $NAM at $MASS"
}

if [[ "$MATCH" == "" ]]; then
    #run $WHAT comb_hgg.root
    #run $WHAT comb_hww.root
    #run $WHAT comb_htt.root
    run $WHAT comb_hzz4l.root
    #run $WHAT comb_hzz2l2nu.root
    #run $WHAT comb_hzz2l2q.root
    #run $WHAT comb.root
else
    run $WHAT $MATCH
fi

