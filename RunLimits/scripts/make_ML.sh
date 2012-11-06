#!/bin/bash

POST=""; POPT="";
if [[ "$1" == "-q" ]]; then POST=".Q"; fi
if [[ "$1" == "-1" ]]; then POST=".1"; fi

FIT="--robustFit=1"; if [[ "$1" == "-q" ]]; then FIT=""; shift; fi;
STRICT=0; if [[ "$1" == "-s" ]]; then STRICT=1; shift; elif [[ "$1" == "-l" ]]; then STRICT=-1; shift; elif [[ "$1" == "-h" ]]; then STRICT=-0.5; shift; fi;
STEP="0.05"; if [[ "$1" == "-t" ]]; then STEP="0.01"; shift; fi;
MINIM=""; if [[ "$1" == "-1" ]]; then MINIM="--minimizerAlgo=Minuit --minimizerAlgoForMinos=SeqMinimizer"; shift; fi;
EXTRA="0"; if [[ "$1" == "-x" ]]; then EXTRA=1; shift; fi;
WHAT="ML"
if [[ "$1" == "--TS" ]]; then WHAT="MLTS"; shift; fi;
if [[ "$1" == "--ES" ]]; then WHAT="MLES"; shift; fi;
if [[ "$1" == "--TB" ]]; then WHAT="MLTB"; shift; fi;

RANGE="auto";
if [[ "$1" == "-Z" ]];   then WHAT="${WHAT}Z";                shift; fi;
if [[ "$1" == "--ZM" ]]; then WHAT="${WHAT}Z"; RANGE="input"; shift; fi;
ZERO="--preFitValue=0.01"; 

MINIM="$MINIM --stepSize=$STEP"

if test -d $1; then MASS=$1; else echo "Usage: $0 mass [what ]"; exit 1; fi; 
cd $MASS; shift;
MATCH=$1;

OPTIONS="$MINIM $ZERO"
if [[ "$STRICT" == 1 ]]; then
   OPTIONS=" --minimizerStrategy=2 --minimizerTolerance=0.0001 $MINIM $ZERO" 
elif [[ "$STRICT" == -1 ]]; then
   OPTIONS=" --minimizerStrategy=0 --minimizerTolerance=0.1 $MINIM $ZERO" 
elif [[ "$STRICT" == -0.5 ]]; then
   OPTIONS=" --minimizerStrategy=0 --minimizerTolerance=0.01 $MINIM $ZERO" 
elif [[ "$POPT" != "" ]]; then
   OPTIONS="$MINIM $ZERO $POPT";
else 
   #OPTIONS="$MINIM $ZERO $FIT --X-rtd FITTER_DYN_STEP  --cminFallbackAlgo Minuit,0.001 " # this seems good at the moment
   OPTIONS="$MINIM $ZERO $FIT --X-rtd FITTER_NEW_CROSSING_ALGO --minimizerAlgoForMinos=Minuit2 --minimizerToleranceForMinos=0.01 --X-rtd FITTER_NEVER_GIVE_UP --minimizerAlgo=Minuit2 --minimizerStrategy=1 --minimizerTolerance=0.1 --cminFallbackAlgo Minuit,0.001 -v 100" # this seems good at the moment
#   OPTIONS="$MINIM $ZERO $FIT --X-rtd FITTER_NEW_CROSSING_ALGO --minimizerAlgoForMinos=Minuit2 --minimizerToleranceForMinos=0.01 --X-rtd FITTER_NEVER_GIVE_UP --X-rtd FITTER_BOUND --minimizerAlgo=Minuit2 --minimizerStrategy=1 --minimizerTolerance=0.001 --cminFallbackAlgo Minuit,0.001 " # this seems good at the moment
   #OPTIONS="$MINIM $ZERO $FIT"
fi
if [[ "$EXTRA" == "0" ]]; then OPTIONS="$OPTIONS --justFit"; fi;
function run {
    WHAT=$1; shift
    NAM=$(echo $1 | sed -e s/comb_*// -e s/.root//   | tr '[a-z]' '[A-Z]')
    OPTS="$OPTIONS"
    if [[ "$RANGE" == "auto" ]]; then
        MIN=0; MAX=3; if echo "$WHAT" | grep -q "Z"; then MIN=-2; fi;
        OPTS="$OPTS  --rMin=$MIN --rMax=$MAX";
    fi;
    if [[ "$MATCH" == "" || "$MATCH" == "$1" ]]; then
        if test -f $1; then
             test -f ${1/.root/.log}.$WHAT$POST && rm ${1/.root/.log}.$WHAT$POST;
             [[ "$COMBINE_NO_LOGFILES" != "1" ]] && DO_LOG="tee -a ${1/.root/.log}.$WHAT$POST" || DO_LOG="dd of=/dev/null" 
             echo "c -M MaxLikelihoodFit $* -n ${NAM}_${WHAT}${POST} -m $MASS $OPTS "     | $DO_LOG; 
	     combine -M MaxLikelihoodFit $* -n ${NAM}_${WHAT}${POST} -m $MASS $OPTS  2>&1 | $DO_LOG; 
        fi;
    fi;
}


run $WHAT $*

