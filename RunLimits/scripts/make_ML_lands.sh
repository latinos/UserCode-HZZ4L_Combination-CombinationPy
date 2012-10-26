#!/bin/bash
STRICT=0; if [[ "$1" == "-s" ]]; then STRICT=1; shift; fi;
MINIM=0; if [[ "$1" == "-1" ]]; then MINIM="--minimizerAlgo=Minuit"; shift; fi;
LD="melaLD"; if [[ "$1" == "--mekd" ]]; then LD="mekd"; shift; fi;
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
    NAM=$(echo $1 | sed -e s/comb_// -e s/.txt//   | tr [a-z] [A-Z])
    if [[ "$MATCH" == "" || "$MATCH" == "$1" ]]; then
        if test -f $1; then
	    #combine -M MaxLikelihoodFit $* -n ${NAM}_${WHAT} -m $MASS $OPTIONS --out . 2>&1 | tee ${1/.root/.log.$WHAT} 
	    $CMSSW_BASE/src/LandS/test/lands.exe -d $* -M MaxLikelihoodFit -m $MASS -L $CMSSW_BASE/lib/*/*so \
		--RebinObservables CMS_zz4l_mass 40 0 0 $LD 20 0 0 -rMin 0 -rMax 10 -n landsCombine${NAM}_${WHAT}.obs 2>&1 | tee ${1/.txt/.log.$WHAT.lands}
	fi;
    fi;
}


if [[ "$MATCH" == "" ]]; then
    run $WHAT comb_hzz4l.root
else
    run $WHAT $MATCH
fi

