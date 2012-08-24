#!/bin/bash

STEPS="200000"; if [[ "$1" == "-i" ]]; then STEPS="$2"; shift; shift; fi;
SEED="-1"; if [[ "$1" == "-s" ]]; then SEED="$2"; shift; shift; fi;

WHAT="BAYES"
POST=""
OPTS=""
WPREFIX=""

if test -d $1; then MASS=$1; else echo "Usage: $0 mass [what ]"; exit 1; fi; 
cd $MASS; shift;

WHO=$1;
NAM=$(echo $1 | sed -e s/comb_*// -e s/.root//   | tr '[a-z]' '[A-Z]' | tr '.' '_')
WORKSPACE="$1"

if test -f $WORKSPACE; then
     test -f ${WORKSPACE/.root/.log}.$WHAT$POST && rm ${WORKSPACE/.root/.log}.$WHAT$POST;
     [[ "$COMBINE_NO_LOGFILES" != "1" ]] && DO_LOG="tee -a ${WORKSPACE/.root/.log}.$WHAT$POST" || DO_LOG="dd of=/dev/null" 
     echo "c -M MarkovChainMC $WORKSPACE --tries 1 --saveChain -i $STEPS -m $MASS -s $SEED -n ${NAM}_${WHAT} $OPTS --propHelperWidthRangeDivisor=20 "    | $DO_LOG; 
     combine -M MarkovChainMC $WORKSPACE --tries 1 --saveChain -i $STEPS -m $MASS -s $SEED -n ${NAM}_${WHAT} $OPTS --propHelperWidthRangeDivisor=10 2>&1 | $DO_LOG;
else 
    echo "Missing workspace $WORKSPACE at $MASS";
fi;
