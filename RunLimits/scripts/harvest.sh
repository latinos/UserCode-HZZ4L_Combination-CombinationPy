#!/bin/bash
VERB=0; if [[ "$1" == "-v" ]]; then VERB=1; shift; fi;
PREFIX=""; REPREF=""; 
if   [[ "$1" == "-4" ]]; then REPREF="$1"; PREFIX="SM4_"; shift; 
elif [[ "$1" == "-F" ]]; then REPREF="$1"; PREFIX="FF_";  shift; fi;
if [[ "$1" != "" ]]; then WHAT=$1; else echo "Usage: $0 what [channel]"; exit 1; fi; 
if [[ "$WHAT" == "A" ]]; then
    shift; 
    for X in ASCLS PLP PLPE MLZ; do 
            bash harvest.sh $REPREF $X $*; 
    done
    exit
elif [[ "$WHAT" == "T" ]]; then
    shift;
    for X in FREQ SMCLS PVAL; do bash harvest.sh $REPREF $X $*; done
    exit
fi;
CHANN=$2;

function run {
    NAM=$(echo $1 | tr '[a-z]' '[A-Z]')
    if [[ "$CHANN" == "" || "$CHANN" == "$1" ]]; then
        FILES=$(ls -l [0-9]*/higgsCombine${NAM}_${WHAT}.*.root 2> /dev/null | awk '{if ($5 > 1000) print}' | wc -l);
        if [[ "$FILES" == "0" ]]; then return; fi;
        if [[ "$VERB" == "1" ]]; then 
            ./hadd2 -f results/higgsCombine${NAM}_${WHAT}.root [0-9]*/higgsCombine${NAM}_${WHAT}.*.root
        else
            ./hadd2 -f results/higgsCombine${NAM}_${WHAT}.root [0-9]*/higgsCombine${NAM}_${WHAT}.*.root > /dev/null 2>&1
        fi;
        echo "higgsCombine${NAM}_${WHAT}.root   ($FILES files)"
    fi;
}

if [[ "$CHANN" == "" ]]; then 
    for W0 in $(cat workspaces.txt); do W=${PREFIX}${W0}; run ${W/comb_/}; done
else
    run $CHANN
fi;
