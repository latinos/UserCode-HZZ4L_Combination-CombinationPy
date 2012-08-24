#!/bin/bash
VERB=0; if [[ "$1" == "-v" ]]; then VERB=1; shift; fi;
if [[ "$1" != "" ]]; then WHAT=$1; else echo "Usage: $0 what [channel]"; exit 1; fi; 
CHANN=$2;

function run {
    NAM=$(echo $1 | tr '[a-z]' '[A-Z]')
    if [[ "$CHANN" == "" || "$CHANN" == "$1" ]]; then
        if [[ "$WHAT" != "BAYES" ]]; then
            if [[ "$VERB" == "1" ]]; then 
                ./hadd2 -f higgsCombine${NAM}_${WHAT}.root [0-9]*/higgsCombine${NAM}_${WHAT}.*.root
            else
                ./hadd2 -f higgsCombine${NAM}_${WHAT}.root [0-9]*/higgsCombine${NAM}_${WHAT}.*root > /dev/null 2>&1
            fi;
        else
            ./hadd2 -f higgsCombine${NAM}_${WHAT}.root  all/crab_0_bayes-*/res/*root
        fi;
        echo higgsCombine${NAM}_${WHAT}.root
    fi;
}

if [[ "$CHANN" == "" ]]; then 
#    run hgg
#    run hww
#    run htt
    run hzz4l
#    run hzz2l2nu
#    run hzz2l2q
#    run comb
else
    run $CHANN
fi;
