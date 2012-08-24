#!/bin/bash
BAYES=0; if [[ "$1" == "-b" ]]; then BAYES=1; shift; fi;
MASS=$1; if [[ "$2" == "" ]]; then echo "Usage: $0 mass workspace"; exit 1; fi;
cd $MASS || exit 1
if test \! -f $2; then echo "Missing $MASS/$2"; exit 2; fi;

if [[ "$BAYES" == "0" ]]; then
    combine $2 -M HybridNew --freq -T 100 --fork 4 --singlePoint 1 --clsAcc 0 --optimizeSim=1 --newToyMC=1 -m $MASS 2>&1 | tee ${2/.root/}.log.TIME

    SECS_PER_TOY=$(tail -n 1 ${2/.root/}.log.TIME  | awk '{print $6*60/25}')
    TOYS_PER_HOUR=$(tail -n 1 ${2/.root/}.log.TIME  | awk '{print 60./($6/25)}')

    echo "  Time per toy:  $SECS_PER_TOY s"
    echo "  Toys per hour: $TOYS_PER_HOUR";
else
    combine $2 -M MarkovChainMC --tries 1 -i 50000 -m $MASS  --optimizeSim=1 --alwaysStepPOI=1 | tee ${2/.root/}.log.MCMC_TIME
    MINS_PER_RUN=$(tail -n 1 ${2/.root/}.log.MCMC_TIME  | awk '{print $6*2}')
    RUNS_PER_HOUR=$(tail -n 1 ${2/.root/}.log.MCMC_TIME  | awk '{print 60/($6*2)}')
    
    echo "  Time per run of 100k iters: $MINS_PER_RUN min (from run of 50k iterations)"
    echo "  Runs per hour (100k iters): $RUNS_PER_HOUR";
fi;
