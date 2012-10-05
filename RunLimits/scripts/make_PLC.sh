#!/bin/bash
STRICT=0; if [[ "$1" == "-s" ]]; then STRICT=1; shift; fi;

WHAT="PLC"
if [[ "$1" == "-S" ]]; then WHAT="PLS"; shift; fi;
if [[ "$1" == "-P" ]]; then WHAT="PLP"; shift; fi;
if [[ "$1" == "--PE" ]]; then WHAT="PLPE"; shift; fi;
if test -d $1; then MASS=$1; else echo "Usage: $0 mass [what ]"; exit 1; fi; 
cd $MASS; shift
MATCH=$1;

OPTIONS=""
if [[ "$STRICT" == 1 ]]; then
    OPTIONS="--minimizerTolerance=0.0001"
fi;

if [[ "$WHAT" == "PLS" ]]; then
    OPTIONS="$OPTIONS --signif"
elif [[ "$WHAT" == "PLP" ]]; then
    OPTIONS="$OPTIONS --signif --pvalue"
elif [[ "$WHAT" == "PLPE" ]]; then
    OPTIONS="$OPTIONS --signif --pvalue --expectSignal=1 -t -1 --toysFreq"
fi;

if [[ "$1" != "" ]] && test -f $1; then
    NAM=$(echo $1 | sed -e s/comb_*// -e s/.root//   | tr '[a-z]' '[A-Z]')
    if [[ "$UPD" == "1" ]]; then test ${1/.root/.log.$WHAT} -nt $1 && return; fi;
    [[ "$COMBINE_NO_LOGFILES" != "1" ]] && DO_LOG="tee -a ${1/.root/.log}.$WHAT$POST" || DO_LOG="dd of=/dev/null"
    combine $* -n ${NAM}_${WHAT}${POST} -m $MASS $OPTIONS 2>&1 | $DO_LOG
    echo "Done $WHAT for $NAM at $MASS"
else
    echo "Missing workspace $1 at mass $MASS"; exit 1; 
fi;

