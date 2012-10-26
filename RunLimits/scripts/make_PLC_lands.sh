#!/bin/bash
STRICT=0; if [[ "$1" == "-s" ]]; then STRICT=1; shift; fi;

WHAT="PLC"
if [[ "$1" == "-S" ]]; then WHAT="PLS"; shift; fi;
if [[ "$1" == "-P" ]]; then WHAT="PLP"; shift; fi;
if [[ "$1" == "--PE" ]]; then WHAT="PLPE"; shift; fi;
if [[ "$1" == "--SE" ]]; then WHAT="PLSE"; shift; fi;
if [[ "$1" == "--mekd" ]]; then LD="mekd"; shift; else LD="melaLD"; fi
if test -d $1; then MASS=$1; else echo "Usage: $0 mass [what ]"; exit 1; fi; 
cd $MASS; shift
MATCH=$1;

OPTIONS=""
if [[ "$STRICT" == 1 ]]; then
    OPTIONS="--minimizerTolerance=0.0001"
fi;

LM="300"
MM="600"
HM="1000"

RMAX="5"
RMAXHM="20"
RMAXMM="10"
RMAXLM="1"

if [[ "$MASS" < "$HM" ]]; then RMAX=${RMAXHM}; fi
if [[ "$MASS" < "$MM" ]]; then RMAX=${RMAXMM}; fi
if [[ "$MASS" < "$LM" ]]; then RMAX=${RMAXLM}; fi

echo "RMAX = $RMAX"

if [[ "$WHAT" == "PLS" ]]; then
    OPTIONS="$OPTIONS --significance 1 -m $MASS -rMin 0 -rMax ${RMAX}"
elif [[ "$WHAT" == "PLSE" ]]; then
    OPTIONS="$OPTIONS -D asimov_sb --significance 1 -m $MASS -rMin 0 -rMax ${RMAX}"
elif [[ "$WHAT" == "PLP" ]]; then
    OPTIONS="$OPTIONS --significance 1 -m $MASS -rMin 0 -rMax ${RMAX}"
elif [[ "$WHAT" == "PLPE" ]]; then
    OPTIONS="$OPTIONS -D asimov_sb --significance 1 -m $MASS -rMin 0 -rMax ${RMAX}"
fi;

if [[ "$1" != "" ]] && test -f $1; then
    NAM=$(echo $1 | sed -e s/comb_*// -e s/.txt//   | tr '[a-z]' '[A-Z]')
    if [[ "$UPD" == "1" ]]; then test ${1/.txt/.log.$WHAT.lands} -nt $1 && return; fi;
    [[ "$COMBINE_NO_LOGFILES" != "1" ]] && DO_LOG="tee -a ${1/.txt/.log}.$WHAT.lands$POST" || DO_LOG="dd of=/dev/null"
    $CMSSW_BASE/src/LandS/test/lands.exe -d $* -M ProfileLikelihood $OPTIONS -L $CMSSW_BASE/lib/*/*so -n landsCombine${NAM}_${WHAT}.obsexp --RebinObservables CMS_zz4l_mass 40 0 0 $LD 20 0 0 2>&1 | $DO_LOG
    echo "Done $WHAT for $NAM at $MASS"
else
    echo "Missing workspace $1 at mass $MASS"; exit 1; 
fi;

