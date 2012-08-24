#!/bin/bash

STEPS="10000"; if [[ "$1" == "-i" ]]; then STEPS="$2"; shift; shift; fi;
JOBS="20"; if [[ "$1" == "-j" ]]; then JOBS="$2"; shift; shift; fi;
INDEX="$1"; STRIDE=1
if [[ "$INDEX" == "" ]]; then echo "Usage: $0 index mass [what ]"; exit 1; fi;
if seq 0 $(( $JOBS - 1 )) | grep -q $INDEX; then
    STRIDE=$(( $STEPS / $JOBS ))
    echo "Will run job $INDEX of $JOBS, processing $STRIDE points.";
    shift; 
else
    echo "Usage: $0 index mass [what ]"; exit 1; 
fi;

WHAT="SCAN"
POST=".$INDEX"
OPT="--algo=grid --points=$STEPS --firstPoint=$(( $INDEX * $STRIDE)) --lastPoint=$(( ($INDEX+1)*$STRIDE - 1))"

if [[ "$1" == "--fast" ]]; then 
    WHAT=${WHAT}_FAST; 
    OPT="${OPT} --fastScan"; 
    shift; 
fi;


if test -d $1; then MASS=$1; else echo "Usage: $0 mass [what ]"; exit 1; fi; 
cd $MASS; shift;

WHO=$1;
NAM=$(echo $1 | sed -e s/comb_*// -e s/.root//   | tr '[a-z]' '[A-Z]' | tr '.' '_')
WORKSPACE="$1"

if test -f $WORKSPACE; then
     test -f ${WORKSPACE/.root/.log}.$WHAT$POST && rm ${WORKSPACE/.root/.log}.$WHAT$POST;
     [[ "$COMBINE_NO_LOGFILES" != "1" ]] && DO_LOG="tee -a ${WORKSPACE/.root/.log}.$WHAT$POST" || DO_LOG="dd of=/dev/null" 
     echo "c -M MultiDimFit $WORKSPACE -m $MASS -n ${NAM}_${WHAT}$POST $OPT "    | $DO_LOG; 
     combine -M MultiDimFit $WORKSPACE -m $MASS -n ${NAM}_${WHAT}$POST $OPT 2>&1 | $DO_LOG;
else 
    echo "Missing workspace $WORKSPACE at $MASS";
fi;
