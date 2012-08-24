#!/bin/bash
GO=0; if [[ "$1" == "--go" ]]; then GO=1; shift; fi;
LSF=0; if [[ "$1" == "--lsf" ]]; then LSF=1; shift; fi;
PRI=0; if [[ "$1" == "--pri" ]]; then PRI=1; shift; fi;
GS=0; if [[ "$1" == "--gs" ]]; then GS=1; shift; fi;
CBL="cern"; if [[ "$1" == "--cbl" ]]; then CBL="$2"; shift; shift; fi;
PFIX=0; 
if   [[ "$1" == "-2" ]]; then PFIX=".2"; shift; 
elif [[ "$1" == "-3" ]]; then PFIX=".3"; shift; 
elif [[ "$1" == "-4" ]]; then PFIX=".4"; shift; 
elif [[ "$1" == "-5" ]]; then PFIX=".5"; shift; 
elif [[ "$1" == "-6" ]]; then PFIX=".6"; shift; 
elif [[ "$1" == "-7" ]]; then PFIX=".7"; shift; 
fi;
NAME="wide"; 
if   [[ "$1" == "--SM"  ]]; then NAME="SM"; shift; 
elif [[ "$1" == "--exp" ]]; then NAME="exp"; shift;
elif [[ "$1" == "--low" ]]; then NAME="low"; shift;
fi;
FGRID=""; if [[ "$1" == "-F" ]]; then FGRID="$2"; shift; shift; fi;
MASS=$1; if [[ "$2" == "" ]]; then echo "Usage: $0 mass workspace"; exit 1; fi;
if test \! -d $MASS; then echo "Mass $MASS not valid"; fi;
cd $MASS || exit 1
if test \! -f $2 ; then echo "Workspace $2 not found for $MASS"; exit 2; fi;


if [[ "$NAME" == "wide" ]]; then
    TRICE_PLC_LIMIT=$(awk '/Limit: r </{print 6*($4  < 1 ? 1 : $4)}'   ${2/.root/.log.PLC} | tail -n 1);
    FRAC_PLC_LIMIT=$(awk '/Limit: r </{print 0.2*$4}' ${2/.root/.log.PLC} | tail -n 1);
    if echo $2 | grep SM4; then
        TRICE_PLC_LIMIT=$(awk '/Limit: r </{print 6*($4  < 0.1 ? 0.1 : $4)}'   ${2/.root/.log.PLC} | tail -n 1);
    fi;
elif [[ "$NAME" == "exp" ]]; then
    TRICE_PLC_LIMIT=$(awk '/Limit: r </{print 4*$4}'   ${2/.root/.log.FREQ_E050} | tail -n 1);
    FRAC_PLC_LIMIT=$(awk '/Limit: r </{print 0.2*$4}' ${2/.root/.log.FREQ_E050} | tail -n 1);
elif [[ "$NAME" == "low" ]]; then
    TRICE_PLC_LIMIT=$(awk '/Limit: r </{print $4}'   ${2/.root/.log.FREQ_E016} | tail -n 1);
    FRAC_PLC_LIMIT=$(awk '/Limit: r </{print 0.3*$4}' ${2/.root/.log.FREQ_E016} | tail -n 1);
elif [[ "$NAME" == "SM" ]]; then
    TRICE_PLC_LIMIT=1; FRAC_PLC_LIMIT=1;
fi;
if [[ "$FGRID" != "" ]]; then
    TRICE_PLC_LIMIT=$(echo $FGRID | awk -F: '{print $2}');
    FRAC_PLC_LIMIT=$(echo $FGRID | awk -F: '{print $1}');
fi;

POINTS=30; INTERLEAVE=1
TOYS=50;   MULTIPLIER=2   # each job does ($POINTS/$INTERLEAVE) * $TOYS * $MULTIPLIER toys
TOTAL_TOYS=2500

if (( $MASS >= 226 )); then
    POINTS=30;
    TOTAL_TOYS=2500; MULTIPLIER=4;
elif (( $MASS >= 141 )); then
    POINTS=30;
    TOTAL_TOYS=2500; MULTIPLIER=4;
fi;

if [[ "$NAME" == "SM" ]]; then
    MULTIPLIER=$(( $MULTIPLIER * 2 ));
    POINTS=$(( $POINTS / 2 ));
elif [[ "$NAME" == "exp" ]]; then
    MULTIPLIER=$(( $MULTIPLIER * 2 ));
    TOTAL_TOYS=5000;
elif [[ "$NAME" == "low" ]]; then
    POINTS=5;
    MULTIPLIER=$(( $MULTIPLIER * 4 ));
    TOTAL_TOYS=5000;
fi;


OPTIONS=""
if [[ "$LSF" == "1" ]]; then
    OPTIONS="$OPTIONS --lsf"
    TOYS=50; MULTIPLIER=1;
    TOTAL_TOYS=2000;
    if (( $MASS >= 226 )); then
        POINTS=10
        TOTAL_TOYS=2500; MULTIPLIER=2;
    elif (( $MASS >= 141 )); then
        POINTS=20;
        TOTAL_TOYS=2500; MULTIPLIER=2;
    fi;
    if [[ "$NAME" == "low" ]]; then
        POINTS=5;
        MULTIPLIER=1;
        TOTAL_TOYS=5000;
        if (( $MASS >= 141 )); then 
            TOTAL_TOYS=10000; MULTIPLIER=8;
        fi;
    fi;
    if [[ "$NAME" == "SM" ]]; then
        TOTAL_TOYS=5000;
        POINTS=5; MULTIPLIER=1;
        if (( $MASS >= 226 )); then
            TOTAL_TOYS=5000; MULTIPLIER=4;
        elif (( $MASS >= 141 )); then
            TOTAL_TOYS=10000; MULTIPLIER=8;
        fi;
    fi;
fi

CRABOPT=""
if [[ "$CBL" != "none" ]]; then
    CRABOPT="$CRABOPT -GRID.ce_black_list=$CBL"
fi;

if [[ "$PRI" == "1" ]]; then OPTIONS="$OPTIONS -P"; fi
if [[ "$GS"  == "1" ]]; then OPTIONS="--glidein --server"; fi;

if [[ "$FGRID" == "" ]]; then
    OPTIONS="$OPTIONS --log"
fi;

JOBS=$(( $TOTAL_TOYS/$TOYS/$MULTIPLIER * $INTERLEAVE))
TPJ=$(( $TOYS*$POINTS*$MULTIPLIER/$INTERLEAVE))
echo "Will run $JOBS jobs each with $POINTS/$INTERLEAVE points, $TOYS*$MULTIPLIER toys each, $TPJ total toys per job, $TOTAL_TOYS total toys per point"
echo "Band to probe: [ $FRAC_PLC_LIMIT , $TRICE_PLC_LIMIT ]"

$CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/test/makeGridUsingCrab.py \
    $2 $FRAC_PLC_LIMIT $TRICE_PLC_LIMIT -n $POINTS -I $INTERLEAVE \
    -T $TOYS -j $JOBS -t $(( $JOBS * $MULTIPLIER )) \
    -o grid-${2/.root/-${NAME}} $OPTIONS \
    -O "--freq --optimizeSim=1 --newToyMC=1 -m $MASS" -r
if [[ "$GO" == "1" ]]; then crab -cfg grid-${2/.root/-${NAME}}.cfg -USER.ui_working_dir=crab_0_${2/.root/_${NAME}}$PFIX $CRABOPT -create -submit; fi

