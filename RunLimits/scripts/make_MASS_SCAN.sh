#!/bin/bash
WHAT="MASS_SCAN"

OPT=""; POST=""
case $1 in
--1d)
    WHAT="${WHAT}_1D"; 
    OPT="${OPT} --algo=grid --points=200 --rMax 10 -P MH --floatOtherPOI=1"; 
    ;;
--2d)
    WHAT="${WHAT}_2D"; 
    OPT="${OPT} --algo=grid --points=10000 --rMax 5"; 
    ;;
--c68)
    WHAT="${WHAT}_2D_C68"; 
    OPT="${OPT} --algo=contour2d --points=40 --rMax 5 --cl 0.68"; 
    ;;
--c95)
    WHAT="${WHAT}_2D_C95"; 
    OPT="${OPT} --algo=contour2d --points=40 --rMax 5 --cl 0.95"; 
    ;;
esac;
shift; 

if [[ "$1" == "--fast" ]]; then 
    WHAT=$(echo $WHAT | sed 's/\(_[12]D\)/&_FAST/'); 
    OPT="${OPT} --fastScan"; 
    shift; 
fi;

STEPPING=200
if [[ "$1" == "--more" ]]; then
    shift;
    if echo $WHAT | grep -q MASS_SCAN_2D; then
        OPT="${OPT/--points=10000/--points=250000}";
        STEPPING=5000;
    elif echo $WHAT | grep -q MASS_SCAN_1D; then
        OPT="${OPT/--points=200/--points=1000}";
    fi;
elif [[ "$1" == "--less" ]]; then
    shift;
    if echo $WHAT | grep -q MASS_SCAN_2D; then
        OPT="${OPT/--points=10000/--points=400}";
        STEPPING=40;
    elif echo $WHAT | grep -q MASS_SCAN_1D; then
        OPT="${OPT/--points=200/--points=50}";
    fi;
fi;

if [[ "$1" == "-i" ]] || [[ "$1" == "-I" ]]; then
    [[ "$1" == "-I" ]] && STEPPING=$(( 5 * $STEPPING ));
    OPT="${OPT} --firstPoint $(( $2*$STEPPING )) --lastPoint $(( ($2+1)*$STEPPING ))"
    POST="${POST}.$2";
    shift; shift;
fi;


if test -d $1; then MASS=$1; else echo "Usage: $0 mass [what ]"; exit 1; fi; 
cd $MASS; shift
MATCH=$1;

WORKSPACE="FloatMass_$1"
if echo $WHAT | grep -q 1D && echo "[$1]" | grep -q comb_hires; then
    WORKSPACE="FloatMassMu_$1";
fi;

MYWHAT=$WHAT${POST}; 
if [[ "$1" != "" ]] && test -f $WORKSPACE; then

    NAM=$(echo $1 | sed -e s/comb_*// -e s/.root//   | tr '[a-z]' '[A-Z]')
    if [[ "$NAM" == "QQH_HGG" ]]; then OPT="${OPT/--rMax 5/--rMax 10}"; fi
    
    echo "Runnining $WHAT for $NAM at $MASS. ";
    [[ "$COMBINE_NO_LOGFILES" != "1" ]] && DO_LOG="tee ${WORKSPACE/.root/.log}.$MYWHAT" || DO_LOG="dd of=/dev/null" 
    ( echo c -M MultiDimFit $WORKSPACE -n ${NAM}_${MYWHAT} -m $MASS $OPT  &&  \
      combine -M MultiDimFit $WORKSPACE -n ${NAM}_${MYWHAT} -m $MASS $OPT 2>&1 ) | $DO_LOG
else
    echo "Missing workspace $WORKSPACE at mass $MASS"; exit 1; 
fi;
