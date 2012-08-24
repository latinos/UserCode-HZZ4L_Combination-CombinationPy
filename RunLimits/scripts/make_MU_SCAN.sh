#!/bin/bash
WHAT="MU_SCAN"

OPT=""; POST=""
OPT="${OPT} --algo=grid -v -1 --points=100"; 

if [[ "$1" == "--fast" ]]; then 
    WHAT="MU_SCAN_FAST"; 
    OPT="${OPT} --fastScan"; 
    shift; 
fi;

if [[ "$1" == "-i" ]]; then
    OPT="${OPT} --firstPoint $(( $2*40 )) --lastPoint $(( $2*40 + 40 ))"
    POST="${POST}.$2";
    shift; shift;
fi;
if [[ "$1" == "-I" ]]; then
    OPT="${OPT} --firstPoint $(( $2*10 )) --lastPoint $(( $2*10 + 10 ))"
    POST="${POST}.$2";
    shift; shift;
fi;


if test -d $1; then MASS=$1; else echo "Usage: $0 mass [what ]"; exit 1; fi; 
cd $MASS; shift
NAM=$(echo $1 | sed -e s/comb_*// -e s/.root//   | tr '[a-z]' '[A-Z]')
WORKSPACE=$1; shift

MIN=0; MAX=2;
BNAM=$(echo $NAM | sed 's/^[78]_//');
case $BNAM in
HBB|VH_HBB|HTT|QQH_HGG|GGH_HTT)
    MIN=-3; MAX=7;;
HZZ|HWW|GGH_HWW|GGH_HGG|HGG|QQH|VH)
    MIN=-1; MAX=3;;   
*TTH*)
    MIN=-5; MAX=5;;
*VH*|*TTH*|*QQH*)
    MIN=-10; MAX=10;;
esac;
if [[ "$2" != "" && "$1" != "" ]]; then
    MIN="$1"; MAX="$2"; shift; shift;
fi;

OPT="${OPT} --rMin $MIN --rMax $MAX"

MYWHAT=$WHAT${POST}; 
if [[ "$WORKSPACE" != "" ]] && test -f $WORKSPACE; then
    echo "Runnining $WHAT for $NAM at $MASS. ";
    [[ "$COMBINE_NO_LOGFILES" != "1" ]] && DO_LOG="tee ${WORKSPACE/.root/.log}.$MYWHAT" || DO_LOG="dd of=/dev/null" 
    ( echo  c -M MultiDimFit $WORKSPACE -n ${NAM}_${MYWHAT} -m $MASS $OPT $*  &&  #\
      combine -M MultiDimFit $WORKSPACE -n ${NAM}_${MYWHAT} -m $MASS $OPT $* 2>&1 ) | $DO_LOG
else
    echo "Missing workspace $WORKSPACE at mass $MASS"; exit 1; 
fi;
