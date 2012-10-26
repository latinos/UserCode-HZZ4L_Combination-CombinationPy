#!/bin/bash
WHAT="ASCLS"

POST=""; POPT="";
if [[ "$1" == "-k"     ]]; then POST=".K";    POPT=" --X-rtd ADDNLL_KAHAN_SUM --X-rtd SIMNLL_KAHAN_SUM "; shift; fi;
if [[ "$1" == "--m2s1" ]]; then POST=".M2S1"; POPT=" --minimizerStrategy=1 --minimizerTolerance=0.01"; shift; fi;
if [[ "$1" == "--m2t2" ]]; then POST=".M2T2"; POPT=" --minimizerTolerance=0.01 "; shift; fi;
if [[ "$1" == "--m2t4" ]]; then POST=".M2T4"; POPT=" --minimizerTolerance=0.0001 "; shift; fi;
if [[ "$1" == "--m1"   ]]; then POST=".M1";   POPT=" --minimizerAlgo=Minuit "; shift; fi;
if [[ "$1" == "--g1"   ]]; then POST=".G1";   POPT=" --minimizerTolerance=0.01 --X-rtd ADDNLL_KAHAN_SUM --X-rtd SIMNLL_KAHAN_SUM "; shift; fi;
if [[ "$1" == "--mekd" ]]; then LD="mekd"; shift; else LD="melaLD"; fi;
CLOPT="";
if [[ "$1" == "--90" ]]; then WHAT="ASCLS90"; CLOPT="--cl 0.90"; shift; fi;
if [[ "$1" == "--99" ]]; then WHAT="ASCLS99"; CLOPT="--cl 0.99"; shift; fi;

STRICT=0; if [[ "$1" == "-s" ]]; then STRICT=1; shift; fi;
if [[ "$1" == "-l" ]]; then STRICT=-1; shift; fi;
UPD=0; if [[ "$1" == "-u" ]]; then UPD=1; shift; fi;
MINIM=""; if [[ "$1" == "-1" ]]; then MINIM="--minimizerAlgo=Minuit"; shift; fi;
TASK="both"; if [[ "$1" == "-O" ]]; then TASK="observed"; shift; elif [[ "$1" == "-E" ]]; then TASK="expected"; shift; fi;
if [[ "$1" == "--cls" ]]; then TASK="singlePoint"; shift; fi;
if test -d $1; then MASS=$1; else echo "Usage: $0 mass [what ]"; exit 1; fi; 
cd $MASS; shift
MATCH=$1;

function run {
    MYWHAT=$1${POST}; shift
    NAM=$(echo $1 | sed -e s/comb_// -e s/.txt//   | tr '[a-z]' '[A-Z]')
    BSP=${1/.root};
    OPTS="$OPTIONS"
#    if [[ "$TASK" == "observed" ]]; then MYWHAT="${MYWHAT}.O"; fi;
#    if [[ "$TASK" == "expected" ]]; then MYWHAT="${MYWHAT}.E"; fi;
#    if [[ "$TASK" == "singlePoint" ]]; then MYWHAT="SM${MYWHAT}"; fi;
    if [[ "$UPD" == "1" ]]; then test ${1/.txt/.log.$MYWHAT.lands} -nt $1 && return; fi;
    if test -f $1; then
         echo "Runnining $TASK asymptotic frequentist limits for $NAM at $MASS. ";
         [[ "$COMBINE_NO_LOGFILES" != "1" ]] && DO_LOG="tee ${1/.txt/.log}.$MYWHAT.lands" || DO_LOG="dd of=/dev/null" 
	 $CMSSW_BASE/src/LandS/test/lands.exe -d $* -M Asymptotic -m $MASS -L $CMSSW_BASE/lib/*/*so \
	     --RebinObservables CMS_zz4l_mass 40 0 0 $LD 20 0 0 -rMin 0 -rMax 10 -n landsCombine${NAM}_${MYWHAT}O.obs  2>&1 | $DO_LOG
	 $CMSSW_BASE/src/LandS/test/lands.exe -d $* -M Asymptotic -m $MASS -L $CMSSW_BASE/lib/*/*so \
	     --RebinObservables CMS_zz4l_mass 40 0 0 $LD 20 0 0 -rMin 0 -rMax 10 -D asimov_b -n landsCombine${NAM}_${MYWHAT}E.exp 2>&1 | $DO_LOG
    fi;
}
if [[ "$MATCH" == "" ]]; then
    for W in $(cat workspaces.txt); do
        run $WHAT $W.root;
    done
else
    run $WHAT $*;
fi;
