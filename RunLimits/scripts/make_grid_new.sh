#!/bin/bash
GO=0; if [[ "$1" == "--go" ]]; then GO=1; shift; fi;
LSF=0; if [[ "$1" == "--lsf" ]]; then LSF=1; shift; fi;
PRI=0; if [[ "$1" == "--pri" ]]; then PRI=1; shift; fi;
GL=0; if [[ "$1" == "--glide" ]]; then GL=1; shift; fi;
GS=0; if [[ "$1" == "--gs" ]]; then GS=1; shift; fi;
CBL="$(perl -npe 'm/^\s*#/ and $_=""; s/\s*#.*//g;' scripts/cms/crab_blacklist.txt | perl -e '@a=<>;chomp(@a);$_=join(",",@a); s/\s+//g; print')";
CWL="none"; if [[ "$1" == "--cwl" ]]; then CWL="$2"; shift; shift; fi;
if [[ "$1" == "--fastcwl" ]]; then 
    CWL=$(perl -npe 'm/^\s*#/ and $_=""; s/\s*#.*//g;' scripts/cms/crab_whitelist.txt | perl -e '@a=<>;chomp(@a);$_=join(",",@a); s/\s+//g; print'); 
    shift;
fi;
PFIX=.0; 
if   [[ "$1" == "-1" ]]; then PFIX=".1"; shift; 
elif [[ "$1" == "-2" ]]; then PFIX=".2"; shift; 
elif [[ "$1" == "-3" ]]; then PFIX=".3"; shift; 
elif [[ "$1" == "-4" ]]; then PFIX=".4"; shift; 
elif [[ "$1" == "-5" ]]; then PFIX=".5"; shift; 
elif [[ "$1" == "-6" ]]; then PFIX=".6"; shift; 
elif [[ "$1" == "-7" ]]; then PFIX=".7"; shift; 
elif [[ "$1" == "-8" ]]; then PFIX=".8"; shift; 
elif [[ "$1" == "-9" ]]; then PFIX=".9"; shift; 
fi;
NAME="wide"; 
if   [[ "$1" == "--SM"  ]]; then NAME="SM"; shift; 
elif [[ "$1" == "--exp" ]]; then NAME="exp"; shift;
elif [[ "$1" == "--low" ]]; then NAME="low"; shift;
elif [[ "$1" == "--up" ]];  then NAME="up"; shift;
elif [[ "$1" == "--sig" ]]; then NAME="sig"; shift;
elif [[ "$1" == "--FC" ]];  then NAME="FC"; shift;
fi;
FGRID=""; FGP=5; 
if [[ "$1" == "-F" ]]; then FGRID="$2"; shift; shift; fi;
if [[ "$1" == "--F2" ]]; then FGRID="$2"; FGP=2; shift; shift; fi;

REMULT=""; if [[ "$1" == "-x" ]]; then REMULT=$2; shift; shift; fi;
ADAPTIVE=""; if [[ "$1" == "-A" ]]; then ADAPTIVE=1; shift; fi;
EXPECT=0; if [[ "$1" == "-E" ]]; then EXPECT=1; shift; fi;

MASS=$1; if [[ "$2" == "" ]]; then echo "Usage: $0 mass workspace"; exit 1; fi;
if test \! -d $MASS; then echo "Mass $MASS not valid"; fi;
cd $MASS || exit 1
if test \! -f $2 ; then echo "Workspace $2 not found for $MASS"; exit 2; fi;

OPTS="--freq -m $MASS"
if [[ "$NAME" == "FC"  ]]; then OPTS="$OPTS --testStat=LHCFC"; fi;
if [[ "$EXPECT"  == "1" ]]; then OPTS="$OPTS --fullBToys"; fi;
if echo $2 | grep -q atlas; then OPTS="$OPTS -w combWS -D combData --guessGenMode"; fi

if echo $MASS | grep -F -q .5; then IMASS=${MASS/.5/}; else IMASS=$MASS; fi

GRID=$(echo grids/${2/.root/.txt} | sed 's/comb_//');
if [[ "$NAME" == "wide" ]]; then
    TRICE_PLC_LIMIT=$(awk "/^$IMASS/{print \$4}" ../$GRID | head -n 1); # put head to cope with half-integer masses if they will arrive
    FRAC_PLC_LIMIT=$(awk "/^$IMASS/{print \$2}" ../$GRID | head -n 1); # put head to cope with half-integer masses if they will arrive
elif [[ "$NAME" == "up" ]]; then
    TRICE_PLC_LIMIT=$(awk "/^$IMASS/{print 2.0*\$4}" ../$GRID | head -n 1); # extend grid higher
    FRAC_PLC_LIMIT=$( awk "/^$IMASS/{print 1.4*\$4}" ../$GRID | head -n 1); # by a factor 2 
elif [[ "$NAME" == "FC" ]]; then
    TRICE_PLC_LIMIT=$(awk "/^$MASS/{print \$4}" ../$GRID | head -n 1); # extend grid higher
    FRAC_PLC_LIMIT=0;
elif [[ "$NAME" == "low" ]]; then
    TRICE_PLC_LIMIT=$(awk '/Limit: (r|mu) </{print $4}'   ${2/.root/.log.FREQ_E016} | tail -n 1);
    FRAC_PLC_LIMIT=$(awk '/Limit: (r|mu) </{print $4/2.5}' ${2/.root/.log.FREQ_E016} | tail -n 1);
elif [[ "$NAME" == "SM" ]]; then
    TRICE_PLC_LIMIT=1; FRAC_PLC_LIMIT=1;
elif [[ "$NAME" == "sig"  ]]; then
    TRICE_PLC_LIMIT=1; FRAC_PLC_LIMIT=1;
    SIG_ATOTAL=$(python -c "
import ROOT; from math import *; from sys import stderr
file = '${GRID}'.replace('grids/', '../grids/pvala_');
for line in open(file,'r'):
    if 'mass' in line or '---' in line: continue
    xs = line.split()
    if float(xs[0]) != $MASS: continue
    p = float(xs[1])
    z = ROOT.ROOT.Math.normal_quantile_c(p,1.0)
    n01z  = int(ceil(6.283 * p * exp(z**2)/(0.1**2)))
    n005z = int(ceil(6.283 * p * exp(z**2)/(0.05**2)))
    n01p  = int(ceil(100/p))
    stderr.write('\tp-value is %g, z = %g\n' % (p,z))
    stderr.write('\ttoys to get 0.1  absolute uncertainty on sigma: %8d\n' % n01z)
    stderr.write('\ttoys to get 0.05 absolute uncertainty on sigma: %8d\n' % n005z)
    stderr.write('\ttoys to get 0.1  relative uncertainty on p-val: %8d\n' % n01p)
    print n005z
" -b)
fi;
if [[ "$FGRID" != "" ]]; then
    TRICE_PLC_LIMIT=$(echo $FGRID | awk -F: '{print $2}');
    FRAC_PLC_LIMIT=$(echo $FGRID | awk -F: '{print $1}');
fi;

POINTS=20; INTERLEAVE=1
TOYS=50;   MULTIPLIER=4   # each job does ($POINTS/$INTERLEAVE) * $TOYS * $MULTIPLIER toys
TOTAL_TOYS=2500
if [[ "$REMULT" != "" ]]; then TOTAL_TOYS=$(( $TOTAL_TOYS * $REMULT )); fi;


if   (( $IMASS <= 130 )); then MULTIPLIER=1; INTERLEAVE=10;
elif (( $IMASS <= 135 )); then MULTIPLIER=1; INTERLEAVE=10;
elif (( $IMASS <= 145 )); then MULTIPLIER=1; INTERLEAVE=10;
elif (( $IMASS <= 150 )); then MULTIPLIER=1; INTERLEAVE=5;
elif (( $IMASS <= 178 )); then MULTIPLIER=1; INTERLEAVE=1; 
elif (( $IMASS <= 199 )); then MULTIPLIER=1; INTERLEAVE=1;
elif (( $IMASS <= 248 )); then MULTIPLIER=1; INTERLEAVE=3; 
else                          MULTIPLIER=1; INTERLEAVE=3; fi;

if echo $2 | grep -q combl.root; then INTERLEAVE=2; TOTAL_TOYS=1000; fi;
if [[ "$NAME" == "FC" ]]; then
    POINTS=$(( $POINTS * 2 ));
    INTERLEAVE=$(( $INTERLEAVE * 2));
fi;

if [[ "$FGRID" != "" ]]; then 
    if (( $MULTIPLIER * $POINTS / $INTERLEAVE / $FGP >= 1)) ; then
        MULTIPLIER=$(( $MULTIPLIER * $POINTS / $INTERLEAVE / $FGP));
        INTERLEAVE=1;
    else
        MULTIPLIER=1;
        INTERLEAVE=$(( $INTERLEAVE * $FGP / $POINTS ));
        if (( $INTERLEAVE == 0 )); then INTERLEAVE=1; fi;
    fi
    POINTS=$FGP; 
    TOTAL_TOYS=5000
    if [[ "$REMULT" != "" ]]; then TOTAL_TOYS=$(( $TOTAL_TOYS * $REMULT )); fi;
fi;

OPTIONS=""
if [[ "$LSF" == "1" ]]; then
    if echo $HOSTNAME | grep -q fnal; then
        OPTIONS="$OPTIONS --condor ";
    else
        OPTIONS="$OPTIONS --lsf -q 2nd";
    fi
    PFIX="${PFIX}_LSF"
fi
if [[ "$NAME" == "up" ]]; then 
    POINTS=2; 
    MULTIPLIER=$(( $MULTIPLIER * 6 / $INTERLEAVE ));
    INTERLEAVE=1; 
elif [[ "$NAME" == "low" ]]; then
    POINTS=$((POINTS / 4));
    if (( $INTERLEAVE > 6 )); then INTERLEAVE=$(( ( $INTERLEAVE + 2 )/ 4 )); else INTERLEAVE=1; fi
    TOTAL_TOYS=$(( $TOTAL_TOYS * 2 ));
fi;


if [[ "$NAME" == "SM" ]] || [[ "$NAME" == "sig" ]]; then
    if (( $INTERLEAVE > 1 )); then
        POINTS=$(( $POINTS / $INTERLEAVE ))
        INTERLEAVE=1;
        TOTAL_TOYS=$(( $TOTAL_TOYS * 8 / $POINTS ));
    else
        echo "Reducing points"
        POINTS=$(( $POINTS / 5 ));
        MULTIPLIER=$(( $MULTIPLIER * 5 ));
        TOTAL_TOYS=$(( $TOTAL_TOYS * 8 / $POINTS ));
    fi;
    if [[ "$NAME" == "sig" ]] && [[ "$ADAPTIVE" == "1" ]]; then
        TOTAL_TOYS=$(( $SIG_ATOTAL  / $POINTS ));
    fi;
fi;

CRABOPT=""
if [[ "$CWL" != "none" ]]; then
    CRABOPT="$CRABOPT -GRID.ce_white_list=$CWL -GRID.remove_default_blacklist=1"
elif [[ "$CBL" != "none" ]]; then
    CRABOPT="$CRABOPT -GRID.ce_black_list=$CBL"
fi;

if [[ "$PRI" == "1" ]]; then OPTIONS="$OPTIONS -P"; fi
if [[ "$GL"  == "1" ]]; then OPTIONS="$OPTIONS --glidein"; fi;
if [[ "$GS"  == "1" ]]; then OPTIONS="$OPTIONS --glidein --server"; fi;

if [[ "$NAME"  == "sig" ]]; then OPTIONS="$OPTIONS --signif"; fi;
if [[ "$FGRID" == "wide" ]]; then
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
    -O "$OPTS" -r
if [[ "$GO" == "1" ]]; then
    if (( $JOBS < 100 )); then
        crab -cfg grid-${2/.root/-${NAME}}.cfg -USER.ui_working_dir=crab_0_${2/.root/_${NAME}}$PFIX $CRABOPT -create -submit;
    else
        crab -cfg grid-${2/.root/-${NAME}}.cfg -USER.ui_working_dir=crab_0_${2/.root/_${NAME}}$PFIX $CRABOPT -create;
        for X in $(seq 1 20 $JOBS); do
            crab -c crab_0_${2/.root/_${NAME}}$PFIX $CRABOPT -submit 20
        done
    fi
fi
