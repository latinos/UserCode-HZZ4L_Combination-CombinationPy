#!/bin/bash
OPTIONS="";

if [[ "$1" == "--blind" ]]; then OPTIONS="-O=-t -O=-1"; shift; fi
MASS="$1"; if [[ "$1" == "" ]]; then echo "Usage: $0 datacards"; exit 1; fi;
shift; test -d $MASS && cd $MASS;

DATACARDS="$*";

if [[ "$DATACARDS" == "" ]]; then DATACARDS=$(ls -1 [Vhv]*txt); fi;
for X in $DATACARDS; do
    echo $X | grep -q clean && continue; # what is clean, thou shalt not clean again
    echo ; echo " ==== ANALYZING NUISANCES OF DATACARD $X (`date`) ==== "
    test -d  $X.d && rm -r $X.d
    python $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/test/sizeUpSystematics.py $OPTIONS $PWD/$X --masses=$MASS --X-keep-global-nuisances --dir $X.d
done
