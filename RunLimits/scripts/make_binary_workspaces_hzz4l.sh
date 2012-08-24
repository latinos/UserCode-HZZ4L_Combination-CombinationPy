#!/bin/bash
IS4=""; if [[ "$1" == "-4" ]]; then IS4="SM4_"; shift; fi;
if test -d $1; then MASS=$1; else echo "Usage: $0 mass [what ]"; exit 1; fi; 
cd $MASS
MATCH=$2;

function t2w {
    if [[ "$MATCH" == "" || "$MATCH" == "$1" ]]; then
        test -f $1 &&  text2workspace.py -b $1
    fi;
}

if [[ "$MATCH" == "" ]]; then
#    t2w ${IS4}comb_hgg.txt
#    t2w ${IS4}comb_hww.txt
#    t2w ${IS4}comb_htt.txt
    t2w ${IS4}comb_hzz4l.txt
#    t2w ${IS4}comb_hzz2l2nu.txt
#    t2w ${IS4}comb_hzz2l2q.txt
#    t2w ${IS4}comb.txt
else
    t2w $MATCH
fi
echo "Done workspaces for $MASS"
