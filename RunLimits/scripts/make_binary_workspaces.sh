#!/bin/bash
PREFIX=""; if [[ "$1" == "-4" ]]; then PREFIX="SM4_"; shift; elif [[ "$1" == "-F" ]]; then PREFIX="FF_"; shift; fi;
UPD=0; if [[ "$1" == "-u" ]]; then UPD=1; shift; fi;
MORPH="shape2"; if [[ "$1" == "-m" ]]; then MORPH=$2; shift; shift; fi;
if test -d $1; then MASS=$1; else echo "Usage: $0 mass [what ]"; exit 1; fi; 
cd $MASS
MATCH=$2;

function t2w {
    if [[ "$MATCH" == "" || "$MATCH" == "$1" ]]; then
        if [[ "$UPD" == "1" ]]; then test ${1/.txt/.root} -nt $1 && return; fi;
        if test -f $1; then
            OUT=${1/.txt/}.root
            test -f $OUT && rm $OUT 
            text2workspace.py -b $1 --default-morphing=$MORPH  -o $OUT
            if test \! -f  $OUT; then
                echo "Failed to convert $1 at mH = $MASS to binary format"
                exit 1;
            fi;
            echo "Converted $1 at mH = $MASS to binary format $OUT "
        fi;
    fi;
}

if [[ "$MATCH" == "" ]]; then
    for W in $(cat ../workspaces.txt); do 
        if echo $W | grep -q 'comb_hgg' && [[ "$PREFIX" != "FF_" ]]; then continue; fi
        t2w ${PREFIX}${W}.txt
    done
else
    t2w ${PREFIX}$MATCH
fi
echo "Done workspaces for $MASS"
