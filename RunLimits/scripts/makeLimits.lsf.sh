#!/bin/bash

TYPE=$1
MASS=$2
OPTION=$3

cd ${LS_SUBCWD}

echo "LSF job running in: " `pwd` with options $TYPE $MASS $OPTION

eval `scram runtime -sh`

if [[ "$TYPE" == "ASCLS" ]]; then
    
    bash make_ASCLS.sh $OPTIONS -l $MASS comb_hzz4l.root

elif [[ "$TYPE" == "PLP" ]]; then
    
    bash make_PLC.sh $OPTIONS -P $MASS comb_hzz4l.root

elif [[ "$TYPE" == "PLPE" ]]; then
    
    bash make_PLC.sh $OPTIONS --PE $MASS comb_hzz4l.root

elif [[ "$TYPE" == "ML" ]]; then
    
    bash make_ML.sh $OPTIONS $MASS comb_hzz4l.root

else 
    echo "Unkown Type: $TYPE"
    exit;

fi



