#!/bin/bash

if [[ -z $1 ]]; then 
    echo "Usage: ./make_parallel_limit.sh <type> <massfile> <options>"
    exit ;
fi


mkdir errFiles
mkdir outFiles
mkdir results

TYPE=$1
MASSFILE=$2
OPTIONS=$3

if [[ "$TYPE" != "ASCLS" && "$TYPE" != "PLP" && "$TYPE" != "PLPE" && "$TYPE" != "ML" ]]; then
    
    echo "Unkown Type: $TYPE"
    echo "Options: ASCLS, PLP, PLPE, ML"
    exit;

fi





for m in $(cat $MASSFILE); 
  do

  bsub -q 1nh -o ${m}/lsflog.txt makeLimits.lsf.sh ${TYPE} ${m} ${OPTIONS}

done

