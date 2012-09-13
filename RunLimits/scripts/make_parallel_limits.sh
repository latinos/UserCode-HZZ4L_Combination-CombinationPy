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

if [[ "$TYPE" != "ASCLS" && "$TYPE" != "PLP" && "$TYPE" != "PLPE" && "$TYPE" != "ML"  && "$TYPE" != "PLS"  && "$TYPE" != "PLSE" ]]; then
    
    echo "Unkown Type: $TYPE"
    echo "Options: ASCLS, PLP, PLPE, PLS, ML"
    exit;

fi





for m in $(cat $MASSFILE); 
  do

  qsub makeLimits.pbs.sh -v MASS=${m},TYPE=${TYPE},OPTIONS=${OPTIONS} -N "${TYPE}_${m}"


done

