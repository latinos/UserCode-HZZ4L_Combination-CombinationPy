#!/bin/bash
WORK=$PWD
SRC=$src
SRC=$CMSSW_BASE/src
QUEUE=1nd
if echo "X$1" | grep -q '^X-[a-zA-Z0-9]'; then 
    QUEUE=$(echo "X$1" | sed 's/^X-//');
    shift;
fi;
bsub -q $QUEUE $WORK/lxbatch_runner.sh $WORK $SRC $*
