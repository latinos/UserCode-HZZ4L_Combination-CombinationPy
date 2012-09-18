#! /bin/bash

### Inputs:
# 1: directory with cards
# 2: first card (signal model1)
# 3: second card (signal model2)
# 4: what to do (default is generate toys)


if [ $# -lt 3 ]
    then
    echo "Need at least three arguments: "
    echo "    1) directory with cards"
    echo "    2) first card (signal model1)"
    echo "    3) second card (signal model2)"
    exit 1
fi


cardDir=$1
card1=$2
card2=$3

cp runSignalSeparation.py submitToLXB_tpl.sh submitToPBS_tpl.csh.pbs tdrstyle.cc $cardDir/

cd $cardDir
outDir="output_LandS/"

if [ -d $outDir ]
    then
    echo "Output directory ${cardDir}/${outDir}/already exisiting. I will not overwrite. Please remove it and try again."
    exit 2
fi

mkdir $outDir


action=0
NJOBS=200  # total number of parallel jobs
NTOYS=20 # toys per parallel job
if [ $# -ge 4 ]
    then
    action=$4
fi

#Step 1: generate toys

if [ $action -eq 0 ]
    then 
    echo "GENERATING TOYS"
    python runSignalSeparation.py -b -m -t "lands" --generateToys --nParallelJobs $NJOBS --toysPerJob $NTOYS -o "$outDir" --card1 "$card1" --card2 "$card2" --TeVStat
#Step 2: fit toys
elif [ $action -eq 1 ]
    then
    echo "FITTING TOYS"
    python runSignalSeparation.py -b -m -t "lands" --fitToys --nParallelJobs $NJOBS  --toysPerJob $NTOYS  -o "$outDir" --card1 "$card1" --card2 "$card2"
#Step 3: plot variables
elif [ $action -eq 2 ]
    then 
    echo "PLOT VARIABLES"
    python runSignalSeparation.py -b -m -t "lands" --plotResults -a --nParallelJobs $NJOBS  --toysPerJob $NTOYS  -o "$outDir" --card1 "$card1" --card2 "$card2"
else
    echo "Requested to perform and unrecognized action: "${action}
    echo "action can be 0:generate toys  ;   1:fit toys   ;   2:plot variables"
    echo "Exiting."
    exit 3
fi

