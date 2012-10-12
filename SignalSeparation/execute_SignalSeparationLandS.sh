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

cp runSignalSeparation.py submitToLXB_tpl.sh submitToPBS_tpl.csh.pbs tdrstyle.cc haddLands.py $cardDir/

cd $cardDir
outDir="output_LandS/"



action=1
NJOBS=40  # total number of parallel jobs
NTOYS=100 # toys per parallel job
MH=125

if [ $# -ge 4 ]
    then
    action=$4
fi

#Step 1: generate toys

if [ $action -eq 1 ]
    then 
    echo "GENERATING TOYS"
    if [ -d $outDir ]
	then
	echo "Output directory ${cardDir}/${outDir}/ already existing. I will not overwrite. Please remove it and try again."
	exit 2
    fi
    
    mkdir $outDir

    python runSignalSeparation.py -b -m -t "lands" --generateToys --nParallelJobs $NJOBS --toysPerJob $NTOYS -o "$outDir" --card1 "$card1" --card2 "$card2"  --mH $MH 
#--TeVStat not commissioned yet

#Step 2: fit toys
elif [ $action -eq 2 ]
    then
    echo "FITTING TOYS"

    if [ ! -d $outDir ]
	then
	echo "Output directory ${cardDir}/${outDir}/ DOES NOT  exist. "
	exit 2
    fi

    python runSignalSeparation.py -b -m -t "lands" --fitToys --nParallelJobs $NJOBS  --toysPerJob $NTOYS  -o "$outDir" --card1 "$card1" --card2 "$card2" --mH $MH

#Step 3: plot variables
elif [ $action -eq 3 ]
    then 
    echo "PLOT VARIABLES"
    if [ ! -d $outDir ]
	then
	echo "Output directory ${cardDir}/${outDir}/ DOES NOT  exist. "
	exit 2
    fi
    mkdir figs/
    python runSignalSeparation.py -b -m -t "lands" --plotResults -a --nParallelJobs $NJOBS  --toysPerJob $NTOYS  -o "$outDir" --card1 "$card1" --card2 "$card2"
elif [ $action -eq -1 ]
    then 
    echo "CLEAN directory with cards"
    rm -f hzz4l_*copy_*txt
    rm -fr LSF*
    rm -fr output_LandS/
    rm -fr figs/
else
    echo "Requested to perform and unrecognized action: "${action}
    echo "action can be 1:generate toys  ;   2:fit toys   ;   3:plot variables  ;  -1:clean up output directory"
    echo "Exiting."
    exit 3
fi

