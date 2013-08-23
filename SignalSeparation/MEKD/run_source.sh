#! /bin/bash

########
DIR=$1
MUTYPE=$2
FITNUIS=$3
NTOYS=$4
MH=$5
SEED=$6
CARD=$7
NITER=$((NTOYS / 1000))
NTOYS=1000
########

#CMSSW
cd $DIR
pwd
#export OSG_APP=/raid/osgpg/pg/app
#export SCRAM_ARCH=slc5_amd64_gcc472
#source $OSG_APP/cmssoft/cms/cmsset_default.sh
eval `scramv1 runtime -sh`


### fixed mu
if [[ "$MUTYPE" == "fixed" ]]
    then 

    if [[ "$FITNUIS" == "TRUE" ]]
	then
	combine -m ${MH} -M HybridNew --testStat=TEV --generateExt=1 --generateNuis=0 --fitNuis=1 ${CARD} --singlePoint 1 --saveHybridResult --seed ${SEED} -T ${NTOYS} --fork 1 -i ${NITER} --clsAcc 0 --fullBToys -n "fixedMu.fitNuis"
    else
	combine -m ${MH} -M HybridNew --testStat=TEV --generateExt=1 --generateNuis=0 --fitNuis=0 ${CARD} --singlePoint 1 --saveHybridResult --seed ${SEED} -T ${NTOYS} --fork 1 -i ${NITER} --clsAcc 0 --fullBToys -n "fixedMu"
    fi
    
### float mu
elif [[ "$MUTYPE" == "float" ]]
    then

    if [[ "$FITNUIS" == "TRUE" ]]
	then
	combine -m $MH -M HybridNew --testStat=TEV --generateExt=1 --generateNuis=0 ${CARD} --singlePoint 1 --saveHybridResult -T ${NTOYS} --fork 1 -i ${NITER} --clsAcc 0 --fullBToys -n "floatMu.fitNuis" --fitNuis=1 --seed ${SEED}
    else
	combine -m $MH -M HybridNew --testStat=TEV --generateExt=1 --generateNuis=0 ${CARD} --singlePoint 1 --saveHybridResult -T ${NTOYS} --fork 1 -i ${NITER} --clsAcc 0 --fullBToys -n "floatMu" --fitNuis=0 --seed ${SEED}
    fi
    
fi