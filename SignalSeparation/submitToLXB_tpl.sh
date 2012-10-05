#!/bin/bash     

if [ $# -lt 1 ]
then
    echo "submitToLXB wants at least one input argument: name of Input directory. Exiting."
    exit 2
fi

InputDir=$1

echo inputdir: $InputDir

cmsswbase=CMSSWBASE
workdir=WORKDIR
#source /scratch1/hep/cms/cmsset_default.csh

export SCRAM_ARCH=slc5_amd64_gcc462
echo SCRAM_ARCH:$SCRAM_ARCH

export PATH="$PATH:${cmsswbase}/src/HiggsAnalysis/LandS/test" 

echo "Setting up $cmsswbase"
cd $cmsswbase/src
eval `scramv1 runtime -sh`
#cmsenv
cd $workdir
echo "Current Directory is $PWD"

COMMAND

# move output file to right directory
mv ONAME ODIR/ONAME

###REMOVE###
