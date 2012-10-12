#!/bin/bash     

if [ $# -lt 1 ]
then
    echo "submitToLXB wants at least one input argument:  name of output directory. Exiting."
    exit 2
fi


OutputDir=$1

#name of Input directory with cards and
InputDir="./"
if [ $# -gt 1 ]
    then
    InputDir=$2
fi

#

curdir=$( pwd )
#cd ${curdir}/$InputDir
#cardStemDir=$( pwd )
#cd -
cmsswbase=CMSSWBASE
workdir=WORKDIR

TOYFILE
TOYCARD
MYSEED
MYSTEM
MYNTOYS

echo inputCardDir: $workdir/$InputDir
echo outputdir: $OutputDir
echo
echo File with toy: $toyFile
echo Dummy card: $copyCard
echo
#source /scratch1/hep/cms/cmsset_default.csh

export SCRAM_ARCH=slc5_amd64_gcc462
echo SCRAM_ARCH:$SCRAM_ARCH

export PATH="$PATH:${cmsswbase}/src/HiggsAnalysis/LandS/test" 

echo "Setting up $cmsswbase"
cd $cmsswbase/src
eval `scramv1 runtime -sh`

#cmsenv
#cd $workdir
#cd $curdir
echo "Moving to $TMPDIR"
cd $TMPDIR

#cp -r ${cardStemDir} .
cp -r ${workdir}/$InputDir .
cp -r $OutputDir .
#cp -r ${InputDir}/../workspaces/ .
#mkdir -p ${OutputDir}

echo "Current Directory is $( pwd )"
echo "Working dir is $workdir"
echo
echo "Files in $PWD :"
ls -lh
echo
#echo "Files in ${workdir}/$InputDir :"
#ls -lh ${workdir}/$InputDir
#echo


COMMAND


echo "Finished main body of the program. List of files in $PWD"
ls -lh
echo
echo
echo

# move output file to right directory

python haddLands.py --seed ${SEED} --toysPerJob ${NTOYS} --stem ${STEM}

cp ONAME ODIR/ONAME
cp ONAME ${workdir}/ODIR/ONAME

RESUBMIT_IF_FAILED=false

if [ ! -f ${workdir}/ODIR/ONAME ]
    then
    echo "ERROR from submitToLXB !!! Output file (ONAME)is not in ${workdir}/ODIR/ <==="
    if $RESUBMIT_IF_FAILED 
	then
	echo "Rerunning the limit calculation"
	COMMAND
	python haddLands.py --seed ${SEED} --toysPerJob ${NTOYS} --stem ${STEM}
	cp ONAME ${workdir}/ODIR/ONAME

    fi
fi

###REMOVE###
