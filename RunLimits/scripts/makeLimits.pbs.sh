#!/bin/bash

#### -------------------------------------------------------------- 
#### Configuration for PBS stuff
#PBS -q submit
#PBS -l walltime=00:30:00
#PBS -l pmem=1gb
#PBS -l nodes=1:ppn=1:phase3
#PBS -N $PBS_JOBNAME
#PBS -o outFiles/$PBS_JOBNAME.out
#PBS -e errFiles/$PBS_JOBNAME.err
#### --------------------------------------------------------------


cd $PBS_O_WORKDIR

export CMSSWVER=CMSSW_5_2_5
export SCRAM_ARCH=slc5_amd64_gcc462
export OSG_APP=/osg/app
export VO_CMS_SW_DIR=${OSG_APP}/cmssoft/cms
export CMS_PATH=${VO_CMS_SW_DIR}
. ${CMS_PATH}/cmsset_default.sh;
eval `scramv1 runtime -sh`
eval `edmPluginRefresh -p ../../lib/$SCRAM_ARCH`



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



