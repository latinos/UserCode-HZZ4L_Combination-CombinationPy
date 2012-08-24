#!/bin/bash

if [ -z $1 ]; then
    echo "Please pass an inputs_<lumi> dir"
    exit;
fi

outDir=SortedDatacards
dcDir=datacardsHCG
wsDir=workspaces

mkdir $outDir

SHAPE="TRUE"
CnC="FALSE"

inputDir=$1

bool_110to140="TRUE"
bool_140to160="TRUE"
bool_160to290="TRUE"
bool_290to350="TRUE"
bool_350to400="TRUE"
bool_400to600="TRUE"

if [ ! -e $inputDir ]; then
	echo "$inputDir does not exist!"
	exit;
fi

elements=6
startMass=( 110.0 140.0 160.0 290.0 350.0 400.0 )
stepSizes=( 0.5 0.5 2.0 5.0 10.0 20.0 )
startVal=( 0 0 0 0 0 0 )
endVal=( 60 40 65 12 5 11 )
### default values
### startVal=( 0 0 0 0 0 0 )
### endVal=( 30 20 65 12 5 11 )

echo "--> ${elements} Different Steps"
a=0
while ((a < ${elements}))
	do
echo "---> Step ${a}"

	c=${startVal[a]}
	while ((c < ${endVal[a]})) 
		do
	
		mStart=${startMass[a]}
		step=${stepSizes[a]}
		area=$(echo "$mStart + ( $step * $c )" | bc)
		massDir=${area%%.0*}
		echo ${massDir}
		
		################################ 
		#A C T I V E   E L E M E N T S #
		mkdir $outDir/$massDir
		if [[ "$SHAPE" == "TRUE" ]]; then
		    cp $inputDir/$wsDir/*hzz4l_4muS.${area}.input.root $outDir/$massDir
                    cp $inputDir/$wsDir/*hzz4l_4eS.${area}.input.root $outDir/$massDir
                    cp $inputDir/$wsDir/*hzz4l_2e2muS.${area}.input.root $outDir/$massDir

		    cp $inputDir/$dcDir/*hzz4l_4muS.${area}.txt $outDir/$massDir
                    cp $inputDir/$dcDir/*hzz4l_4eS.${area}.txt $outDir/$massDir
                    cp $inputDir/$dcDir/*hzz4l_2e2muS.${area}.txt $outDir/$massDir

		    #cp $inputDir/${dcDir}_NoThSys/hzz4l_4muS.${area}.txt $outDir/$massDir/XS_hzz4l_4muScorr.${area}.txt
                    #cp $inputDir/${dcDir}_NoThSys/hzz4l_4eS.${area}.txt $outDir/$massDir/XS_hzz4l_4eScorr.${area}.txt
                    #cp $inputDir/${dcDir}_NoThSys/hzz4l_2e2muS.${area}.txt $outDir/$massDir/XS_hzz4l_2e2muScorr.${area}.txt

		fi

		if [[ "$CnC" == "TRUE" ]]; then
		    cp $inputDir/hzz4l_4muC.${area}.txt $outDir/$massDir
                    cp $inputDir/hzz4l_4eC.${area}.txt $outDir/$massDir
                    cp $inputDir/hzz4l_2e2muC.${area}.txt $outDir/$massDir

		fi		
		
		################################ 		

		let c=$c+1
		
	done

	let a=$a+1
done

