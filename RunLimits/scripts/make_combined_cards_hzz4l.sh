#!/bin/bash

if [[ "$1" == "" ]]; then echo "Usage: $0 mass sqrts"; exit 1; fi;
if test -d $1; then MASS=$1; else echo "Usage: $0 mass sqrts"; exit 1; fi; 
if [[ "$1" != "" ]]; then sqrts=$2; else echo "Usage: $0 mass sqrts"; exit 1; fi; 

cd $MASS


COMB=""
SM4_COMB=""
XS_COMB=""



if [[ "$2" == "0" ]]; then COMBINED4L="TRUE"; else COMBINED4L="FALSE"; fi
if [[ "$2" == "-1" ]]; then COMBINED2l2tau="TRUE"; else COMBINED2l2tau="FALSE"; fi
if [[ "$2" == "-2" ]]; then COMBINED4L_JETS="TRUE"; else COMBINED4L_JETS="FALSE"; fi


#if expr index $MASS .5 > /dev/null; then MASSD="$MASS"; else MASSD="$MASS.0"; fi
if echo $MASS | grep -F -q .5; then MASSD="$MASS"; else MASSD="$MASS.0"; fi

if [[ "$COMBINED4L" == "TRUE" ]]; then
    HZZ4L="hzz4l_2e2mu_7=hzz4l_2e2muS_7TeV.txt hzz4l_4e_7=hzz4l_4eS_7TeV.txt hzz4l_4mu_7=hzz4l_4muS_7TeV.txt hzz4l_2e2mu_8=hzz4l_2e2muS_8TeV.txt hzz4l_4e_8=hzz4l_4eS_8TeV.txt hzz4l_4mu_8=hzz4l_4muS_8TeV.txt"

elif [[ "$COMBINED4L_JETS" == "TRUE" ]]; then
    HZZ4L="hzz4l_2e2mu_7_0=hzz4l_2e2muS_7TeV_0.txt hzz4l_4e_7_0=hzz4l_4eS_7TeV_0.txt hzz4l_4mu_7_0=hzz4l_4muS_7TeV_0.txt hzz4l_2e2mu_8_0=hzz4l_2e2muS_8TeV_0.txt hzz4l_4e_8_0=hzz4l_4eS_8TeV_0.txt hzz4l_4mu_8_0=hzz4l_4muS_8TeV_0.txt hzz4l_2e2mu_7_1=hzz4l_2e2muS_7TeV_1.txt hzz4l_4e_7_1=hzz4l_4eS_7TeV_1.txt hzz4l_4mu_7_1=hzz4l_4muS_7TeV_1.txt hzz4l_2e2mu_8_1=hzz4l_2e2muS_8TeV_1.txt hzz4l_4e_8_1=hzz4l_4eS_8TeV_1.txt hzz4l_4mu_8_1=hzz4l_4muS_8TeV_1.txt"

elif [[ "$COMBINED2l2tau" == "TRUE" ]]; then

    if (( $MASS > 178 ))
	then  HZZ4L="hzz4l_*.txt hzz2l2t_*.txt"
    else
	HZZ4L="hzz4l_*.txt"
    fi
else 
    HZZ4L="hzz4l_2e2mu=hzz4l_2e2muS_${sqrts}TeV.txt hzz4l_4e=hzz4l_4eS_${sqrts}TeV.txt hzz4l_4mu=hzz4l_4muS_${sqrts}TeV.txt"
fi
SM4_HZZ4L="hzz4l_2e2mu=SM4_hzz4l_2e2muS_${sqrts}TeV.txt hzz4l_4e=SM4_hzz4l_4eS_${sqrts}TeV.txt hzz4l_4mu=SM4_hzz4l_4muS_${sqrts}TeV.txt"
XS_HZZ4L="hzz4l_2e2mu=XS_hzz4l_2e2muS_${sqrts}TeV.txt hzz4l_4e=XS_hzz4l_4eS_${sqrts}TeV.txt hzz4l_4mu=XS_hzz4l_4muS_${sqrts}TeV.txt"

## ZZ 4L
COMB="$COMB $HZZ4L"
combineCards.py -S $HZZ4L > comb_hzz4l.txt
#combineCards.py -S $HZZ4LU > comb_hzz4l_uncorrected.txt

#COMBU=$(echo $COMB | sed "s/Scorr\.$MASS\.0/S\.$MASS\.0/g");
#combineCards.py -S $COMB > comb.txt
#combineCards.py -S $COMBU > comb_uncorrected.txt

echo "Done cards for mass $MASS"

###### SM4

## ZZ 4L
#SM4_COMB="$SM4_COMB $SM4_HZZ4L"
#combineCards.py -S $SM4_HZZ4L > SM4_comb_hzz4l.txt

#echo "Done SM4 cards for mass $MASS"

###### XS

## ZZ 4L
#XS_COMB="$XS_COMB $XS_HZZ4L"
#combineCards.py -S $XS_HZZ4L > XS_comb_hzz4l.txt

#echo "Done XS cards for mass $MASS"
