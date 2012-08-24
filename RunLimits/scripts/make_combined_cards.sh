#!/bin/bash
if [[ "$1" == "" ]]; then echo "Usage: $0 mass"; exit 1; fi;
if test -d $1; then MASS=$1; else echo "Usage: $0 mass"; exit 1; fi; 
cd $MASS
WHAT=$2

if echo $MASS | grep -F -q .5; then 
    FRAC=1; MASSG=$MASS; MASSD="$MASS"; MASS=${MASS/.5/}; MASSP1=$((MASS+1));  
else 
    FRAC=0; MASSG=$MASS; MASSD="$MASS.0"; 
fi

VHBB="vhbb_DC_ALL_BDT.${MASS}.0.txt" # ${MASSD}
HWW2L_01J="hww_0jsf_shape=hwwsf_0j_shape.txt  hww_0jof_shape=hwwof_0j_shape.txt hww_1jsf_shape=hwwsf_1j_shape.txt  hww_1jof_shape=hwwof_1j_shape.txt"
HWW2L_2J="hww_2j_cut=hww_2j_cut.txt"
HWW2L="$HWW2L_01J $HWW2L_2J"
HWWWL="vhww_3l=vh3l_cut.txt"
HWWLJ="" # MISSING
HZZ4L="hzz4l_2e2mu=hzz4l_2e2muS.${MASSD}.txt hzz4l_4e=hzz4l_4eS.${MASSD}.txt hzz4l_4mu=hzz4l_4muS.${MASSD}.txt"
HZZ2Q="hzz2l2q_ee0b=hzz2l2q_ee0b.${MASSG}.txt  hzz2l2q_ee1b=hzz2l2q_ee1b.${MASSG}.txt  hzz2l2q_ee2b=hzz2l2q_ee2b.${MASSG}.txt "
HZZ2Q="hzz2l2q_mm0b=hzz2l2q_mm0b.${MASSG}.txt  hzz2l2q_mm1b=hzz2l2q_mm1b.${MASSG}.txt  hzz2l2q_mm2b=hzz2l2q_mm2b.${MASSG}.txt $HZZ2Q"
HZZ2N="hzz2l2nu_ee=hzz2l2nu_ee_${MASS}.txt  hzz2l2nu_mm=hzz2l2nu_mm_${MASS}.txt"
HZZ2T="hzz2l2t_${MASS}.txt"
HTT0="htt_et_0=htt_et_0.txt  htt_et_1=htt_et_1.txt htt_et_2=htt_et_2.txt"
HTT0="htt_mt_0=htt_mt_0.txt  htt_mt_1=htt_mt_1.txt htt_mt_2=htt_mt_2.txt $HTT0"
HTT0="htt_em_0=htt_em_0.txt  htt_em_1=htt_em_1.txt htt_em_2=htt_em_2.txt $HTT0";
HTTM="htt_mm_0=htt_mm_0.txt  htt_mm_1=htt_mm_1.txt htt_mm_2=htt_mm_2.txt"
HTTV="vhtt=vhtt.txt" 

COMB=""
COMBP=""
COMBL=""
COMBGGH="";
COMBVBF="";   
COMBVH="";

## VH, H->BB 
if (( $MASS < 135 )) || (( $MASS == 135 && $FRAC == 0 )); then 
    [[ "$FRAC" == "1" ]] && ( test -f vhbb_DC_ALL_BDT.${MASS}.0.txt || cp -v -s ../$MASS/vhbb_DC_ALL_BDT.${MASS}.0.txt . -v )
    test -L comb_vhbb.txt || cp -s $VHBB comb_vhbb.txt; 
    COMB="$COMB vhbb=$VHBB";  
    COMBL="$COMBL vhbb=$VHBB";  
    COMBVH="$COMBVH vhbb=$VHBB";  
fi

## Gamma Gamma
if (( $MASS <= 150 )); then 
    ## first produce the following under /common/
    ##    text2workspace.py hggmva_databinned.txt
    test -f comb_hgg.txt  || cp -s ../common/hggmva_databinned.txt  comb_hgg.txt
    test -f comb_hgg.root || cp -s ../common/hggmva_databinned.root comb_hgg.root
    #test -f comb_hgg_vbf.txt  || cp -s ../common/comb_hgg_vbf.txt  comb_hgg_vbf.txt
    #test -f comb_hgg_vbf.root || cp -s ../common/comb_hgg_vbf.root comb_hgg_vbf.root
    #test -f comb_hgg_novbf.txt  || cp -s ../common/comb_hgg_novbf.txt  comb_hgg_novbf.txt
    #test -f comb_hgg_novbf.root || cp -s ../common/comb_hgg_novbf.root comb_hgg_novbf.root
    #COMB="$COMB .=comb_hgg_vbf.txt .=comb_hgg_novbf.txt"; 
    #COMBGGH="$COMBGGH .=comb_hgg_novbf.txt"
    #COMBVBF="$COMBVBF .=comb_hgg_vbf.txt"
    COMB="$COMB hgg=comb_hgg.txt"; 
    COMBP="$COMBP hgg=comb_hgg.txt"; 
fi

## ======== W W ======================
COMBWW="";
##### Good old H->WW->lvlv
[[ "$FRAC" == "1" ]] && ( test -f hwwsf_0j_shape.txt || cp -v -s ../$MASSP1/hww[so]*{txt,root} . -v )
[[ "$FRAC" == "1" ]] && ( test -f hww_2j_cut.txt || cp -v -s ../$MASSP1/hww_2j_cut.txt . -v )
combineCards.py -S $HWW2L > comb_hww2l.txt
COMBWW="$COMBWW $HWW2L"
COMBGGH="$COMBGGH $HWW2L_01J";
COMBVBF="$COMBVBF $HWW2L_2J";

##### Brand new WH->WWW->lvlvlv
if (( $MASS <= 200)); then
    [[ "$FRAC" == "1" ]] && ( test -f vh3l_cut.txt || cp -v -s ../$MASSP1/vh3l*txt . -v )
    combineCards.py -S $HWWWL > comb_vhww3l.txt
    COMBWW="$COMBWW $HWWWL"
    COMBVH="$COMBVH $HWWWL"
fi;

#### WW->lvqq: Mr PhD, where art thou?

#### All H->WW modes combined
combineCards.py -S $COMBWW > comb_hww.txt
COMB="$COMB $COMBWW"
COMBL="$COMBL $COMBWW"
## =====================================



## ======== Z Z ======================
COMBZZ=""

#### ZZ 4L
COMBZZ="$COMBZZ $HZZ4L"
combineCards.py -S $HZZ4L > comb_hzz4l.txt

#### ZZ 2Q
if (( $MASS >= 200)); then 
   COMBZZ="$COMBZZ $HZZ2Q";
   combineCards.py -S $HZZ2Q > comb_hzz2l2q.txt
elif (( $MASS >= 130 && $MASS <= 164 )); then
   if [[ "$FRAC" == "1" ]]  ; then
        echo "Missing ZZ2L2Q at $MASSG, will take from $MASSP1";
        ( test -f hzz2l2q_ee0b.${MASSG}.txt || (cp -v -s ../$MASSP1/hzz2l2q*txt  . -v ; rename .$MASSP1. .$MASSG. hzz2l2q*txt) )
        ( test -f hzz2l2q_ee0b.input.root   || (cp -v -s ../$MASSP1/hzz2l2q*root . -v ; rename .$MASSP1. .$MASSG. hzz2l2q*root) )
   fi;
   COMBZZ="$COMBZZ $HZZ2Q";
   combineCards.py -S $HZZ2Q > comb_hzz2l2q.txt
fi;

## ZZ 2N
if (( $MASS >= 250)); then 
   COMBZZ="$COMBZZ $HZZ2N";
   combineCards.py -S $HZZ2N > comb_hzz2l2nu.txt
fi;

## ZZ 2T
if (( $MASS >= 190)); then 
    COMBZZ="$COMBZZ hzz2l2t=$HZZ2T";
    test -L comb_hzz2l2t.txt || cp -s $HZZ2T comb_hzz2l2t.txt; 
fi;

combineCards.py -S $COMBZZ > comb_hzz.txt

COMBGGH="$COMBGGH $COMBZZ";
COMB="$COMB $COMBZZ";
COMBP="$COMBP $COMBZZ"
## =====================================


## ======== Tau Tau ======================
COMBTT=""
if (( $MASS < 140 )) || (( $MASS == 140 && $FRAC == 0 )); then 
   [[ "$FRAC" == "1" ]] && ( test -f htt_mt_0.txt || cp -v -s ../$MASS/htt*{txt,root} . -v )
   [[ "$FRAC" == "1" ]] && ( test -f vhtt.txt         || cp -v -s ../$MASS/vhtt*txt         . -v )
   [[ "$FRAC" == "1" ]] && ( test -f vhtt_shapes.root || cp -v -s ../$MASS/vhtt_shapes.root . -v )
   combineCards.py -S $HTT0 > comb_htt0.txt
   combineCards.py -S $HTTM > comb_httm.txt
   combineCards.py -S $HTTV > comb_vhtt.txt
   COMBTT="$HTT0 $HTTM $HTTV"

   combineCards.py -S $COMBTT > comb_htt.txt
   COMB="$COMB $COMBTT"
   COMBL="$COMBL $COMBTT"
elif (( $MASS < 145 )) || (( $MASS == 145 && $FRAC == 0 )); then 
   [[ "$FRAC" == "1" ]] && ( test -f htt_mt_0.txt || cp -v -s ../$MASS/htt*{txt,root} . -v )
   combineCards.py -S $HTT0 > comb_htt0.txt
   combineCards.py -S $HTTM > comb_httm.txt
   COMBTT="$HTT0 $HTTM"

   combineCards.py -S $COMBTT > comb_htt.txt
   COMB="$COMB $COMBTT"
   COMBL="$COMBL $COMBTT"
fi
## =======================================

combineCards.py -S $COMB > comb.txt

if (( $MASS < 145 )) || (( $MASS == 145 && $FRAC == 0 )); then
    combineCards.py -S $COMBL > combl.txt
fi
if (( $MASS <= 150)); then
    combineCards.py -S $COMBP > combp.txt
fi

echo "Done cards for mass $MASSD"
