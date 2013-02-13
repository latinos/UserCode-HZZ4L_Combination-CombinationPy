
#!/bin/bash
#############################################################
#
# Driver script for creating Hybrid or Frequentist grids
#
# author: Giovanni Petrucciani, UCSD
#         from a similar script by Luca Lista, INFN
#
##############################################################

i="$1"
if [ "$i" = "" ]; then
  echo "Error: missing job index"
  exit 1;
fi
echo "max events from CRAB: $MaxEvents"
n="$MaxEvents"
if [ "$n" = "" ]; then
  n="$2"
fi
if [ "$n" = "" ]; then
  echo "Error: missing number of experiments"
  exit 2;
fi

## Save memory on batch systems by avoinding a redundant fork when only one child will be ever spawned
nchild=1;
if  [[ "$nchild" == "1" && "$n" == "1" ]]; then
    nchild=0;
fi;

NTOYS=200
echo "## Starting at $(date); $NTOYS toys"

(( ($i + 1) % 1 == 0 )) &&  ./combine floatMu.root -M HybridNew --testStat=TEV --generateExt=1 --generateNuis=0 --singlePoint 1 --saveHybridResult  --clsAcc 0 --fullBToys -m 126.0 --fork $nchild -T $NTOYS -v 0 -n hzz4l_SpinCP_CRAB --saveToys -s $((10000 + $i)) -i $n



hadd hzz4l_SpinCP_CRAB.root higgsCombine*.root
echo "## Done at $(date)"
