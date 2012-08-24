#!/bin/bash
BIN=1; if [[ "$1" == "-n" ]]; then BIN=0; shift; fi;
MASS=$1

for WSP in $*; do
    TXT="$WSP.txt"; 
    test -f $MASS/$TXT.gz && TXT="$TXT.gz"
    test -f $MASS/$TXT || continue;
    [[ "$BIN" == "1" ]] && bash make_binary_workspaces.sh $MASS $TXT
    bash make_PLC.sh    -P $MASS ${WSP}.root
    bash make_PLC.sh    --PE $MASS ${WSP}.root
    bash make_ASCLS.sh  -O $MASS ${WSP}.root
    bash make_ASCLS.sh  -E $MASS ${WSP}.root
    #bash make_ML.sh     -Z $MASS ${WSP}.root
done;
