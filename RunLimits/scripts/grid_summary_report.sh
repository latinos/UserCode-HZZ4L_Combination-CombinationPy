#!/bin/bash
WSP="comb"; if [[ "$1" != "" ]]; then WSP="$1"; fi;
echo "Status of $WSP at $(date)";
echo -e "$M\tRETR\tPERC%\tRUN\tPERC%\tALL";
XALL=0; XRETR=0; XRUN=0;
for M in $(cat masses.txt); do 
    ls $M/status.crab_0_${WSP}_[weS]* > /dev/null 2>&1 || continue
    RETRIEVED=$(cat $M/status.crab_0_${WSP}_[weS]* | grep -c '^[0-9].*\(Retrieved\|Done\|DONE\)')
    RUNNING=$(cat $M/status.crab_0_${WSP}_[weS]* | grep -c '^[0-9].*Run')
    ALL=$(cat $M/status.crab_0_${WSP}_[weS]* | grep -c '^[0-9].*')
    if [[ "$ALL" == "0" ]] ; then continue; fi;
    PERC=$(python -c "print '%4.0f' % ($RETRIEVED*100./$ALL);");
    RPERC=$(python -c "print '%4.0f' % ($RUNNING*100./$ALL);");
    XALL="$(( $XALL + $ALL ))";
    XRETR="$(( $XRETR + $RETRIEVED ))";
    XRUN="$(( $XRUN + $RUNNING ))";
    echo -e "$M\t$RETRIEVED\t$PERC%\t$RUNNING\t$RPERC%\t$ALL";
done
XPERC=$(python -c "print '%4.0f' % ($XRETR*100./$XALL);");
XRPERC=$(python -c "print '%4.0f' % ($XRUN*100./$XALL);");
echo -e "TOTAL\t$XRETR\t$XPERC%\t$XRUN\t$XRPERC%\t$XALL";
