#!/bin/bash
FORCE=0; if [[ "$1" == "-f" ]]; then FORCE=1; shift; fi;
MASS=$1; if [[ "$1" == "" ]]; then echo "Usage: $0 mass workspace"; exit 1; fi;
cd $MASS || exit 1
ls crab_0_* > /dev/null 2>&1 || exit 0;

for d in crab_0_*; do 
    test -L $d && continue;
    if echo $d | grep -q LSF && echo $HOSTNAME | grep -v -q lxplus; then continue; fi
    if test -f status.$d; then
        if grep -q 'Waiting\|Submitted\|Ready\|Running\|RUN\|DONE\|Aborted\|Done\|Created' status.$d; then
            crab -c $d -status 2>&1 | tee status.$d;    
        elif [[ "$FORCE" == "1" ]]; then
            crab -c $d -status 2>&1 | tee status.$d;    
        elif grep -q 'Retrieved' status.$d; then
            echo "Job $d is completed, nothing to do for it"
        else
            crab -c $d -status 2>&1 | tee status.$d;    
        fi;
    else
        crab -c $d -status 2>&1 | tee status.$d;
    fi 
    if grep -q 'Created' status.$d; then crab -c $d -submit; fi;
    if grep -q 'Done\|DONE' status.$d; then crab -c $d -get; fi;
    if grep -q 'Aborted' status.$d; then 
        crab -c $d $(perl -e '
    %bl=();@a=(); 
    while(<>){ 
        ($j,$e,$s) = m/^(\d+)\s+\S+\s+(?:Aborted|Done)\s+Aborted(\s+\S+?\.(\S+))?/ or next; 
        push @a,$j; 
        $bl{$s}=1 if $e;
    } 
    $blc = "";
    if (scalar(keys(%bl))) { $blc = " -GRID.ce_black_list  " . join(",", keys(%bl))  }
    print "-resubmit " . join(",",@a) . " $blc\n"' status.$d);
    fi;
done
