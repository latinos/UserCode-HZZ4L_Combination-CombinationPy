#!/bin/bash
ls [0-9]*/status.crab* > all-tasks
grep '^[0-9]' [0-9]*/status.crab* | grep -v 'Terminated\|Retrieved' > ongoing-jobs
for S in $(cat all-tasks); do 
    D=${S/status./};
    if test -d $D; then
        if grep -F -q  "$S:" ongoing-jobs; then 
            echo "$D ongoing"; else echo "$D done"; 
        fi; 
    else
        echo "Task $D does not exist. Removing status file." 1>&2;
        rm $S;
    fi;
done | awk '/done/{print $1}'
