#!/bin/bash
count=0

checkStatus(){
    while [ '1' == '1' ]; do
        mesg=`squeue -u ACMEtany | grep ACMEtany`
        if [ "$mesg" == ''  ]; then
            break
        fi
        sleep 1m
    done
}

for i in `ls mc.*.sh | sort -t '.' -nk 3 `; do
    # skip finished jobs
    flog=`grep '#SBATCH -o' $i | awk '{print $3}'`

    sbatch $i
#    let count++
#    rem=$(( $count % 2 ))
#    if [ $rem -eq 0 ]; then 
#	checkStatus
#    fi
done

