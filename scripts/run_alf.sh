#!/bin/bash

MAX_TRIALS=10

if [ $# -lt 1 ]; then
    echo -e "\tusage: $0 <DRW FILE 1> ... <DRW FILE N>"
    exit 1
fi

for i in "$@"; do 
    if [ ! -f $i ]; then
        echo -e "\tusage: $0 <DRW FILE 1> ... <DRW FILE N>"
        exit 1
    fi
    rm -f *.db.map *.db.tree
    EXIT_STATUS=1
    OUTDIR=$(grep 'mname' $i | cut -f2 -d\');
    m=0
    while [ $EXIT_STATUS -ne 0 -a $m -le $MAX_TRIALS ]; do
        echo Running ALFSIM \(iteration $m\)
        x=0
        while [ -e ${OUTDIR}_$((x+1)) ]; do
            x=$((x+1))
        done
        if [ $x -ne 0 ]; then
            OUTDIR=${OUTDIR}_$x
        fi
        if [ $m -ne 0 -a -e $OUTDIR ]; then
            echo removing $OUTDIR
            rm -rf $OUTDIR
            sleep 1;
        fi
        alfsim $i
        EXIT_STATUS=$?
        m=$((m+1))
    done
done

