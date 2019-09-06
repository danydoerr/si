#!/bin/bash

if [ $# -lt 1 -o ! -f $1 ]; then
    echo -e "\tusage: $0 <LOGFILE 1> ... <LOGFILE N>"
    exit 1
fi

calc(){ awk "BEGIN { print $* }"; }

echo -e ORGANISM\\tINVERSION\\tTRANSPOSITION\\tDUPLICATION\\tLOSS\\tHGT

for i in $(seq 15); do
    TRANSLOC=0
    INVERSION=0
    DUPL=0
    LOSS=0
    HGT=0
    for f in $@; do
        TRANSLOC_i=$(grep -ce "translocation of size [0-9]\+ in organism $i" $f)
        INVERSION_i=$(grep -ce "inversion in organism $i" $f)
        DUPL_i=$(grep -ce "gene duplication in organism $i with gene" $f)
        LOSS_i=$(grep -ce "gene loss in organism $i with gene" $f)
        HGT_i=$(grep -ce "lgt from organism $i" $f)
        TRANSLOC=$(($TRANSLOC + $TRANSLOC_i))
        INVERSION=$(($INVERSION + $INVERSION_i))
        DUPL=$(($DUPL + $DUPL_i))
        LOSS=$(($LOSS + $LOSS_i))
        HGT=$(($HGT + $HGT_i))
    done
    TRANSLOC=$(calc $TRANSLOC./$#)
    INVERSION=$(calc $INVERSION/$#)
    DUPL=$(calc $DUPL/$#)
    LOSS=$(calc $LOSS/$#)
    HGT=$(calc $HGT/$#)
    echo -e $i\\t$INVERSION\\t$TRANSLOC\\t$DUPL\\t$LOSS\\t$HGT
done

