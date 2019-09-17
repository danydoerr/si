#!/bin/sh

NOT_SUCCESS=1
OUT=""
TRIAL=1
echo "trial" 1>&2
while [ "$NOT_SUCCESS" != "0" ]; do
    echo "$TRIAL" 1>&2
    #OUT=$(/prj/si/Transposition-Median/grappa -t3 -T3 -l -m -b 10000 -f $1 2> /dev/null)
    OUT=$(/prj/si/Transposition-Median/grappa -t3 -T3 -l -m -b 300 -J -f $1 2> /dev/null)
    NOT_SUCCESS=$?
    TRIAL=$(($TRIAL + 1))
done
echo "... success!" 1>&2
echo "$OUT"
