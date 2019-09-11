#!/bin/bash

TAXA_NO="$1"
EXTEDG="$2"
CHECK_INT1=`echo "$TAXA_NO" | egrep "^[0-9]+$"`
CHECK_FLOAT1=`echo "$EXTEDG" | egrep "^[0-9]+(\.[0-9]*)?$"`

if [ -z "$CHECK_INT1" -o -z "$CHECK_FLOAT1" ]; then
    echo "usage: $0 <NUMBER OF TAXA> <OUTER EDGE LENGTH>"
    exit 1
fi

function appendZero {
    FCHAR=`head -c1 <<< "$1"`

    if [ "$FCHAR" = "." ]; then
        echo "0$1"
    elif [ -z `echo "$1" | grep '\.'` ]; then
	echo "$1.0"
    else
        echo "$1"
    fi
}

NEWICK_TREE="1:\$EXTEDG"
STP=`bc <<< "$TAXA_NO - 2"`
for i in `seq 2 $STP`;
do
    if [ $i -eq $STP ]; then
        NEWICK_TREE="\($NEWICK_TREE, $STP:\$EXTEDG\)"
        SNDLST=`bc <<< "$TAXA_NO -1"`
        NEWICK_TREE="\(\($SNDLST:\$EXTEDG, $TAXA_NO:\$EXTEDG\):\$HLFEDG,$NEWICK_TREE:\$HLFEDG\)\;"
    else
        NEWICK_TREE="\($NEWICK_TREE, $i:\$EXTEDG\):\$INTEDG"
    fi 
done
INTEDG=`bc <<< "scale=4; $EXTEDG / 5"`
INTEDG=`appendZero $INTEDG`
HLFEDG=`bc <<< "scale=4; $INTEDG / 2"`
HLFEDG=`appendZero $HLFEDG`
GEN_TREE=`eval echo $NEWICK_TREE`
echo $GEN_TREE
