#!/bin/bash

TAXA_NO="$1"
EXTEDG="$2"
CHECK_INT1=`echo "$TAXA_NO" | egrep "^[0-9]+$"`
CHECK_FLOAT1=`echo "$EXTEDG" | egrep "^[0-9]+(\.[0-9]*)?$"`

if [ -z "$CHECK_INT1" -o -z "$CHECK_FLOAT1" ]; then
    echo "usage: $0 <NUMBER OF TAXA> <OUTER EDGE LENGTH>"
    exit 1
fi

# XXX PARAMETERS

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


# recursive tree reconstruction algorithm
function rec_construct_tree {
    # parameters: 
    #   1: taxa index
    #   2: current index
    #   3: count index (needs to be 1)

    COUNT=`bc <<< "$3 + 1"`
    EDGLEN=\$INTEDG
    if [ "$3" -eq "-1" ]; then
        # reset count
        COUNT=2
        EDGLEN=\$HLFEDG
    fi

    # there are 4 states in tree construction
    case $1 in
        1)
            TREE="$2"
            if [ $3 -eq 1 ]; then
                TREE="$TREE\;"
            else
                TREE="$TREE:\$EXTEDG"
            fi
            echo $TREE
            exit
            ;;
        2)
            FST=`cut -d',' -f1 <<< $2`
            SND=`cut -d',' -f2 <<< $2`
            TREE="\($FST:\$EXTEDG, $SND:\$EXTEDG\)"
            if [ $3 -eq 1 ]; then
                TREE="$TREE\;"
            else
                TREE="$TREE:$EDGLEN"
            fi
            echo $TREE
            exit
            ;;
        4)
            FST=`cut -d',' -f1 <<< $2`
            SND=`cut -d',' -f2 <<< $2`
            TRD=`cut -d',' -f3 <<< $2`
            FTH=`cut -d',' -f4 <<< $2`
            if [ "$3" -eq "1" ]; then
                TREE="\(\($FST:\$EXTEDG, $SND:\$EXTEDG\):\$HLFEDG, \($TRD:\$EXTEDG,$FTH:\$EXTEDG\):\$HLFEDG\)\;"
            else
                TREE="\(\($FST:\$EXTEDG, $SND:\$EXTEDG\):$EDGLEN,\($TRD:\$EXTEDG,$FTH:\$EXTEDG\):$EDGLEN\):$EDGLEN"
            fi
            echo $TREE
            exit
            ;;
        *)
            TXIND=`bc <<< "$1 / 3"`
            SEQ=`seq -s',' 1 $TXIND`
            CURTX=`cut -d',' -f$SEQ <<< $2`
            if [ "$3" -eq 1 ]; then
                FST=`rec_construct_tree $TXIND $CURTX "-1"`
            else
                FST=`rec_construct_tree $TXIND $CURTX $COUNT`
            fi
            STRT=`bc <<< "$TXIND + 1"`
            STOP=`bc <<< "$TXIND * 2 + $1 % 3"`
            SEQ=`seq -s',' $STRT $STOP`
            CURTX=`cut -d',' -f$SEQ <<< $2`
            MDL=`bc <<< "$TXIND + $1 % 3"`
            SND=`rec_construct_tree $MDL $CURTX $COUNT`
            STRT=`bc <<< "$STOP + 1"`
            SEQ=`seq -s',' $STRT $1`
            CURTX=`cut -d',' -f$SEQ <<< $2`
            TRD=`rec_construct_tree $TXIND $CURTX $COUNT`
            if [ "$3" -eq 1 ]; then
                echo "\($FST, \($SND, $TRD\):\$HLFEDG\)\;"
            else
                echo "\($FST, \($SND, $TRD\):\$INTEDG\):$EDGLEN"
            fi
            ;;
    esac
}

TAXA=`seq -s ',' 1 $TAXA_NO`
NEWICK_TREE=`rec_construct_tree $TAXA_NO "$TAXA" 1` 
# give a gradient of different edge lengths
INTEDG=`bc <<< "scale=4; $EXTEDG / 5"`
INTEDG=`appendZero $INTEDG`
HLFEDG=`bc <<< "scale=4; $INTEDG / 2"`
HLFEDG=`appendZero $HLFEDG`
GEN_TREE=`eval echo $NEWICK_TREE`
EXTEDG=`bc <<< "scale=4; $EXTEDG"`
echo $GEN_TREE

