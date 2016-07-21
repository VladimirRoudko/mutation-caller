#!/bin/bash
FILE=$1
if ! [[ -f $FILE ]]; then
    echo "file $FILE does not exist..." 1>&2
    exit
fi


