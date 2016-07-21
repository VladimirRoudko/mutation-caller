#!/bin/bash
ID_LIST_PATH=$1
CheckSamples() {
path=$1
count=0
pairs=0
pair1=""
while read -r line
do
    read id flag pair <<< $line
    if [[ $pair1 == $pair ]]; then
    pairs=$((pairs+1))
    else
        count=$((count+1))
    fi
    pair1=$pair
done < $path;
if [[ $count = $pairs ]]; then
    echo $count
else
    echo "amount of unique files is different from amount of tumor-normal pairs" >> main.log
fi
}
CheckSamples $ID_LIST_PATH




