#!/bin/bash
### script copies input folder to the ../data folder
INPUT_FOLDER_PATH=$1
INPUT_FOLDER_NAME=$2
mkdir data
cp -avr $INPUT_FOLDER_PATH data/$INPUT_FOLDER_NAME & wait;
echo "printed the input folder to the data folder ... done" >> main.log
