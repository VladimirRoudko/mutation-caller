#!/bin/bash

################################# Creates output folder:
INPUT_FOLDER_NAME=$1
mkdir output_$INPUT_FOLDER_NAME
echo "created output folder output_$INPUT_FOLDER_NAME ... done" >> main.log



if ! [[ -d "output_$INPUT_FOLDER_NAME"  ]] ; then
echo "error: output_$INPUT_FOLDER_NAME is not created" >> main.log
exit
fi
