#!/bin/bash
## This version designed to run with qsub module with grid array job definition.
## Use amount of jobs according to amount of samples to analyse.
#$ -S /bin/bash
#$ -pe smp 2
#$ -l mem=40G
#$ -cwd
#$ -N mutation-calling
#$ -o /dev/null
#$ -e /dev/null


mkdir -p log
exec 1>log/mutation.$SGE_TASK_ID.out
exec 2>log/mutation.$SGE_TASK_ID.err




INPUT_FOLDER_NAME=$1
ID_LIST=$2
sh varscan.calling.sh $INPUT_FOLDER_NAME $ID_LIST gatk_preprocessing varscan_output $SGE_TASK_ID
echo "varscan.caller is done" 1>&2
sh strelka.calling.sh $INPUT_FOLDER_NAME $ID_LIST gatk_preprocessing strelka_output $SGE_TASK_ID
echo "strelka.caller is done" 1>&2
#sh somatic-sniper.calling.sh $INPUT_FOLDER_NAME $ID_LIST gatk_preprocessing somatic-sniper_output varscan_output $SGE_TASK_ID
#echo "somatic-sniper is done" 1>&2
sh mutect.calling.sh $INPUT_FOLDER_NAME $ID_LIST gatk_preprocessing mutect_output $SGE_TASK_ID
echo "mutect is done" 1>&2


