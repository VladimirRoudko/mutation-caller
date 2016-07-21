#!/bin/sh
#$ -S /bin/bash
#$ -pe smp 1
#$ -l mem=5G
#$ -cwd
#$ -N main.script
#$ -o /dev/null
#$ -e /dev/null
## This version designed to run with qsub module with grid array job definition.
## Use amount of jobs according to amount of samples to analyse.
## also divide pipeline on several files depending on job they do.
## run the command:
## qsub main.sh -i /work/DATA2/vova/data/samples_2 -l samples.txt -t bam -m PE -n DNAseq -rg RG -sn SN -ln LN -pu PU -p illumina -sc SC &     


source /etc/profile.d/sge.sh
mkdir -p log
exec 1>log/main.out
exec 2>log/main.err
help1="-h"
help2="--help"



if [[ -z $1 ]]
then
        echo "for help type -h or --help"
        exit 113
elif [[ $1 = $help1 || $1 = $help2 ]]
then
        echo "Required parameters to run the pipeline: "
        echo "-i <path_to_folder>       path to the folder with the input files"
        echo "-l <file_list_of_id>      Three-column, tab-delimited .txt file with the list of tumor/normal datasets names."
        echo "                                                  First column: name of read filenames without extensions."
        echo "                                                  Second column: tumor/normal flag (T/N)."
        echo "                                                  Third column: Number, matching pairs of datasets (1/1, 2/2 and etc)"
        echo "                                                  important to keep the order: first file is Normal, second - Tumor"
        echo "-t <fastq/bam>            type of input files: either fastq or bam"
        echo "-m <SE/PE>                mode of sequencing: single-end or pair-end respectively"
        echo "-n <RNAseq/DNAseq>        type of sequencing: RNAseq or DNAseq"
        echo "-rg <read_group_name>     provide the read group name"
        echo "-sn <sample_name>         sample name"
        echo "-ln <library_name>        provide library name"
        echo "-pu <platform_unit>       provide the platform unit"
        echo "-p <platform>             provide the platform type: illumina/solexa/helios and etc"
        echo "-sc <sequence_center>     name of sequencing center"
        echo
        exit 113
fi



INPUT_FOLDER_PATH=$2
INPUT_FOLDER_NAME=$(echo $INPUT_FOLDER_PATH | awk -F'[/]' '{print $NF}')
if [[ -z $INPUT_FOLDER_NAME ]]; then
    INPUT_FOLDER_NAME=$(echo $INPUT_FOLDER_PATH | awk -F'[/]' '{print $(NF-1)}')
    len=${#INPUT_FOLDER_PATH}
    INPUT_FOLDER_PATH=${INPUT_FOLDER_PATH::len-1}
fi
ID_LIST=$4
TYPE=$6
MODE=$8
SEQ=${10}
RG=${12}
SN=${14}
LN=${16}
PU=${18}
P=${20}
SC=${22}


sh copy.files.sh $INPUT_FOLDER_PATH $INPUT_FOLDER_NAME


sh create.output.sh $INPUT_FOLDER_NAME



declare -i jobs
jobs=$(sh check.sample.list.sh data/$INPUT_FOLDER_NAME/$ID_LIST)
echo "jobs: $jobs" >> main.log
re='^[0-9]+$'
if ! [[ $jobs =~ $re ]] ; then
   echo "error: $jobs" >> main.log
   exit
fi


echo "running preprocessing script pipeline.v3.sh using input bam files ... " >> main.log
qsub -t 1-$jobs gatk.processing.sh $TYPE $MODE $SEQ $RG $SN $LN $PU $P $SC $INPUT_FOLDER_NAME $ID_LIST gatk_preprocessing
echo "running preprocessing script pipeline.v3.sh using input bam files ... done" >> main.log




check="gatk is finished for gatk.$jobs"
if [[ -f log/gatk.$jobs.err ]]; then
   key=$(tail -n 1 log/gatk.$jobs.err)
until [[ $key = $check ]];
do
    "waiting for the last bam file to be processed by gatk ... " >> main.log
    sleep 600
done
echo "running mutation callers with gatk-processed files ... " >> main.log
qsub -t 1-$jobs mutation-caller.sh $INPUT_FOLDER_NAME $ID_LIST
fi


check2="mutect2 calling ... done"
if [[-f log/mutation.$jobs.err ]]; then
mkey=$(tail -n 1 log/mutation.$jobs.err)


until [[ $mkey = $check2 ]];
do
    "waiting for the last caller to finish its job ... " >> main.log
    sleep 600
done
echo "running mutation callers with gatk-processed files ... done" >> main.log




echo "filtering robust mutations, called by all 4 callers ... " >> main.log

echo "copying the annotated mutation files to the common folder ... " >> main.log
if ! [[ -d output_$INPUT_FOLDER_NAME/annotated_mutations
 mkdir output_$INPUT_FOLDER_NAME/annotated_mutations
fi


cp output_$INPUT_FOLDER_NAME/varscan_output/*.snpeff_annotation.* output_$INPUT_FOLDER_NAME/annotated_mutations/
cp output_$INPUT_FOLDER_NAME/stelka_output/*.analysis/results/*.snpeff_annotation.* output_$INPUT_FOLDER_NAME/annotated_mutations/
cp output_$INPUT_FOLDER_NAME/somatic-sniper_output/*.snpeff_annotation.* output_$INPUT_FOLDER_NAME/annotated_mutations/
cp output_$INPUT_FOLDER_NAME/mutect_output/*.snpeff_annotation.* output_$INPUT_FOLDER_NAME/annotated_mutations/


############################## add common tags for snps and indels at the end of each file ###########################


: '

THIS PART NEEDS TO BE DEFINED.

############################### select common mutations among snps and indels separately: #######################################
file_count=$(ls -lrt file? |wc -l)
file1=$(ls -1 *.vcf | head -1)
sort -u file1 > temp;cat temp > file1;rm temp
while read i
do
result_count=$(grep -lw "$i" file? | wc -l)
if [ $result_count -eq $file_count ]; then
 echo $i >> common.mutations.vcf
fi
done < file1




'


