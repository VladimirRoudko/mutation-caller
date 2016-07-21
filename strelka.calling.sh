#!/bin/bash
########################### strelka version: 2.0.13 ###################
INPUT_FOLDER_NAME=$1
ID_LIST=$2
GATK_FOLDER_NAME=$3
OUTPUT_FOLDER_NAME=$4
TASK_ID=$5
echo "Calling mutations with strelka" 1>&2




if ! [[ -d output_$INPUT_FOLDER_NAME/$OUTPUT_FOLDER_NAME ]]; then
    mkdir output_$INPUT_FOLDER_NAME/$OUTPUT_FOLDER_NAME
    echo "made folder for strelka step" 1>&2
fi
OUTPUT_FOLDER=output_$INPUT_FOLDER_NAME/$OUTPUT_FOLDER_NAME




input=data/$INPUT_FOLDER_NAME/$ID_LIST
filenames=()
while read -r line
    do
    read id flag pair <<< $line
if [[ $pair = $TASK_ID ]]; then
filenames=("${filenames[@]}" $id)
fi
done < $input;




echo "calling with strelka ... " 1>&2
cp output_$INPUT_FOLDER_NAME/$GATK_FOLDER_NAME/${filenames[0]}.markadapt.piped.markdupl.indelrealign.recal.bai output_$INPUT_FOLDER_NAME/$GATK_FOLDER_NAME/${filenames[0]}.markadapt.piped.markdupl.indelrealign.recal.bam.bai
cp output_$INPUT_FOLDER_NAME/$GATK_FOLDER_NAME/${filenames[1]}.markadapt.piped.markdupl.indelrealign.recal.bai output_$INPUT_FOLDER_NAME/$GATK_FOLDER_NAME/${filenames[1]}.markadapt.piped.markdupl.indelrealign.recal.bam.bai



STRELKA_INSTALL_DIR=/work/software/strelka
if ! [[ -f $OUTPUT_FOLDER/config.ini ]]; then
cp $STRELKA_INSTALL_DIR/etc/strelka_config_bwa_default.ini $OUTPUT_FOLDER/config.ini
fi
$STRELKA_INSTALL_DIR/bin/configureStrelkaWorkflow.pl \
--normal=output_$INPUT_FOLDER_NAME/$GATK_FOLDER_NAME/${filenames[0]}.markadapt.piped.markdupl.indelrealign.recal.bam \
--tumor=output_$INPUT_FOLDER_NAME/$GATK_FOLDER_NAME/${filenames[1]}.markadapt.piped.markdupl.indelrealign.recal.bam \
--ref=/work/scratch/vladimir/annotations/b37/human_g1k_v37_decoy.fasta \
--config=$OUTPUT_FOLDER/config.ini --output-dir=$OUTPUT_FOLDER/${filenames[0]}.${filenames[1]}.$TASK_ID.analysis & wait;




cd $OUTPUT_FOLDER/${filenames[0]}.${filenames[1]}.$TASK_ID.analysis
make -j 8 & wait;
cd ../../../
echo "calling with strelka ... done " 1>&2



rm $OUTPUT_FOLDER/${filenames[0]}.${filenames[1]}.$TASK_ID.analysis/results/all.*

for i in $OUTPUT_FOLDER/${filenames[0]}.${filenames[1]}.$TASK_ID.analysis/results/*.vcf
do
mv $i ${filenames[0]}.${filenames[1]}.strelka.$i
done;

for i in $OUTPUT_FOLDER/${filenames[0]}.${filenames[1]}.$TASK_ID.analysis/results/*.vcf
do
    sh snpeff.annotation.sh $i $TASK_ID &
done; wait

echo "calling is done" 1>&2




