#!/bin/bash
INPUT_FOLDER_NAME=$1
ID_LIST=$2
GATK_FOLDER_NAME=$3
OUTPUT_FOLDER_NAME=$4
TASK_ID=$5
echo "Calling mutations with mutect2" 1>&2



if ! [[ -d output_$INPUT_FOLDER_NAME/$OUTPUT_FOLDER_NAME ]]; then
    mkdir output_$INPUT_FOLDER_NAME/$OUTPUT_FOLDER_NAME
    echo "made folder for mutect2 step" 1>&2
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



echo "mutect2 calling ... " 1>&2
gatk -T MuTect2 \
    -R /work/scratch/vladimir/annotations/b37/human_g1k_v37_decoy.fasta -nct 10 \
    -I:tumor output_$INPUT_FOLDER_NAME/$GATK_FOLDER_NAME/${filenames[1]}.markadapt.piped.markdupl.indelrealign.recal.bam \
    -I:normal output_$INPUT_FOLDER_NAME/$GATK_FOLDER_NAME/${filenames[0]}.markadapt.piped.markdupl.indelrealign.recal.bam \
    --dbsnp /work/scratch/vladimir/annotations/b37/dbsnp_138.b37.vcf
    -o $OUTPUT_FOLDER/${filenames[0]}.${filenames[1]}.mutect2.vcf & wait;


grep "#" log/mutation.$TASK_ID.out >> $OUTPUT_FOLDER/${filenames[0]}.${filenames[1]}.mutect2.called.vcf
grep "PASS" log/mutation.$TASK_ID.out >> $OUTPUT_FOLDER/${filenames[0]}.${filenames[1]}.mutect2.called.vcf
cp log/mutation.$TASK_ID.out $OUTPUT_FOLDER/${filenames[0]}.${filenames[1]}.mutect2.all.vcf

echo "mutect2 calling ... done" 1>&2



sh snpeff.annotation.sh $OUTPUT_FOLDER/${filenames[0]}.${filenames[1]}.mutect2.called.vcf $TASK_ID & wait;

echo "calling is done" 1>&2
echo "mutect2 calling ... done" 1>&2

