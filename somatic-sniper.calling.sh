#!/bin/bash

### Somatic sniper is highly dependent on the output of varscan.  It uses varscan filtering process extensively
INPUT_FOLDER_NAME=$1
ID_LIST=$2
GATK_FOLDER_NAME=$3
OUTPUT_FOLDER_NAME=$4
VARSCAN_FOLDER_NAME=$5
TASK_ID=$6
echo "Calling mutations with somatic-sniper" 1>&2


if ! [[ -d output_$INPUT_FOLDER_NAME/$OUTPUT_FOLDER_NAME ]]; then
    mkdir output_$INPUT_FOLDER_NAME/$OUTPUT_FOLDER_NAME
    echo "made folder for somatic-sniper step" 1>&2
fi
OUTPUT_FOLDER=output_$INPUT_FOLDER_NAME/$OUTPUT_FOLDER_NAME
OUTPUT_FOLDER_VARSCAN=output_$INPUT_FOLDER_NAME/$VARSCAN_FOLDER_NAME


input=data/$INPUT_FOLDER_NAME/$ID_LIST
filenames=()
while read -r line
    do
    read id flag pair <<< $line
if [[ $pair = $TASK_ID ]]; then
filenames=("${filenames[@]}" $id)
fi
done < $input;


bam-somaticsniper -q 1 -Q 40 -G -L -F vcf \
    -f /work/scratch/vladimir/annotations/b37/human_g1k_v37_decoy.fasta \
    output_$INPUT_FOLDER_NAME/$GATK_FOLDER_NAME/${filenames[1]}.markadapt.piped.markdupl.indelrealign.recal.bam \
    output_$INPUT_FOLDER_NAME/$GATK_FOLDER_NAME/${filenames[0]}.markadapt.piped.markdupl.indelrealign.recal.bam \
    $OUTPUT_FOLDER/${filenames[0]}.${filenames[1]}.somaticsniper.snv.vcf & wait;


samtools mpileup -A -B -g -u \
    -f /work/scratch/vladimir/annotations/b37/human_g1k_v37_decoy.fasta \
    output_$INPUT_FOLDER_NAME/$GATK_FOLDER_NAME/${filenames[1]}.markadapt.piped.markdupl.indelrealign.recal.bam | \
    bcftools view -c /dev/stdin | \
    vcfutils.pl varFilter -Q 20 | \
    awk 'NR > 55 {print}' > $OUTPUT_FOLDER/${filenames[1]}.indels.pileup &


perl /work/software/somatic-sniper/src/scripts/snpfilter.pl \
    --snp-file $OUTPUT_FOLDER/${filenames[0]}.${filenames[1]}.somaticsniper.snv.vcf \
    --indel-file ${filenames[1]}.indels.pileup \
    --out-file ${filenames[0]}.${filenames[1]}.somaticsniper.snv.vcf.filter

until [[ -f $OUTPUT_FOLDER_VARSCAN/${filenames[0]}.${filenames[1]}.*.fpfilter ]];
do
    wait 300;
done

perl /work/software/somatic-sniper/src/scripts/fpfilter.pl \
    -snp-file $OUTPUT_FOLDER/${filenames[0]}.${filenames[1]}.somaticsniper.snv.vcf.filter \
    -readcount-file $OUTPUT_FOLDER_VARSCAN/${filenames[0]}.${filenames[1]}.readcount & wait;

perl /work/software/somatic-sniper/src/scripts/highconfidence.pl \
    -snp-file $OUTPUT_FOLDER/${filenames[0]}.${filenames[1]}.somaticsniper.snv.vcf.filter.fp_pass \
    --out-file $OUTPUT_FOLDER/${filenames[0]}.${filenames[1]}.somaticsniper.snv.vcf.filter.fp_pass.hc & wait;
    

sh snpeff.annotation.sh \
    $OUTPUT_FOLDER/${filenames[0]}.${filenames[1]}.somaticsniper.snv.vcf.filter.fp_pass.hc $TASK_ID & wait;
echo "calling is done" 1>&2

