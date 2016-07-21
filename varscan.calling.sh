#!/bin/bash

################################  MUTATION CALLING  ##################################################################
################################   VARSCAN.v2.4.2   ##################################################################
INPUT_FOLDER_NAME=$1
ID_LIST=$2
GATK_FOLDER_NAME=$3
OUTPUT_FOLDER_NAME=$4
TASK_ID=$5
echo "Calling mutations with varscan" 1>&2


if ! [[ -d output_$INPUT_FOLDER_NAME/$OUTPUT_FOLDER_NAME ]]; then
    mkdir output_$INPUT_FOLDER_NAME/$OUTPUT_FOLDER_NAME
    echo "made folder for varscan step" 1>&2
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



echo "Converting bam to pileup format with samtools ... " 1>&2
samtools mpileup \
    -B -q 1 \
    -f /work/scratch/vladimir/annotations/b37/human_g1k_v37_decoy.fasta \
    output_$INPUT_FOLDER_NAME/$GATK_FOLDER_NAME/${filenames[0]}.markadapt.piped.markdupl.indelrealign.recal.bam \
    output_$INPUT_FOLDER_NAME/$GATK_FOLDER_NAME/${filenames[1]}.markadapt.piped.markdupl.indelrealign.recal.bam > \
    $OUTPUT_FOLDER/${filenames[0]}.${filenames[1]}.markadapt.piped.markdupl.indelrealign.recal.bam.mpileup & wait
echo "Converting bam to pileup format with samtools ... done" 1>&2

echo "identifing somatic mutations ... " 1>&2
sh exit.code.sh $OUTPUT_FOLDER/${filenames[0]}. markadapt.piped.markdupl.indelrealign.recal.bam.pileup
varscan somatic \
$OUTPUT_FOLDER/${filenames[0]}.${filenames[1]}.markadapt.piped.markdupl.indelrealign.recal.bam.mpileup \
$OUTPUT_FOLDER/${filenames[0]}.${filenames[1]}.varscan --mpileup 1 --min-coverage 10 --min-var-freq 0.08 --somatic-p-value 0.05 --output-vcf 1 & wait;
echo "identifing somatic mutations ... done" 1>&2

echo "processing somatic mutation files ... " 1>&2
for i in $OUTPUT_FOLDER/${filenames[0]}.${filenames[1]}.*.vcf
do
    sh exit.code.sh $i
    varscan processSomatic $i &
done; wait
echo "processing somatic mutation files ... done" 1>&2




echo "filtering somatic mutations using somaticFilter ... " 1>&2
for i in $OUTPUT_FOLDER/${filenames[0]}.${filenames[1]}.*.Somatic.hc.vcf
do
    sh exit.code.sh $i
    varscan somaticFilter $i \
    -indel-file $OUTPUT_FOLDER/${filenames[0]}.${filenames[1]}.indel.vcf \
    -output-file $i.filter &
done; wait
echo "filtering somatic mutations using somaticFilter ... done" 1>&2




echo "filtering falsepositives ... " 1>&2
echo "obtaining read counts for the variants ... " 1>&2
for i in $OUTPUT_FOLDER/${filenames[0]}.${filenames[1]}.*.filter
do
    sh exit.code.sh $i
    awk '{var1=$1; var2=$2; print var1,var2,var2}' OFS='\t' $i > $i.mod
    tail -n +2 $i.mod > $i.mod.tmp
    mv $i.mod.tmp $i.mod
    rm $i.mod.tmp
    bam-readcount -q 1 -b 20 \
    -f /work/scratch/vladimir/annotations/b37/human_g1k_v37_decoy.fasta \
    -l $i.mod \
    output_$INPUT_FOLDER_NAME/$GATK_FOLDER_NAME/${filenames[1]}.markadapt.piped.markdupl.indelrealign.recal.bam > \
    $i.readcount &
done; wait
echo "obtaining read counts for the variants ... done" 1>&2




echo "running fpfilter ... " 1>&2
for i in $OUTPUT_FOLDER/${filenames[0]}.${filenames[1]}.*.filter
do
    sh exit.code.sh $i.readcount
    perl /work/software/varscan/fpfilter.pl $i $i.readcount --output-basename $i.fpfilter &
done; wait
echo "running fpfilter ... done" 1>&2
echo "filtering falsepositives ... done" 1>&2



for i in $OUTPUT_FOLDER/${filenames[0]}.${filenames[1]}.*.pass
do
sh snpeff.annotation.sh $i $TASK_ID &
done; wait




echo "detecting copynumber variation ... " 1>&2
        varscan copynumber \
        $OUTPUT_FOLDER/${filenames[0]}.${filenames[1]}.markadapt.piped.markdupl.indelrealign.recal.bam.mpileup \
        $OUTPUT_FOLDER/${filenames[0]}.${filenames[1]} \
        -min-coverage 10 --data-ratio 1.0 --min-segment-size 20 --max-segment-size 100 --mpileup 1 & wait;
echo "detecting copynumber variation ... done" 1>&2



echo "calling copynumber variation ... " 1>&2
    sh exit.code.sh $OUTPUT_FOLDER/${filenames[0]}.${filenames[1]}.copynumber
    varscan copyCaller $OUTPUT_FOLDER/${filenames[0]}.${filenames[1]}.copynumber \
    --output-file $OUTPUT_FOLDER/${filenames[0]}.${filenames[1]}.copynumber.called \
    --homdel-file $OUTPUT_FOLDER/${filenames[0]}.${filenames[1]}.copynumber.homdel & wait;
echo "calling copynumber variation ... done" 1>&2




echo "filtering true copynumber variation ... " 1>&2
    sh exit.code.sh $OUTPUT_FOLDER/${filenames[0]}.${filenames[1]}.copynumber.called
    Rscript varscan.R $OUTPUT_FOLDER/${filenames[0]}.${filenames[1]}.copynumber.called $OUTPUT_FOLDER/${filenames[0]}.${filenames[1]}.copynumber.called.segments.p_value & wait;
echo "filtering true copynumber variation ... done" 1>&2



rm $OUTPUT_FOLDER/${filenames[0]}.${filenames[1]}.markadapt.piped.markdupl.indelrealign.recal.bam.mpileup

echo "calling is done" 1>&2
