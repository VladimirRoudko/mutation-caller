#!/bin/sh
#$ -S /bin/bash
#$ -pe smp 2
#$ -l mem=80G
#$ -cwd
#$ -N gatk_processing
#$ -o /dev/null
#$ -e /dev/null


mkdir -p log
exec 1>log/gatk.$SGE_TASK_ID.out
exec 2>log/gatk.$SGE_TASK_ID.err



bam="bam"
fastq="fastq"
SE="SE"
PE="PE"
TYPE=$1
MODE=$2
SEQ=$3
RG=$4
SN=$5
LN=$6
PU=$7
P=$8
SC=$9
INPUT_FOLDER_NAME=${10}
ID_LIST=${11}
OUTPUT_FOLDER_NAME=${12}
DNAseq="DNAseq"


if ! [[ -d output_$INPUT_FOLDER_NAME/$OUTPUT_FOLDER_NAME ]]; then
    mkdir output_$INPUT_FOLDER_NAME/$OUTPUT_FOLDER_NAME
    echo "made folder for preprocessing step" 1>&2
fi
OUTPUT_FOLDER=output_$INPUT_FOLDER_NAME/$OUTPUT_FOLDER_NAME
input=data/$INPUT_FOLDER_NAME/$ID_LIST
filenames=()
while read -r line
    do
    read id flag pair <<< $line
if [[ $pair = $SGE_TASK_ID ]]; then
filenames=("${filenames[@]}" $id)
fi
done < $input;


if [[ $TYPE = $bam && $MODE = $PE ]]; then
    echo "raw reads are in .bam files" 1>&2
    echo "converting BAM to unaligned BAM (uBAM) and adding read group information ... " 1>&2
    for i in "${filenames[@]}"
    do
    sh exit.code.sh data/$INPUT_FOLDER_NAME/$i.bam
    picard RevertSam I=data/$INPUT_FOLDER_NAME/$i.bam \
                    O=$OUTPUT_FOLDER/$i.tmp.ubam \
                    SANITIZE=true \
                    ATTRIBUTE_TO_CLEAR=XT \
                    ATTRIBUTE_TO_CLEAR=XN \
                    ATTRIBUTE_TO_CLEAR=AS \
                    ATTRIBUTE_TO_CLEAR=OC \
                    ATTRIBUTE_TO_CLEAR=OP \
                    SORT_ORDER=queryname \
                    RESTORE_ORIGINAL_QUALITIES=true \
                    REMOVE_DUPLICATE_INFORMATION=true \
                    REMOVE_ALIGNMENT_INFORMATION=true \
                    TMP_DIR=$OUTPUT_FOLDER/tmp &
    done; wait
    for i in "${filenames[@]}"
    do
    sh exit.code.sh $OUTPUT_FOLDER/$i.tmp.ubam
    picard AddOrReplaceReadGroups \
        I=$OUTPUT_FOLDER/$i.tmp.ubam \
        O=$OUTPUT_FOLDER/$i.ubam \
        RGID=$RG \
        RGLB=$LN \
        RGPL=$P \
        RGPU=$PU \
        RGSM=$SN \
        TMP_DIR=$OUTPUT_FOLDER/tmp &
    done; wait
fi
rm $OUTPUT_FOLDER/*.tmp.ubam
echo "converting inupt data to unaligned BAM (uBAM) and adding read group information ... done" 1>&2


echo "marking adapter sequences ... " 1>&2
for i in "${filenames[@]}"
do
        sh exit.code.sh $OUTPUT_FOLDER/$i.ubam
        picard MarkIlluminaAdapters I=$OUTPUT_FOLDER/$i.ubam \
                                    O=$OUTPUT_FOLDER/$i.markadapt.ubam \
                                    M=$OUTPUT_FOLDER/$i.markadapt.metrics.txt \
                                    TMP_DIR=$OUTPUT_FOLDER/tmp &
done; wait
echo "marking adapter sequences ... done" 1>&2



set -o pipefail
echo "aligning reads with BWA-MEM, merging with uBAM and generating clean, aligned BAM files ..." 1>&2
sh exit.code.sh /work/scratch/vladimir/annotations/b37/human_g1k_v37_decoy.fasta
for i in "${filenames[@]}"
do
sh exit.code.sh $OUTPUT_FOLDER/$i.markadapt.ubam
picard SamToFastq I=$OUTPUT_FOLDER/$i.markadapt.ubam \
                        FASTQ=/dev/stdout \
                        CLIPPING_ATTRIBUTE=XT \
                        CLIPPING_ACTION=2 \
                        INTERLEAVE=true \
                        NON_PF=true \
                        TMP_DIR=$OUTPUT_FOLDER/tmp | \
bwa mem -M -t 10 -p /work/scratch/vladimir/annotations/b37/human_g1k_v37_decoy.fasta /dev/stdin | \
picard MergeBamAlignment ALIGNED=/dev/stdin \
                                UNMAPPED=$OUTPUT_FOLDER/$i.markadapt.ubam \
                                O=$OUTPUT_FOLDER/$i.markadapt.piped.bam \
                                R=/work/scratch/vladimir/annotations/b37/human_g1k_v37_decoy.fasta \
                                CREATE_INDEX=true \
                                ADD_MATE_CIGAR=true \
                                CLIP_ADAPTERS=false \
                                CLIP_OVERLAPPING_READS=true \
                                INCLUDE_SECONDARY_ALIGNMENTS=true \
                                MAX_INSERTIONS_OR_DELETIONS=-1 \
                                PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
                                ATTRIBUTES_TO_RETAIN=XS \
                                TMP_DIR=$OUTPUT_FOLDER/tmp &
done; wait
echo "aligning reads with BWA-MEM, merging with uBAM and generating clean, aligned BAM files ... done" 1>&2
rm $OUTPUT_FOLDER/*.ubam



echo "marking duplicates ... " 1>&2
if [[ $SEQ = $DNAseq ]]; then
        for i in "${filenames[@]}"
        do
        sh exit.code.sh $OUTPUT_FOLDER/$i.markadapt.piped.bam
        picard MarkDuplicatesWithMateCigar \
                        INPUT=$OUTPUT_FOLDER/$i.markadapt.piped.bam \
                        OUTPUT=$OUTPUT_FOLDER/$i.markadapt.piped.markdupl.bam \
                        METRICS_FILE=$OUTPUT_FOLDER/$i.markadapt.piped.markdupl.metrics.txt \
                        OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
                        CREATE_INDEX=true \
                        TMP_DIR=$OUTPUT_FOLDER/tmp &
        done; wait
fi



rm $OUTPUT_FOLDER/*.markadapt.piped.bam
rm $OUTPUT_FOLDER/*.markadapt.piped.bai
echo "marking duplicates ... done" 1>&2



echo "performing local realignment around indels ... " 1>&2



for i in "${filenames[@]}"
do
        sh exit.code.sh $OUTPUT_FOLDER/$i.markadapt.piped.markdupl.bam
        gatk -T RealignerTargetCreator \
                -R /work/scratch/vladimir/annotations/b37/human_g1k_v37_decoy.fasta -nt 24 \
                -known /work/scratch/vladimir/annotations/b37/1000G_phase1.indels.b37.vcf \
                -known /work/scratch/vladimir/annotations/b37/Mills_and_1000G_gold_standard.indels.b37.vcf \
                -I $OUTPUT_FOLDER/$i.markadapt.piped.markdupl.bam \
                -o $OUTPUT_FOLDER/$i.realignertargetcreator.intervals &
done; wait



for i in "${filenames[@]}"
do
        gatk -T IndelRealigner \
                -R /work/scratch/vladimir/annotations/b37/human_g1k_v37_decoy.fasta \
                -targetIntervals $OUTPUT_FOLDER/$id.realignertargetcreator.intervals \
                -known /work/scratch/vladimir/annotations/b37/1000G_phase1.indels.b37.vcf \
                -known /work/scratch/vladimir/annotations/b37/Mills_and_1000G_gold_standard.indels.b37.vcf \
                -I $OUTPUT_FOLDER/$i.markadapt.piped.markdupl.bam \
                -o $OUTPUT_FOLDER/$i.markadapt.piped.markdupl.indelrealign.bam &
done; wait



rm $OUTPUT_FOLDER/*.markadapt.piped.markdupl.bam
rm $OUTPUT_FOLDER/*.markadapt.piped.markdupl.bai
echo "performing local realignment around indels ... done" 1>&2




echo "analyzing patterns of covariation in the sequence dataset ... " 1>&2
for i in "${filenames[@]}"
do
        sh exit.code.sh $OUTPUT_FOLDER/$i.markadapt.piped.markdupl.indelrealign.bam &
        gatk -T BaseRecalibrator \
                -R /work/scratch/vladimir/annotations/b37/human_g1k_v37_decoy.fasta -nct 5 \
                -I $OUTPUT_FOLDER/$i.markadapt.piped.markdupl.indelrealign.bam \
                -knownSites /work/scratch/vladimir/annotations/b37/1000G_phase1.indels.b37.vcf \
                -knownSites /work/scratch/vladimir/annotations/b37/Mills_and_1000G_gold_standard.indels.b37.vcf \
                -knownSites /work/scratch/vladimir/annotations/b37/dbsnp_138.b37.excluding_sites_after_129.vcf \
                -knownSites /work/scratch/vladimir/annotations/b37/dbsnp_138.b37.vcf \
                -o $OUTPUT_FOLDER/$i.recal_data.table &
done; wait
echo "analyzing patterns of covariation in the sequence dataset ... done" 1>&2




echo "doing a second pass to analyze covariation remaining after recalibration ... " 1>&2
for i in "${filenames[@]}"
do
        sh exit.code.sh $OUTPUT_FOLDER/$i.recal_data.table &
        gatk -T BaseRecalibrator \
                -R /work/scratch/vladimir/annotations/b37/human_g1k_v37_decoy.fasta -nct 5 \
                -I $OUTPUT_FOLDER/$i.markadapt.piped.markdupl.indelrealign.bam \
                -knownSites /work/scratch/vladimir/annotations/b37/1000G_phase1.indels.b37.vcf \
                -knownSites /work/scratch/vladimir/annotations/b37/Mills_and_1000G_gold_standard.indels.b37.vcf \
                -knownSites /work/scratch/vladimir/annotations/b37/dbsnp_138.b37.excluding_sites_after_129.vcf \
                -knownSites /work/scratch/vladimir/annotations/b37/dbsnp_138.b37.vcf \
                -BQSR $OUTPUT_FOLDER/$i.recal_data.table \
                -o $OUTPUT_FOLDER/$i.post_recal_data.table &
done; wait
echo "doing a second pass to analyze covariation remaining after recalibration ... done" 1>&2




echo "generating before/after plots ... " 1>&2
for i in "${filenames[@]}"
do
        gatk -T AnalyzeCovariates \
                -R /work/scratch/vladimir/annotations/b37/human_g1k_v37_decoy.fasta \
                -before $OUTPUT_FOLDER/$i.recal_data.table \
                -after $OUTPUT_FOLDER/$i.post_recal_data.table \
                -plots $OUTPUT_FOLDER/$i.recalibration_plots.pdf &
done; wait
echo "generating before/after plots ... done" 1>&2




echo "applying recalibration to the sequence data ..." 1>&2
for i in "${filenames[@]}"
do
        sh exit.code.sh $OUTPUT_FOLDER/$i.recal_data.table &
        gatk -T PrintReads \
                -R /work/scratch/vladimir/annotations/b37/human_g1k_v37_decoy.fasta -nct 5\
                -I $OUTPUT_FOLDER/$i.markadapt.piped.markdupl.indelrealign.bam \
                -BQSR $OUTPUT_FOLDER/$i.recal_data.table \
                -o $OUTPUT_FOLDER/$i.markadapt.piped.markdupl.indelrealign.recal.bam &
done; wait
echo "applying recalibration to the sequence data ... done" 1>&2
rm $OUTPUT_FOLDER/*.markadapt.piped.markdupl.indelrealign.bam
rm $OUTPUT_FOLDER/*.markadapt.piped.markdupl.indelrealign.bai
rm data/$INPUT_FOLDER_NAME/*.bam
echo "gatk is finished for gatk.$SGE_TASK_ID" 1>&2



####  this part works. Bugs fixed.
############# Generated cleaned, re-aligned, re-calibrated bam files ready for mutation calling #####################

