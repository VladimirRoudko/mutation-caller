#!/bin/bash

INPUT_VCF=$1
pair=$2
echo "annotating resulted vcf file with snpEff ... " 1>&2

snpeff -c /work/software/snpeff/snpEff.config -v GRCh37.75 $INPUT_VCF > $INPUT_VCF.snpeff_annotation.$pair.vcf & wait;

echo "annotating resulted vcf file with snpEff ... done" 1>&2



