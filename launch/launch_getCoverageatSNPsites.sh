#!/bin/bash

###########################################################
#   get  coverage at specific SNP sites   	  #
###########################################################

# Copyright (c) 2024 Claire Redin @medigenome, geneva, Switzerland.
# Contact: Claire Redin <claire.redin@medigenome.ch>


# checks for appropriate input
if [ $# -eq 2 ]; then
 sample_bams=$1 #full path to bam files. col1: sample ID, col2: paht to bam
 SNP_bed=$2 #full path to SNP bed file

else
 echo -e "\n\nGet depth of coverage from bam file at spectific SNP sites\n\nAuthor: Claire Redin (claire.redin@medigenome.ch)\n\n"
 echo "Usage:"
 echo "launch_getCoverageatSNPsites.sh [sample_bams] [SNP_bed]"
 echo "sample_bams: 2 columns files. col1: sample ID, col2: path to bam files with bai files in same folder"
 echo "SNP_bed: path to bed file ( /!\ same genome version as bam file)"
 exit 1
fi

while read ID bam; do
samtools depth -a ${bam} \
	-b ${SNP_bed} -H -Q 0 \
	-o ${ID}.cov 
done < "$sample_bams"

