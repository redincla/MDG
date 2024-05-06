#!/bin/bash

###############################################################################
#   Copy .bam and .bai files locally for tools that run in containers   	  #
###############################################################################

# Copyright (c) 2024 Claire Redin @medigenome, geneva, Switzerland.
# Contact: Claire Redin <claire.redin@medigenome.ch>


# checks for appropriate input
if [ $# -eq 2 ]; then
 saphetor_bam_map=$1 #full path to file: 2-col map with: col1: ID, col2: path-to-bam
 WRKDIR=$2
else
 echo -e "\n\nCopy .bam and .bai files from saphetor folder location locally for tools that need a local copy. \n DO NOT FORGET TO THEN REMOVE FILES WHEN DONE, OTHERWISE IT WILL EAT UP SPACE\n\nAuthor: Claire Redin (claire.redin@medigenome.ch)\n\n"
 echo "Usage:"
 echo "copy_bam_from_saphetor.sh [sample_list] [WRKDIR]"
 echo "saphetor_bam_map: full path to file: 2-col map with: col1: ID, col2: path-to-bam"
 echo "full path where to drop the bam files"
 exit 1
fi


# Need to change the header to have MDG sample ID instead of a random Saphetor ID that would appear in the output plots from Somalier
while read ID bam_file; do
samtools view -H ${bam_file} | sed "s/^\(@RG\t.*\tSM:\).*\(\t.*\)$/\1${ID}\2/" | samtools reheader - ${bam_file} > ${WRKDIR}/${ID}.recalibrated.renamed.bam
samtools index -b ${WRKDIR}/${ID}.recalibrated.renamed.bam
done < "$saphetor_bam_map"