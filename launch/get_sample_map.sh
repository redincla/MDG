#!/bin/bash

###########################################################
#   get  map location in saphetor for MDG samples   	  #
###########################################################

# Copyright (c) 2024 Claire Redin @medigenome, geneva, Switzerland.
# Contact: Claire Redin <claire.redin@medigenome.ch>


# checks for appropriate input
if [ $# -eq 1 ]; then
 sample_list=$1 #full path to file: list of samples to search Bam files from

else
 echo -e "\n\nFind bam folder location in Saphetor mount for a list of MDG samples\n\nAuthor: Claire Redin (claire.redin@medigenome.ch)\n\n"
 echo "Usage:"
 echo "get_sample_map.sh [sample_list]"
 echo "sample_list: list of samples IDs"
 exit 1
fi


WRKDIR=$(pwd)
cd /mnt/analyses/

while read ID; do
folder=$(find . -name "*-${ID}_*" -print -quit | cut -f 1-4 -d/)
printf "%s\t%s/bam/recalibrated.bam\ttest\n" "${ID}" "${folder}" >> ${WRKDIR}/saphetor_map
done < ${WRKDIR}/$1
