#!/bin/bash

###########################################################
#   launch somalier on a batch of bam samples		   	  #
###########################################################

# Copyright (c) 2024 Claire Redin @medigenome, geneva, Switzerland.
# Contact: Claire Redin <claire.redin@medigenome.ch>
# somalier: https://github.com/brentp/somalier

# checks for appropriate input
if [ $# -eq 4 ]; then
	WRKDIR=$1
	sample_list_name=$2
	ped_file_name=$3
	genome=$4

	else
 echo -e "\n\nPerform relatedness, sex and ancestry QC-check a list of bam samples\nAuthor: Claire Redin (claire.redin@medigenome.ch)\n"
 echo -e "Usage:\n"
 echo -e "launch-somalier.sh [WRKDIR] [sample_list_name] [ped_file_name] [genome]\n\n"
 echo -e "- WRKDIR: folder of path to be mounted in docker image. Should contain all necessary files for somalier to run (SNP sites, genome .fa file, sample_list file). Output results will be dropped here as well.\n"
 echo -e "- sample_list_name: name of file containing sample list copied in WRKDIR folder\n"
 echo -e "- ped_file_name: name of 6-columns ped file for cohort copied in WRKDIR folder\n"
 echo -e "- genome version: version of the genome hg19/hg38 \n"
 exit 1
fi

###extract SNP sites from bam to perform relatedness upon
cd ${WRKDIR}
if ! [ -e ${WRKDIR}/somalier_output ]; then
	mkdir ${WRKDIR}/somalier_output
fi

### cp necessary somalier files temporaly within WRKDIR and check genome version
if [[ "${genome}" == "hg19" || "${genome}" == "19" ]]; then
	cp /mnt/raw-data/references/genome/hg19/hg19.fa.gz genome.fa.gz
	cp /mnt/raw-data/references/genome/hg19/hg19.fa.gz.fai genome.fa.gz.fai
	cp /mnt/raw-data/references/somalier_sites/sites.hg19.vcf.gz sites.vcf.gz
	cp -r /mnt/raw-data/references/somalier_sites/1kg-somalier .

elif [[ "${genome}" == "hg38" || "${genome}" == "38" ]]; then
	cp /mnt/raw-data/references/genome/hg38/hg38.fa.gz genome.fa.gz
	cp /mnt/raw-data/references/genome/hg38/hg38.fa.gz.fai genome.fa.gz.fai
	cp /mnt/raw-data/references/somalier_sites/sites.hg38.vcf.gz sites.vcf.gz
	cp -r /mnt/raw-data/references/somalier_sites/1kg-somalier .

else
    echo -e "Genome version incorrect: accepted entries are [hg19/19] or [hg38/38]"
    exit 1

fi


### run somalier
while read bam; do
docker run -v ${WRKDIR}:/mnt/files brentp/somalier:latest somalier extract -d /mnt/files/somalier_output --sites /mnt/files/sites.vcf.gz -f /mnt/files/genome.fa.gz ${bam}
done < ${WRKDIR}/${sample_list_name}

#### perform relatedness analyse
docker run -v ${WRKDIR}:/mnt/files -w /mnt/files brentp/somalier:latest somalier relate --ped /mnt/files/${ped_file_name} /mnt/files/somalier_output/*somalier


#### assign ancestry based on ancestry from 1KG data
docker run -v ${WRKDIR}:/mnt/files -w /mnt/files brentp/somalier:latest somalier ancestry --labels ../../ancestry-labels-1kg.tsv /mnt/files/1kg-somalier/*.somalier ++ /mnt/files/somalier_output/*.somalier
