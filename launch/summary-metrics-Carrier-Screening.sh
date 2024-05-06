#!/bin/bash

########################
#   gather WGS metrics   #
########################

# Copyright (c) 2024 Claire Redin @medigenome, geneva, Switzerland.
# Contact: Claire Redin <claire.redin@medigenome.ch>

# checks for appropriate input
if [ $# -eq 2 ]; then
 sample_list=$1 #full path to file: list of samples to gather NGS metrics on
 WRKDIR=$2 #full directory where QC files are

else
 echo -e "\n\nGather NGS summary metrics by sample\n\nAuthor: Claire Redin (claire.redin@medigenome.ch)\n\n"
 echo "Usage:"
 echo "summary-metrics-Carrier-Screening.sh [sample_list] [working directory] "
 echo "sample_list: full path to file: list of samples to gather WGS metrics on"
 echo "working directory: full path"
 exit 1
fi

cd $WRKDIR

while read sample_ID; do
        num_TOTAL_READS=$(sed -n 2p *${sample_ID}_Sample-Metrix-QC.txt  | cut -f 2)
        num_READ_PAIRS=$(echo "scale=0;$num_TOTAL_READS / 2" | bc -l)
		pct_READS_ABOVE_Q30=$(sed -n 4p *${sample_ID}_Sample-Metrix-QC.txt  | cut -f 2)
		MEAN_QUALITY=$(sed -n 5p *${sample_ID}_Sample-Metrix-QC.txt  | cut -f 2)
		MEDIAN_INSERT_SIZE=$(sed -n 6p *${sample_ID}_Sample-Metrix-QC.txt  | cut -f 2)		
		pct_UNIQUE_READS=$(sed -n 9p *${sample_ID}_Sample-Metrix-QC.txt  | cut -f 2)
		pct_MAPPED_READS=$(sed -n 10p *${sample_ID}_Sample-Metrix-QC.txt  | cut -f 2)
		pct_MAPPED_OFF_TARGET=$(sed -n 13p *${sample_ID}_Sample-Metrix-QC.txt  | cut -f 2)
		pct_MAPPED_ON_TARGET=$(echo "scale=2;100-$pct_MAPPED_OFF_TARGET" | bc -l)
		pct_TARGET_1X=$(sed -n 14p *${sample_ID}_Sample-Metrix-QC.txt  | cut -f 3)
		pct_TARGET_10X=$(sed -n 15p *${sample_ID}_Sample-Metrix-QC.txt  | cut -f 3)
		pct_TARGET_20X=$(sed -n 16p *${sample_ID}_Sample-Metrix-QC.txt  | cut -f 3)
		pct_TARGET_30X=$(sed -n 17p *${sample_ID}_Sample-Metrix-QC.txt  | cut -f 3)
		pct_TARGET_50X=$(sed -n 18p *${sample_ID}_Sample-Metrix-QC.txt  | cut -f 3)		
		MEDIAN_TARGET_COV=$(sed -n 20p *${sample_ID}_Sample-Metrix-QC.txt | cut -f 2)	
		Ti_Tv=$(sed -n 22p *${sample_ID}_Sample-Metrix-QC.txt | cut -f 2)	
		CONTAMINATION=$(sed -n 23p *${sample_ID}_Sample-Metrix-QC.txt | cut -f 2)			
		Het_Hom_ratio=$(sed -n 24p *${sample_ID}_Sample-Metrix-QC.txt | cut -f 2)	
		num_SNP_RAW=$(sed -n 25p *${sample_ID}_Sample-Metrix-QC.txt | cut -f 2)			
		num_SNP_PF=$(sed -n 26p *${sample_ID}_Sample-Metrix-QC.txt | cut -f 2)	
		SNP_Het_Hom_ratio=$(sed -n 27p *${sample_ID}_Sample-Metrix-QC.txt | cut -f 2)	
		num_InDel_RAW=$(sed -n 28p *${sample_ID}_Sample-Metrix-QC.txt | cut -f 2)			
		num_InDel_PF=$(sed -n 29p *${sample_ID}_Sample-Metrix-QC.txt | cut -f 2)	
		InDel_Het_Hom_ratio=$(sed -n 30p *${sample_ID}_Sample-Metrix-QC.txt | cut -f 2)

        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"  \
        "${sample_ID}" "${num_TOTAL_READS}" "${num_READ_PAIRS}" "${pct_READS_ABOVE_Q30}" "${MEAN_QUALITY}" "${MEDIAN_INSERT_SIZE}" "${pct_UNIQUE_READS}" "${pct_MAPPED_READS}" "${pct_MAPPED_OFF_TARGET}" "${pct_MAPPED_ON_TARGET}" \
		"${pct_TARGET_1X}" "${pct_TARGET_10X}" "${pct_TARGET_20X}" "${pct_TARGET_30X}" "${pct_TARGET_50X}" "${MEDIAN_TARGET_COV}" \
		"${Ti_Tv}" "${CONTAMINATION}" "${Het_Hom_ratio}" "${num_SNP_RAW}" "${num_SNP_PF}" "${SNP_Het_Hom_ratio}" "${num_InDel_RAW}" "${num_InDel_PF}" "${InDel_Het_Hom_ratio}"	>> ${WRKDIR}/Agregate_WGS_metrics.tsv

done < ${sample_list}
