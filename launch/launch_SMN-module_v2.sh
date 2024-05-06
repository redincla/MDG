#!/bin/bash

#############################
#  wrapper for SMN module   #
#############################

# checks for appropriate input
if [ $# -eq 4 ]; then
 sample_list=$1 #Col1: list of samples to process, Col2: full path to full-map
 cohort_name=$2
 SNP_list=$3
 WRKDIR=$4 #full directory of processed results to be dropped

else
	echo -e "\n\nLaunching SMN CN modules on cohort\n\nAuthor: Claire Redin (claire.redin@medigenome.ch)\n\n"
	echo "Usage:"
	echo "launch_SMN-module.sh [sample_list] [cohort_name] [SNP_list] [output directory]"
	echo "sample_list: Column 1= sample IDs, Column 2 = full path to bam containing folder, Column 3 = [TEST/PON]"
	echo "cohort name: string for cohort ID"
	echo "SNP_list: bedfile of SNPs to use for module in hg38"
	echo "output directory: full path"
	exit 1
fi


### Set local parameters
export WRKDIR

if ! [ -e $WRKDIR ]; then
	mkdir $WRKDIR
fi

cd ${WRKDIR}


### Prepare folder structure
if ! [ -e ${WRKDIR}/json ]; then
	    mkdir ${WRKDIR}/json
fi

if ! [ -e ${WRKDIR}/SMNmodule ]; then
	    mkdir ${WRKDIR}/SMNmodule
fi

# Write input json
touch ${WRKDIR}/json/${cohort_name}.SMNCNcalling.input.json

cut -f 1 $1 > ${WRKDIR}/SMNmodule/${cohort_name}_sample_IDs
cut -f 2 $1 > ${WRKDIR}/SMNmodule/${cohort_name}_sample_bams

grep "PON" $1 | cut -f 1 > ${WRKDIR}/SMNmodule/${cohort_name}_PON_IDs
grep "PON" $1 | cut -f 2 > ${WRKDIR}/SMNmodule/${cohort_name}_PON_bams

sample_bams=$(<${WRKDIR}/SMNmodule/${cohort_name}_sample_bams tr '\n' ',') 
sample_bams=$(echo "${sample_bams}" | sed 's/,/","/g' | sed 's/^/"/' | sed 's/,"$//')

sample_IDs=$(<${WRKDIR}/SMNmodule/${cohort_name}_sample_IDs tr '\n' ',') 
sample_IDs=$(echo "${sample_IDs}" | sed 's/,/","/g' | sed 's/^/"/' | sed 's/,"$//')

PON_bams=$(<${WRKDIR}/SMNmodule/${cohort_name}_PON_bams tr '\n' ',') 
PON_bams=$(echo "${PON_bams}" | sed 's/,/","/g' | sed 's/^/"/' | sed 's/,"$//')

PON_IDs=$(<${WRKDIR}/SMNmodule/${cohort_name}_PON_IDs tr '\n' ',') 
PON_IDs=$(echo "${PON_IDs}" | sed 's/,/","/g' | sed 's/^/"/' | sed 's/,"$//')

    cat <<EOF > ${WRKDIR}/json/${cohort_name}.SMNCNcalling.input.json
{	
  "SMNSampleCN.input_bams_samples": [ cat $sample_bams ],
  "SMNSampleCN.input_bams_PON": [ cat $PON_bams ],
  "SMNSampleCN.sample_list": [ cat $sample_IDs ],
  "SMNSampleCN.PON_list": [ cat $PON_IDs ],
  "SMNSampleCN.SNP_list": "$3",
  "SMNSampleCN.output_file_basename_PON": "PON",
  "SMNSampleCN.SMN_CNplot_Rscript": "/mnt/raw-data/scripts/wdls/tasks/SMN-CN-Kmeans_v2.R",
  "SMNSampleCN.output_file_basename_samples": "${cohort_name}"
}
EOF

sed -i 's/cat //g' ${WRKDIR}/json/${cohort_name}.SMNCNcalling.input.json

# Write cromwell options with final output directory
  touch ${WRKDIR}/json/${cohort_name}.SMNCNcalling.options.json
  cat <<EOF > ${WRKDIR}/json/${cohort_name}.SMNCNcalling.options.json
{
    "final_workflow_outputs_dir": "${WRKDIR}/SMNmodule",
    "use_relative_output_paths": true
}
EOF
	
	
### Launch SMN module script
java -Dconfig.file=/home/credin/.cromwell_config.file -jar /mnt/raw-data/tools/cromwell/cromwell-86.jar run /mnt/raw-data/scripts/wdls/workflows/SMN-modules_v2.wdl -i ${WRKDIR}/json/${cohort_name}.SMNCNcalling.input.json -o ${WRKDIR}/json/${cohort_name}.SMNCNcalling.options.json

