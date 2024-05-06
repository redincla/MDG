version 1.0

## Copyrigght Medigenome, 2024
##
## This WDL defines tasks used for SMN CN calling based on targeted sequencing data.
##

#################################################################
# WORKFLOW DEFINITION
#################################################################
workflow SMNSampleCN {
input {
	Array[File] input_bams_samples
	File input_bams_PON
	File SNP_list
	Array[File] sample_list
	File SMN_CNplot_Rscript
	String output_file_basename_samples
	String output_file_basename_PON
}


call GetCoverageAtSNPsitesCohort as GetCoveragePON {
	input:
		input_bams = input_bams_PON,
		SNP_list = SNP_list,
		output_file_basename = output_file_basename_PON
}

 
call CoveragePON {
	input:
		coverage_PON=GetCoveragePON.output_cov,
		output_file_basename = output_file_basename_PON
}

scatter (i in range(length(sample_list))){

	call GetCoverageAtSNPsitesSingleSample as GetCoverageSamples {
		input:
			input_bam = input_bams_samples[i],
			SNP_list = SNP_list,
			sample_ID = sample_list[i]
}


	call GetSMNCN {
		input:
			sample_ID = sample_list[i],
			coverage_sample = GetCoverageSamples.output_cov,
			coverage_PON = CoveragePON.output_cov,
			sites_file = SNP_list
	}
}

call GatherSMNCNcohort {
	input:
		SMN_CN_samples = GetSMNCN.output_CN,
		output_file_basename = output_file_basename_samples
}

call PlotSMNCNcohort {
	input:
		cohort_SMN_CN = GatherSMNCNcohort.cohort_SMN_CN,
		SMN_CNplot_Rscript = SMN_CNplot_Rscript
}
	
output {
	File PON_output_cov = CoveragePON.output_cov
	Array[File] samples_output_cov = GetCoverageSamples.output_cov
	File cohort_SMN_CN = GatherSMNCNcohort.cohort_SMN_CN
	File cohort_SMN_CN_plot = PlotSMNCNcohort.cohort_SMN_CN_plot
	File cohort_SMN_CN_ratios = PlotSMNCNcohort.cohort_SMN_CN_ratios
	File cohort_SMN_CN_kmeans = PlotSMNCNcohort.cohort_SMN_CN_kmeans
}
}


#################################################################
# TASK DEFINITION
#################################################################

############
### Computes read depth at specific positions for a single sample
############

task GetCoverageAtSNPsitesSingleSample {
input {
	File input_bam
	File SNP_list
	String sample_ID
}

command {
	samtools depth -a \
	~{input_bam} \
	-b ~{SNP_list} -H -Q 0 \
	-o ~{sample_ID}.cov \
}

output {
	File output_cov = "~{sample_ID}.cov"
}
}

############
### Computes read depth at specific positions for multiple samples
############

task GetCoverageAtSNPsitesCohort {
input {
	File input_bams
	File SNP_list
	String output_file_basename
}

command {
	samtools depth -a \
	-b ~{SNP_list} -H \
	-f ~{input_bams} -Q 0 \
	-o ~{output_file_basename}.cov \
}

output {
	File output_cov = "~{output_file_basename}.cov"
}
}

############
### Compute mean coverage from PON coverage file - minimum size: 8
############ #coverage file output from GetCoverageAtSNPsites for PON

task CoveragePON {
input {
	File coverage_PON
	String output_file_basename
 }

command <<<
	paste <(cut -f 1,2 ~{coverage_PON}) <(awk '{sum=0; for (i=3;i<=NF;i++) sum+=$i; sum=sum/(NF-2); print sum}' ~{coverage_PON}) > ~{output_file_basename}.mean
>>>

output {
	File output_cov = "~{output_file_basename}.mean"
}
}


############
### Compute SMN CNs
############

task GetSMNCN {
input {
	String sample_ID
	File coverage_sample #coverage file output from GetCoverageAtSNPsites for specific sample
	File coverage_PON  #mean coverage file for PON from CoveragePON
	File sites_file  #bed file with SNP positions to run module on. string "norm" required for retrieving normalization sites
 }

command <<<
	norm_coordinates=$(grep -P "norm" ~{sites_file} | cut -f 3)
	norm_coverage_sample=$(grep -P "$norm_coordinates" ~{coverage_sample} | cut -f 3)
	norm_coverage_PON=$(grep -P "$norm_coordinates" ~{coverage_PON} | cut -f 3)
	paste <(awk -v norm1="$norm_coverage_sample" '{print $1,$2,($3/norm1)}' ~{coverage_sample}) \
	<(awk -v norm2="$norm_coverage_PON" '{print $1,$2,($3/norm2)}' ~{coverage_PON}) | \
	awk 'NR>1{print $1,$2,$3/$6}' >  ~{sample_ID}.CN.tmp
	echo -e "\t~{sample_ID}" >> ~{sample_ID}.SMN.CN
	while read chr1 start1 end1 site; do
		while read chr2	start2 ratio; do	
			if [ "${end1}" == "${start2}" ]; then
			rounded_ratio=$(printf "%.2f" "${ratio}")
			echo -e "${site}\t${rounded_ratio}" >> ~{sample_ID}.SMN.CN
			fi
		done < ~{sample_ID}.CN.tmp
	done < ~{sites_file}
>>>

output {
	File output_CN = "~{sample_ID}.SMN.CN"
}
}

############
### Gather SMN calls for a cohort
############

task GatherSMNCNcohort {
input {
	Array[File] SMN_CN_samples
	String output_file_basename
 }

command {
	paste ~{sep=" " SMN_CN_samples} > "~{output_file_basename}.SMN.CN"

}

output {
	File cohort_SMN_CN = "~{output_file_basename}.SMN.CN"
}
}

############
### Plot SMN calls for a cohort
############

task PlotSMNCNcohort {
input {
	File cohort_SMN_CN
	File SMN_CNplot_Rscript
 }

command <<<
	/usr/bin/Rscript ~{SMN_CNplot_Rscript} ~{cohort_SMN_CN}
>>>

output {
	File cohort_SMN_CN_plot = "SMN_PCA_plot.png"
	File cohort_SMN_CN_ratios = "ratios_table.txt"
	File cohort_SMN_CN_kmeans = "kmeans_result_table.txt"
}
}