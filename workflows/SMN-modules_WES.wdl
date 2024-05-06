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
	Array[File] input_bams_PON
	File SNP_list
	Array[File] sample_list
	Array[File] PON_list
	File SMN_CNplot_Rscript
	File MedianDepthPON_Rscript
	String output_file_basename_samples
	String output_file_basename_PON
}

scatter (i in range(length(input_bams_PON))){

	call GetCoverageAtSNPsitesSingleSample as GetCoveragePONs {
		input:
			input_bam = input_bams_PON[i],
			SNP_list = SNP_list,
			sample_ID = PON_list[i]
	}

	call NormalizesDepth as NormalizesDepthPONs {
		input:
			sample_ID = PON_list[i],
			SNP_list = SNP_list,
			coverage_sample = GetCoveragePONs.output_cov
	}
}

call MedianDepthPON {
		input:
			normalized_cov = NormalizesDepthPONs.output_cov_matrix_normalized,
			MedianDepthPON_Rscript = MedianDepthPON_Rscript,
			output_file_basename = output_file_basename_samples
}

scatter (i in range(length(sample_list))){

	call GetCoverageAtSNPsitesSingleSample as GetCoverageSamples {
		input:
			input_bam = input_bams_samples[i],
			SNP_list = SNP_list,
			sample_ID = sample_list[i]
}

	call NormalizesDepth as NormalizesDepthSamples {
		input:
			sample_ID = sample_list[i],
			SNP_list = SNP_list,
			coverage_sample = GetCoverageSamples.output_cov
}
	
	call GetSMNCN {
		input:
			sample_ID = sample_list[i],
			normalized_cov_sample = NormalizesDepthSamples.output_cov_matrix_normalized,
			normalized_cov_PON = MedianDepthPON.PON_median_cov
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
		output_file_basename = output_file_basename_samples,
		SMN_CNplot_Rscript = SMN_CNplot_Rscript
}
	
output {
	File cohort_SMN_CN = GatherSMNCNcohort.cohort_SMN_CN
	File cohort_SMN_CN_ratios = PlotSMNCNcohort.cohort_SMN_CN_ratios
	File cohort_SMN_plot_i7 = PlotSMNCNcohort.cohort_SMN_plot_i7
	File interactive_cohort_SMN_plot_i7 = PlotSMNCNcohort.interactive_cohort_SMN_plot_i7
	File cohort_SMN_plot_E8 = PlotSMNCNcohort.cohort_SMN_plot_E8
	File interactive_cohort_SMN_plot_E8 = PlotSMNCNcohort.interactive_cohort_SMN_plot_E8
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
### Compute median normalized coverage from PON coverage file - minimum size: 8
############ 

task MedianDepthPON {
input {
	Array[File] normalized_cov
	File MedianDepthPON_Rscript
	String output_file_basename
 }

command <<<
	/usr/bin/Rscript ~{MedianDepthPON_Rscript} ~{sep=" " normalized_cov}
	mv PON_median_matrix.tsv ~{output_file_basename}.PON_median_matrix.tsv
>>>

output {
	File PON_median_cov = "~{output_file_basename}.PON_median_matrix.tsv"
}
}


############
### Compute normalized depth intra sample
############

task NormalizesDepth {
input {
	String sample_ID
	File coverage_sample #coverage file output from GetCoverageAtSNPsites for specific sample
	File SNP_list  #bed file with SNP positions to run module on. string "norm" required for retrieving normalization sites
 }

command <<<
	grep -P "norm" ~{SNP_list} > norm_sites 
	i=1
	while read line; do
		norm_coordinates=$(echo $line | awk '{print $3}')
		norm_coverage=$(grep -P "$norm_coordinates" ~{coverage_sample} | cut -f 3)
		awk -v norm1="$norm_coverage" 'NR>1 {print ($3/norm1)}' ~{coverage_sample} > ~{sample_ID}.${i}.cov.normalized
		i=$((i+1))
	done < norm_sites
	
	paste ""~{sample_ID}.*.cov.normalized"" >> normalized.matrix
	
	R --no-save << RSCRIPT
	library(data.table)
	data <- as.matrix(read.table("normalized.matrix", header=TRUE, row.names = NULL))
	data <- apply(data, 2, as.numeric)
	median <- apply(data, 1, median)
	write.table(median, file = "normalized.cov.median", sep = "\t", quote = FALSE, row.names = FALSE)
	RSCRIPT
	
	mv normalized.cov.median ~{sample_ID}.cov.normalized
	mv normalized.matrix ~{sample_ID}.normalized.matrix
>>>

output {
	File output_cov_normalized = "~{sample_ID}.cov.normalized"
	File output_cov_matrix_normalized = "~{sample_ID}.normalized.matrix"
}
}

############
### Compute SMN CNs
############

task GetSMNCN {
input {
	String sample_ID
	File normalized_cov_sample #coverage matrix output from GetCoverageAtSNPsites for specific sample
	File normalized_cov_PON  #median normalized coverage matrix for PONs
 }

command <<<
	awk 'NR==FNR {m1[NR]=$0; next} {split(m1[FNR], row1); split($0, row2); for (i=1; i<=NF; i++) printf("%.6f ", row1[i]/row2[i]); print ""}' ~{normalized_cov_sample} ~{normalized_cov_PON} >> matrix.SMN.CN
	
	R --no-save << RSCRIPT
	library(data.table)
	data <- as.matrix(read.table("matrix.SMN.CN", header=FALSE))
	medians <- apply(data, 1, median)
	write.table(medians, file = "median.SMN.CN", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
	RSCRIPT
	cat <(echo "~{sample_ID}") median.SMN.CN > "~{sample_ID}.SMN.CN"
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
	String output_file_basename
 }

command <<<
	/usr/bin/Rscript ~{SMN_CNplot_Rscript} ~{cohort_SMN_CN}
	mv SMN_plot_i7.png ~{output_file_basename}.SMN_plot_i7.png
	mv SMN_plot_i7.html ~{output_file_basename}.SMN_plot_i7.html
	mv SMN_plot_E8.png ~{output_file_basename}.SMN_plot_E8.png
	mv SMN_plot_E8.html ~{output_file_basename}.SMN_plot_E8.html
	mv ratios_table.txt ~{output_file_basename}.ratios_table.txt
>>>

output {
	File cohort_SMN_plot_i7 = "~{output_file_basename}.SMN_plot_i7.png"
	File interactive_cohort_SMN_plot_i7 = "~{output_file_basename}.SMN_plot_i7.html"
	File cohort_SMN_plot_E8 = "~{output_file_basename}.SMN_plot_E8.png"
	File interactive_cohort_SMN_plot_E8 = "~{output_file_basename}.SMN_plot_E8.html"
	File cohort_SMN_CN_ratios = "~{output_file_basename}.ratios_table.txt"
}
}