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
	File MedianSMNRatios_Rscript
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

	scatter (j in range(length(PON_list))){
	
		call GetSMNratios {
			input:
				sample_ID = sample_list[i],
				normalized_cov_sample = NormalizesDepthSamples.output_cov_matrix_normalized,
				normalized_cov_PON = NormalizesDepthPONs.output_cov_matrix_normalized[j]
		}
	}
	
	call GetMedianSMNratios {
		input:
			sample_ID = sample_list[i],
			SMN_ratios = GetSMNratios.output_ratios,
			MedianSMNRatios_Rscript = MedianSMNRatios_Rscript,
	}
}

call GatherSMNCNcohort {
	input:
		SMN_CN_samples = GetMedianSMNratios.median_ratios,
		output_file_basename = output_file_basename_samples
}

call PlotSMNCNcohort {
	input:
		cohort_SMN_CN = GatherSMNCNcohort.cohort_SMN_CN,
		output_file_basename = output_file_basename_samples,
		SNP_list = SNP_list,
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
### Compute SMN CNs ratios using each SNP normalization site: one test sample compared with one PON normal
############

task GetSMNratios {
input {
	String sample_ID
	File normalized_cov_sample #normalized coverage matrix output from GetCoverageAtSNPsites for specific sample
	File normalized_cov_PON  #normalized coverage matrix from GetCoverageAtSNPsites for a single PON sample
 }

command <<<	
	R --no-save << RSCRIPT
	library(data.table)
	sample <- as.matrix(read.table("~{normalized_cov_sample}", header=FALSE))
	PON <- as.matrix(read.table("~{normalized_cov_PON}", header=FALSE))
	PON[PON == 0] <- NA
	ratios <- sample / PON
	write.table(ratios, file = "ratios.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
	RSCRIPT
	mv ratios.txt "~{sample_ID}.SMN.ratios"
>>>

output {
	File output_ratios = "~{sample_ID}.SMN.ratios"
}
}


############
### Compute median SMNratios compared to the full PON set
############ 

task GetMedianSMNratios {
input {
	Array[File] SMN_ratios
	File MedianSMNRatios_Rscript
	String sample_ID
 }

command <<<
	/usr/bin/Rscript ~{MedianSMNRatios_Rscript} ~{sep=" " SMN_ratios}
	cat <(echo "~{sample_ID}") median_ratios.tsv > ~{sample_ID}.median_ratios.tsv
>>>

output {
	File median_ratios = "~{sample_ID}.median_ratios.tsv"
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
	File SNP_list
	String output_file_basename
 }

command <<<
	header=$(head -n 1 ~{cohort_SMN_CN})
	first_line=$(grep -n SMN ~{SNP_list} | head -n 1 | cut -f 1 -d: | awk '{print $1 + 1}' | bc)
	last_line=$(grep -n SMN ~{SNP_list} | tail -n 1 | cut -f 1 -d: | awk '{print $1 + 1}' | bc)
	sed -n "${first_line},${last_line}p" ~{cohort_SMN_CN} > tmp
	cat <(echo "$header") tmp > SMN_CN
	rm tmp
	/usr/bin/Rscript ~{SMN_CNplot_Rscript} SMN_CN
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