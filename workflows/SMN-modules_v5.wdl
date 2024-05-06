version 1.0

## Copyrigght Medigenome, 2024
##
## This WDL defines tasks used for SMN CN calling based on targeted sequencing data.
## version 4: addig HBA + GBA callers in the workflow wrapper
## version 5: adapting script on newest design (more sites covered)
#################################################################
# WORKFLOW DEFINITION 
#################################################################
workflow CarrierScreeningSampleCN {
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
	
	call GetSiteCN {
		input:
			sample_ID = sample_list[i],
			normalized_cov_sample = NormalizesDepthSamples.output_cov_matrix_normalized,
			normalized_cov_PON = MedianDepthPON.PON_median_cov
	}
}

call GatherSitesCNcohort {
	input:
		sites_CN_samples = GetSiteCN.output_CN,
		output_file_basename = output_file_basename_samples
}

call PlotSMNCNcohort {
	input:
		cohort_sites_CN = GatherSitesCNcohort.cohort_sites_CN,
		output_file_basename = output_file_basename_samples,
		SNP_list = SNP_list,
		SMN_CNplot_Rscript = SMN_CNplot_Rscript
}
	
	
call PlotHBACN {
	input:
		cohort_sites_CN = GatherSitesCNcohort.cohort_sites_CN,
		output_file_basename = output_file_basename_samples,
		SNP_list = SNP_list
}

call PlotGBACN {
	input:
		cohort_sites_CN = GatherSitesCNcohort.cohort_sites_CN,
		output_file_basename = output_file_basename_samples,
		SNP_list = SNP_list
}

output {
	File cohort_sites_CN = GatherSitesCNcohort.cohort_sites_CN
	File cohort_SMN_plot_i7 = PlotSMNCNcohort.cohort_SMN_plot_i7
	File interactive_cohort_SMN_plot_i7 = PlotSMNCNcohort.interactive_cohort_SMN_plot_i7
	File cohort_SMN_plot_E8 = PlotSMNCNcohort.cohort_SMN_plot_E8
	File interactive_cohort_SMN_plot_E8 = PlotSMNCNcohort.interactive_cohort_SMN_plot_E8
	File cohort_SMN_plot_i8a = PlotSMNCNcohort.cohort_SMN_plot_i8a
	File interactive_cohort_SMN_plot_i8a = PlotSMNCNcohort.interactive_cohort_SMN_plot_i8a
	File cohort_SMN_plot_i8b = PlotSMNCNcohort.cohort_SMN_plot_i8b
	File interactive_cohort_SMN_plot_i8b = PlotSMNCNcohort.interactive_cohort_SMN_plot_i8b
	File cohort_HBA_plot = PlotHBACN.cohort_HBA_plot
	File cohort_HBA_CN = PlotHBACN.cohort_HBA_CN
	File cohort_GBA_plot = PlotGBACN.cohort_GBA_plot
	File cohort_GBA_CN = PlotGBACN.cohort_GBA_CN
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
	-o ~{sample_ID}.cov 

echo "locus" >> colnames
while read chr2 start2 depth; do
	while read chr1 start1 end1 locus; do
		if [ "$start2" = "$end1" ]; then
			echo "$locus" >> colnames
		fi
	done < "~{SNP_list}"
done < "~{sample_ID}.cov"
paste colnames "~{sample_ID}.cov" > tmp
mv tmp "~{sample_ID}.cov"

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
	-o ~{output_file_basename}.cov 
}

output {
	File output_cov = "~{output_file_basename}.cov"
}
}

############
### Compute median normalized coverage from PON coverage file - recommended min. size: 8  Remove 1st column with loci before launching R script else will fail
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
		norm_coverage=$(grep -P "$norm_coordinates" ~{coverage_sample} | cut -f 4)
		awk -v norm1="$norm_coverage" 'NR>1 {print ($4/norm1)}' ~{coverage_sample} > ~{sample_ID}.${i}.cov.normalized
		i=$((i+1))
	done < norm_sites
	
	paste ""~{sample_ID}.*.cov.normalized"" >> normalized.matrix
	
	
#	R --no-save << RSCRIPT  # normalization vector not used elsewhere in script
#	library(data.table)
#	data <- as.matrix(read.table("normalized.matrix", header=TRUE, row.names = NULL))
#	data <- apply(data, 2, as.numeric)
#	median <- apply(data, 1, median)
#	write.table(median, file = "normalized.cov.median", sep = "\t", quote = FALSE, row.names = FALSE)
#	RSCRIPT
	
	paste <(cut -f 4 ~{SNP_list}) normalized.matrix > ~{sample_ID}.normalized.matrix
#	mv normalized.cov.median ~{sample_ID}.cov.normalized # not same number of lines as matrix ???? to double check
 
>>>

output {
#	File output_cov_normalized = "~{sample_ID}.cov.normalized"
	File output_cov_matrix_normalized = "~{sample_ID}.normalized.matrix"
}
}

############
### Compute Copy Number status at each SNP site
############

task GetSiteCN {
input {
	String sample_ID
	File normalized_cov_sample #coverage matrix output from GetCoverageAtSNPsites for specific sample
	File normalized_cov_PON  #median normalized coverage matrix for PONs
 }

command <<<
	awk 'NR==FNR {m1[NR]=$0; next} {split(m1[FNR], row1); split($0, row2); printf "%s ", $1; for (i=2; i<=NF; i++) printf("%.6f ", row1[i]/row2[i]); print ""}' ~{normalized_cov_sample} ~{normalized_cov_PON} >> matrix.sites.CN
	
	R --no-save << RSCRIPT
	library(data.table)
	data <- as.matrix(read.table("matrix.sites.CN", header=FALSE))
	medians <- apply(data[, -1], 1, median)  # Exclude the first column when computing median
	write.table(medians, file = "median.sites.CN", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
	RSCRIPT
	cat <(echo "~{sample_ID}") median.sites.CN > tmp
	cat <(echo "locus") <(awk '{print $1}'  matrix.sites.CN) > colnames
	paste colnames tmp > "~{sample_ID}.sites.CN"
>>>

output {
	File output_CN = "~{sample_ID}.sites.CN"
}
}

############
### Compute and plot HBA CNs
############

task PlotHBACN {
input {
	File cohort_sites_CN
	File SNP_list
	String output_file_basename
 }

command <<<
	header=$(head -n 1 ~{cohort_sites_CN})
	first_line=$(grep -n HBA ~{SNP_list} | head -n 1 | cut -f 1 -d: | awk '{print $1 + 1}' | bc)
	last_line=$(grep -n HBA ~{SNP_list} | tail -n 1 | cut -f 1 -d: | awk '{print $1 + 1}' | bc)
	sed -n "${first_line},${last_line}p" ~{cohort_sites_CN} > tmp
	cat <(echo "$header") tmp > HBA_CN
	rm tmp
	
	R --no-save << RSCRIPT
	library(data.table)
	data <- as.matrix(read.table("HBA_CN", sep = "\t", header = TRUE))
	HBA_copies<-2*(as.data.frame[data, -1])
	rownames(HBA_copies)<-c("upstream_a","upstream_b","upstream_c","upstream_d","upstream_e","HBA2_E1","HBA2_i2a","HBA2_i2b","HBA2_i2c","HBA2_i2d","HBA2_E3","HBA2_3UTR","intergenic_a","intergenic_b","intergenic_c","intergenic_d","intergenic_e","intergenic_f","HBA1_E2","HBA1_i2a","HBA1_i2b","HBA1_i2c","HBA1_3UTR","downstream_a","downstream_b","downstream_c")
	colors <- c("blue", "red", "green", "purple", "orange", "cyan", "magenta", "black","darkgreen")
	line_types <- c(1, 2, 3, 4)
	p<-matplot(HBA_copies, type="l", xlab="", ylab="copy number", main="HBA1/HBA2 copies", col=colors, lty=line_types, lwd=2, ylim=c(0, 4),xaxt="n")
	legend(legend=colnames(HBA_copies), col=colors, lty=line_types, lwd=2,x="topright",cex=0.6, ncol=3)
	axis(1, at = 1:nrow(HBA_copies), labels = rownames(HBA_copies),cex.axis=0.6,las=2)
	p
	dev.off()
	write.table(HBA_copies, file = "HBA_copies.txt", sep = "\t", quote = FALSE)
	RSCRIPT
	
	mv Rplots.pdf ~{output_file_basename}.HBA.plot.pdf
	mv HBA_copies.txt ~{output_file_basename}.HBA.CN.tsv
>>>

output {
	File cohort_HBA_plot = "~{output_file_basename}.HBA.plot.pdf"
	File cohort_HBA_CN = "~{output_file_basename}.HBA.CN.tsv"
}
}

############
### Compute and plot GBA CNs
############

task PlotGBACN {
input {
	File cohort_sites_CN
	File SNP_list
	String output_file_basename
 }

command <<<
	header=$(head -n 1 ~{cohort_sites_CN})
	first_line=$(grep -n GBA ~{SNP_list} | head -n 1 | cut -f 1 -d: | awk '{print $1 + 1}' | bc)
	last_line=$(grep -n GBA ~{SNP_list} | tail -n 1 | cut -f 1 -d: | awk '{print $1 + 1}' | bc)
	sed -n "${first_line},${last_line}p" ~{cohort_sites_CN} > tmp
	cat <(echo "$header") tmp > GBA_CN
	rm tmp
	
	R --no-save << RSCRIPT
	library(data.table)
	data <- as.matrix(read.table("GBA_CN", sep = "\t", header = TRUE))
	GBA_copies<-2*(as.data.frame[data, -1])
	colors <- c("blue", "red", "green", "purple", "orange", "cyan", "magenta", "black","darkgreen")
	line_types <- c(1, 2, 3, 4)
	p<-matplot(GBA_copies, type="l", xlab="", ylab="copy number", main="GBA/GBAP1 copies", col=colors, lty=line_types, lwd=2, ylim=c(0, 5),xaxt="n")
	legend(legend=colnames(GBA_copies),col=colors, lty=line_types, lwd=2,x="topright",cex=0.6, ncol=3)
	axis(1, at = 1:nrow(GBA_copies), labels = c(rep("GBAP1",33),rep("GBAdel",3),rep("intergenic",46),rep("GBA1",37)),cex.axis=0.6,las=2)
	p
	dev.off()
	write.table(GBA_copies, file = "GBA_copies.txt", sep = "\t", quote = FALSE)
	RSCRIPT
	
	mv Rplots.pdf ~{output_file_basename}.GBA.plot.pdf
	mv GBA_copies.txt ~{output_file_basename}.GBA.CN.tsv
>>>

output {
	File cohort_GBA_plot = "~{output_file_basename}.GBA.plot.pdf"
	File cohort_GBA_CN = "~{output_file_basename}.GBA.CN.tsv"
}
}

############
### Gather SNP sites CN calls for a cohort
############

task GatherSitesCNcohort {
input {
	Array[File] sites_CN_samples
	String output_file_basename
 }

command {
#	paste <(cut -f 2- -d ' ' ~{sep=" " sites_CN_samples}) > "~{output_file_basename}.sites.CN"
	paste $(for file in ~{sep=" " sites_CN_samples}; do cut -f 2- -d ' ' $file; done) > "~{output_file_basename}.sites.CN"
}

output {
	File cohort_sites_CN = "~{output_file_basename}.sites.CN"
}
}

############
### Plot SMN calls for a cohort
############

task PlotSMNCNcohort {
input {
	File cohort_sites_CN
	File SMN_CNplot_Rscript
	File SNP_list
	String output_file_basename
 }

command <<<
	header=$(head -n 1 ~{cohort_sites_CN})
	first_line=$(grep -n SMN ~{SNP_list} | head -n 1 | cut -f 1 -d: | awk '{print $1 + 1}' | bc)
	last_line=$(grep -n SMN ~{SNP_list} | tail -n 1 | cut -f 1 -d: | awk '{print $1 + 1}' | bc)
	sed -n "${first_line},${last_line}p" ~{cohort_sites_CN} > tmp
	cat <(echo "$header") tmp > SMN_CN
	rm tmp
	/usr/bin/Rscript ~{SMN_CNplot_Rscript} SMN_CN
	
	mv SMN_plot_i7.png ~{output_file_basename}.SMN_plot_i7.png
	mv SMN_plot_i7.html ~{output_file_basename}.SMN_plot_i7.html
	mv SMN_plot_E8.png ~{output_file_basename}.SMN_plot_E8.png
	mv SMN_plot_E8.html ~{output_file_basename}.SMN_plot_E8.html
	mv SMN_plot_i8a.png ~{output_file_basename}.SMN_plot_i8a.png
	mv SMN_plot_i8a.html ~{output_file_basename}.SMN_plot_i8a.html
	mv SMN_plot_i8b.png ~{output_file_basename}.SMN_plot_i8b.png
	mv SMN_plot_i8b.html ~{output_file_basename}.SMN_plot_i8b.html
	mv ratios_table.txt ~{output_file_basename}.ratios_table.txt
>>>

output {
	File cohort_SMN_plot_i7 = "~{output_file_basename}.SMN_plot_i7.png"
	File interactive_cohort_SMN_plot_i7 = "~{output_file_basename}.SMN_plot_i7.html"
	File cohort_SMN_plot_E8 = "~{output_file_basename}.SMN_plot_E8.png"
	File interactive_cohort_SMN_plot_E8 = "~{output_file_basename}.SMN_plot_E8.html"
	File cohort_SMN_plot_i8a = "~{output_file_basename}.SMN_plot_i8a.png"
	File interactive_cohort_SMN_plot_i8a = "~{output_file_basename}.SMN_plot_i8a.html"
	File cohort_SMN_plot_i8b = "~{output_file_basename}.SMN_plot_i8b.png"
	File interactive_cohort_SMN_plot_i8b = "~{output_file_basename}.SMN_plot_i8b.html"
	File cohort_SMN_CN_ratios = "~{output_file_basename}.ratios_table.txt"
}
}