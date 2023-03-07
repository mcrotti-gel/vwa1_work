#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/* --------------------------
  Set up variables
 ----------------------------*/

workflow {
	// single participant vcfs
	input_vcfs_ch = Channel
		.fromPath(params.input_vcfs)
		.ifEmpty { exit 1, "Cannot find input file : ${params.input_vcfs}" }
		.splitCsv(skip:1)
		.map {file_name, vcf, vcf_idx -> [ file_name, file(vcf), file(vcf_idx) ] }


	// regions file
	regions_ch = Channel
		.fromPath(params.regions)
		.ifEmpty { exit 1, "Cannot find input file : ${params.regions}" }

	// consequence file
	consequence_ch = Channel
		.fromPath(params.consequence)

	
	if (params.vcf_type == "standard") {
		query_standard_vcf(input_vcfs_ch, regions_ch, consequence_ch)
	} else if (params.vcf_type == "structural") {
		query_structural_vcf(input_vcfs_ch, regions_ch)
	} else {
		println "Error: wrong vcf type specified!"
	}


}


// Run query

process query_standard_vcf {

	tag "$file_name"
	publishDir "${params.outdir}/standard_vcf_out", mode: 'copy'

	input:
	tuple val(file_name), file(vcf), file(vcf_idx)
	each file(regions)
	each file(consequence)

	output:
	file("${file_name}_VWA1_results.txt")

	script:

	"""
	tabix -h ${vcf} -R ${regions} | bcftools norm -m -any > ${file_name}_VWA1_results.txt
	
	"""

}

process query_structural_vcf {

	tag "$file_name"
	publishDir "${params.outdir}/structural_vcf_out", mode: 'copy'

	input:
	tuple val(file_name), file(vcf), file(vcf_idx)
	each file(regions)

	output:
	file("VWA1_SV_results.txt")

	script:

	"""
	tabix -h ${vcf} -R ${regions} | \
	bcftools norm -m -any | \
	bcftools query -f '[%SAMPLE=%GT,%GQ]\t%SVTYPE\t%CHROM\t%POS\t%END\t%REF\t%ALT\t%FILTER\t%QUAL\n' >> VWA1_SV_results.txt
	
	sed -i '1s/^[%SAMPLE=%GT,%GQ]\t%SVTYPE\t%CHROM\t%POS\t%END\t%REF\t%ALT\t%FILTER\t%QUAL\n/' VWA1_SV_results.txt

	"""

}