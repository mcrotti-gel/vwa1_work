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

	// Processes
	query_vcf(input_vcfs_ch, regions_ch, consequence_ch)
	


}


// Run query

process query_vcf {

	tag "$file_name"
	publishDir "${params.outdir}/query_vcf_out", mode: 'copy'

	input:
	tuple val(file_name), file(vcf), file(vcf_idx)
	each file(regions)
	each file(consequence)

	output:
	file("${file_name}_VWA1_results.txt")

	shell:

	if( params.vcf_type == "standard" )
		"""
		tabix -h !{vcf} -R !{regions} | bcftools norm -m -any | \
		bcftools view -f PASS -i '(MIN(FMT/DP)>10 & MIN(FMT/GQ)>15) | (MIN(FMT/DPI)>10 & MIN(FMT/GQ)>15)' | \
		bcftools query -f '[%SAMPLE]\t%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t[%GT]\t[%GQ]\t([%DP]|[%DPI])\tINFO/CSQT\t%INFO/AF1000G\t%INFO\n' > !{file_name}_out.txt
		grep -f !{consequence} !{file_name}_out.txt || true > !{file_name}_VWA1_results.txt
		"""

	else if( params.vcf_type == "structural" )

		"""
		tabix -h ${vcf} -R ${regions} | \
		bcftools norm -m -any | \
		bcftools query -f '[%SAMPLE=%GT,%GQ]\t%SVTYPE\t%CHROM\t%POS\t%END\t%REF\t%ALT\t%FILTER\t%QUAL\n' >> VWA1_SV_results.txt
		"""

	else 
		error "Error: wrong vcf type specified!"

}
