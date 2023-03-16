#!/usr/bin/env nextflow
nextflow.enable.dsl=2


// Run query

process QUERY_VCF {

	tag "$file_name"

	input:
	tuple val(file_name), path(vcf), path(vcf_idx)
	each path(regions)
	each path(consequence)

	output:
	path("${file_name}_VWA1_results.txt")

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
		bcftools query -f '[%SAMPLE=%GT,%GQ]\t%SVTYPE\t%CHROM\t%POS\t%END\t%REF\t%ALT\t%FILTER\t%QUAL\n' >> !{file_name}_VWA1_results.txt
		"""

	else 
		error "Error: wrong vcf type specified!"

}

// Concatenate outfiles

process CONCATENATE {

	publishDir "${params.outdir}/final_out", mode: 'copy'

	input:
	path(sample_output)

	output:
	path("VWA1_results.txt")

	shell:

	if ( params.vcf_type == "standard" )
	'''
	cat !{sample_output} >> VWA1_results.txt
	sed -i '1i SAMPLE\tCHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tGT\tGQ\tDP_DPI\tCSQT\tAF1000G' VWA1_results.txt
	'''
	else if( params.vcf_type == "structural" )

	'''
	cat !{sample_output} >> VWA1_results.txt
	sed -i '1i [SAMPLE=GT,GQ]\tSVTYPE\tCHROM\tPOS\tEND\tREF\tALT\tFILTER\tQUAL' VWA1_results.txt
	'''


}

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
	sample_output_ch = QUERY_VCF( input_vcfs_ch, regions_ch, consequence_ch )
	CONCATENATE( sample_output_ch.collect() )

}