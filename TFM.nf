#!/usr/bin/env nextflow

/*
#################################################
## Riboseq analysis pipeline for lncRNA detection
#################################################
*/

// Channel for reading raw fastq files

Channel
	.fromFilePairs(params.raw_fastq, size: 1)
        .ifEmpty { exit 1, "Cannot find any fastq files matching: ${params.raw_fastq}" }
	.into { fastqc_ch; trimgalore_ch }

// Starting the pipeline

process Fastqc{

        tag "$sampleID-quality-control"
        publishDir "$params.outdir/quality-control", mode: params.publish_dir_mode, pattern: "*.html"
        container 'ubuntu-tfm:v1'

        input:
         set sampleID, file(raw_fastq) from fastqc_ch

	output:
	file '*_fastqc.{zip,html}' into fastqc_results

	script:
	"""
	 /opt/bin/FastQC/fastqc ${params.fastqc_cmd_args} $params.basedocker/data/$raw_fastq -o .    
        """
}

process TrimGalore{

        tag "$sampleID-trimming"
        publishDir "$params.outdir/trim", mode: params.publish_dir_mode
        container 'ubuntu-tfm:v1'

        input:
         set sampleID, file(raw_fastq) from trimgalore_ch

        output:
         tuple val(sampleID), file ('*_trimmed.fq.gz') into fastq_trimmed, fastq_trimmed2, fastq_trimmed3, fastq_trimmed4
	 file '*.fastq.gz_trimming_report.txt' into trim_report
	
	script:
	"""
        /opt/bin/trim_galore ${params.trimgalore_cmd_args} $params.basedocker/data/$raw_fastq -o .
        """
}

// File objects requiring PRICE

genome_file = file(params.genome)
annotation_file = file(params.annot)

// Check if files are in the indicated directory

if ( !genome_file.exists() ) exit 1, "genome file not found: ${params.genome}"
if ( !annotation_file.exists() ) exit 1, "annotation file not found: ${params.annot}"

process Index_Price{

	tag "index-price"
        publishDir "$params.outdir/price", mode: params.publish_dir_mode
        container 'ubuntu-tfm:v1'

        input:
         file(genome) from genome_file
         file(annotation) from annotation_file

	output:
	file '*.oml' into oml_price

	script:
	"""
	/opt/bin/gedi ${params.indexprice_cmd_args} -s $params.basedocker/reference/$genome -a $params.basedocker/reference/$annotation -f $params.basedocker/reference -o ./hs_38.oml
	"""
}

process Bowtie_Price{

        tag "$sampleID-bowtie-price"
        publishDir "$params.outdir/price", mode: params.publish_dir_mode
        container 'ubuntu-tfm:v1'

        input:
        tuple val(sampleID), file (fastq) from fastq_trimmed
	
        output:
	 file '*.txt' into bowtie_report 
	 tuple val(sampleID), file ('*.sorted.bam') into sorted_bam
         file '*.sorted.bam.bai' into index_bam

        script:
        """
        bowtie -p $params.threads $params.basedocker/reference/gencode.v37.annotation.genomic -q $params.basedocker/results/trim/$fastq -S >${sampleID}.sam 2> ${sampleID}_bowtie_mapping_report.txt
        /opt/bin/samtools sort -@ $params.threads ./${sampleID}.sam > ${sampleID}.sorted.bam
        /opt/bin/samtools index -@ $params.threads ./${sampleID}.sorted.bam >${sampleID}.sorted.bam.bai
	 """
}

process Price{

        tag "price"
        publishDir "$params.outdir/price", mode: params.publish_dir_mode 
        container 'ubuntu-tfm:v1'

        input:
         tuple val(sampleID), file (bam) from sorted_bam
	 file (hs38) from oml_price
	 
	output: 
	 file "*-hsv.orfs.tsv" into price_results

        script:
        """
        /opt/bin/gedi ${params.price_cmd_args} -reads $bam -genomic $hs38 -prefix ${sampleID}-hsv  
        """
}


// File objects requiring RiboCode

rRNA_file = file(params.rRNA)
annot_file = file(params.annot)
genom_file = file(params.genome)

// Check if files are in the indicated directory

if ( !rRNA_file.exists() ) exit 1, "rRNA file not found: ${params.rRNA}"
if ( !annot_file.exists() ) exit 1, "annot file not found: ${params.annot}"
if ( !genom_file.exists() ) exit 1, "genom file not found: ${params.genome}"

process Index_RiboCode{

        tag "index-ribocode"
        publishDir "$params.outdir/RiboCode", mode: params.publish_dir_mode 
        container 'ubuntu-tfm:v1'

        input:
         file (rRNA) from rRNA_file

        script:
        """
        bowtie-build $params.basedocker/reference/$rRNA $params.basedocker/reference/rRNA
        """
}


process Bowtie_RiboCode{

        tag "bowtie-ribocode"
        publishDir "$params.outdir/RiboCode", mode: params.publish_dir_mode
        container 'ubuntu-tfm:v1'

        input:
         tuple val(sampleID), file (fastq) from fastq_trimmed2
	
        output:
         tuple val(sampleID), file ("*_un_aligned.fastq") into unaligned_ribocode
	 file "*_rRNA.align" into aligned_ribocode

        script:
        """
         bowtie ${params.bowtie_ribocode_cmd_args} $params.basedocker/reference/rRNA --un ./${sampleID}_un_aligned.fastq -q $params.basedocker/results/trim/$fastq -S ./${sampleID}_rRNA.align
        """
}


process Genome_STAR{

	tag "genome-STAR"
	publishDir "$params.outdir/RiboCode", mode: params.publish_dir_mode
        container 'ubuntu-tfm:v1'

	input:
         file (geno) from genom_file
	 file (annot) from annot_file

        script:
	if( "${params.starIndex}" == 'true' )
	        """
	       /opt/bin/STAR ${params.genomeSTAR_cmd_args} --genomeDir $params.basedocker/reference --genomeFastaFiles $params.basedocker/reference/$geno --sjdbGTFfile $params.basedocker/reference/$annot 
	        """
	else
		"""
		echo "no star index"
		"""
}


process STAR{

        tag "STAR"
        publishDir "$params.outdir/RiboCode", mode: params.publish_dir_mode
        container 'ubuntu-tfm:v1'

        input:
         tuple val(sampleID), file (unaligned) from unaligned_ribocode

        output:
         tuple val(sampleID), file ("*_star.bamAligned.toTranscriptome.out.bam") into STAR_results, STAR_results2, STAR_results3 
        
	script:
        """
        /opt/bin/STAR ${params.STAR_cmd_args} --genomeDir $params.basedocker/reference --readFilesIn ./$unaligned --outFileNamePrefix ./${sampleID}_star.bam
        """
}

// File objects requiring RiboCode

gtf_file = file(params.annot)
reference_file = file(params.genome)

// Check if files are in the indicated directory

if ( !gtf_file.exists() ) exit 1, "gtf file not found: ${params.annot}"
if ( !reference_file.exists() ) exit 1, "reference file not found: ${params.genome}"

process Transcripts{

        tag "transcript"
        publishDir "$params.outdir/RiboCode/", mode: params.publish_dir_mode
        container 'ubuntu-tfm:v1'

        input:
         file (gtf) from gtf_file
	 file (reference) from reference_file
	 tuple val(sampleID), file (transcript) from STAR_results2
	
	output:
	 file "*_ORFs_result_collapsed_out_out.txt" into ribocode
		
	script:
        """
        prepare_transcripts -f $params.basedocker/reference/$reference -g $params.basedocker/reference/$gtf -o $params.basedocker/results/RiboCode/Annot
	metaplots -a $params.basedocker/results/RiboCode/Annot -r $params.basedocker/results/RiboCode/$transcript -o $params.basedocker/results/RiboCode/${sampleID}
        RiboCode -a $params.basedocker/results/RiboCode/Annot -c $params.basedocker/results/RiboCode/${sampleID}_pre_config.txt -l no -g -o $params.basedocker/results/RiboCode/${sampleID}_ORFs_result
	cat $params.basedocker/results/RiboCode/${sampleID}_ORFs_result_collapsed.txt >./${sampleID}_ORFs_result_collapsed_out_out.txt
	"""
}

process getTranscriptsRiboCode{
	tag "getTranscripts"
	publishDir "$params.outdir/RiboCode", mode: params.publish_dir_mode

	input:
	 path x from ribocode

	output:
	 file '*ribocode_out.txt' into ribocodeOut

	script:
  	 """
	 cat ${x} | awk 'BEGIN{FS="\\t"}{if(\$2=="novel" && \$28<=0.000005){print \$2,\$3,\$28}}' >${x}_ribocode_out.txt
	 """	

 }

process getTranscriptsPrice{
        tag "getTranscriptsPrice"
        publishDir "$params.outdir/price", mode: params.publish_dir_mode

        input:
         path x from price_results

        output:
         file '*price_out.txt' into priceOut

        script:
         """
	 cat ${x} | awk 'BEGIN{FS="\\t"}{if(\$6=="ncRNA" && \$9<=0.05){print \$6,\$2,\$9}}' >${x}_price_out.txt
         """

 }

Channel
        .fromFilePairs(params.raw_fastq, size: 1)
        .ifEmpty { exit 1, "Cannot find any fastq files matching: ${params.raw_fastq}" }
        .set { fastqcFiles }

process getTranscriptFasta{
	tag "transcriptsFasta"
	publishDir "$params.outdir/fastas", mode: params.publish_dir_mode

	input:
	set sampleID, file(raw_fastq) from fastqcFiles

	output:
	file '*transcript.fasta' into fasta1, fasta2

	script:
	"""
        cat ${params.outdir}/RiboCode/${sampleID}_ORFs_result_collapsed_out_out.txt_ribocode_out.txt ${params.outdir}/price/${sampleID}-hsv.orfs.tsv_price_out.txt |  awk '{split(\$2,a,".");print a[1]}' >listFiles
        while read lines; do cat ${params.lncfasta} | awk -v var="\${lines}" 'BEGIN{RS=">"}; \$0 ~ var {print ">"\$0}';done<listFiles | sed -r '/^\s*\$/d'>${sampleID}_transcript.fasta
	"""
}

Channel
        .fromFilePairs(params.raw_fastq, size: 1)
        .ifEmpty { exit 1, "Cannot find any fastq files matching: ${params.raw_fastq}" }
        .set { fastqcFiles2 }

process DeepCPP{

        tag "deepcpp"
        publishDir "$params.outdir/DeepCPP", mode: params.publish_dir_mode
        container 'ubuntu-tfm:v1'
	memory='16 GB'
	
	
        input:
	set sampleID, file(raw_fastq) from fastqcFiles2
	path x from fasta1	

	output:
	file '*_predict_results.csv' into deepcpp

        script:
        """
        /opt/bin/DeepCPP.py -i /opt/bin/input_files/ -f $params.basedocker/results/fastas/${sampleID}_transcript.fasta -o ${sampleID}_ -s human -t sorf
	"""
}

Channel
        .fromFilePairs(params.raw_fastq, size: 1)
        .ifEmpty { exit 1, "Cannot find any fastq files matching: ${params.raw_fastq}" }
        .set { fastqcFiles3 }

process Mipepid{
	
        tag "mipepid"
        publishDir "$params.outdir/Mipepid", mode: params.publish_dir_mode
        container 'ubuntu-tfm:v1'

        input:
         tuple val(sampleID), file (fastq) from fastqcFiles3
	 path x from fasta2

	output:
	 file '*_mipepid_result' into mipepid
	 
        script:
        """
        /opt/bin/mipepid.py $params.basedocker/results/fastas/${sampleID}_transcript.fasta ${sampleID}_mipepid_result
        """
}


// Pipeline execution characteristics

workflow.onComplete {
    println "Pipeline completed!"
    println "Started at $workflow.start"
    println "Finished at $workflow.complete"
    println "Time elapsed: $workflow.duration"
   }


