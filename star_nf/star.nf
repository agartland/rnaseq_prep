#!/usr/bin/env nextflow

// nextflow star.nf -profile awscloud -w s3://fh-pi-gilbert-p/agartlan/rnaseq_working --url s3://fh-pi-gilbert-p/agartlan/TBVPX203-RNA
// nextflow star.nf -profile local --max_libs 1

params.genome = 's3://fh-pi-gilbert-p/agartlan/Ensembl_hg38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz'
params.gtf = 's3://fh-pi-gilbert-p/agartlan/Ensembl_hg38/Homo_sapiens.GRCh38.94.gtf.gz'

// params.url = 's3://fh-pi-gilbert-p/agartlan/HVTN602-RNA'
// params.url = 's3://fh-pi-gilbert-p/agartlan/TBVPX203-RNA'
params.max_libs = -1

Channel.fromPath("$params.url" + '/quant_keys.csv')
        .splitCsv(header: true)
        .take(params.max_libs)
        .set { infiles_ch }

genome_file = file(params.genome)
gtf_file = file(params.gtf)

process star_quant {

    label 'multicore'
    
    publishDir path:"$prefix", mode:'copy'
    
    tag "star on $sampid at $prefix"
    
    input:
    set sampid, prefix from infiles_ch

    file genome from genome_file
    file gtf from gtf_file
 
    output:
    file("star_quant_${sampid}") into quant_folder_ch

    script:
    """
    aws s3 cp "${prefix}/" ./ --recursive --exclude "*" --include "*.fastq.gz" --quiet &

    mkfifo gen_fifo gtf_fifo
    zcat $genome > gen_fifo &
    zcat $gtf > gtf_fifo &

    mkdir index
    star --runThreadN $task.cpus --runMode genomeGenerate \
         --genomeDir ./index \
         --genomeFastaFiles gen_fifo \
         --sjdbGTFfile gtf_fifo --sjdbOverhang 49

    wait

    mkdir star_quant_${sampid}
    mkfifo R1 R2
    zcat *_R1_*.fastq.gz > R1 &
    zcat *_R2_*.fastq.gz > R2 &
    star --runThreadN $task.cpus \
         --genomeDir ./index \
         --readFilesIn R1 R2 \
         --quantMode GeneCounts \
         --outFileNamePrefix ./star_quant_${sampid}
    """
}

process multiqc {

    tag "Post-star multiqc report generation"

    label 'singlecore'

    publishDir path:params.url, mode:'copy'
       
    input:
    file('*') from quant_folder_ch.collect()
    
    output:
    file('post_star_report.html')  
     
    script:
    """
    multiqc -o ./ -n post_star_report.html . 
    """
}

workflow.onComplete { 
    println ( workflow.success ? "\nSTAR quantification report saved to $params.url\n" : "ERROR: quantification not completed" )
}