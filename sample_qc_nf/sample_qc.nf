#!/usr/bin/env nextflow

// nextflow sample_qc.nf -profile awscloud -w s3://fh-pi-gilbert-p/agartlan/rnaseq_working --url s3://fh-pi-gilbert-p/agartlan/TBVPX203-RNA
// nextflow sample_qc.nf -profile local --max_libs 1

// params.url = 's3://fh-pi-gilbert-p/agartlan/HVTN602-RNA'
// params.url = 's3://fh-pi-gilbert-p/agartlan/TBVPX203-RNA'
params.max_libs = -1

Channel.fromPath("$params.url" + '/quant_keys.csv')
        .splitCsv(header: true)
        .take(params.max_libs)
        .set { infiles_ch }

process sample_fastqc {
    label 'multicore'
    
    //publishDir path:"$prefix", mode:'copy'
    
    tag "Sample-level QC of $sampid"
    
    input:
    set sampid, prefix from infiles_ch
 
    output:
    file("sample_fastqc_${sampid}_logs") into fastqc_folder_ch

    script:
    """
    aws s3 cp "${prefix}/" ./ --recursive --exclude "*" --include "*.fastq.gz" --quiet

    mkdir sample_fastqc_${sampid}_logs
    zcat *.fastq.gz | fastqc -t $task.cpus -o sample_fastqc_${sampid}_logs -f fastq -q stdin:${sampid}
    """
}

process sample_multiqc {

    tag "Sample-level fastqc multiqc report generation"

    label 'singlecore'

    publishDir path:params.url, mode:'copy'
       
    input:
    file('*') from fastqc_folder_ch.collect()
    
    output:
    file('sample_fastqc_report.html')  
     
    script:
    """
    multiqc -o ./ -n sample_fastqc_report.html . 
    """
}

workflow.onComplete { 
    println ( workflow.success ? "\nSample-level fastqc report saved to $params.url\n" : "ERROR: fastqc not completed" )
}