#!/usr/bin/env nextflow

// nextflow grouped_qc_trim.nf -profile awscloud -w s3://fh-pi-gilbert-p/agartlan/rnaseq_working
// nextflow grouped_qc_trim.nf -profile local --max_rows 2 -w /home/agartlan/gitrepo/rnaseq_prep/grouped_qc_trim/work

params.url = 's3://fh-pi-gilbert-p/agartlan/TBVPX203-RNA'
params.max_rows = -1

Channel.fromPath("$params.url" + '/grouped_qc_keys.csv')
        .splitCsv(header: true)
        .take(params.max_rows)
        .set { infiles_ch }

process qc_and_trim {

    //publishDir path:"$prefix", mode:'copy', pattern:'*.fastq.gz'

    time '5h'

    label 'multicore'
    
    tag "Grouped QC trim of $lane_read_id"

    input:
    set lane_read_id, prefix, file_list from infiles_ch

    output:
    file("pre_fastqc_${lane_read_id}_logs") into pre_fastqc_ch

    script: 
    """
    for FILE in $file_list; do
        aws s3 cp \$FILE ./ --quiet
    done
    wait

    mkdir pre_fastqc_${lane_read_id}_logs
    
    # zcat \$(ls *.fastq.gz) | head -n 2000 | gzip > ${lane_read_id}.fastq.gz
    # zcat \$(ls *.fastq.gz) | gzip > ${lane_read_id}.fastq.gz

    zcat *.fastq.gz | fastqc -t $task.cpus -o pre_fastqc_${lane_read_id}_logs -f fastq -q stdin:${lane_read_id}
    """
}

process pre_multiqc {

    tag "Pre-trimming multiqc report generation"

    label 'singlecore'
    
    publishDir path:params.url, mode:'copy'
       
    input:
    file('*') from pre_fastqc_ch.collect()
    
    output:
    file('pre_multiqc_report.html')  
     
    script:
    """
    multiqc -o ./ -n pre_multiqc_report.html . 
    """
}

workflow.onComplete { 
    println ( workflow.success ? "\nReports saved to: $params.url\n" : "ERROR: reports not generated" )
}