#!/usr/bin/env nextflow

// nextflow salmon.nf -profile awscloud -w s3://fh-pi-gilbert-p/agartlan/rnaseq_working --max_libs 2
// nextflow salmon.nf -profile local --max_libs 1

// params.transcriptome = 's3://fh-pi-gilbert-p/agartlan/hg38/refMrna_reduced.fa.gz'
params.transcriptome = 's3://fh-pi-gilbert-p/agartlan/hg38/refMrna.fa.gz'
params.url = 's3://fh-pi-gilbert-p/agartlan/HVTN602-RNA'
params.max_libs = -1

Channel.fromPath("$params.url" + '/quant_keys.csv')
        .splitCsv(header: true)
        .take(params.max_libs)
        .set { infiles_ch }

transcriptome_file = file(params.transcriptome)
process index {

    label 'multicore'

    tag "Generating salmon index from $transcriptome"

    input:
    file transcriptome from transcriptome_file  
 
    output:
    file 'index' into index_ch

    script:
    """
    salmon index --threads $task.cpus -t $transcriptome -i index
    """
}

process salmon_quant {

    label 'multicore'
    
    publishDir path:"$prefix", mode:'copy'
    
    tag "salmon on $sampid at $prefix"
    
    input:
    file index from index_ch
    set sampid, prefix from infiles_ch
 
    output:
    file("squant_filtered_${sampid}_logs") into quant_folder_ch

    script:
    """
    aws s3 cp "${prefix}/" ./ --recursive --exclude "*" --include "filtered*.fastq.gz" --quiet

    # zcat \$(ls filtered_P1_*.fastq.gz) | head -n 1000 | gzip > red_filtered_P1_red.fastq.gz
    # zcat \$(ls filtered_P2_*.fastq.gz) | head -n 1000 | gzip > red_filtered_P2_red.fastq.gz

    salmon quant --threads $task.cpus --gcBias\
            -i $index -l A \
            -1 \$(ls filtered_P1_*.fastq.gz) \
            -2 \$(ls filtered_P2_*.fastq.gz) \
            -o squant_filtered_${sampid}

    mkdir tmp_${sampid}
    mv squant_filtered_${sampid}/quant.sf tmp_${sampid}/quant.sf
    cp -r squant_filtered_${sampid} squant_filtered_${sampid}_logs
    mv tmp_${sampid}/quant.sf squant_filtered_${sampid}/quant.sf
    rm -r tmp_${sampid}
    """
}

process multiqc {

    tag "Post-salmon multiqc report generation"

    label 'singlecore'

    publishDir path:params.url, mode:'copy'
       
    input:
    file('*') from quant_folder_ch.collect()
    
    output:
    file('post_salmon_report.html')  
     
    script:
    """
    multiqc -o ./ -n post_salmon_report.html . 
    """
}

workflow.onComplete { 
    println ( workflow.success ? "\nQuantification report saved to $params.url\n" : "ERROR: quantification not completed" )
}