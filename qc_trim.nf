#!/usr/bin/env nextflow

// nextflow qc_trim.nf -profile awscloud -w s3://fh-pi-gilbert-p/agartlan/rnaseq_working
// nextflow qc_trim.nf -profile local --max_samples 1

//params.output = 's3://fh-pi-gilbert-p/agartlan/rnaseq-test/'
//params.input = 's3://fh-pi-gilbert-p/agartlan/rnaseq-test/9*_R{1,2}_001.fastq.gz'

/*params.input = '/home/agartlan/fast/LamarAndrew/RawFASTQexamples/9*_R{1,2}_001.fastq.gz'
params.output = '.'
params.index = '/home/agartlan/fast/LamarAndrew/UCSC_h38/h38_refMrna_index'*/

params.url = 's3://fh-pi-gilbert-p/agartlan/HVTN602-RNA'
params.max_samples = -1

/*params.parser = 's3://fh-pi-gilbert-p/agartlan/HVTN602-RNA/parse_keys.py'
// params.parser = 'parse_keys.py'

Channel.fromPath(params.parser).set{ parser_ch }
//Channel.fromFilePairs(params.input, flat:true).into{infiles1_ch; infiles2_ch}

process parse_bucket {
    executor 'local'

    label 'singlecore'

    tag "Parsing with $parser ($params.url)"

    input:
    file parser from parser_ch

    output:
    //stdout into parsed_csv_ch
    file(csv_out) into parsed_csv_ch

    script:
    """
    python $parser $params.url > csv_out
    """
}
*/
Channel.fromPath("$params.url" + '/qc_keys.csv')
        .splitCsv(header: true)
        .take(params.max_samples)
        .set { infiles_ch }

// tmp.subscribe { println "value: $it" }

process qc_and_trim {

    publishDir path:"$prefix", mode:'copy', pattern:'*.fastq.gz'

    label 'multicore'
    
    tag "QC trim of $seqid on $forfile and $revfile"

    input:
    set seqid, prefix, forfile, revfile from infiles_ch

    output:
    file("pre_fastqc_${seqid}_R1_logs") into pre_P1_fastqc_ch
    file("pre_fastqc_${seqid}_R2_logs") into pre_P2_fastqc_ch

    file("post_fastqc_${seqid}_R1_logs") into post_P1_fastqc_ch
    file("post_fastqc_${seqid}_R2_logs") into post_P2_fastqc_ch

    file("trimmomatic_output/${seqid}.log") into trimmomatic_logs_ch

    file("filtered_P1_${seqid}.fastq.gz")
    file("filtered_P2_${seqid}.fastq.gz")

    script: 
    """
    aws s3 cp "$forfile" ./${seqid}_R1.fastq.gz &
    aws s3 cp "$revfile" ./${seqid}_R2.fastq.gz &
    wait

    mkdir pre_fastqc_${seqid}_R1_logs
    mkdir pre_fastqc_${seqid}_R2_logs
    fastqc -t $task.cpus -o pre_fastqc_${seqid}_R1_logs -f fastq -q ${seqid}_R1.fastq.gz
    fastqc -t $task.cpus -o pre_fastqc_${seqid}_R2_logs -f fastq -q ${seqid}_R2.fastq.gz

    mkdir trimmomatic_output
    TrimmomaticPE -threads $task.cpus -trimlog trimmomatic_output/${seqid}_dump.log \
        ${seqid}_R1.fastq.gz ${seqid}_R2.fastq.gz\
        filtered_P1_${seqid}.fastq.gz filtered_UP1_dump.fastq.gz \
        filtered_P2_${seqid}.fastq.gz filtered_UP2_dump.fastq.gz \
        SLIDINGWINDOW:10:20 MINLEN:30 \
        2> trimmomatic_output/${seqid}.log

    mkdir post_fastqc_${seqid}_R1_logs
    mkdir post_fastqc_${seqid}_R2_logs
    fastqc -t $task.cpus -o post_fastqc_${seqid}_R1_logs -f fastq -q filtered_P1_${seqid}.fastq.gz
    fastqc -t $task.cpus -o post_fastqc_${seqid}_R2_logs -f fastq -q filtered_P2_${seqid}.fastq.gz
    """
}

process post_multiqc {

    tag "Post-trimming multiqc report generation"

    label 'singlecore'

    publishDir path:params.url, mode:'copy'
       
    input:
    file('*') from trimmomatic_logs_ch.mix(post_P1_fastqc_ch).mix(post_P2_fastqc_ch).collect()
    
    output:
    file('post_multiqc_report.html')  
     
    script:
    """
    multiqc -o ./ -n post_multiqc_report.html . 
    """
}

process pre_multiqc {

    tag "Pre-trimming multiqc report generation"

    label 'singlecore'
    
    publishDir path:params.url, mode:'copy'
       
    input:
    file('*') from pre_P1_fastqc_ch.mix(pre_P2_fastqc_ch).collect()
    
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