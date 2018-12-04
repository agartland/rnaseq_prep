#!/usr/bin/env nextflow

// nextflow qc_trim.nf --profile awsbatch -w s3://fh-pi-gilbert-p/agartlan/rnaseq_working 
params.output = 's3://fh-pi-gilbert-p/agartlan/rnaseq-test/'
params.input = 's3://fh-pi-gilbert-p/agartlan/rnaseq-test/9*_R{1,2}_001.fastq.gz'

/*params.input = '/home/agartlan/fast/LamarAndrew/RawFASTQexamples/9*_R{1,2}_001.fastq.gz'
params.output = '.'
params.index = '/home/agartlan/fast/LamarAndrew/UCSC_h38/h38_refMrna_index'*/

Channel.fromFilePairs(params.input, flat:true).into{infiles1_ch; infiles2_ch}

process qc_and_trim {
    executor 'awsbatch'
    container 'quay.io/afioregartland/rnaseq_prep'
    queue 'mixed'
    cpus 8
    memory '8GB'
    
    tag "QC trim of $seqid on $forfile and $revfile"

    input:
    set seqid, file(forfile), file(revfile) from infiles1_ch

    output:
    file("pre_fastqc_${forfile}_logs") into pre_P1_fastqc_ch
    file("pre_fastqc_${revfile}_logs") into pre_P2_fastqc_ch

    file("post_fastqc_${forfile}_logs") into post_P1_fastqc_ch
    file("post_fastqc_${revfile}_logs") into post_P2_fastqc_ch

    file("trimmomatic_output/${seqid}.log") into trimmomatic_logs_ch


    script: 
    """
    mkdir pre_fastqc_${forfile}_logs
    mkdir pre_fastqc_${revfile}_logs
    fastqc -t $task.cpus -o pre_fastqc_${forfile}_logs -f fastq -q $forfile
    fastqc -t $task.cpus -o pre_fastqc_${revfile}_logs -f fastq -q $revfile

    mkdir trimmomatic_output
    TrimmomaticPE -threads $task.cpus -trimlog trimmomatic_output/${seqid}_dump.log \
        $forfile $revfile\
        filtered_P1_${seqid}.fastq.gz filtered_UP1_dump.fastq.gz \
        filtered_P2_${seqid}.fastq.gz filtered_UP2_dump.fastq.gz \
        SLIDINGWINDOW:10:20 MINLEN:30 \
        2> trimmomatic_output/${seqid}.log

    mkdir post_fastqc_${forfile}_logs
    mkdir post_fastqc_${revfile}_logs
    fastqc -t $task.cpus -o post_fastqc_${forfile}_logs -f fastq -q filtered_P1_${seqid}.fastq.gz
    fastqc -t $task.cpus -o post_fastqc_${revfile}_logs -f fastq -q filtered_P2_${seqid}.fastq.gz

    """
}
/*
process pre_fastqc {
    container = 'quay.io/afioregartland/rnaseq_prep'
    
    tag "Pre-FASTQC of $seqid on $forfile and $revfile"

    input:
    set seqid, file(forfile), file(revfile) from infiles1_ch

    output:
    file("pre_fastqc_${forfile}_logs") into pre_P1_fastqc_ch
    file("pre_fastqc_${revfile}_logs") into pre_P2_fastqc_ch

    script: 
    """
    mkdir pre_fastqc_${forfile}_logs
    mkdir pre_fastqc_${revfile}_logs
    fastqc -t $task.cpus -o pre_fastqc_${forfile}_logs -f fastq -q $forfile
    fastqc -t $task.cpus -o pre_fastqc_${revfile}_logs -f fastq -q $revfile
    """
}

process trim {
    container = 'quay.io/afioregartland/rnaseq_prep'
    
    tag "trimmomatic on $seqid with $forfile and $revfile"

    input:
    set seqid, file(forfile), file(revfile) from infiles2_ch

    output:
    file("trimmomatic_output/${seqid}.log") into trimmomatic_logs_ch
    set file("filtered_P1_${seqid}.fastq.gz"), file("filtered_P2_${seqid}.fastq.gz") into trimmed_ch
    //set file("filtered_P1_${seqid}.fastq.gz") into trimmed1_ch
    //set file("filtered_P2_${seqid}.fastq.gz") into trimmed2_ch

    script: 
    """
    mkdir trimmomatic_output
    TrimmomaticPE -threads $task.cpus -trimlog trimmomatic_output/${seqid}_dump.log \
        $forfile $revfile\
        filtered_P1_${seqid}.fastq.gz filtered_UP1_dump.fastq.gz \
        filtered_P2_${seqid}.fastq.gz filtered_UP2_dump.fastq.gz \
        SLIDINGWINDOW:10:20 MINLEN:30 \
        2> trimmomatic_output/${seqid}.log
    """
}
*/
/*trimmed_ch.into{trimmed_fastqc_ch; trimmed_squant_ch}
trimmed_squant_ch.collect().into{collected_trimmed_squant_ch}*/

/*trimmed1_ch.into{trimmed1_fastqc_ch; trimmed1_squant_ch}
trimmed2_ch.into{trimmed2_fastqc_ch; trimmed2_squant_ch}
trimmed1_fastqc_ch.mix(trimmed2_fastqc_ch).into{trimmed_fastqc_ch}
trimmed1_squant_ch.collect().into{collected1_trimmed_squant_ch}
trimmed2_squant_ch.collect().into{collected2_trimmed_squant_ch}*/

//tmp1_ch.merge(tmp2_ch).into{trimmed_fastqc_ch; trimmed_squant_ch}
//tmp1.into{trimmomatic_P1_filtered_ch; trimmed_P1_ch}
//tmp2.into{trimmomatic_P2_filtered_ch; trimmed_P2_ch}
/*
process post_fastqc {
    container = 'quay.io/afioregartland/rnaseq_prep'
    
    tag "Post-FASTQC on $forfile and $revfile"

    input:
    set file(forfile),file(revfile) from trimmed_ch

    output:
    file("post_fastqc_${forfile}_logs") into post_P1_fastqc_ch
    file("post_fastqc_${revfile}_logs") into post_P2_fastqc_ch

    script: 
    """
    mkdir post_fastqc_${forfile}_logs
    mkdir post_fastqc_${revfile}_logs
    fastqc -t $task.cpus -o post_fastqc_${forfile}_logs -f fastq -q $forfile
    fastqc -t $task.cpus -o post_fastqc_${revfile}_logs -f fastq -q $revfile
    """
}*/

process multiqc {
    container 'quay.io/afioregartland/rnaseq_prep'
    publishDir params.output, mode:'copy'
       
    input:
    file('*') from pre_P1_fastqc_ch.mix(pre_P2_fastqc_ch).mix(trimmomatic_logs_ch).mix(post_P1_fastqc_ch).mix(post_P2_fastqc_ch).collect()
    
    output:
    file('multiqc_report.html')  
     
    script:
    """
    multiqc . 
    """
}

workflow.onComplete { 
    println ( workflow.success ? "\nDone! Open the following report in your browser --> $params.output/multiqc_report.html\n" : "Oops .. something went wrong" )
}