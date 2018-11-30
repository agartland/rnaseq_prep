#!/usr/bin/env nextflow

//params.output = 's3://fh-pi-gilbert-p/agartlan/rnaseq-test/'
//params.input = 's3://fh-pi-gilbert-p/agartlan/rnaseq-test/*.fastq.gz'

params.input = '/home/agartlan/fast/LamarAndrew/RawFASTQexamples/9*_R{1,2}_001.fastq.gz'
params.output = '.'

//Channel.fromPath(params.input).into{infiles1_ch; infiles2_ch}
Channel.fromFilePairs(params.input).into{infiles1_ch; infiles2_ch}
//infiles_ch = Channel.fromPath('s3://fh-pi-gilbert-p/agartlan/rnaseq-test/shorty.fastq.gz')
//outfiles_ch = Channel.fromPath('s3://fh-pi-gilbert-p/agartlan/rnaseq-test/')
//myFile = file('s3://fh-pi-gilbert-p/agartlan/rnaseq-test/HVTN602-M/421400319V03BC006-50055/421400319V03BC006_S6_L001_R2_001.fastq.gz')
//println myFile.text

process pre_fastqc {
    
    tag "Pre-FASTQC of $seqid on $infiles"

    input:
    set seqid,file(infiles) from infiles1_ch

    output:
    file("pre_fastqc_${seqid}_logs") into pre_fastqc_ch

    script: 
    """
    mkdir pre_fastqc_${seqid}_logs
    fastqc -t $task.cpus -o pre_fastqc_${seqid}_logs -f fastq -q $infiles
    """
}

process trim {
    
    tag "trimmomatic on $seqid"

    input:
    set seqid,file(infiles) from infiles2_ch
    //file input_fasta from files_ch

    output:
    file("trimmomatic_output/${seqid}.log") into trimmomatic_logs_ch
    file("filtered_P1_${seqid}") into trimmomatic_P1_filtered_ch
    file("filtered_P2_${seqid}") into trimmomatic_P2_filtered_ch
    //file("filtered_UP_${forfile}") into trimmomatic_UP1_filtered_ch
    //file("filtered_UP_${revfile}") into trimmomatic_UP2_filtered_ch

    script: 
    """
    mkdir trimmomatic_output
    TrimmomaticPE -threads $task.cpus -trimlog trimmomatic_output/${seqid}_dump.log \
        $infiles \
        filtered_P1_${seqid} filtered_UP1_dump.fastq.gz \
        filtered_P2_${seqid} filtered_UP2_dump.fastq.gz \
        SLIDINGWINDOW:10:20 MINLEN:30 \
        2> trimmomatic_output/${seqid}.log
    """
}

process post_P1_fastqc {
    
    tag "Post-FASTQC on P1 $infile"

    input:
    file infile from trimmomatic_P1_filtered_ch

    output:
    file("post_fastqc_${infile}_logs") into post_P1_fastqc_ch

    script: 
    """
    mkdir post_fastqc_${infile}_logs
    fastqc -t $task.cpus -o post_fastqc_${infile}_logs -f fastq -q $infile
    """
}

process post_P2_fastqc {
    
    tag "Post-FASTQC on P2 $infile"

    input:
    file infile from trimmomatic_P2_filtered_ch

    output:
    file("post_fastqc_${infile}_logs") into post_P2_fastqc_ch

    script: 
    """
    mkdir post_fastqc_${infile}_logs
    fastqc -t $task.cpus -o post_fastqc_${infile}_logs -f fastq -q $infile
    """
}



process multiqc {
    
    publishDir params.output, mode:'copy'
       
    input:
    file('*') from pre_fastqc_ch.mix(trimmomatic_logs_ch).mix(post_P1_fastqc_ch).mix(post_P2_fastqc_ch).collect()
    
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