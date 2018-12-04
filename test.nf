#!/usr/bin/env nextflow

// nextflow test.nf -w s3://fh-pi-gilbert-p/agartlan/rnaseq_working
//params.output = 's3://fh-pi-gilbert-p/agartlan/rnaseq-test/'
//params.input = 's3://fh-pi-gilbert-p/agartlan/rnaseq-test/*.fastq.gz'

params.input = '/home/agartlan/fast/LamarAndrew/RawFASTQexamples/9*_R{1,2}_001.fastq.gz'
params.output = '.'
params.index = '/home/agartlan/fast/LamarAndrew/UCSC_h38/h38_refMrna_index'
index_ch = Channel.fromPath(params.index).collect()

//Channel.fromPath(params.input).into{infiles1_ch; infiles2_ch}
Channel.fromFilePairs(params.input).into{infiles1_ch; infiles2_ch}
//infiles_ch = Channel.fromPath('s3://fh-pi-gilbert-p/agartlan/rnaseq-test/shorty.fastq.gz')
//outfiles_ch = Channel.fromPath('s3://fh-pi-gilbert-p/agartlan/rnaseq-test/')
//myFile = file('s3://fh-pi-gilbert-p/agartlan/rnaseq-test/HVTN602-M/421400319V03BC006-50055/421400319V03BC006_S6_L001_R2_001.fastq.gz')
//println myFile.text

process pre_fastqc {
    container = 'quay.io/afioregartland/rnaseq_prep'
    
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
    container = 'quay.io/afioregartland/rnaseq_prep'
    
    tag "trimmomatic on $seqid"

    input:
    set seqid,file(infiles) from infiles2_ch
    //file input_fasta from files_ch

    output:
    file("trimmomatic_output/${seqid}.log") into trimmomatic_logs_ch
    // set file("filtered_P1_${seqid}.fastq.gz"), file("filtered_P2_${seqid}.fastq.gz") into trimmed_ch
    set file("filtered_P1_${seqid}.fastq.gz") into trimmed1_ch
    set file("filtered_P2_${seqid}.fastq.gz") into trimmed2_ch

    script: 
    """
    mkdir trimmomatic_output
    TrimmomaticPE -threads $task.cpus -trimlog trimmomatic_output/${seqid}_dump.log \
        $infiles \
        filtered_P1_${seqid}.fastq.gz filtered_UP1_dump.fastq.gz \
        filtered_P2_${seqid}.fastq.gz filtered_UP2_dump.fastq.gz \
        SLIDINGWINDOW:10:20 MINLEN:30 \
        2> trimmomatic_output/${seqid}.log
    """
}

/*trimmed_ch.into{trimmed_fastqc_ch; trimmed_squant_ch}
trimmed_squant_ch.collect().into{collected_trimmed_squant_ch}*/

trimmed1_ch.into{trimmed1_fastqc_ch; trimmed1_squant_ch}
trimmed2_ch.into{trimmed2_fastqc_ch; trimmed2_squant_ch}
trimmed1_fastqc_ch.mix(trimmed2_fastqc_ch).into{trimmed_fastqc_ch}
trimmed1_squant_ch.collect().into{collected1_trimmed_squant_ch}
trimmed2_squant_ch.collect().into{collected2_trimmed_squant_ch}

//tmp1_ch.merge(tmp2_ch).into{trimmed_fastqc_ch; trimmed_squant_ch}
//tmp1.into{trimmomatic_P1_filtered_ch; trimmed_P1_ch}
//tmp2.into{trimmomatic_P2_filtered_ch; trimmed_P2_ch}

process post_fastqc {
    container = 'quay.io/afioregartland/rnaseq_prep'
    
    tag "Post-FASTQC on $forfile and $revfile"

    input:
    set file(forfile),file(revfile) from trimmed_fastqc_ch

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
}

process salmon_quant {
    container = 'quay.io/afioregartland/rnaseq_prep'

    publishDir params.output, mode:'copy'
    
    //tag "salmon on $trimmed11, $trimmed21, $trimmed12, $trimmed22"
    tag "salmon $trimmed1 and $trimmed2"
    
    input:
    file index from index_ch
    //set file(trimmed11), file(trimmed21), file(trimmed12), file(trimmed22) from collected_trimmed_squant_ch
    set file(trimmed1) from collected1_trimmed_squant_ch
    set file(trimmed2) from collected2_trimmed_squant_ch
 
    output:
    file("squant_all")

    script:
    """
    salmon quant --threads $task.cpus --gcBias\
            -i $index -l A -1 ${trimmed1} -2 ${trimmed2}\
            -o squant_all
    """
}

process multiqc {
    container = 'quay.io/afioregartland/rnaseq_prep'
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