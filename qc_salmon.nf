#!/usr/bin/env nextflow
 
// params.transcriptome = "ftp://ftp.ensembl.org/pub/release-90/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"
params.transcriptome = 'http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/refMrna.fa.gz'
params.metadata = "$baseDir/testmetadata2.csv"
params.output = "."
params.max_samples = 2

metadata_file = file(params.metadata)
// transcriptome_file = file(params.transcriptome)

if( !metadata_file.exists() ) error "Metadata file does not exist: $metadata_file"
// if( !transcriptome_file.exists() ) error "Transcriptome file does not exist: $transcriptome_file" 

/*
process createIndex {
    
    tag "$transcriptome"

    input:
    file transcriptome from transcriptome_file  
 
    output:
    file 'index' into index_ch

    script:
    """
    salmon index --threads $task.cpus -t $transcriptome -i index
    """
}
 */

process parseMetadata {

    tag "$params.metadata"

    input:
    file(metadata) from metadata_file

    output:
    stdout into url_csv_ch

    """
    #!/usr/bin/env python3
    import pandas as pd
    import sys
    import re
    
    sra = []
    df = pd.read_csv("$metadata")
    for i,row in df.iterrows():
        sample,strand = re.match('.*(4214[\\w_]+)_(R[12])', row['filename']).groups()
        print(",".join([sample, strand, row['filename']]))
        # print(row['filename'])
    """
}

/* 
 * Parse the CSV file 
 * Takes only the first `max_samples` entries
 * Creates two separate channels for PE reads
 */

url_csv_ch
    .splitCsv()
    .take(params.max_samples)
    .set {files_ch}

//files_ch = Channel.fromPath(url_csv_ch)
Channel.fromPath('s3://fh-pi-gilbert-p/agartlan/rnaseq-test/HVTN602-M/421400319V03BC006-50055/421400319V03BC006_S6_L001_R2_001.fastq.gz').println() 
//myFile = file('s3://fh-pi-gilbert-p/agartlan/rnaseq-test/HVTN602-M/421400319V03BC006-50055/421400319V03BC006_S6_L001_R2_001.fastq.gz')
//println myFile.text


process pre_fastqc {
    
    tag "FASTQC on $url"

    input:
    set sample,strand, file(url) from files_ch
    //file input_fasta from files_ch

    output:
    file("pre_fastqc_${sample}_${strand}") into fastqc_ch

    script: 
    // note -- fastq skips quietly any missing read files, add an explicit check to stop the task if any input is missing
    """
    #!/usr/bin/env python3

    print('$url')
    print('$sample', '$strand')
    with open('$url', 'r') as fh:
        print(fh.readline())
    """

} 

/*
    """
    mkdir /scratch/logs/pre_fastqc_${sample}_${strand}
    fastqc -t $task.cpus -o pre_fastqc_${sample}_${strand} -f fastq -q $url
    """  
*/

/*
process quant {
    
    tag "$dbxref"
    
    input:
    file index from index_ch
    set dbxref,sample_type,strand_specific,url from encode_files_ch1
 
    output:
    file("${sample_type}-${dbxref}") into quant_ch

    script:
    def libType = strand_specific == "True" ? "SF" : "U"
    """
    wget -q ${url}/${dbxref}_1.fastq.gz &
    wget -q ${url}/${dbxref}_2.fastq.gz &
    wait 
    salmon quant --threads $task.cpus --libType=${libType} -i index -1 ${dbxref}_1.fastq.gz -2 ${dbxref}_2.fastq.gz -o ${sample_type}-${dbxref}
    """
}
*/ 

/* 
process multiqc {
    
    publishDir params.output, mode:'copy'
       
    input:
    file('*') from quant_ch.mix(fastqc_ch).collect()
    
    output:
    file('multiqc_report.html')  
     
    script:
    """
    multiqc . 
    """
}
*/ 

/*
workflow.onComplete { 
    println ( workflow.success ? "\nDone! Open the following report in your browser --> $params.output/multiqc_report.html\n" : "Oops .. something went wrong" )
}
*/