#!/usr/bin/env python3

"""
https://raw.githubusercontent.com/agartland/rnaseq_prep/master/qc_salmon.py

Usage examples:
python3 qc_salmon.py --r1 s3://afg-rnaseq-test/421400319V05BC019_S18_L001_R1_001.fastq.gz --r2 s3://afg-rnaseq-test/421400319V05BC019_S18_L001_R2_001.fastq.gz



Basic docker
docker info
docker stats
docker ps

docker run --name rna -tdi -v D:/data/Docker:/shared/data -m="3g" quay.io/afioregartland/rnaseq_prep
docker attach rna

LAPTOP
docker run --name rna -tdi -v C:/Andrew/data:/shared/data -m="3g" quay.io/afioregartland/rnaseq_prep
docker attach rna

Salmon index making:
salmon index -t /shared/data/UCSC_h38/refMrna.fa.gz -i /shared/data/UCSC_h38/h38_refMrna_index

salmon index -t /shared/data/UCSC_h38/mrna.fa.gz -i /shared/data/UCSC_h38/h38_mrna_index (too many short transcripts in mrna file)

ps aux --sort --rss

Test FASTQ file:
/shared/data/Raw\ FASTQ\ file\ examples/HVTN602-M/421400319V03BC006-50055/421400319V03BC006_S6_L001_R1_001.fastq.gz

FASTQC example (GC content and quality filtering):
./fastqc /shared/data/Raw\ FASTQ\ file\ examples/HVTN602-M/421400319V03BC006-50055/421400319V03BC006_S6_L001_R1_001.fastq.gz -O /shared/data/qc_out

fastqc /shared/data/Docker/filtered_1P.fq.gz -O /shared/data/Docker

Trimmomatic for trimming bases off the 3' end that are low quality
TrimmomaticPE -threads 2 -trimlog /shared/data/Docker/trimlog.log -basein /shared/data/Docker/161019_A02_1_C02_TAAGGCGA-GCGTAAGA_L001_R1_002.fastq.gz -baseout /shared/data/Docker/filtered.fq.gz SLIDINGWINDOW:10:20 MINLEN:30

https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-0956-2

Salmon for transcript quantifcation
salmon quant -i /shared/data/Docker/h38_refMrna_index -l A -1 /shared/data/Docker/filtered_1P.fq.gz -2 /shared/data/Docker/filtered_2P.fq.gz -o filtered_transcripts_quant

Tximport for making gene-level count matrix


deseq2 for testing for differential gene expression
"""

import argparse
import subprocess
import re
import boto3
import botocore

def parse_s3_url(url):
    bucket, key = [s for s in re.match(r's3://([-\w]+)/([-\w\.]+)', url).groups()]
    return bucket,key

def grabfile(url, localfile):
    bucket, key = parse_s3_url(url)
    s3 = boto3.resource('s3')

    try:
        s3.Bucket(bucket).download_file(key, localfile)
    except botocore.exceptions.ClientError as e:
        if e.response['Error']['Code'] == "404":
            print("The object does not exist.")
        raise

def putfile(url, localfile):
    bucket, key = parse_s3_url(url)
    s3 = boto3.client('s3')
    s3.upload_file(localfile, bucket, key)

def runFASTQC(localfile):
    cmd = ['fastqc',
            localfile,
            '--noextract',
            '-q']
    subprocess.check_call(' '.join(cmd), shell=True)
    outfile = localfile.replace('.gz', '') + '_fastqc.html'
    return outfile

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Fetch PE reads from S3 for: QC, trimming, salmon mapping and gene-level quantification.')
    parser.add_argument('--r1', type=str,
                        help='Location of R1 reads file (S3 bucket or filesystem)')
    parser.add_argument('--r2', type=str,
                        help='Location of R2 reads file (S3 bucket or filesystem)')
    args = parser.parse_args()

    r1Fn = 'R1.fa.gz'
    r2Fn = 'R2.fa.gz'

    grabfile(args.r1, r1Fn)
    grabfile(args.r2, r2Fn)

    outfile1 = runFASTQC(r1Fn)
    outfile2 = runFASTQC(r2Fn)

    putfile(args.r1 + '_fastqc.html', outfile1)
    putfile(args.r1 + '_fastqc.html', outfile2)