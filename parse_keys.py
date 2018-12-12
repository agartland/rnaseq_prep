#!/usr/bin/env python

import boto3
import re
import pandas as pd
from os.path import join as opj
import sys
import io

# bucketName = 'fh-pi-gilbert-p'
# prefix = 'agartlan/HVTN602-RNA'
# url = 's3://fh-pi-gilbert-p/agartlan/HVTN602-RNA'
#  python parse_keys.py s3://fh-pi-gilbert-p/agartlan/HVTN602-RNA

def parse_s3_url(url):
    bucket, key = [s for s in re.match(r's3://([^/]+)/(.*)', url).groups()]
    return bucket, key

def s3ls(bucketName, prefix):
    s3 = boto3.resource('s3')
    bucket = s3.Bucket(name=bucketName)
    paths = [obj.key for obj in b.objects.filter(Prefix=prefix)]
    return paths

def parseFilename(p, pattern, cols):
    try:
        t = re.match(pattern, p.split('/')[-1]).groups()
    except AttributeError:
        print(p)
        raise
    o = {c:e for c,e in zip(cols, t)}
    o['key'] = p
    return o

s3p = 's3://'

url = sys.argv[1]
bucketName, prefix = parse_s3_url(url)

s3 = boto3.resource('s3')
bucket = s3.Bucket(name=bucketName)
paths = [obj.key for obj in bucket.objects.filter(Prefix=prefix) if '.fastq.gz' in obj.key and not 'filtered' in obj.key]

pattern = r'(42\d\d\d\d\d\d\d)[Vv](\d+)([^_]*)_(S\d+)_(L\d+)_(R[12])_001.fastq.gz'
cols = ['ptid', 'visit', 'sample_code', 'sample_num', 'lane', 'read']
dlist = [parseFilename(p, pattern, cols) for p in paths]


plist = [{'full_path':opj(s3p, bucketName, p),
           'filename':p.split('/')[-1],
           'run_folder':p.split('/')[-2],
           'prefix':opj(s3p, bucketName, '/'.join(p.split('/')[:-1]))} for p in paths]

df = pd.concat((pd.DataFrame(dlist), pd.DataFrame(plist)), axis=1)

df.loc[:, 'visit'] = df['visit'].map(lambda s: 'V%d' % int(s))
df.loc[:, 'seqid'] = df.apply(lambda r: '_'.join(r[['ptid', 'visit', 'sample_num', 'run_folder', 'lane']]), axis=1)
df.loc[:, 'sampid'] = df.apply(lambda r: '_'.join(r[['ptid', 'visit', 'sample_num', 'run_folder']]), axis=1)

"""Only duplicates are for 421400175V05BC021-40052: Removed from raw data 12-5-2018 on email from Lamar"""
# df = df.loc[~df['key'].str.contains('flowG_175v5')]

ind = ['seqid', 'prefix', 'read']
assert df.loc[df.duplicated(subset=ind, keep=False)].sort_values('filename')['full_path'].shape[0] == 0

"""Emit pairs of read files"""
out = df.set_index(ind)['full_path'].unstack('read').reset_index()

sio = io.StringIO()
outcols = ['seqid', 'prefix', 'R1', 'R2']
out[outcols].to_csv(sio, header=True, index=False)
sio.seek(0)
s3.Bucket(bucketName).put_object(Key=opj(prefix, 'qc_keys.csv'), Body=sio.getvalue())
print('Wrote QC keys to %s' % opj(prefix, 'qc_keys.csv'))

"""Emit the seqid and prefix only"""
out = df.drop_duplicates(['sampid', 'prefix'])

sio = io.StringIO()
outcols = ['sampid', 'prefix']
out[outcols].to_csv(sio, header=True, index=False)
sio.seek(0)
s3.Bucket(bucketName).put_object(Key=opj(prefix, 'quant_keys.csv'), Body=sio.getvalue())
print('Wrote quant keys to %s' % opj(prefix, 'quant_keys.csv'))