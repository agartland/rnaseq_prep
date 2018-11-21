#!/usr/bin/env python

"""
Usage examples:
python /home/agartlan/gitrepo/rnaseq_prep/qc_salmon.py --r1 --r2
"""

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Fetch PE reads from S3 for: QC, trimming, salmon mapping and gene-level quantification.')
    parser.add_argument('--r1', type=str,
                        help='Location of R1 reads file (S3 bucket or filesystem)')
    parser.add_argument('--r2', type=str,
                        help='Location of R2 reads file (S3 bucket or filesystem)')
    args = parser.parse_args()

    import itertools
    import pandas as pd
    import numpy as np
    from os.path import join as opj
    import os
    from functools import partial
    import time
    import sys
    import feather
    
    """Make sure the utils are on path before importing"""
    sys.path.append(args.utils)

