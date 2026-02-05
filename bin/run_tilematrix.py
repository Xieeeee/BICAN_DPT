#!/projects/ps-renlab/y2xie/anaconda3/envs/seurat/bin/python

import argparse

parser = argparse.ArgumentParser(description='generate gene matrix and TSS matrix ')
parser.add_argument('--infile', type=str, dest="infile", help='input h5ad snapatac2 file')
parser.add_argument('--outprfx', type=str, dest="outprfx", help='output prefix')
parser.add_argument('--res', type=str, dest="res", help='resolution to calculate tile matrix')

args = parser.parse_args()

import numpy as np
import pandas as pd
from scipy import stats as ss
import snapatac2 as snap
import scanpy as sc
import anndata
import os
from time import perf_counter as pc

def run():
    start_time = pc()
    ### read args
    inf = str(args.infile)
    outp = str(args.outprfx)
    ress = str(args.res)    
    DPT_process(inf, outp, ress)
    end_time = pc()
    print('Used (secs): ', end_time - start_time)

def format_genomic_size(size_bp):
    """
    Convert genomic size (in base pairs) to bp, kb, Mb, or Gb
    without decimals.
    """
    size_bp = int(size_bp)
    if size_bp < 1_000:
        return f"{int(size_bp)}bp"
    elif size_bp < 1_000_000:
        return f"{int(size_bp // 1_000)}kb"
    elif size_bp < 1_000_000_000:
        return f"{int(size_bp // 1_000_000)}Mb"
    else:
        return f"{int(size_bp // 1_000_000_000)}Gb"

def DPT_process(inf, outp, ress):
    data = snap.read(inf, backed=None)
    prefix=format_genomic_size(ress)
    ### peak matrix
    peak_matrix = snap.pp.add_tile_matrix(data, bin_size=int(ress), counting_strategy = 'fragment', max_frag_size = 5000, min_frag_size = 10, n_jobs = 2, inplace=False)
    peak_matrix.write(f'{outp}_{prefix}.h5ad', compression='gzip')


if __name__ == "__main__":
    run()


