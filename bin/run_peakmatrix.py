#!/projects/ps-renlab/y2xie/anaconda3/envs/seurat/bin/python

import argparse

parser = argparse.ArgumentParser(description='generate gene matrix and TSS matrix ')
parser.add_argument('--infile', type=str, dest="infile", help='input h5ad snapatac2 file')
parser.add_argument('--outprfx', type=str, dest="outprfx", help='output prefix')
parser.add_argument('--peak', type=str, dest="peakf", help='peak bed file to calculate')

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
    peakf = str(args.peakf)    
    DPT_process(inf, outp, peakf)
    end_time = pc()
    print('Used (secs): ', end_time - start_time)

def DPT_process(inf, outp, peakf):
    data = snap.read(inf, backed=None)
    ### peak matrix
    peak_matrix = snap.pp.make_peak_matrix(data, peak_file=peakf,  max_frag_size = 5000, min_frag_size = 10)
    peak_matrix.write(outp + '_peak_matrix.h5ad', compression='gzip')


if __name__ == "__main__":
    run()


