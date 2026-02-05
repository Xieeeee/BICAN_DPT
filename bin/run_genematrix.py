#!/projects/ps-renlab/y2xie/anaconda3/envs/seurat/bin/python

import argparse

parser = argparse.ArgumentParser(description='generate gene matrix and TSS matrix ')
parser.add_argument('--infile', type=str, dest="infile", help='input h5ad snapatac2 file')
parser.add_argument('--outprfx', type=str, dest="outprfx", help='output prefix')

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
    DPT_process(inf, outp)
    end_time = pc()
    print('Used (secs): ', end_time - start_time)

def DPT_process(inf, outp):
    data = snap.read(inf, backed=None)
    ### gene matrix
    gene_matrix = snap.pp.make_gene_matrix(data, snap.genome.hg38, max_frag_size = 5000, min_frag_size = 10)
    gene_matrix.write(outp + '_gene_matrix_raw.h5ad', compression='gzip')
    # sc.pp.filter_genes(gene_matrix, min_cells= 5)
    # sc.pp.normalize_total(gene_matrix)
    # sc.pp.log1p(gene_matrix)
    # sc.external.pp.magic(gene_matrix, solver="approximate")
    # gene_matrix.write(outp + '_gene_matrix.h5ad', compression='gzip')
    ### TSS matrix
    peak_matrix = snap.pp.make_peak_matrix(data, peak_file="/projects/ps-renlab2/y2xie/ps-renlab/y2xie/projects/genome_ref/refdata-cellranger-GRCh38.p13_v43/GRCh38/genes/tss_1500bp.bed",  max_frag_size = 5000, min_frag_size = 10)
    peak_matrix.write(outp + '_TSS_matrix.h5ad', compression='gzip')


if __name__ == "__main__":
    run()


