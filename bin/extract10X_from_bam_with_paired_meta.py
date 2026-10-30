#!/bin/python

import argparse

parser = argparse.ArgumentParser(description='filter bam based on barcodes field')
parser.add_argument('--cluster', type=str, dest="cluster", help='two columns meta: preindex, Target')
parser.add_argument('--inbam', type=str, dest="inbam", help='directory path of input bam')
parser.add_argument('--outprfx', type=str, dest="outprfx", help='output prefix. Output will be outprfx_[Target].sam')

args = parser.parse_args()

import numpy as np
import pandas as pd
import pysam
from time import perf_counter as pc

def run():
    start_time = pc()

    clusterf = args.cluster
    bamf = str(args.inbam)
    outPrefix = args.outprfx

    print("filtering bam files by preindex at readname...")
    generate_bams(clusterf, bamf, outPrefix)
    end_time = pc()
    print('Used (secs): ', end_time - start_time)

def generate_bams(clusterf, bamf, outPrefix):
    clust_dat = pd.read_csv(clusterf, sep="\t", names = ["barcode", "cluster"]) ### four columns meta: sample.bam, library, barcode, cluster
    clust_dat_dict = pd.Series(clust_dat['cluster'].values, index=clust_dat['barcode']).to_dict()
    outf_dict = dict()
    templateBamF = pysam.AlignmentFile(bamf)
    for cls_id in pd.unique(clust_dat["cluster"]):
        outsamfname = "".join((outPrefix, "_", str(cls_id), ".sam"))
        opysam = pysam.AlignmentFile(outsamfname, "wh", template=templateBamF)
        outf_dict[cls_id] = opysam
    generate_bam_worker(clust_dat_dict, outf_dict, bamf, outPrefix)

def generate_bam_worker(clust_dat_dict, outf_dict, bamf, outPrefix):
    print("extract reads from", bamf)
    bamF = pysam.AlignmentFile(bamf, "rb")
    for b in bamF.fetch(until_eof=True):
        try:
            preindex = b.query_name.split(':')[-1]
        except KeyError:
            continue
        if preindex in clust_dat_dict:
            cls = clust_dat_dict[preindex]
            osam = outf_dict[cls]
            osam.write(b)
        else:
            continue
    bamF.close()

if __name__ == "__main__":
    """filter bam based on QNAMES"""
    run()


