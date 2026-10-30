#!/bin/python

import argparse

parser = argparse.ArgumentParser(description='counts fragments and unique fragments associacted with each barcode name')
parser.add_argument('--input', type=str, dest="iput", help='input fragment files.')
parser.add_argument('--output', type=str, dest="oput", help='output name')
parser.add_argument('--barcode_field', type=int, dest="bf", default=4, help='nth field containing barcode. Default 4.')
parser.add_argument('--count_field', type=int, dest="cf", default=5, help='th field containing fragments counts. Default 5.')

args = parser.parse_args()

import os
import numpy as np
import pandas as pd
from time import perf_counter as pc
import gzip

def run():
    """ Run standard NMF on rank """
    start_time = pc()
    """ init input files """
    iiput = str(args.iput)
    ooput = str(args.oput)
    bf = args.bf
    cf = args.cf
    print("count contacts per cell...")
    countFrag(iiput, ooput, bf, cf)
    end_time = pc()
    print('Used (secs): ', end_time - start_time)

### chr1    10007   10025   GCTTGCGTCCGCAATC-1      1
def countFrag(iiput, ooput, bf, cf):
    cid = {}
    with oppf(iiput, 'r') as infile:
        for line in infile:
            try:
                dline = line.decode('utf8')
            except AttributeError:
                dline = line
                pass
            if dline[:1] == "#":
                continue
            narr = dline.strip("\n").split("\t")
            cell1 = narr[bf - 1]
            count = narr[cf - 1]
            if cell1 in cid:
                cid[cell1]['total'] += int(count)
                cid[cell1]['unique'] += 1
            else: ### initialize
                cid[cell1] = {'total': int(count), 'unique': 1}
    cidd = pd.DataFrame.from_dict(cid)
    cidd.T.to_csv(ooput, sep='\t')

def oppf(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode) 
    else:
        return open(filename, mode)

if __name__ == "__main__":
    run()


