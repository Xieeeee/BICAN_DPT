#!/bin/bash
#SBATCH -J cts_partition
#SBATCH -p hotel
#SBATCH -q hotel
#SBATCH -N 1
#SBATCH -A htl195
#SBATCH -n 4
#SBATCH -c 2
#SBATCH --mem=4G
#SBATCH -t 01:00:00

source /tscc/nfs/home/y2xie/anaconda3/etc/profile.d/conda.sh
conda activate ldsc

sname=`basename $sumstats .sumstats.gz`;

cd /tscc//projects/ps-renlab2/y2xie/projects/BICAN/BICAN_BGC/ldsc
python /tscc//projects/ps-renlab2/y2xie/scripts/git/ldsc/ldsc.py --h2-cts ${sumstats} --ref-ld-chr /tscc/projects/ps-renlab2/y2xie/projects/genome_ref/ldsc/1000G_EUR_Phase3_baseline/baseline. --out partition/${cts_name}/${sname} --ref-ld-chr-cts ${cts_name}.allpeak.ldcts --w-ld-chr /tscc/projects/ps-renlab2/y2xie/projects/genome_ref/ldsc/weights_hm3_no_hla/weights. 
