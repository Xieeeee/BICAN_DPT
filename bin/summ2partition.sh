#!/bin/bash
#SBATCH -J partition
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -N 1
#SBATCH -A csd788
#SBATCH -n 4
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 02:00:00

source /tscc/nfs/home/y2xie/anaconda3/etc/profile.d/conda.sh
conda activate ldsc

### pass in:
# state: state (E1 - E9)
# sumstats: sumstats.gz file
# ipath: working directory location
cd ${ipath}
sname=`basename $sumstats .sumstats.gz`

python /tscc/projects/ps-renlab2/y2xie/scripts/git/ldsc/ldsc.py --h2-cts ${sumstats} --ref-ld-chr /tscc/projects/ps-renlab2/y2xie/projects/genome_ref/ldsc/1000G_EUR_Phase3_baseline/baseline. --out partition/${state}/${sname} --ref-ld-chr-cts partition/${state}.allpeak.ldcts --w-ld-chr /tscc/projects/ps-renlab2/y2xie/projects/genome_ref/ldsc/weights_hm3_no_hla/weights. 
