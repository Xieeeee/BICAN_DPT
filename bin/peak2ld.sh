#!/bin/bash
#SBATCH -J peak2ld
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
# peak file: bed
# chrom id: r
# outdir: outdir
# bname: need to pass in # =`basename $bed .txt.hg19`

python /tscc/projects/ps-renlab2/y2xie/scripts/git/ldsc/make_annot.py --bed-file ${bed} --bimfile /tscc/projects/ps-renlab2/y2xie/projects/genome_ref/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.${r}.bim --annot-file ${outdir}/${bname}.${r}.annot.gz
python /tscc/projects/ps-renlab2/y2xie/scripts/git/ldsc/ldsc.py --l2 --bfile /tscc/projects/ps-renlab2/y2xie/projects/genome_ref/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.${r} --ld-wind-cm 1 --annot ${outdir}/${bname}.${r}.annot.gz --thin-annot --out ${outdir}/${bname}.${r} --print-snps /tscc/projects/ps-renlab2/y2xie/projects/genome_ref/ldsc/hapmap3_snps/hm.${r}.snp
