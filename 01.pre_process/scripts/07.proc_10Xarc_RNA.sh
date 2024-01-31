#!/bin/bash
#PBS -q hotel
#PBS -N cellranger_RNA
#PBS -l nodes=1:ppn=16
#PBS -l walltime=72:00:00

### need to submit: start, end and meta
###############################
current="/projects/ps-renlab2/y2xie/projects/BICAN/exp/97.IGM_sequencing_231109"
fastq_dir="${current}/01.rawdata/"
# meta_dir="${current}/01.rawdata/metadata/"
map_dir="${current}/03.mapping/"
mtx_dir="${current}/04.matrices/"
script_dir="${current}/scripts/"


cellranger_mm10="/projects/ps-renlab/y2xie/projects/genome_ref/mm10"
cellranger_hg38="/projects/ps-renlab/y2xie/projects/genome_ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0"
cellranger_mix="/projects/ps-renlab/y2xie/projects/genome_ref/GRCh38_and_mm10"
################################

cd ${current}
mkdir -p ${map_dir} ${mtx_dir}

### create the sub-meta for running
cd ${fastq_dir}
# awk -v start=$start -v end=$end 'NR>=start && NR<end {print}' ${meta_dir}/${meta} > ${script_dir}/Multiome_RNA_ref_${start}_${end}.txt

while read s genome c
do
	if [[ "${genome}" == "mm10" ]]
	then 
		cellranger_ref=${cellranger_mm10}
	elif [[ "${genome}" == "hg38" ]]
	then
		cellranger_ref=${cellranger_hg38}
	elif [[ "${genome}" == "mix" ]]
        then
                cellranger_ref=${cellranger_mix}
	fi			
	/projects/ps-renlab/y2xie/packages/cellranger-6.1.2/cellranger count --id=${s} --transcriptome=${cellranger_ref} --fastqs=${fastq_dir} --sample=${s} --include-introns --chemistry=ARC-v1
done <${script_dir}/Multiome_RNA_ref.txt

/projects/ps-renlab2/pbs/qstat.sh
