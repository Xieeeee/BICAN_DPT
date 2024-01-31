#!/bin/bash
#PBS -q hotel
#PBS -N cellranger_arc_DNA_${start}_${end}
#PBS -l nodes=1:ppn=8
#PBS -l walltime=128:00:00

### 231127: add the temporary working directory: submit_dir, since when generating multiple files with same sample names, cellranger will fail.
###############################
current="/projects/ps-renlab2/y2xie/projects/BICAN/exp/95.IGM_sequencing_230124"
fastq_dir="${current}/01.rawdata/"
submit_dir="${fastq_dir}/submit/"
meta_dir="${current}/scripts/"
map_dir="${current}/03.mapping/10X/"
bw_dir="${current}/06.bw/"
macs2_dir="${current}/09.macs2/"
script_dir="${current}/scripts/"
mm10_ref="/projects/ps-renlab/y2xie/projects/genome_ref/refdata-cellranger-arc-mm10-2020-A-2.0.0"
hg38_ref="/projects/ps-renlab/y2xie/projects/genome_ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0"
mix_ref="/projects/ps-renlab2/y2xie/projects/genome_ref/refdata-cellranger-atac-GRCh38-and-mm10-2020-A-2.0.0"
################################

mkdir -p ${submit_dir} ${map_dir} ${bw_dir} ${macs2_dir}

awk -v start=$start -v end=$end 'NR>=start && NR<end {print}' ${meta_dir}/${meta} > ${script_dir}/Multiome_DNA_ref_${start}_${end}.txt

cd ${submit_dir}
### run atac processing

while read s genome type c #name, genome, peak type, note
do
	if [[ "${genome}" == "mm10" ]]
	then 
		cellranger_ref=${mm10_ref}
		bl=/projects/ps-renlab/y2xie/projects/genome_ref/mm10-blacklist.v2.bed
		ref_peak=/projects/ps-renlab/y2xie/projects/genome_ref/gencode.vM25.annotation.gtf
		macs2_genome="mm"
	elif [[ "${genome}" == "hg38" ]]
	then
		cellranger_ref=${hg38_ref}
		bl=/projects/ps-renlab/y2xie/projects/genome_ref/hg38-blacklist.v2.bed
		ref_peak=/projects/ps-renlab/y2xie/projects/genome_ref/gencode.vH35.annotation.gtf
		macs2_genome="hs"
	elif [[ "${genome}" == "mix" ]]
	then
		cellranger_ref=${mix_ref}
		macs2_genome=4.57e9
		bl=/projects/ps-renlab/y2xie/projects/genome_ref/mix-blacklist.bed
		ref_peak=/projects/ps-renlab/y2xie/projects/genome_ref/mix_arc.annotation.gtf
	fi	
	
	### trim if necessary: always to 24bp
	### there are many R2...
	for nseq in ${fastq_dir}/${s}*R2*fastq.gz
	do
		nseq_name=`basename $nseq .fastq.gz`
		ln -s ${nseq} ${nseq_name}.fq.gz
		seqq=$(zcat ${nseq_name}.fq.gz | head -2 | tail -1)
		seq_length=$(expr length $seqq)
		if [[ $seq_length -gt 24 ]]
		then
			/projects/ps-renlab/y2xie/projects/scifi-multiome/upstoools/upstools trimfq ${nseq_name}.fq.gz 1 24 ### ${fastq_dir}/${nseq_name}_trim.fq.gz
			# unlink ${nseq}
			mv ${nseq_name}_trim.fq.gz ${nseq_name}.fastq.gz
			unlink ${nseq_name}.fq.gz
		else
			unlink ${nseq_name}.fq.gz
			ln -s ${nseq} ./
		fi
		### link other files
		for nseq in ${fastq_dir}/${s}*R{1,3}*fastq.gz
		do
			ln -s ${nseq} ./
		done
	done

	### cellranger	
	/projects/ps-renlab/y2xie/packages/cellranger-atac-2.0.0/cellranger-atac count --id=${s} --reference=${cellranger_ref} --fastqs=${submit_dir}/ --sample=${s} --chemistry=ARC-v1
	mv ${submit_dir}/${s} ${current}/
	
	### rmdup
	python /projects/ps-renlab2/y2xie/scripts/scifi/scifi.CB_to_BB.py --in ${current}/${s}/outs/possorted_bam.bam
	java -Xmx8G -XX:ParallelGCThreads=16 -jar /projects/ps-renlab/y2xie/packages/picard.jar MarkDuplicates INPUT=${current}/${s}/outs/possorted_bam.bam.BB.bam TMP_DIR=${map_dir} METRICS_FILE=${map_dir}/${s}_dup.qc VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=TRUE OUTPUT=${map_dir}/${s}_rmdup.bam BARCODE_TAG=BB REMOVE_DUPLICATES=TRUE

	### bamCoverage
	samtools index ${map_dir}/${s}_rmdup.bam
	bamCoverage -b ${map_dir}/${s}_rmdup.bam -o ${bw_dir}/${s}_rmdup.bw -p max --normalizeUsing RPKM -bl ${bl}
	computeMatrix reference-point --referencePoint TSS -b 2000 -a 2000 -R ${ref_peak} -S ${bw_dir}/${s}_rmdup.bw --skipZeros -o ${bw_dir}/${s}_rmdup.mtx.gz -p max
	plotProfile -m ${bw_dir}/${s}_rmdup.mtx.gz --refPointLabel "TSS" -o ${bw_dir}/${s}_rmdup_profile.pdf --outFileNameData ${bw_dir}/${s}_rmdup_profile.txt
	plotHeatmap -m ${bw_dir}/${s}_rmdup.mtx.gz --refPointLabel "TSS" -o ${bw_dir}/${s}_rmdup_heatmap.pdf --outFileNameMatrix {bw_dir}/${s}_rmdup_heatmap.mtx.gz

	### peak calling
	if [[ "${type}" == "narrow" ]]
	then
		macs2 callpeak -t ${map_dir}/${s}_rmdup.bam -g ${macs2_genome} -n ${s} -f BAMPE --outdir ${macs2_dir} -q 0.05
	elif [[ "${type}" == "broad" ]]
	then
		### estimated extend size
		size=$(/projects/ps-renlab/y2xie/anaconda3/bin/python /projects/ps-renlab/y2xie/scripts/random/getSize.py --bam ${map_dir}/${s}_rmdup.bam)
		size=${size##* } ### median size
		macs2 callpeak -t ${map_dir}/${s}_rmdup.bam -g ${macs2_genome} -n ${s} -f BAMPE --outdir ${macs2_dir} -q 0.05 --nomodel --extsize ${size} --nolambda --broad-cutoff 0.1 --broad
	fi

	### FRiP
	zcat ${current}/${s}/outs/fragments.tsv.gz | egrep -v '#' | bedtools intersect -wa -u -a stdin -b ${macs2_dir}/${s}_peaks.${type}Peak > ${macs2_dir}/${s}_FRiP_fragments.tsv
	python /projects/ps-renlab2/y2xie/scripts/scifi/10XcountFrag.py --input ${current}/${s}/outs/fragments.tsv.gz --output ${macs2_dir}/${s}_clean_fragments.tsv.gz_Count.xls
	python /projects/ps-renlab2/y2xie/scripts/scifi/10XcountFrag.py --input ${macs2_dir}/${s}_FRiP_fragments.tsv --output ${macs2_dir}/${s}_FRiP_fragments.tsv_Count.xls
	rm ${macs2_dir}/${s}_peaks.${type}Peak.tmp ${macs2_dir}/${s}_FRiP_fragments.tsv

done <${script_dir}/Multiome_DNA_ref_${start}_${end}.txt

/projects/ps-renlab2/pbs/qstat.sh
