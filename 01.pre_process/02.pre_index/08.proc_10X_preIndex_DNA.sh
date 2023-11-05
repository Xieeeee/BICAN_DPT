#!/bin/bash
#PBS -q hotel
#PBS -N cellranger_arc_DNA_preindex
#PBS -l nodes=1:ppn=4
#PBS -l walltime=72:00:00

### This is for pre-indexed Droplet Paired-Tag

###############################
# bcl2="230705_VH00454_179_AACVFMWM5"
# current="/projects/ps-renlab2/y2xie/projects/DPT/45.human_brain_DPT_230705"
bcl2_dir="${current}/${bcl2}/"
fastq_dir="${current}/01.rawdata/"
map_dir="${current}/03.mapping/10X/"
bw_dir="${current}/06.bw/"
macs2_dir="${current}/09.macs2/"
script_dir="${current}/scripts/"
mm10_ref="/projects/ps-renlab/y2xie/projects/genome_ref/refdata-cellranger-arc-mm10-2020-A-2.0.0"
hg38_ref="/projects/ps-renlab/y2xie/projects/genome_ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0"
mix_ref="/projects/ps-renlab2/y2xie/projects/genome_ref/refdata-cellranger-atac-GRCh38-and-mm10-2020-A-2.0.0"
################################
mkdir -p ${fastq_dir}/preindex ${map_dir} ${bw_dir} ${macs2_dir}

### demultiplex
cd ${current}
# cellranger-arc mkfastq --run ${bcl2_dir} --csv ${script_dir}/SampleSheet.atac.csv --output-dir ${fastq_dir} #--use-bases-mask Y*,I8,Y*,Y*

### run atac processing
cd ${fastq_dir}/
while read s genome type c #name, genome, peak type, note
do
	if [[ "${genome}" == "mm10" ]]
	then 
		cellranger_ref=${mm10_ref}
		bl=/projects/ps-renlab/y2xie/projects/genome_ref/mm10-blacklist.v2.bed
		ref_peak=/projects/ps-renlab/y2xie/projects/genome_ref/gencode.vM25.annotation.gtf
		macs2_genome="mm"
		use_chrom="(?i)^chr"
	elif [[ "${genome}" == "hg38" ]]
	then
		cellranger_ref=${hg38_ref}
		bl=/projects/ps-renlab/y2xie/projects/genome_ref/hg38-blacklist.v2.bed
		ref_peak=/projects/ps-renlab/y2xie/projects/genome_ref/gencode.vH35.annotation.gtf
		macs2_genome="hs"
		use_chrom="(?i)^chr"
	elif [[ "${genome}" == "mix" ]]
	then
		cellranger_ref=${mix_ref}
		macs2_genome=4.57e9
		bl=/projects/ps-renlab/y2xie/projects/genome_ref/mix-blacklist.bed
		ref_peak=/projects/ps-renlab/y2xie/projects/genome_ref/mix_arc.annotation.gtf
		use_chrom="^(mm10_chr|GRCh38_chr)"
	fi	
	
	### get preindex
	for f in ${fastq_dir}/${bcl2##*_}/${s}/${s}*L001*_001.fastq.gz; do fname=`basename $f _001.fastq.gz`; cat ${fastq_dir}/${bcl2##*_}/${s}/${fname%%_*}*${fname##*_}*.fastq.gz > ${fastq_dir}/${fname%%_*}_${fname##*_}.fq.gz; done
	/projects/ps-renlab/y2xie/projects/scifi-multiome/upstoools/upstools sepType_DPT ${fastq_dir}/${s} ${nlen} ###output: ${fastq_dir}/${s}_??_prefix.fq.gz
	cd ${fastq_dir}/preindex
	for f in ${fastq_dir}/${s}*_prefix.fq.gz; do fname=`basename $f _prefix.fq.gz`; ln -s $f ${fname%%_*}_S1_${fname##*_}_001.fastq.gz; done
	cd ${fastq_dir}/

	## cellranger	
	cellranger-atac count --id=${s} --project=${bcl2##*_} --reference=${cellranger_ref} --fastqs=${fastq_dir}/preindex --sample=${s} --chemistry=ARC-v1
	mv ${fastq_dir}/${s} ${current}/
	
	## rmdup, considering pre-index barcode
	python /projects/ps-renlab2/y2xie/scripts/scifi/scifi.preindex_CB_to_BB.py --in ${current}/${s}/outs/possorted_bam.bam
	java -Xmx8G -XX:ParallelGCThreads=16 -jar /projects/ps-renlab/y2xie/packages/picard.jar MarkDuplicates INPUT=${current}/${s}/outs/possorted_bam.bam.BB.bam TMP_DIR=${map_dir} METRICS_FILE=${map_dir}/${s}_dup.qc VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=TRUE OUTPUT=${map_dir}/${s}_rmdup.bam BARCODE_TAG=BB REMOVE_DUPLICATES=TRUE

	## Also need to re-generate fragments file
	## generate fragment file (CB field)
	mv ${current}/${s}/outs/fragments.tsv.gz ${current}/${s}/outs/Undemultiplex_fragments.tsv.gz

	## snapatac2 fails, version issue
	/projects/ps-renlab/y2xie/anaconda3/envs/snapatac2/bin/python /projects/ps-renlab2/y2xie/scripts/scifi/snapatac_pp_fragment.py --input ${current}/${s}/outs/possorted_bam.bam.BB.bam --output ${current}/${s}/outs/fragments.tsv --field BB
	samtools index ${current}/${s}/outs/possorted_bam.bam.BB.bam
	/projects/ps-renlab/y2xie/anaconda3/envs/seurat/bin/sinto fragments -b ${current}/${s}/outs/possorted_bam.bam.BB.bam -f ${current}/${s}/outs/fragments.tsv -t BB --use_chrom ${use_chrom}

	### sort fragments
	sort -k1,1 -k2,2n ${current}/${s}/outs/fragments.tsv > ${current}/${s}/outs/fragments.tsv.tmp
	mv ${current}/${s}/outs/fragments.tsv.tmp ${current}/${s}/outs/fragments.tsv
	bgzip -f ${current}/${s}/outs/fragments.tsv
	tabix -p bed ${current}/${s}/outs/fragments.tsv.gz

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
	zcat ${current}/${s}/outs/fragments.tsv.gz |  egrep -v '#' | bedtools intersect -wa -u -a stdin -b ${macs2_dir}/${s}_peaks.${type}Peak > ${macs2_dir}/${s}_FRiP_fragments.tsv
	python /projects/ps-renlab2/y2xie/scripts/scifi/10XcountFrag.py --input ${current}/${s}/outs/fragments.tsv.gz --output ${macs2_dir}/${s}_clean_fragments.tsv.gz_Count.xls
	python /projects/ps-renlab2/y2xie/scripts/scifi/10XcountFrag.py --input ${macs2_dir}/${s}_FRiP_fragments.tsv --output ${macs2_dir}/${s}_FRiP_fragments.tsv_Count.xls
	rm ${macs2_dir}/${s}_peaks.${type}Peak.tmp ${macs2_dir}/${s}_FRiP_fragments.tsv
done <${script_dir}/Multiome_DNA_ref.txt

/projects/ps-renlab2/pbs/qstat.sh
