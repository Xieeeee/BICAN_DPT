#!/bin/bash
#PBS -q hotel
#PBS -N cellranger_arc_DNA_preindex_${nlen}_${start}_${end}
#PBS -l nodes=1:ppn=8
#PBS -l walltime=128:00:00

### This is for pre-indexed Droplet Paired-Tag
### v2: split pre-index before generate bigwig and calling peaks 
### pass: start, end, nlen, meta, pair_meta
### pair_meta: two columns file, with first column being preindex and second column being Target

###############################
current="/projects/ps-renlab2/y2xie/projects/BICAN/exp/97.IGM_sequencing_231109"
fastq_dir="${current}/01.rawdata"
submit_dir="${fastq_dir}/submit/"
meta_dir="${current}/scripts"
map_dir="${current}/03.mapping/10X/"
bw_dir="${current}/06.bw/"
macs2_dir="${current}/09.macs2/"
script_dir="${current}/scripts/"

pair_meta="/projects/ps-renlab2/y2xie/scripts/DPT/preindex_Target.xls"
mm10_ref="/projects/ps-renlab/y2xie/projects/genome_ref/refdata-cellranger-arc-mm10-2020-A-2.0.0"
hg38_ref="/projects/ps-renlab/y2xie/projects/genome_ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0"
mix_ref="/projects/ps-renlab2/y2xie/projects/genome_ref/refdata-cellranger-atac-GRCh38-and-mm10-2020-A-2.0.0"
################################

mkdir -p ${submit_dir} ${fastq_dir}/preindex ${map_dir} ${bw_dir} ${macs2_dir}

awk -v start=$start -v end=$end 'NR>=start && NR<end {print}' ${meta_dir}/${meta} > ${script_dir}/Multiome_DNA_ref_${nlen}_${start}_${end}.txt

### run atac processing
cd ${fastq_dir}
while read s genome c #name, genome, peak type, note
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
	
	## merge fastq
	## fix error 231206: there is not always L001 there
	for f in ${fastq_dir}/${s}*L00?*_001.fastq.gz
	do 
		fname=`basename $f _001.fastq.gz`
		### files checking
		if [[ ! -f "${fastq_dir}/${fname%%_*}_${fname##*_}.fq.gz" ]]
		then
			cat ${fastq_dir}/${fname%%_*}*${fname##*_}*.fastq.gz > ${fastq_dir}/${fname%%_*}_${fname##*_}.fq.gz
		fi
	done

	/projects/ps-renlab/y2xie/projects/scifi-multiome/upstoools/upstools sepType_DPT ${fastq_dir}/${s} ${nlen} ###output: ${fastq_dir}/${s}_??_prefix.fq.gz

	### trim R2 if necessary
	### if length still > 24bp, trim the last n bp
	cd ${fastq_dir}/preindex
	seqq=$(zcat ${fastq_dir}/${s}_R2_prefix.fq.gz | head -2 | tail -1)
	seq_length=$(expr length $seqq)
	if [[ $seq_length -gt 24 ]]
	then
		/projects/ps-renlab/y2xie/projects/scifi-multiome/upstoools/upstools trimfq ${fastq_dir}/${s}_R2_prefix.fq.gz 1 24 ### ${fastq_dir}/${s}_R2_prefix_trim.fq.gz
		ln -nsf ${fastq_dir}/${s}_R2_prefix_trim.fq.gz ${s}_S1_R2_001.fastq.gz
	fi

	### link other files
	for f in ${fastq_dir}/${s}*_prefix.fq.gz; do fname=`basename $f _prefix.fq.gz`; ln -s $f ${fname%%_*}_S1_${fname##*_}_001.fastq.gz; done

	## cellranger	
	cd ${fastq_dir}/
	cellranger-atac count --id=${s} --reference=${cellranger_ref} --fastqs=${fastq_dir}/preindex --sample=${s} --chemistry=ARC-v1
	mv ${fastq_dir}/${s} ${current}/
	
	## rmdup, considering pre-index barcode
	python /projects/ps-renlab2/y2xie/scripts/scifi/scifi.preindex_CB_to_BB.py --in ${current}/${s}/outs/possorted_bam.bam
	java -Xmx8G -XX:ParallelGCThreads=16 -jar /projects/ps-renlab/y2xie/packages/picard.jar MarkDuplicates INPUT=${current}/${s}/outs/possorted_bam.bam.BB.bam TMP_DIR=${map_dir} METRICS_FILE=${map_dir}/${s}_dup.qc VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=TRUE OUTPUT=${map_dir}/${s}_rmdup.bam BARCODE_TAG=BB REMOVE_DUPLICATES=TRUE

	## Also need to re-generate fragments file
	## generate fragment file (CB field)
	mv ${current}/${s}/outs/fragments.tsv.gz ${current}/${s}/outs/Undemultiplex_fragments.tsv.gz

	## snapatac2 fails, version issue
	# /projects/ps-renlab/y2xie/anaconda3/envs/snapatac2/bin/python /projects/ps-renlab2/y2xie/scripts/scifi/snapatac_pp_fragment.py --input ${current}/${s}/outs/possorted_bam.bam.BB.bam --output ${current}/${s}/outs/fragments.tsv --field BB
	samtools index ${current}/${s}/outs/possorted_bam.bam.BB.bam
	/projects/ps-renlab/y2xie/anaconda3/envs/seurat/bin/sinto fragments -b ${current}/${s}/outs/possorted_bam.bam.BB.bam -f ${current}/${s}/outs/fragments.tsv -t BB --use_chrom ${use_chrom}

	### sort fragments
	sort -k1,1 -k2,2n ${current}/${s}/outs/fragments.tsv > ${current}/${s}/outs/fragments.tsv.tmp
	mv ${current}/${s}/outs/fragments.tsv.tmp ${current}/${s}/outs/fragments.tsv
	bgzip -f ${current}/${s}/outs/fragments.tsv
	tabix -p bed ${current}/${s}/outs/fragments.tsv.gz
	python /projects/ps-renlab2/y2xie/scripts/scifi/10XcountFrag.py --input ${current}/${s}/outs/fragments.tsv.gz --output ${macs2_dir}/${s}_clean_fragments.tsv.gz_Count.xls

	### in v2: split bam files by preindex table!
	/projects/ps-renlab/y2xie/anaconda3/envs/seurat/bin/python /projects/ps-renlab2/y2xie/scripts/DPT/extract10X_from_bam_with_paired_meta.py --cluster ${pair_meta} --inbam ${map_dir}/${s}_rmdup.bam --outprfx ${map_dir}/${s}

	for target in ${map_dir}/${s}_*.sam;
	do
		targetf=`basename $target .sam`
		samtools sort -o ${map_dir}/${targetf}.bam ${target}
		samtools index ${map_dir}/${targetf}.bam
		bamCoverage -b ${map_dir}/${targetf}.bam -o ${bw_dir}/${targetf}.bw -p max --normalizeUsing RPKM -bl ${bl}

		### peak calling
		if [[ "${targetf}" == *"H3K27ac"* ]]
		then
			macs2 callpeak -t ${map_dir}/${targetf}.bam -g ${macs2_genome} -n ${targetf} -f BAMPE --outdir ${macs2_dir} -q 0.05
			type="narrow"
		else
			### estimated extend size
			size=$(/projects/ps-renlab/y2xie/anaconda3/bin/python /projects/ps-renlab/y2xie/scripts/random/getSize.py --bam ${map_dir}/${targetf}.bam)
			size=${size##* } ### median size
			macs2 callpeak -t ${map_dir}/${targetf}.bam -g ${macs2_genome} -n ${targetf} -f BAMPE --outdir ${macs2_dir} -q 0.05 --nomodel --extsize ${size} --nolambda --broad-cutoff 0.1 --broad
			type="broad"
		fi

		### FRiP
		zcat ${current}/${s}/outs/fragments.tsv.gz |  egrep -v '#' | bedtools intersect -wa -u -a stdin -b ${macs2_dir}/${targetf}_peaks.${type}Peak > ${macs2_dir}/${targetf}_FRiP_fragments.tsv
		python /projects/ps-renlab2/y2xie/scripts/scifi/10XcountFrag.py --input ${macs2_dir}/${targetf}_FRiP_fragments.tsv --output ${macs2_dir}/${targetf}_FRiP_fragments.tsv_Count.xls
		rm ${macs2_dir}/${targetf}_peaks.${type}Peak.tmp ${macs2_dir}/${targetf}_FRiP_fragments.tsv
	done
	
	### summarize all target files from the same library
	computeMatrix reference-point --referencePoint TSS -b 2000 -a 2000 -R ${ref_peak} -S ${bw_dir}/${s}*.bw --skipZeros -o ${bw_dir}/${s}_rmdup.mtx.gz -p max
	plotProfile -m ${bw_dir}/${s}_rmdup.mtx.gz --refPointLabel "TSS" -o ${bw_dir}/${s}_rmdup_profile.pdf --outFileNameData ${bw_dir}/${s}_rmdup_profile.txt
	plotHeatmap -m ${bw_dir}/${s}_rmdup.mtx.gz --refPointLabel "TSS" -o ${bw_dir}/${s}_rmdup_heatmap.pdf --outFileNameMatrix {bw_dir}/${s}_rmdup_heatmap.mtx.gz

done <${script_dir}/Multiome_DNA_ref_${nlen}_tmp.txt  #Multiome_DNA_ref_${nlen}_${start}_${end}.txt

/projects/ps-renlab2/pbs/qstat.sh
