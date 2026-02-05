## ChromHMM analysis with histone modifications, ATAc and DNA methylation

### link file
```
while read a b; do ln -s ../../78.bw/MiniAtlas_merged_dual_filt_clean.cluster.${a}_H3K27ac.bw ./; echo -e ${a}"\t""H3K27ac""\t"MiniAtlas_merged_dual_filt_clean.cluster.${a}_H3K27ac.bw >> config; ln -s ../../78.bw/MiniAtlas_merged_dual_filt_clean.cluster.${a}_H3K27me3.bw ./; echo -e ${a}"\t""H3K27me3""\t"MiniAtlas_merged_dual_filt_clean.cluster.${a}_H3K27me3.bw >> config; ln -s ../../../ref/hba_ATAC/bigwig/subclass/${b}.bw ./; echo -e ${a}"\t""ATAC""\t"${b}.bw >> config; done <../../06.integration/snATAC_DPT_RNA_250722.pred_subclass_match.simplify.txt
```

### filter bigwig
```
for f in *bw; do fname=`basename $f`; nohup bash ../filter_bigwig.sh ${fname} & done
```

### Prepare input from bigwig for binarizeSignal
```
conda activate py2example
### Try multiple resolution to pick the best. Default is 200
for bsize in {200,500,1000}
do
	mkdir -p ../binned/${bsize}
	nohup python ~/scripts/DPT/chromhmm_signalize.py rebin --binsize ${bsize} --binned_dir ../binned/${bsize} --chromsizes_file /projects/ps-renlab2/y2xie/projects/genome_ref/hg38.main.chrom.sizes config hg38 &
done

for bsize in {200,500,1000}
do
	mkdir -p ../combine/${bsize}
	nohup python ~/scripts/DPT/chromhmm_signalize.py combine --binned_dir ../binned/${bsize} --binsize ${bsize} --combined_dir ../combine/${bsize} config &
done
```

### BinarizeSignal
```
current="/projects/ps-renlab2/y2xie/projects/BICAN/analysis/10.chromhmm"
signaldir=${current}/combine
binneddir=${current}/binned
meta=${bamdir}/cellmarkfiletable
bindir=${current}/binary
outdir=${current}/output

for bsize in {200,500,1000}
do
	mkdir -p ${bindir}/${bsize}
	nohup java -mx12800M -jar /projects/ps-renlab2/y2xie/packages/ChromHMM/ChromHMM.jar BinarizeSignal -gzip ${signaldir}/${bsize} ${bindir}/${bsize} &
done
```

### additionally "binarize" DNA hyperCG score using 0.9 as cutoff
```
for bsize in {200,500,1000}
do
	mkdir -p ${bindir}/${bsize}/DNAme ${signaldir}/${bsize}/DNAme
	python ~/scripts/DPT/chromhmm_signalize.py combine --binned_dir ${binneddir}/${bsize}/DNAme --binsize ${bsize} --combined_dir ${signaldir}/${bsize}/DNAme config
done
```

### convert value to prevent chromhmm error
```
for bsize in {200,500,1000}
do
	mkdir ${bindir}/${bsize}/DNAme
	for file in ${signaldir}/${bsize}/DNAme/*
	do
		fname=`basename $file _signal.txt`
		nohup python convert_value.py ${file} ${bindir}/${bsize}/DNAme/${fname}_binary.txt &
	done
done
```

### concat histone binary files with DNAme
```
for bsize in {200,500,1000}
do
	for f2 in ${bindir}/noDNAme/${bsize}/*.gz; do
		fname=`basename $f2`
		dirname=`dirname $f2`
		f1="${dirname}/DNAme/${fname}"
		out=${bindir}/${bsize}/${fname}

		if [[ ! -f "$f1" ]]; then
		echo "Missing DNAme file for: $f2" >&2
		continue
		fi

		if ! diff -q <(gzip -dc "$f1" | head -n1) <(gzip -dc "$f2" | head -n1) >/dev/null; then
		echo "Header mismatch: $f1 vs $f2 — skipping" >&2
		continue
		fi

		gzip -dc "$f1" | \
		awk -v OFS='\t' 'NR==FNR { if (FNR>1) a[FNR]=$0; next }
		                FNR==1 { print; next }
		                { print $0, a[FNR] }' \
		/dev/fd/0 <(gzip -dc "$f2") | gzip > "$out"

		echo "Wrote: $out"
	done
done
```

### LearnModel
```
for bsize in {200,500,1000}
do
	for r in {3..15}
	do
        mkdir -p ${outdir}/${bsize}_${r}
        nohup java -mx9600M -jar /projects/ps-renlab2/y2xie/packages/ChromHMM/ChromHMM.jar LearnModel -b ${bsize} -noautoopen -gzip -r 1000 -p 16 -printposterior -f ${bindir}/${bsize}/valid.file -l /projects/ps-renlab2/y2xie/projects/genome_ref/hg38.main.chrom.sizes ${bindir}/${bsize} ${outdir}/${bsize}_${r}/ ${r} hg38 &
	done
	sleep 40m
done
```
