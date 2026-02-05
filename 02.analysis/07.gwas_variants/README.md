## GWAS enrichment analysis with LDSC

### Prepare background peak sets
```
ipath="/tscc/projects/ps-renlab2/y2xie/projects/BICAN/analysis/14.ldsc"
mkdir ${ipath}/background;
for f in ${ipath}/peaks/AST/*txt; do state=`basename $f .txt`
cat ${ipath}/peaks/*/${state}.txt | sortBed | mergeBed > ${ipath}/background/${state}.allpeak.bed
done
```

### liftOver files
```
for ct in $(ls ${ipath}/peaks); do for bed in /tscc/projects/ps-renlab2/y2xie/projects/BICAN/analysis/14.ldsc/peaks/${ct}/*txt; do nohup bash lifOver2ldsc.sh ${bed} & done; done

for bed in ${ipath}/background/*bed; do nohup bash lifOver2ldsc.sh ${bed} & done
```

### Make annotation and calculate LD score (use hg19 files)
```
MAX_JOBS=1000
for ct in $(ls ../peaks)
do
for bed in /tscc/projects/ps-renlab2/y2xie/projects/BICAN/analysis/14.ldsc/peaks/${ct}/*txt.hg19
do
bname=`basename $bed .txt.hg19`
mkdir -p ../anno/${ct}_${bname};
for r in $(seq 1 22)
do while [ $(squeue -u y2xie | wc -l) -ge $((MAX_JOBS + 1)) ]; do echo "Job limit ($MAX_JOBS) reached. Waiting..."; sleep 60   # check every minute
done
sbatch --export=bed=${bed},bname=${bname},r=$r,outdir=/tscc/projects/ps-renlab2/y2xie/projects/BICAN/analysis/14.ldsc/anno/${ct}_${bname} ../peak2ld.sh
done; done; done
```

### Also process background peak
```
for bed in ${ipath}/background/*bed.hg19; 
do
bname=`basename $bed .bed.hg19`; mkdir -p ${ipath}/anno/${bname};
for r in $(seq 1 22); do while [ $(squeue -u y2xie | wc -l) -ge $((MAX_JOBS + 1)) ]; do echo "Job limit ($MAX_JOBS) reached. Waiting..."; sleep 60 
done;
sbatch --export=bed=${bed},bname=${bname},r=$r,outdir=/tscc/projects/ps-renlab2/y2xie/projects/BICAN/analysis/14.ldsc/anno/${bname} ${ipath}/peak2ld.sh;
done; done
```

### prepare ldcts file
```
### ldcts format: 

${cts_name}.ldcts:
Astrocyte       Cahoy_1000Gv3_ldscores/Cahoy.1.,Cahoy_1000Gv3_ldscores/Cahoy.control.
Oligodendrocyte Cahoy_1000Gv3_ldscores/Cahoy.2.,Cahoy_1000Gv3_ldscores/Cahoy.control.
Neuron  Cahoy_1000Gv3_ldscores/Cahoy.3.,Cahoy_1000Gv3_ldscores/Cahoy.control.
```

```
for f in ${ipath}/peaks/AST/*txt; do state=`basename $f .txt`; mkdir -p ${ipath}/partition/${state}; for ct in $(ls ${ipath}/peaks/); do echo -e ${ct}"\t"anno/${ct}_${state}/${state}.","anno/${state}.allpeak/${state}.allpeak. >> ${ipath}/partition/${state}.allpeak.ldcts; done; done
```

### run cts partition
```
for f in ${ipath}/peaks/AST/*txt; 
do 
state=`basename $f .txt`; 
for sumstats in /tscc/projects/ps-renlab2/y2xie/projects/genome_ref/GWAStraits/*.sumstats.gz;
do
while [ $(squeue -u y2xie | wc -l) -ge $((MAX_JOBS + 1)) ]; do echo "Job limit ($MAX_JOBS) reached. Waiting..."; sleep 60
done
sbatch --export=state=$state,sumstats=${sumstats},ipath=/tscc/projects/ps-renlab2/y2xie/projects/BICAN/analysis/14.ldsc ../summ2partition.sh
done
done
```

