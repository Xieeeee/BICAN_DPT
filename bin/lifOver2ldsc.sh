bed=$1
liftOver ${bed} /tscc/projects/ps-renlab2/y2xie/ps-renlab/y2xie/projects/genome_ref/liftOver/hg38ToHg19.over.chain.gz ${bed}.hg19 ${bed}.unMapped -minMatch=0.95
liftOver ${bed}.hg19 /tscc/projects/ps-renlab2/y2xie/ps-renlab/y2xie/projects/genome_ref/liftOver/hg19ToHg38.over.chain ${bed}.hg38 ${bed}.unMapped -minMatch=0.95
bedtools intersect -u -a ${bed} -b ${bed}.hg38 > ${bed}.reciprocal
liftOver ${bed}.reciprocal /tscc/projects/ps-renlab2/y2xie/ps-renlab/y2xie/projects/genome_ref/liftOver/hg38ToHg19.over.chain.gz ${bed}.hg19.reciprocal ${bed}.unMapped -minMatch=0.95

awk '($3-$2) <= 1000' ${bed}.hg19.reciprocal > ${bed}.hg19