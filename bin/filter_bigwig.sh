### retain only main scaffold in bigwig

bigwig=$1
bname=`basename $bigwig .bw`;
bigWigToBedGraph $bigwig ${bname}.bedGraph
grep -E 'chr[1-9]|chr1[0-9]|chr2[0-2]|chrX|chrY' ${bname}.bedGraph > ${bname}.tmp.bedGraph
bedGraphToBigWig ${bname}.tmp.bedGraph /projects/ps-renlab2/y2xie/projects/genome_ref/hg38.main.chrom.sizes ${bigwig}.main.bw
rm ${bname}.bedGraph ${bname}.tmp.bedGraph

