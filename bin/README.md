## Scripts
This directory contains in-house scripts and functions to process and analyze Droplet Paired-Tag data.

### Pre-processing pipeline
| File | Description |
|---|---|
| `07.proc_10Xarc_RNA.sh` | Run `cellranger count` on RNA libraries |
| `08.proc_10Xarc_DNA.sh` | Process DNA libraries: trim R2, `cellranger-atac count`, deduplicate, bigWig, MACS2 peaks, and FRiP |
| `08.proc_10X_preIndex_DNA.sh` | Same as above, for pre-indexed DNA libraries |
| `08.proc_10X_preIndex_DNA_v2.sh` | Same as above, but splits each library by pre-index before making bigWig files and calling peaks |
| `scifi.CB_to_BB.py` | Copy the cell barcode to the `BB` tag of a BAM file (libraries without pre-index) |
| `scifi.preindex_CB_to_BB.py` | Add pre-index + cell barcode to the `BB` tag of a BAM file |
| `extract10X_from_bam_with_paired_meta.py` | Split a BAM file by pre-index into one file per target |
| `preindex_Target.xls` | Pre-index sequence and its target (H3K27ac or H3K27me3) |
| `10XcountFrag.py` | Count total and unique fragments per barcode from a fragment file |
| `getSize.py` | Print mean and median fragment size of a BAM file (used as MACS2 `--extsize` for broad peaks) |
| `upstoools/` | C++ tool for pre-index extraction and fastq trimming (see below) |

### Matrices and RNA processing
| File | Description |
|---|---|
| `run_genematrix.py` | Make gene and TSS matrices from a SnapATAC2 `.h5ad` file |
| `run_peakmatrix.py` | Make a peak matrix from a SnapATAC2 `.h5ad` file and a peak BED file |
| `run_tilematrix.py` | Make a tile (bin) matrix at a given resolution from a SnapATAC2 `.h5ad` file |
| `BICAN.ProcessRNA_hicat.R` | Make a Seurat object from RNA cells that pass DNA and RNA filtering |
| `BICAN.ProcessRNA_multiplex_hicat.R` | Same as above, for pre-indexed (multiplexed) libraries |
| `seurat_find_doublet.R` | Score doublets in a Seurat object with a stricter cutoff |

### Downstream analysis
| File | Description |
|---|---|
| `chromhmm_signalize.py` | Bin bigWig or bedGraph signal into ChromHMM `BinarizeSignal` input |
| `convert_value.py` | Binarize DNA methylation signal for ChromHMM |
| `filter_bigwig.sh` | Keep main chromosomes in a bigWig file |
| `fast_diff.py` | SnapATAC2 differential peak test with a downsampled, stratified background |
| `JSS_cts.py` | Call cell-type-specific peaks by Jensen-Shannon specificity (JSS) |
| `JSS_cts_fdr.py` | Same as above, with empirical p-values and FDR |
| `region_lisi_asw.R` | Calculate LISI and silhouette width of brain regions in an embedding |
| `run_permanova.R` | PERMANOVA of brain regions in an embedding from an `.h5ad` file |
| `run_DM.R` | Diffusion map and Palantir pseudotime |
| `region_edgeR.R` | edgeR test of each cortical region vs. all other regions, on pseudoreplicate pseudobulk |
| `lifOver2ldsc.sh` | Lift over a BED file from hg38 to hg19 for LDSC |
| `peak2ld.sh` | Make LDSC annotation and LD scores for one chromosome (SLURM) |
| `summ2partition.sh`, `cts_partition.sh` | Run LDSC cell-type-specific partitioned heritability (SLURM) |

### R helper functions (`R_helper/`)
| File | Description |
|---|---|
| `basics.R` | Color palettes, plot themes, and general functions |
| `DPT_help.R` | Droplet Paired-Tag functions: FRiP and barcode filtering, RNA and DNA clustering, CPM/RPKM |
| `BICAN_help.R` | RNA processing for this project: SoupX, doublet removal, clustering |
| `spatial_help.R` | MERFISH label transfer, spatial prediction, and diffusion maps |
| `hicat_vg.R` | Variable gene selection with Brennecke's method (`find_vg`), from scrattch.hicat |
| `hicat_doubletfinder_mo.R` | Doublet detection with artificial nearest neighbors (`doubletfinder`), adapted from DoubletFinder |

### upstoools
To handle the 1st round barcode (pre-index) in the cellranger-atac pipeline, we extract the barcode from the start of the R2 fastq sequence and append it to the read name, which is kept during alignment. After alignment, the pre-index is added back to the cell barcode in a new BAM tag, `BB`. We use this modified barcode to identify explicit collisions and for downstream analysis.

To build `upstools`, run:
```
cd upstoools
sh make.sh
```
