## Analysis notebooks
This directory contains notebooks to analyze the data and reproduce figures in our manuscript.

`00.qc`: library QC
- `01.library_QC_summary.ipynb`: TSS enrichment, library complexity, cell recovery, and correlation between samples

`01.clustering`: clustering and annotation of RNA and DNA data with Seurat and SnapATAC2
- `01.RNA_Seurat_clustering.ipynb`: pool RNA libraries, remove doublets, cluster, and annotate from class to cluster
- `01.MapMyCell.ipynb`: map RNA cells to reference taxonomies of matched brain regions with MapMyCells
- `02.DNA_snapatac2_clustering.ipynb`: filter DNA cells, split H3K27ac and H3K27me3, and cluster
- `03.Clustering_summary.ipynb`: figures for UMAPs, annotation hierarchy, marker genes, subclass QC, and integration

`02.integration`: CCA integration of RNA with scRNA-seq references and with snATAC-seq
- `00.siletti_process_1.ipynb`: subset the Siletti et al. 2023 reference to matched dissection regions
- `00.siletti_process_2.ipynb`: process and subset the Siletti reference by cluster label
- `00.smartseq_process_1.ipynb`: process the Allen SMART-seq reference (multiple cortical areas)
- `01.RNA_integration.ipynb`: integrate RNA with the references and calculate overlap scores
- `02.snatac_process.ipynb`: merge snATAC-seq data, subset by subclass, and make gene activity and cCRE matrices
- `03.ATAC_integration.ipynb`: integrate snATAC-seq with RNA, separately for neurons and non-neurons

`03.ccre`: classification and characterization of cCREs by chromatin state
- `01.ccre_chmm_annotation.ipynb`: classify H3K27ac+ and H3K27me3+ cCREs; genomic features, conservation, motifs, VISTA enhancers, and cell-type-specific and epi-variable cCREs

`04.chromhmm`: chromatin states from H3K27ac, H3K27me3, ATAC, and DNA methylation with ChromHMM
- `README.md`: commands to make ChromHMM inputs and train models
- `10.chromhmm.ipynb`: select bin size and state number, annotate states, compare with Roadmap Epigenomics, and compare states between subclasses

`05.abc`: enhancer-promoter links with the activity-by-contact (ABC) model
- `README.md`: link to the ABC model used
- `01.abc_analysis.ipynb`: filter H3K27ac+ ABC links; examples, GO terms, motifs, and comparison with eQTL studies
- `02.abc_enhancer_score.ipynb`: per-cell H3K27ac and H3K27me3 enhancer scores

`06.grn`: active and repressive gene regulatory networks with SCENIC+
- `00.scenic_plus_input_process.ipynb`: make RNA, H3K27ac, and H3K27me3 inputs for SCENIC+
- `01.scenic_plus_output_check.ipynb`: regulon specificity scores by cell type
- `02.GRN_analysis.ipynb`: figures for active and repressive networks, including interneuron (MGE/CGE) networks

`07.gwas_variants`: disease risk variant enrichment and effect prediction
- `README.md`: commands for LDSC partitioned heritability
- `01.ldsc_summary.ipynb`: LDSC enrichment by chromatin state and subclass
- `02.finemapped_ovlp.ipynb`: overlap of fine-mapped variants with ABC links; motifbreakR and deep learning predictions

`08.heterogenity`: structural and regional heterogeneity in RNA expression and epigenome
- `01.structure_heterogenity.ipynb`: variance explained by brain structure; forebrain vs. hindbrain glia
- `02.cortical_heterogenity.ipynb`: regional differences in cortex and anterior-posterior axis analysis

`09.spatial`: MERFISH analysis and spatial imputation
- `01.MERFISH_integration.ipynb`: integrate MERFISH with RNA, assign cortical layers, and predict spatial location of Droplet Paired-Tag cells
- `02.MERFISH_analysis.ipynb`: spatial prediction figures, SST laminar distribution, and differential tests along cortical depth

`10.hic`: annotation of short-range and long-range loops called from snm3C-seq
- `01.loop_annotation.ipynb`: annotate loops with histone peaks; loop strength, anchor gene expression, and support for ABC links
- `02.SLR_annotation.ipynb`: super long-range (SLR) interactions: APA, histone marks at anchors, and anchor genes
- `03.cooler.ipynb`: contact matrices by cell type; call SLRs with Peakachu and FitHiC

`11.conservation`: conserved and divergent TF programs between human and mouse
- `13.cross_species_integration.ipynb`: integrate human and mouse cells and compare active and repressive TF activity
