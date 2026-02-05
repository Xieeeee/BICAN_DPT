
suppressPackageStartupMessages(library(dplyr))
library(ggplot2)
suppressPackageStartupMessages(library(Seurat))
source("/projects/ps-renlab2/y2xie/scripts/basics.R")
source("/projects/ps-renlab2/y2xie/scripts/DPT_help.R")
source("/projects/ps-renlab2/y2xie/scripts/hicat_vg.R")
source("/projects/ps-renlab2/y2xie/scripts/hicat_doubletfinder_mo.R")


args <- commandArgs(trailingOnly = TRUE)
rpath <- args[1] ### path to 10X generated RNA matrix
valid1 <- args[2] ### path to valid cells dataframe, containing colnames "rna" and "atac" indicating valid barcodes
valid2 <- args[3] ### path to FRiP PF cells dataframe, default first column is the DNA valid barcode
prefix <- args[4] ### saved RDS prefix

valid <- read.table(valid1, header = T)
dna_qc <- read.table(valid2, header = F)
valid <- valid[valid$atac %in% dna_qc$V1, ]
cat(paste0("Number of cells pass dual modalities filtering: ", nrow(valid)))
if(nrow(valid) == 0){
	stop("No cells pass filtering. Check the valid barcode files supplied!")
}

print("Start preprocessing...")
mmg <- ProcessRNA_hicat(rpath, valid)
saveRDS(mmg, paste0(prefix, ".rds"))