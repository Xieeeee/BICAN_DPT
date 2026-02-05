
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
dna <- args[3] ### DNA library id.
valid2 <- args[4] ### path to FRiP PF cells dataframe. will search for all files associated with the DNA library id
prefix <- args[5] ### saved RDS prefix

valid <- read.table(valid1, header = T)
files <- list.files(valid2, pattern = paste0(dna, ".*_PF_cells.txt"))
dna_qc <- list()
i = 1
for (file in files){
    dna_qc[[i]] <- read.table(paste0(valid2, "/", file), header = F)
    i = i + 1
}
dna_qc <- do.call(rbind, dna_qc)

valid <- valid[valid$atac_bc %in% dna_qc$V1, ] %>% distinct
cat(paste0("Number of cells pass dual modalities filtering: ", nrow(valid)))
if(nrow(valid) == 0){
	stop("No cells pass filtering. Check the valid barcode files supplied!")
}

print("Start preprocessing...")
mmg <- ProcessRNA_hicat(rpath, valid)
saveRDS(mmg, paste0(prefix, ".rds"))
