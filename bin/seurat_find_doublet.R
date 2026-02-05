### re-identify doublets with stricter cutoff ###

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(Seurat))
library(ggplot2)
source("/projects/ps-renlab2/y2xie/scripts/basics.R")
source("/projects/ps-renlab2/y2xie/scripts/DPT_help.R")
source("/projects/ps-renlab2/y2xie/scripts/hicat_vg.R")
source("/projects/ps-renlab2/y2xie/scripts/hicat_doubletfinder_mo.R")

args <- commandArgs(trailingOnly = TRUE)
obj <- args[1] ### path to processed RNA seurat object

cat(paste0("file being processed: ", obj, "\n"))
mmg <- readRDS(obj)
doublet <- doubletfinder(mmg, VariableFeatures(mmg), plot = F, proportion.artificial = 0.2)
cat(paste0("calculated doublet rate: ", 100*(table(doublet[[1]] > 0.25)/length(doublet[[1]]))[2], "%\n"))
mmg$new_doublet_score <- doublet[[1]]
mmg$new_doublet_ident <- ifelse(doublet[[1]] > 0.25, "Doublet", "Singlet")

saveRDS(mmg, obj)
cat("complete.\n")
