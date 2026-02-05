
suppressPackageStartupMessages(library(dplyr))
library(SeuratExtend)
library(ggplot2)
library(data.table)
library(future)
library(future.apply)

source("/projects/ps-renlab2/y2xie/scripts/basics.R")
source("/projects/ps-renlab2/y2xie/scripts/DPT_help.R")

# SeuratExtend::activate_python("/home/y2xie/miniconda3/envs/seurat/bin/python")
library(reticulate)
use_python("/home/y2xie/miniconda3/envs/seurat/bin/python", required = TRUE)

args <- commandArgs(trailingOnly = TRUE)
rpath <- args[1] ### path to input
outprfx <- args[2] ### outprefix

obj <- qs::qread(rpath)
obj <- Palantir.RunDM(obj, reduction = "harmony")

### Find out start cell
scell <- obj@meta.data %>% filter(group == 1) %>% rownames
dlist <- list()
for (f in colnames(obj@reductions$ms@cell.embeddings)){
	dlist[[f]] <- data.frame(pc = f, scc = cor(obj$group, obj@reductions$ms@cell.embeddings[, f], method = "spearman"))
}
dlist <- do.call(rbind, dlist) 
qemb <- dlist %>% arrange(desc(abs(scc))) %>% head(1) %>% rownames
qsign <- dlist %>% arrange(desc(abs(scc))) %>% head(1) %>% select(scc) %>% unlist %>% sign

scell1 <- as.data.frame(obj@reductions$ms@cell.embeddings[scell, qemb, drop = F]) %>% setNames("ddim") %>% slice_min(order_by = ddim, with_ties = F) %>% rownames
scell2 <- as.data.frame(obj@reductions$ms@cell.embeddings[scell, qemb, drop = F]) %>% setNames("ddim") %>% slice_max(order_by = ddim, with_ties = F) %>% rownames
if (qsign == -1){scell = scell2}else{scell = scell1}

# obj@meta.data[scell, c("orig.ident", "region", "group")]
obj <- Palantir.Pseudotime(obj, start_cell = scell)
ps <- obj@misc$Palantir$Pseudotime
obj$Pseudotime <- ps$Pseudotime
qs::qsave(obj, paste0(outprfx, ".qs"))
write.table(obj@meta.data[,c("region", "group", "Pseudotime")], paste0(outprfx, ".pseudotime.txt"), sep = "\t", quote = F, row.names = T, col.names = T)
write.table(dlist, paste0(outprfx, ".ms.scc.txt"), sep = "\t", quote = F, row.names = T, col.names = T)
