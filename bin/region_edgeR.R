
suppressPackageStartupMessages(library(dplyr))
library(ggplot2)
suppressPackageStartupMessages(library(Seurat))
source("/projects/ps-renlab2/y2xie/scripts/basics.R")
source("/projects/ps-renlab2/y2xie/scripts/DPT_help.R")
source("/projects/ps-renlab2/y2xie/scripts/DPT/BICAN_help.R")
library(harmony)

args <- commandArgs(trailingOnly = TRUE)
rpath <- args[1] ### path to input
outprfx <- args[2] ### outprefix

set.seed(921)

mmg <- qs::qread(rpath)
rmeta <- read.csv("/projects/ps-renlab2/y2xie/projects/BICAN/analysis/04.clustering/02.RNA/03.l3_clustering/MiniAtlas_RNA_merged_dual_filt_clean_corrected_250722_metadata.txt", sep = '\t')
good_cells <- rownames(rmeta[rmeta$keep == 1,])
good_cells <- intersect(good_cells, colnames(mmg))

### 251110: only slice cortical regions for comparsion
cortx <- rownames(rmeta[rmeta$structure == "Cortex",])
good_cells <- intersect(good_cells, cortx)
rmeta <- rmeta[good_cells, ]

### split into two pseudoreplicates
rmeta <- rmeta %>% mutate(donor = sample(c("A","B"), n(), replace = TRUE))

# collect cells to keep (balanced by donor)
regions <- sort(unique(rmeta$region))
dsel <- unique(rmeta$donor)
cells_per_donor_cap <- 500 ### try to balance donor 

cells_keep <- unlist(lapply(regions, function(r) {
    unlist(lapply(dsel, function(dn) {
    ids <- rownames(rmeta)[rmeta$region == r & rmeta$donor == dn]
    if (is.finite(cells_per_donor_cap) && length(ids) > cells_per_donor_cap) {
      sample(ids, cells_per_donor_cap, replace = F)
    } else ids
  }))
}), use.names = FALSE)

mmg <- subset(mmg, cells = cells_keep)
mmg$donor <- rmeta[colnames(mmg),]$donor

exprr <- AggregateExpression(mmg, group.by = c("donor", "region"))
designs = stringr::str_split_fixed(colnames(exprr$RNA), n = 2, pattern = "_") %>% as.data.frame %>% setNames(c("sample", "group"))
rownames(designs) <- colnames(exprr$RNA)
edgeR_ovr_pairwise(exprr$RNA, sample = designs, group = "group", pairwise = FALSE, 
          prefix = paste0(outprfx, "_"))



