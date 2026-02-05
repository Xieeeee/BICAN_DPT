
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(Seurat))
library(ggplot2)
library(lisi)
library(cluster)
library(RColorBrewer)
library(harmony)
source("/projects/ps-renlab2/y2xie/scripts/basics.R")
source("/projects/ps-renlab2/y2xie/scripts/DPT_help.R")
source("/projects/ps-renlab2/y2xie/scripts/DPT/BICAN_help.R")

args <- commandArgs(trailingOnly = TRUE)
rpath <- args[1] ### path to input qs file
outprfx <- args[2] ### outprefix

compute_asw_region <- function(emb, labels) {
  # Average silhouette width using Euclidean distance in embedding
  if(length(unique(labels)) < 2) return(NA_real_)
  d <- dist(emb, method="euclidean")
  sw <- silhouette(as.integer(factor(labels)), d)
  return(sw) #(sw[, "sil_width"], na.rm=TRUE)
}

mmg <- qs::qread(rpath)

# qfeature <- VariableFeatures(mmg)
qregion <- c('A46','A38','A25','Idg','MTG','ITG','M1C','Ig','A5-A7','TH-TL','Pro','V2','A19','V1C')

### to calculate lisi, try control the number of cells per region per donor to be similar
set.seed(921)
qcell <- list()
for (f in qregion){
    for (d in unique(mmg$donor)){
        tmp <- mmg@meta.data %>% filter(region == f, donor == d) %>% rownames
        if (length(tmp) >= 200){
            tmp <- sample(tmp, size = 200, replace = F)
        }
        qcell[[paste0(f, ":", d)]] <- tmp
    }
}

tmp <- subset(mmg, cells = unlist(qcell))
tmp[["RNA"]] <- as(tmp[["RNA"]], "Assay")
tmp <- seurat_onestep_clust(tmp, batch.label = "donor", res = 0.3)

X1 <- tmp@reductions$harmony@cell.embeddings
tmeta <- tmp@meta.data[,c("region", 'donor','target')]
lisi1 <- compute_lisi(X1, tmeta, colnames(tmeta))
write.table(lisi1, paste0(outprfx, ".lisi"), sep = "\t", quote = F, row.names = T, col.names = T)

asw <- compute_asw_region(X1, tmeta[,"region"])
write.table(asw, paste0(outprfx, ".asw"), sep = "\t", quote = F, row.names = T, col.names = T)
qs::qsave(tmp, paste0(outprfx, ".sample.qs"))
