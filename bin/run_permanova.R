#!/usr/bin/env Rscript

suppressPackageStartupMessages({
library(argparse)
library(reticulate)
library(anndata)
library(vegan)     # adonis2
library(cluster)   # silhouette
})


parser <- ArgumentParser(
  description = "Quantify region separability in an embedding from an h5ad file"
)

parser$add_argument(
  "-i", "--h5ad", required = TRUE,
  help = "Input .h5ad file"
)
parser$add_argument(
  "-r", "--region_col", required = TRUE,
  help = "obs column name indicating region labels"
)
parser$add_argument(
  "-e", "--embedding", required = TRUE,
  help = "obsm key for embedding (e.g. X_pca, X_lsi)"
)
parser$add_argument(
  "-o", "--out_prefix", required = TRUE,
  help = "Output prefix for results TSV"
)
parser$add_argument(
  "-n", "--max_cells_per_region", type = "integer", default = 500,
  help = "Maximum number of cells per region (for downsampling) [default: 500]"
)
parser$add_argument(
  "-d", "--ndim", type = "integer", default = 30,
  help = "Number of embedding dimensions to use [default: 30]"
)
parser$add_argument(
  "-s", "--seed", type = "integer", default = 1,
  help = "Random seed for downsampling [default: 1]"
)

args <- parser$parse_args()

if (!file.exists(args$h5ad)) {
  stop("h5ad not found: ", args$h5ad)
}

set.seed(921)
adata <- read_h5ad(args$h5ad)
obs <- adata$obs
region <- obs[[args$region_col]]

obsm_r <- adata$obsm
emb <- obsm_r[[args$embedding]]
emb <- as.matrix(emb)

valid <- !is.na(region)
emb <- emb[valid, , drop = FALSE]
region <- droplevels(as.factor(region[valid]))

if (nlevels(region) < 2L) {
  stop("Fewer than 2 non-NA region levels after filtering; cannot compute separability.")
}

ndim <- min(args$ndim, ncol(emb))
emb <- emb[, seq_len(ndim), drop = FALSE]
idx_by_region <- split(seq_along(region), region)

keep_idx <- unlist(lapply(idx_by_region, function(ix) {
  if (length(ix) > args$max_cells_per_region) {
    sample(ix, args$max_cells_per_region)
  } else {
    ix
  }
}), use.names = FALSE)

keep_idx <- sort(keep_idx)
emb_sub <- emb[keep_idx, , drop = FALSE]
region_sub <- droplevels(region[keep_idx])

n_cells   <- length(region_sub)
n_regions <- nlevels(region_sub)

cat("[INFO] Computing distance matrix (Euclidean)...\n")
dist_mat <- dist(emb_sub)

meta_sub <- data.frame(region = region_sub)

cat("[INFO] Running PERMANOVA (adonis2)...\n")
adon <- adonis2(dist_mat ~ region, data = meta_sub)

region_row <- adon[1, ]
R2_region  <- as.numeric(region_row[["R2"]])
p_region   <- as.numeric(region_row[["Pr(>F)"]])

cat("[INFO] Computing mean silhouette by region...\n")
# numeric labels for silhouette
lab_int <- as.integer(region_sub)
sil     <- silhouette(lab_int, dist_mat)
mean_sil <- mean(sil[, "sil_width"])

# output --------------------------------------------------------------

out_df <- data.frame(
  h5ad                 = basename(args$h5ad),
  embedding            = args$embedding,
  region_col           = args$region_col,
  n_cells              = n_cells,
  n_regions            = n_regions,
  max_cells_per_region = args$max_cells_per_region,
  ndim                 = ndim,
  R2_region            = R2_region,
  p_region             = p_region,
  mean_silhouette      = mean_sil,
  stringsAsFactors     = FALSE
)

out_file <- paste0(args$out_prefix, "_region_separability.tsv")
cat("[INFO] Writing results to: ", out_file, "\n", sep = "")
write.table(
  out_df,
  file      = out_file,
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

cat("[INFO] Done.\n")

