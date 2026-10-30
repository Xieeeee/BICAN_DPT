#' Doublet detection in single-cell RNA sequencing data
#'
#' Adopted from https://www.cell.com/cell-systems/fulltext/S2405-4712(19)30073-0 
#' (https://github.com/chris-mcginnis-ucsf/DoubletFinder) 
#'
#' This function generates artificial nearest neighbors from existing single-cell RNA
#' sequencing data. First, real and artificial data are merged. Second, dimension reduction
#' is performed on the merged real-artificial dataset using PCA. Third, the proportion of
#' artificial nearest neighbors is defined for each real cell. Finally, real cells are rank-
#' ordered and predicted doublets are defined via thresholding based on the expected number
#' of doublets.
 
doubletfinder <- function(mmg, select.genes, proportion.artificial = 0.2,
                          k = NULL, plot = FALSE) {
  
    if(is.null(k)) {
      k <- round(pmin(100, ncol(mmg) * 0.01))
      k <- round(pmax(10, k))
    }
      
  
  library(RANN)
  
  ## Step 1: Generate artificial doublets from Seurat object input
  print("Creating artificial doublets...")
  ### seurat v5 to v3
  mmg[["RNA"]] <- as(object = mmg[["RNA"]], Class = "Assay")
  data <- mmg[["RNA"]]@counts
  real.cells <- colnames(data)
  
  n_real.cells <- length(real.cells)
  n_doublets <- round((n_real.cells/(1-proportion.artificial))-n_real.cells)
  real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
  real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
  doublets <- (data[ , real.cells1] + data[ , real.cells2])
  colnames(doublets) <- paste("X", 1:n_doublets, sep="")
  
  data_wdoublets <- cbind(data, doublets)
  data_wdoublets <- CreateSeuratObject(data_wdoublets) 
  data_wdoublets[["RNA"]] <- as(object = data_wdoublets[["RNA"]], Class = "Assay")
  data_wdoublets <- NormalizeData(data_wdoublets, normalization.method = "LogNormalize", scale.factor = 1e6)
  VariableFeatures(data_wdoublets) <- select.genes
  norm.dat <- data_wdoublets[["RNA"]]@data
  data_wdoublets <- ScaleData(data_wdoublets)
  
  print("Running PCA")
  data_wdoublets <- RunPCA(data_wdoublets, features = select.genes, npcs = 50, verbose = F)
  npc <- significant_pc_test(data_wdoublets@reductions$pca@cell.embeddings, min_pc = 15)
  pca.coord <- data_wdoublets@reductions$pca@cell.embeddings[ , 1:npc]

  ## Compute pANN
  print("Initialize pANN structure")
  knn.result = RANN::nn2(pca.coord, k=k)
  knn.idx = knn.result[[1]]
  knn.dist = knn.result[[2]]
  
  num = ncol(data)
  knn.result1 = RANN::nn2(pca.coord[1:num,], pca.coord[(num+1):nrow(pca.coord),], k = 10)
  knn.dist1 = knn.result1[[2]]
  
  dist.th = mean(as.vector(knn.dist1)) + 1.64 * sd(as.vector(knn.dist1))
  
  doublet.freq = knn.idx  > ncol(data) & knn.dist < dist.th
  doublet.score = pmax(rowMeans(doublet.freq),rowMeans(doublet.freq[,1:ceiling(k/2)]))
  names(doublet.score) = row.names(pca.coord)
  artificial.doublet.score=doublet.score[colnames(doublets)]
  doublet.score=doublet.score[colnames(data)]
    
  if (plot == TRUE) {
    print("plotting")
    ds = pmax(rowMeans(doublet.freq),rowMeans(doublet.freq[,1:ceiling(k/2)])) 
    ds=as.data.frame(ds)
    ds$sample <- colnames(data_wdoublets)
    
    ds$group <- ""
    idx <- startsWith(ds$sample,"X")
    ds[idx, "group"] <- "artifical doublets"
    idx <- !startsWith(ds$sample,"X")
    ds[idx, "group"] <- "samples"
    
    plot.title <- gsub("^.*?-","",ds[1,2])
    
    p <- ggplot2::ggplot(ds, aes(x = ds, fill=group, color=group)) +
    geom_density(alpha=0.4)+scale_color_manual(values=c("#F9627D","#2F3061"))+
    scale_fill_manual(values=c("#F9627D","#2F3061")) + 
    scale_x_continuous(breaks = seq(0, 1, 0.1)) + 
    labs(title=plot.title)
    
    return(list(doublet.score,artificial.doublet.score, p))
    
  } else {  return(list(doublet.score,artificial.doublet.score)) }
  
}