
## 2024-04-16
## Yang Xie (y2xie@ucsd.edu)
## Functions for Droplet Paired-Tag analysis in R used for BICAN project

library(Matrix)
suppressPackageStartupMessages(library(dplyr)) 
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SoupX)) 
suppressPackageStartupMessages(library(harmony)) 

source("/projects/ps-renlab2/y2xie/scripts/basics.R")
source("/projects/ps-renlab2/y2xie/scripts/DPT_help.R")
source("/projects/ps-renlab2/y2xie/scripts/hicat_vg.R")
source("/projects/ps-renlab2/y2xie/scripts/hicat_doubletfinder_mo.R")

ProcessRNA_hicat <- function(rna_path, valid, cpu = 4, subset = TRUE, cluster = TRUE) {
	cat("Loading and validating input...\n")
	mmg <- load_and_validate_rna(rna_path, valid, cpu = cpu)
	mmg <- run_qc_and_clustering(mmg, subset = subset, cluster = cluster)
	if (cluster){
		mmg <- run_soupX(mmg, rna_path)
	}
	return(mmg)
}

load_and_validate_rna <- function(rna_path, valid, cpu = 4) {
	mmg <- Read10X(rna_path)
	if (!("rna" %in% colnames(valid))) {
		stop("rna barcode information is not found in valid barcodes file. Make sure column with name 'rna' exists.")
	}
	if ("atac_bc" %in% colnames(valid)){
        cat("The provided barcode file seems to be a multiplexed reaction. Determine multiplet based on the observed collision rate...\n")
        nondup <- valid$rna[!duplicated(valid$rna) & !duplicated(valid$rna, fromLast = TRUE)] ### RNA barcodes appear once
        dup <- valid$rna[duplicated(valid$rna) | duplicated(valid$rna, fromLast = TRUE)] %>% unique ### RNA barcodes appear multiple times
        cat("Number of RNA barcodes matched with multiple ATAC barcodes: ", length(unique(dup)), "\n")
        observed_collide <- 100 * (length(unique(valid$rna)) - length(nondup))/length(unique(valid$rna))
        cat(paste0("observed collision rate:", round(observed_collide, 2), "%\n"))
        # estimated_collide <- 100 * (length(unique(valid$rna)) * 0.08)/10000
        # collide <- estimated_collide - observed_collide
        preindex <- substr(valid$atac_bc, 1, (nchar(valid$atac_bc[1]) - 16))
        collide_ratio <- sum(table(preindex)^2)/(sum(table(preindex))^2 - sum(table(preindex)^2))
        collide <- observed_collide * collide_ratio
        if(is.na(collide)){
            collide <- 0
        }
        part <- max(0.2, collide/100) ### update 240709: change to default
        cat(paste0("Estimated doublets rate after removing explict multiplet: ", round(collide, 2), "%\n"))
    }else if (length(unique(valid$rna)) != nrow(valid)){
        stop("Duplicated RNA barcodes detected, but the valid barcode file is not for multiplexed sample.")
    }else{
        nondup <- unique(valid$rna)
        collide <- 100 * (length(nondup) * 0.08)/10000
        part <- max(0.20, collide/100)
    }
    mmg <- mmg[, nondup]
    mmg <- CreateSeuratObject(mmg)

    ### to accomodate seurat5...
    mmg[["RNA"]] <- as(object = mmg[["RNA"]], Class = "Assay")
    mmg <- NormalizeData(mmg, normalization.method = "LogNormalize", scale.factor = 1e6)
	mmg <- FindVariableFeatures(mmg, nfeatures = 3000)
    mmg <- run_doublet_detection(mmg, part = part, cpu = cpu)
    return(mmg)
}

run_soupX <- function(mmg, rna_path){
    tmp <- as(object = mmg[["RNA"]], Class = "Assay")
	toc <- tmp@counts
	tod <- Read10X(rna_path)
	sc <- SoupChannel(tod, toc)
	sc <- setDR(sc, as.data.frame(mmg@reductions$umap@cell.embeddings))
	sc <- setClusters(sc, mmg$seurat_clusters)
	sc <- autoEstCont(sc)
	out <- adjustCounts(sc, roundToInt=TRUE)
	background <- sc$soupProfile[order(sc$soupProfile$est,decreasing=TRUE),]
	mmg[["RNA.soup"]] <- CreateAssayObject(counts = out)
	Misc(object = mmg, slot = "soupX_bg") <- background
	return(mmg)
}

run_doublet_detection <- function(mmg, part = 0.05, cpu = cpu) {
	library(foreach)
	library(doParallel)
	cl <- makeCluster(cpu)
	clusterEvalQ(cl, {
		library(Seurat)
		source("/tscc/projects/ps-renlab2/y2xie/scripts/hicat_doubletfinder_mo.R")
		source("/tscc/projects/ps-renlab2/y2xie/scripts/hicat_vg.R")
		source("/tscc/projects/ps-renlab2/y2xie/scripts/DPT/BICAN_help.R")
	})
	# clusterExport(cl, varlist = c("mmg", "part"))
	registerDoParallel(cl)

	dscores <- foreach(i = 1:50, .combine = 'c') %dopar% {
		doubletfinder(mmg, VariableFeatures(mmg), proportion.artificial = part, plot = FALSE)
	}
	stopCluster(cl)

	even_indices <- seq(1, length(dscores), 2)
	dconstant <- Reduce(intersect, lapply(dscores[even_indices], function(x) names(x)[which(x > 0.2)]))

	avscore <- do.call(cbind, lapply(dscores[even_indices], function(x) as.data.frame(x))) %>%
	setNames(as.character(1:10)) %>%
	rowMeans() %>% as.data.frame() %>%
	setNames("avg_dscore") %>%
	tibble::rownames_to_column("bc") %>%
	mutate(group = ifelse(bc %in% dconstant, "Constant", "Varied")) %>%
	mutate(group = factor(group, levels = c("Varied", "Constant")))

	mmg$doublet_score <- avscore$avg_dscore
	mmg$doublet_ident <- ifelse(colnames(mmg) %in% dconstant, "Doublet", "Singlet")
	return(mmg)
}

run_qc_and_clustering <- function(mmg, subset = TRUE, cluster = TRUE) {
	mmg$percent.mt <- PercentageFeatureSet(mmg, pattern = "^MT-")
	mmg$percent.ribo <- PercentageFeatureSet(mmg, pattern = "(^RPL|^RPS|^MRP)")

	gene_hcf <- median(log10(mmg$nFeature_RNA)) + 3 * mad(log10(mmg$nFeature_RNA))
	count_hcf <- median(log10(mmg$nCount_RNA)) + 3 * mad(log10(mmg$nCount_RNA))
  
	cat(paste0(
	"\ngene highest cutoff: ", round(10^gene_hcf), 
	"\nreads highest cutoff: ", min(50000, round(10^count_hcf)), "\n"))
	cat(paste0("Cells number: ", ncol(mmg), "\n"))
	if (subset){
		mmg <- subset(mmg, subset = nCount_RNA < min(50000, 10^count_hcf) & nFeature_RNA < 10^gene_hcf & percent.mt < 10)
		cat(paste0("Cells number after quality segment: ", ncol(mmg), "\n"))
	}
	if (cluster){
		mmg <- seurat_onestep_clust(mmg)
	}
	return(mmg)
}


### implement of onestep_clust() from hicat
seurat_onestep_clust <- function(mmg, reduction = "pca", var = "none", batch.label = "none", npc = "none", res = "none"){
    cat("Perform log CPM normalization...\n")
    mmg <- NormalizeData(mmg, normalization.method = "LogNormalize", scale.factor = 1e6)
    norm.dat <- mmg[["RNA"]]@data
    
    cat("Perform variable features selection using Brennecke's method...\n")
    filt_genes <- row.names(norm.dat)[which(Matrix::rowSums(norm.dat > 1) >= 10)]
    norm.count <- norm.dat
    norm.count@x = 2^(norm.count@x) - 1
    vg <- find_vg(norm.count[filt_genes,], rescaled = F, return_type = "data", verbose = F)
    sg <- vg[which(vg$loess.padj < 0.5 | vg$dispersion > 3), "gene"]
    cat(paste0("number of variable genes calculated by hicat method: ", length(sg), "\n"))
    
    if (length(sg) > 3000){
        sg <- head(sg[order(vg[sg, "loess.padj"], -vg[sg, "z"])], 3000)
    }
    Misc(object = mmg, slot = "hicat_vg") <- vg
    VariableFeatures(mmg) <- sg

    if(length(sg) <= 50){
        cat("too few varibale genes found...")
        mmg <- FindVariableFeatures(mmg, nfeatures = 3000)
    }
    
    cat("Perform data scaling and PCA...\n")
    if (var != "none") {
        mmg <- ScaleData(mmg, vars.to.regress = var)
    }else {
        mmg <- ScaleData(mmg)
    }
    mmg <- RunPCA(mmg, features = VariableFeatures(object = mmg), npcs = 50, verbose = F)
    
    ### select top pc
    if (npc == "none") {
        npc <- significant_pc_test(mmg@reductions$pca@cell.embeddings, min_pc = 15)
    }
    cat(paste0("Number of PC selected for downstream analysis: ", npc, "\n")) 
    
    if (batch.label != "none") {
        library(harmony)
        mmg <- RunHarmony(mmg, group.by.vars = batch.label, dims.use = 1:npc)
        reduction <- "harmony"
    }
    
    if(ncol(mmg) >= 500000){
    k <- 50
    }else if(ncol(mmg) >= 200000 & ncol(mmg) < 500000){
        k <- 25
    }else{
        k <- 15
    }
    cat(paste0("Number of neighbors selected for downstream analysis: ", k, "\n"))
    
    cat("Perform UMAP embedding...\n")
    mmg <- RunUMAP(mmg, seed.use = 131, reduction = reduction, dims = 1:npc,
                  n.neighbors = k, min.dist = 0.1, n.components = 2L, 
                  umap.method = "uwot")

    if (res == "none") {
        res <- 0.3
    }
    cat(paste0("resolution used for leiden clustering: ", res, "\n")) 

    mmg <- FindNeighbors(object = mmg, k.param = k, reduction = reduction, dims = 1:npc, 
                     nn.method = "annoy", annoy.metric = "cosine", verbose = F)
    
    ### in hicat the resolution for one step clustering is set to be 0.01. 
    ### update 240206: Extremely slow when using leidenalg. Use leidenbase instead for faster clustering.
    suppressWarnings(mmg <- FindClusters(object = mmg, algorithm = 4, method = "igraph", weights = T, resolution = res, verbose = F))

    ### update 240208: do not have an efficient way of dealing with singleton. From now on export RNA_snn and use leidenalg for clustering.
#     cat("Convering SNN graph to igraph...\n")
#     graph <- mmg@graphs$RNA_snn
#     adj.matrix <- Matrix::sparseMatrix(i = graph@i+1, p = graph@p, x = graph@x, dims = graph@Dim, dimnames = graph@Dimnames)
#     sgraph <- igraph::graph_from_adjacency_matrix(adjmatrix = adj.matrix, weighted = TRUE, diag = FALSE)
#     sgraph <- igraph::simplify(sgraph, remove.multiple = T)
    
#     cat("Performing leiden clustering with leidenbase...\n")
#     lres <- leidenbase::leiden_find_partition(sgraph, resolution_parameter = res, num_iter = 10, partition_type = "RBConfigurationVertexPartition", seed = 921)
#     cat(paste0("Modulaity: ", lres$modularity, "\n"))
#     mmg$leiden <- lres$membership
    
    ### to accomodate seurat5...
    mmg[["RNA"]] <- as(object = mmg[["RNA"]], Class = "Assay5")
    return(mmg)
}

### calculate number of PC used for clustering; this is an implementation of ALLCools significant_pc_test function
significant_pc_test <- function(pc, p_cutoff = 0.1, min_pc = 15, downsample = 25000){
    if(nrow(pc) > downsample){
        pc <- pc[sample(nrow(pc), size = downsample), ]
        print(paste0("Downsample PC matrix to ", downsample, " cells to calculate significant PC components"))
    }
    for (i in 1:(ncol(pc)-1)){
        cur_pc = pc[, i]
        next_pc = pc[, i + 1]
        test_result <- ks.test(cur_pc, next_pc)
        if (test_result$p.value > p_cutoff){
            break
        }
    }
    n_components <- min(i + 1, ncol(pc))
    min_pc <- min(15, ncol(pc))
    if (n_components < min_pc){
        print(paste0('only ', n_components, ' passed the P cutoff. In order to proceed following analysis, will use first ', min_pc, ' PCs'))
        n_components <- min_pc
    } 
    return(n_components)
}
