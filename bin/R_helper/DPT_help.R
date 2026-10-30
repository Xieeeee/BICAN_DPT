
## 2022-12-19
## Yang Xie (y2xie@ucsd.edu)
## Functions for Droplet Paired-Tag analysis in R

### set color ###
library(Matrix)
library(dichromat)
library(viridis)

tori <- c("#477a96", "#c3533d", "#85212b", "#bcc9d1", "#6a9d84")
colfunc2 <- colorRampPalette(c(tori))

### mouse to human gene conversion or vice versa
convert_mouse_to_human <- function(gene_list, direction = "m2h"){
    mouse_human_genes = read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt", sep = "\t")
    output = c()
    if(direction == "m2h"){
        qry <- "mouse, laboratory"
        ref <- "human"
    }else if(direction == "h2m"){
        ref <- "mouse, laboratory"
        qry <- "human"
    }else{
        stop("direction can only be m2h or h2m.")
    }
    for(gene in gene_list){
        class_key = (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name==qry))[['DB.Class.Key']]
        if(!identical(class_key, integer(0)) ){
            human_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name==ref))[,"Symbol"]
            for(human_gene in human_genes){
                output = append(output,human_gene)
            }
        }
    }
    return(output)
}

### Process homer annotatePeak.pl results
sum_anno <- function(df){
	    t1 <- table(df$Annotation) %>% as.data.frame %>%
		        mutate(Var1 = gsub("^([^\\(]+).*", "\\1", Var1)) %>%
			    group_by(Var1) %>% summarise(Freq = sum(Freq))

		        t1a1 <- t1[grep("\\|", t1$Var1, invert = T), ] %>%
				    mutate(Var1 = gsub("^(CpG)-.*$", "\\1", Var1)) %>%
				        group_by(Var1) %>% summarise(Freq = sum(Freq)) %>% setNames(c("Var1", "Freq"))

				    t1a2 <- t1[grep("\\|", t1$Var1), ] %>% 
					        mutate(Var1 = gsub("^([^|]+)\\|([^|]+)\\|([^|]+)$", "\\2", Var1)) %>%
						    mutate(Var1 = gsub("\\?", "", Var1)) %>%
						        mutate(simple = ifelse(Var1 %in% c("LINE", "LTR", "SINE", "Simple_repeat"), Var1, "Other repeats")) %>%
							    group_by(simple) %>% summarise(Freq = sum(Freq)) %>% setNames(c("Var1", "Freq"))

						        anno <- rbind(t1a1, t1a2)
							    return(anno)
}

sum_detailed_anno <- function(df){
	    t1 <- table(df$Detailed.Annotation) %>% as.data.frame %>%
		        mutate(Var1 = gsub("^([^\\(]+).*", "\\1", Var1)) %>%
			    group_by(Var1) %>% summarise(Freq = sum(Freq))

		        t1a1 <- t1[grep("\\|", t1$Var1, invert = T), ] %>%
				    mutate(Var1 = gsub("^(CpG)-.*$", "\\1", Var1)) %>%
				        group_by(Var1) %>% summarise(Freq = sum(Freq)) %>% setNames(c("Var1", "Freq"))

				    t1a2 <- t1[grep("\\|", t1$Var1), ] %>% 
					        mutate(Var1 = gsub("^([^|]+)\\|([^|]+)\\|([^|]+)$", "\\2", Var1)) %>%
						    mutate(Var1 = gsub("\\?", "", Var1)) %>%
						        mutate(simple = ifelse(Var1 %in% c("LINE", "LTR", "SINE", "Simple_repeat"), Var1, "Other repeats")) %>%
							    group_by(simple) %>% summarise(Freq = sum(Freq)) %>% setNames(c("Var1", "Freq"))

						        anno <- rbind(t1a1, t1a2)
							    return(anno)
}

### read in fragment counts to calculate frip ###
ImportArcFRiP <- function(raw_count, frip_count){
    tmp1 <- read.table(paste0(raw_count), header = T, row.names = 1)
    tmp2 <- read.table(paste0(frip_count), header = T, row.names = 1)
    frip <- merge(x = tmp1, y = tmp2, by = 0) %>% setNames(c("bc", "counts_total", "counts_unique", "counts_in_peaks", "unique_in_peaks"))
    frip$frip <- frip$unique_in_peaks / frip$counts_unique
    return(frip)
}

### Find FRiP cutoff ###
CutoffArcFRiP <- function(frip, xcut_low = 100, ycut_low = 0.05, ycut_high = 0.8, prefix){

}

### Plot calculated frip ###
PlotArcFRiP <- function(frip, xcut_low = 100, ycut_low = 0.05, ycut_high = 0.8, prefix){
    valid <- frip[frip$frip > ycut_low & frip$frip < ycut_high & frip$counts_unique > xcut_low, ]
    label <- data.frame(anno = paste0("Unique reads cutoff: ", xcut_low, "\n",
                                      "FRiP cutoff: ", ycut_low, "-", ycut_high, "\n", 
                                      "PF cells: ", nrow(valid), "\n")) 
                                      # "Median total reads: ", as.integer(median(valid$counts_total))))
    t1 <- frip %>%
    ggplot(aes(x = counts_unique, y = 100*frip)) + 
    geom_point(size = 0.5, alpha = 0.1, color = "grey") + 
    geom_vline(xintercept = xcut_low, color = colfunc2(1), linetype="dashed") + 
    geom_hline(yintercept = 100*ycut_low, color = colfunc2(1), linetype="dashed") + 
    geom_hline(yintercept = 100*ycut_high, color = colfunc2(1), linetype="dashed") +
    geom_point(data = valid, aes(x = counts_unique, y = 100*frip),
               size = 0.1, alpha = 0.5, color = colfunc2(3)[[2]]) + 
    geom_text(data = label, aes(x = max(frip$counts_unique), y = max(100*as.numeric(frip$frip)), 
                                hjust = 1, vjust = 1, label = anno), size = 2.5) +  
    theme_Publication() + xlab("log10(Fragments)") + ylab("% FRiP") + 
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 10)) +   
    scale_x_log10()
    
    write.table(valid, file = paste0(prefix, "_PF_cells.txt"), row.names = F, col.names = F, sep = "\t", quote = F)
    ggsave(t1, filename = paste0(prefix, "_FRiP_macs2.png"), height = 6, width = 6, dpi = 300)
    return(valid)
}

### Plot calculated frip, color by local jitter density for better visualization ###
PlotArcFRiP_d <- function(frip, xcut_low = 100, ycut_low = 0.05, ycut_high = 0.8, prefix){
    filt_data <- frip
    ### assign color density 
    x <- densCols(log10(filt_data$counts_unique), 100*filt_data$frip, colramp=colorRampPalette(viridis(12)))
    frip$dens <- col2rgb(x)[1,] + 1L

    ### calculate density again to exclude low counts population
    filt_data <- filt_data %>% filter(counts_unique > 100)
    x <- densCols(log10(filt_data$counts_unique), 100*filt_data$frip, colramp=colorRampPalette(viridis(12)))
    frip[frip$counts_unique > 100, ]$dens <- (col2rgb(x)[1,] + 1L)

    valid <- frip[frip$frip > ycut_low & frip$frip < ycut_high & frip$counts_unique > xcut_low, ]
    label <- data.frame(anno = paste0("Unique reads cutoff: ", xcut_low, "\n",
                                      "FRiP cutoff: ", ycut_low, "-", ycut_high, "\n", 
                                      "PF cells: ", nrow(valid), "\n")) 
                                      # "Median total reads: ", as.integer(median(valid$counts_total))))
    t1 <- frip %>%
    ggplot(aes(x = counts_unique, y = 100*frip)) + 
    geom_point(size = 0.5, alpha = 0.1, aes(color = dens)) + 
    geom_vline(xintercept = xcut_low, color = colfunc2(1), linetype="dashed") + 
    geom_hline(yintercept = 100*ycut_low, color = colfunc2(1), linetype="dashed") + 
    geom_hline(yintercept = 100*ycut_high, color = colfunc2(1), linetype="dashed") +
    geom_text(data = label, aes(x = max(frip$counts_unique), y = max(100*as.numeric(frip$frip)), 
                                hjust = 1, vjust = 1, label = anno), size = 2.5) +  
    theme_Publication() + xlab("log10(Fragments)") + ylab("% FRiP") + 
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 10)) +   
    scale_x_log10() + 
    theme(legend.position = "") + 
    scale_color_viridis(trans = "log")
    
    write.table(valid, file = paste0(prefix, "_PF_cells.txt"), row.names = F, col.names = F, sep = "\t", quote = F)
    ggsave(t1, filename = paste0(prefix, "_FRiP_macs2.png"), height = 6, width = 6, dpi = 300)
    return(valid)
}

### Plot CUT&TAG library fragment length pattern ###
PlotArcFragment <- function(frag_path, prefix){
    suppressMessages(library(data.table))
    frag <- fread(frag_path)
    frag$len <- frag$V3 - frag$V2
    frags <- frag[sample(nrow(frag), 10000, replace = F), ]
    frags <- frags[, c("len", "V5")]
    frags_pv <- frags[rep(1:nrow(frags), frags[["V5"]]), ]
    t6 <- frags_pv %>% dplyr::filter(len < 1000) %>%
    ggplot(aes(x = len)) + 
    geom_histogram(color = "white", bins = 100) +
    theme_Publication() + 
    scale_x_continuous(breaks = seq(0,1000,100)) + 
    xlab("Fragment length")
    ggsave(t6, filename = paste0(prefix, "_fragments.pdf"), height = 6, width = 6, dpi = 300)
    return(t6)
}

### After in-house processing of 10X Multiome, merge both modalities for plotting ###
PairArc <- function(dmat, rmat, names = "atac"){
    translate <- read.table("/projects/ps-renlab2/y2xie/ps-renlab/y2xie/projects/genome_ref/arc_bc-translation.txt", header = T, row.names = 1, sep = "\t")
    translate$atac <- paste0(translate$atac, "-1")
    translate$rna <- paste0(translate$rna, "-1")
    # rownames(translate) <- paste0(translate$rna, "-1")
    if (names == "atac"){ ### cells name in dna
        idx <- match(colnames(rmat), translate$rna)
        colnames(rmat) <- translate[idx, names]
        colnames(dmat) <- colnames(dmat)
    }else if (names == "rna"){ ### cells name in rna
        idx <- match(colnames(dmat), translate$atac)
        colnames(dmat) <- translate[idx, names]
        colnames(rmat) <- colnames(rmat)
    }
    dna <- colSums(dmat) %>% as.data.frame() %>% setNames("dna") %>% tibble::rownames_to_column("bc")
    rna <- colSums(rmat) %>% as.data.frame() %>% setNames("rna") %>% tibble::rownames_to_column("bc")
    merge1 <- merge(x = dna, y = rna, by = "bc")
    return(merge1)
}

### Plot reads count for each pair of barcode. Save plot and segmented cells with name prefix ###
PlotArcPair <- function(pair, dcutoff = 100, rcutoff = 100, prefix, names = "atac"){
    # if(!identical(colnames(pair), c("bc", "dna_count", "rna_count"))){
    if((!all(c("bc", "dna_count", "rna_count") %in% colnames(pair)))){
        stop('colnames of input dataframe should be c("bc", "dna_count", "rna_count")')
    }
    translate <- read.table("/projects/ps-renlab2/y2xie/ps-renlab/y2xie/projects/genome_ref/arc_bc-translation.txt", header = T, row.names = 1, sep = "\t")
    translate$atac <- paste0(translate$atac, "-1")
    translate$rna <- paste0(translate$rna, "-1")
    
    pair_pf <- pair[pair$dna_count > dcutoff & pair$rna_count > rcutoff, ]
    label <- data.frame(anno = paste0("DNA cutoff: ", min(pair_pf$dna_count), "\n", "DNA PF cells: ", nrow(pair[pair$dna_count > dcutoff, ]), "\n", 
                                      "RNA cutoff: ", min(pair_pf$rna_count), "\n",  "RNA PF cells: ", nrow(pair[pair$rna_count > rcutoff, ]), "\n", 
                                      "Both PF cells: ", nrow(pair_pf)))
    t1 <- pair %>% dplyr::filter(dna_count > 5 & rna_count > 5) %>%
    ggplot(aes(x = dna_count, y = rna_count)) + 
    geom_point(size = 0.1, color = "grey") + 
    geom_vline(xintercept = dcutoff, color = colfunc2(1), linetype="dashed") + 
    geom_hline(yintercept = rcutoff, color = colfunc2(1), linetype="dashed") + 
    geom_point(data = pair_pf, aes(x = dna_count, y = rna_count), size = 0.1, color = colfunc2(3)[[2]]) + 
    theme_Publication() + 
    geom_text(data = label, aes(x = max(pair$dna_count), y = min(pair$rna_count), hjust = 1, vjust = -0.5, label = anno), size = 2.5) + 
    ylab("Unique UMI (RNA)") + xlab("Unique fragemnts (DNA)") + 
    scale_x_log10() + 
    scale_y_log10()
    ggsave(t1, filename = paste0(prefix, "_valid_cells.png"), dpi = 300, height = 6, width = 6)
    
    ### export names
    idx <- match(pair_pf$bc, translate[, names])
    valid_bc <- translate[idx, c("atac", "rna")]
    valid_bc <- merge(x = valid_bc, y = pair_pf, by.x = names, by.y = "bc") ### atac, rna, bc, dna_count, rna_count
    write.table(valid_bc, file = paste0(prefix, "_valid_cells.xls"), row.names = F, col.names = T, sep = "\t", quote = F)
    return(valid_bc)
}

### Plot reads count for each pair of barcode, color by scatter density. Save plot and segmented cells with name prefix ###
PlotArcPair_d <- function(pair, dcutoff = 100, rcutoff = 100, prefix, names = "atac"){
    if((!all(c("bc", "dna_count", "rna_count") %in% colnames(pair)))){
        stop('colnames of input dataframe should be c("bc", "dna_count", "rna_count")')
    }
    translate <- read.table("/projects/ps-renlab2/y2xie/ps-renlab/y2xie/projects/genome_ref/arc_bc-translation.txt", header = T, row.names = 1, sep = "\t")
    translate$atac <- paste0(translate$atac, "-1")
    translate$rna <- paste0(translate$rna, "-1")
    
    pair_pf <- pair[pair$dna_count > dcutoff & pair$rna_count > rcutoff, ]
    label <- data.frame(anno = paste0("DNA cutoff: ", min(pair_pf$dna_count), "\n", "DNA PF cells: ", nrow(pair[pair$dna_count > dcutoff, ]), "\n", 
                                      "RNA cutoff: ", min(pair_pf$rna_count), "\n",  "RNA PF cells: ", nrow(pair[pair$rna_count > rcutoff, ]), "\n", 
                                      "Both PF cells: ", nrow(pair_pf)))

    filt_data <- pair %>% dplyr::filter(dna_count > 30 & rna_count > 30)
    ### assign color density 
    x <- densCols(log10(filt_data$dna_count), log10(filt_data$rna_count), colramp=colorRampPalette(viridis(12)))
    pair$dens <- min(col2rgb(x)[1,] + 1L)
    pair[pair$dna_count > 30 & pair$rna_count > 30, ]$dens <- col2rgb(x)[1,] + 1L

    t1 <- pair %>% dplyr::filter(dna_count > 5 & rna_count > 5) %>%
    ggplot(aes(x = dna_count, y = rna_count)) + 
    geom_point(size = 0.1, aes(color = dens)) + 
    geom_vline(xintercept = dcutoff, color = colfunc2(1), linetype="dashed") + 
    geom_hline(yintercept = rcutoff, color = colfunc2(1), linetype="dashed") + 
    theme_Publication() + 
    geom_text(data = label, aes(x = max(pair$dna_count), y = min(pair$rna_count), hjust = 1, vjust = -0.5, label = anno), size = 2.5) + 
    ylab("Unique UMI (RNA)") + xlab("Unique fragemnts (DNA)") + 
    scale_x_log10() + 
    scale_y_log10() + 
    theme(legend.position = "") + 
    scale_color_viridis(trans = "log")

    ggsave(t1, filename = paste0(prefix, "_valid_cells.png"), dpi = 300, height = 6, width = 6)
    
    ### export names
    idx <- match(pair_pf$bc, translate[, names])
    valid_bc <- translate[idx, c("atac", "rna")]
    valid_bc <- merge(x = valid_bc, y = pair_pf, by.x = names, by.y = "bc") ### atac, rna, bc, dna_count, rna_count
    write.table(valid_bc, file = paste0(prefix, "_valid_cells.xls"), row.names = F, col.names = T, sep = "\t", quote = F)
    return(valid_bc)
}

### Plot reads count for each pair of barcode by providing valid barcode ###
PlotArcPair2 <- function(pair, dvalid, rvalid, prefix, names = "atac"){
    if(!identical(colnames(pair), c("bc", "dna_count", "rna_count"))){
        stop('colnames of input dataframe should be c("bc", "dna_count", "rna_count")')
    }
    translate <- read.table("/projects/ps-renlab2/y2xie/ps-renlab/y2xie/projects/genome_ref/arc_bc-translation.txt", header = T, row.names = 1, sep = "\t")
    translate$atac <- paste0(translate$atac, "-1")
    translate$rna <- paste0(translate$rna, "-1")
    pair$rna_bc <- translate[match(pair$bc, translate$atac), "rna"]
    validd <- which(pair$bc %in% dvalid)
    validr <- which(pair$rna_bc %in% rvalid)
    pair_pf <- pair[intersect(validd, validr), ]
    dcutoff <- min(pair_pf$dna_count)
    rcutoff <- min(pair_pf$rna_count)
    label <- data.frame(anno = paste0("DNA cutoff: ", dcutoff, "\n", "DNA PF cells: ", length(dvalid), "\n", 
                                      "RNA cutoff: ", rcutoff, "\n",  "RNA PF cells: ", length(rvalid), "\n", 
                                      "Both PF cells: ", nrow(pair_pf)))
    t1 <- pair %>% dplyr::filter(dna_count > 5 & rna_count > 5) %>%
    ggplot(aes(x = dna_count, y = rna_count)) + 
    geom_point(size = 0.1, color = "grey") + 
    geom_vline(xintercept = dcutoff, color = colfunc2(1), linetype="dashed") + 
    geom_hline(yintercept = rcutoff, color = colfunc2(1), linetype="dashed") + 
    geom_point(data = pair_pf, aes(x = dna_count, y = rna_count), size = 0.1, color = colfunc2(3)[[2]]) + 
    theme_Publication() + 
    geom_text(data = label, aes(x = max(pair$dna_count), y = min(pair$rna_count), hjust = 1, vjust = -0.5, label = anno), size = 2.5) + 
    ylab("Unique UMI (RNA)") + xlab("Unique fragemnts (DNA)") + 
    scale_x_log10() + 
    scale_y_log10()
    ggsave(t1, filename = paste0(prefix, "_valid_cells.png"), dpi = 300, height = 6, width = 6)
    
    valid_bc <- pair_pf[, c("bc", "rna_bc", "dna_count", "rna_count")] %>% setNames(c("atac", "rna", "dna_count", "rna_count"))
    # write.table(valid_bc, file = paste0(prefix, "_valid_cells.xls"), row.names = F, col.names = T, sep = "\t", quote = F)
    return(valid_bc)
}

### Fast run Seurat ###
RunRNA <- function(obj, reduction = "pca", var = "none", batch.label = "none"){
    obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2500)
    if (var != "none") {
        obj <- ScaleData(obj, vars.to.regress = var)
    } else {
        obj <- ScaleData(obj)
    }
    obj <- RunPCA(obj, features = VariableFeatures(object = obj), verbose = F)
    if (batch.label != "none") {
        obj <- RunHarmony(obj, group.by.vars = batch.label)
        reduction <- "harmony"
    }
    obj <- RunUMAP(obj, dims = 1:30, min.dist = 0.1, seed.use=131, reduction = reduction,
                   n.components = 2L, umap.method = "uwot", n.neighbors = 25, uwot.sgd=TRUE, verbose=F)
    return(obj)
}

RunRNA2 <- function(obj, reduction = "pca", var = "none", batch.label = "none", k = 15) 
{
    if (var != "none") {
        obj <- SCTransform(obj, vars.to.regress = var, verbose = FALSE)
    }
    else {
        obj <- SCTransform(obj, verbose = FALSE)
    }
    obj <- RunPCA(obj, verbose = F)
    if (batch.label != "none") {
        obj <- RunHarmony(obj, group.by.vars = batch.label)
        reduction <- "harmony"
    }
    obj <- RunUMAP(obj, dims = 1:50, reduction = reduction, 
                   return.model = TRUE, seed.use = 921, metric = "cosine",
                   n.neighbors = k, min.dist = 0.1, n.components = 2L, umap.method = "uwot", 
                   uwot.sgd = TRUE, verbose = F)
    return(obj)
}

### Fast run Signac ###
RunDNA <- function(obj){
    obj <- RunTFIDF(obj, method = 3)
    obj <- FindTopFeatures(obj, min.cutoff = 'q15')
    obj <- RunSVD(obj)
    obj <- RunUMAP(object = obj, seed.use = 2022, reduction = 'lsi', dims = 2:30, 
                   min.dist = 0.01, n.neighbors = 15, verbose = F)
    # obj <- FindNeighbors(object = obj, reduction = 'lsi', dims = 2:30, verbose=F)
    # obj <- FindClusters(object = obj, resolution = 0.3)
    return(obj)
}

### sparse matrix: collapse columns with same colname, return sparse still ###
OP2 <- function (x) {
    nms <- colnames(x)
    uniquenms <- unique(nms)
    sparseMatrix(i = x@i + 1, j = match(nms, uniquenms)[x@j + 
        1], x = x@x, dimnames = list(rownames(x), uniquenms), 
	dims = c(nrow(x), length(uniquenms)),
        repr = "C")
}

### calculate cpm/rpkm by cluster (aggreagated cells for calculation) ###
    ### for histone, nCount_histone is used. for RNA, nCount_RNA is used. ###
ArcXPM <- function (obj_mtx, meta, method = "CPM", group.by, gname = "gene_id", gene_length = "/projects/ps-renlab2/y2xie/ps-renlab/y2xie/projects/genome_ref/mm10-2020-A_build/mm10_gene_allowlist.bed") {
    colnames(obj_mtx) <- meta[colnames(obj_mtx), group.by]
    # obj_mtx <- as(obj_mtx, "dgTMatrix")
    obj_mtx <- as(obj_mtx, "TsparseMatrix")
    obj_mtx_collapse <- OP2(obj_mtx)
    spars <- length(obj_mtx_collapse@x)/obj_mtx_collapse@Dim[1]/obj_mtx_collapse@Dim[2]
    cat(paste0("sparsity: ", spars, "\n"))
    if (spars > 0.2) {
        obj_mtx_collapse <- as(obj_mtx_collapse, "matrix")
        cat("coarse dgTMatrix into Matrix.\n")
    }
    if (method == "CPM") {
        readSums <- aggregate(meta$nCount_histone, list(meta[, group.by]), sum)
        colnames(readSums) <- c("tmp", "sums")
        readSums <- setNames(as.numeric(readSums$sums), readSums$tmp)
        cat("check readSums: ", length(names(readSums)), "\n")
        cat("check obj_mtx_collapse: ", length(colnames(obj_mtx_collapse)), 
            "\n")
        readSums <- readSums[order(match(names(readSums), colnames(obj_mtx_collapse)))]
        obj_collapse_XPM <- t(t(obj_mtx_collapse) * 10^6/readSums)
    }
    else if (method == "RPKM") {
        ### mm10: /projects/ps-renlab/y2xie/projects/genome_ref/mm10-2020-A_build/mm10_gene_allowlist.bed
        ### gname needs to be either ensembl or gene_id
        length <- read.table(file = gene_length, header = F)
        length <- length %>% setNames(c("chr", "start", "end", "strand", "ensembl", "gene_id"))
        length$length <- length$end - length$start
        len_mtx <- as.data.frame(length[match(rownames(obj_mtx_collapse), length[,gname]), ])
        readSums <- aggregate(meta$nCount_RNA, list(meta[, group.by]), sum)
        colnames(readSums) <- c("tmp", "sums")
        readSums <- setNames(as.numeric(readSums$sums), readSums$tmp)
        cat("check readSums: ", length(names(readSums)), "\n")
        cat("check obj_mtx_collapse: ", length(colnames(obj_mtx_collapse)), 
            "\n")
        readSums <- readSums[order(match(names(readSums), colnames(obj_mtx_collapse)))]
        obj_collapse_XPM <- t(t(obj_mtx_collapse) * 10^9/readSums)/len_mtx$length
    }
    return(obj_collapse_XPM)
}

ArcXPM2 <- function (obj_mtx, meta, method = "CPM", group.by, gname = "gene_id", gene_length = "/projects/ps-renlab2/y2xie/ps-renlab/y2xie/projects/genome_ref/mm10-2020-A_build/mm10_gene_allowlist.bed") {
    colnames(obj_mtx) <- meta[colnames(obj_mtx), group.by]
    obj_mtx <- as(obj_mtx, "dgTMatrix")
    uniqnms <- unique(colnames(obj_mtx))
    obj_mtx_collapse <- list()
    for (col_name in uniqnms) {
      cta <- which(colnames(obj_mtx) == col_name)
      obj_mtx_collapse[[col_name]] <- rowSums(obj_mtx[, cta])
    }
    obj_mtx_collapse <- do.call(cbind, obj_mtx_collapse)
    spars <- 1 - sum(obj_mtx_collapse != 0) / prod(dim(obj_mtx_collapse))
    cat(paste0("sparsity: ", spars, "\n"))
    if (method == "CPM") {
        readSums <- aggregate(meta$nCount_histone, list(meta[, group.by]), sum)
        colnames(readSums) <- c("tmp", "sums")
        readSums <- setNames(as.numeric(readSums$sums), readSums$tmp)
        cat("check readSums: ", length(names(readSums)), "\n")
        cat("check obj_mtx_collapse: ", length(colnames(obj_mtx_collapse)), 
            "\n")
        readSums <- readSums[order(match(names(readSums), colnames(obj_mtx_collapse)))]
        obj_collapse_XPM <- t(t(obj_mtx_collapse) * 10^6/readSums)
    }
    else if (method == "RPKM") {
        ### mm10: /projects/ps-renlab/y2xie/projects/genome_ref/mm10-2020-A_build/mm10_gene_allowlist.bed
        ### gname needs to be either ensembl or gene_id
        length <- read.table(file = gene_length, header = F)
        length <- length %>% setNames(c("chr", "start", "end", "strand", "ensembl", "gene_id"))
        length$length <- length$end - length$start
        len_mtx <- as.data.frame(length[match(rownames(obj_mtx_collapse), length[,gname]), ])
        readSums <- aggregate(meta$nCount_RNA, list(meta[, group.by]), sum)
        colnames(readSums) <- c("tmp", "sums")
        readSums <- setNames(as.numeric(readSums$sums), readSums$tmp)
        cat("check readSums: ", length(names(readSums)), "\n")
        cat("check obj_mtx_collapse: ", length(colnames(obj_mtx_collapse)), 
            "\n")
        readSums <- readSums[order(match(names(readSums), colnames(obj_mtx_collapse)))]
        obj_collapse_XPM <- t(t(obj_mtx_collapse) * 10^9/readSums)/len_mtx$length
    }
    return(obj_collapse_XPM)
}

### calculate big cor: propagate::bigcor ### 
bigcor <- function(x, y, chunks = 2000, method = "pearson"){
    cor <- list()
    nchunks <- ceiling(seq_along(1:ncol(x))/chunks)
    range_chunks <- split(seq_len(ncol(x)), nchunks)
    for(chunk in names(range_chunks)){
        x_qry <- x[,range_chunks[[chunk]]]
        cor[[chunk]] <- cor(x_qry, y, method = method)
    }
    cor <- do.call(rbind, cor)
    return(cor)
}

### calculate cluster overlap score after integration; from Yang Li ### 
cal_ovlpScore <- function(t1, t2){ 
  t1.table <- table(t1)
  t2.table <- table(t2)
  t1.pct <- apply(t1.table, 2, function(x){x/sum(x)})
  t2.pct <- apply(t2.table, 2, function(x){x/sum(x)})
  t1.labels <- colnames(t1.pct)
  t2.labels <- colnames(t2.pct)
  ovlpScore.df <- data.frame(anno1=as.character(), anno2=as.character(), ovlpScore=as.numeric())
  for(t1.label in t1.labels){
    for(t2.label in t2.labels){
      t1.pct.df <- data.frame(t1.pct[,t1.label])
      colnames(t1.pct.df) <- "t1"
      t1.pct.df$ident <- rownames(t1.pct.df)
      t2.pct.df <- data.frame(t2.pct[,t2.label])
      colnames(t2.pct.df) <- "t2"
      t2.pct.df$ident <- rownames(t2.pct.df)
      comp.df <- dplyr::full_join(t1.pct.df, t2.pct.df, by="ident")
      comp.df[is.na(comp.df)] <- 0
      comp.df$ident <- NULL
      comp.df <- t(comp.df)
      ovlpScore <- sum(apply(comp.df, 2, min))
      out <- data.frame(anno1=t1.label, anno2=t2.label, ovlpScore=ovlpScore)
      ovlpScore.df <- rbind(ovlpScore.df, out)
    }
  }
  return(ovlpScore.df)
}


### calculate enrichment of gene with aggregated cell by group matrix
### can be used to calculate QC score with specific set of gene
ssgsea <- function(X, gene_sets, alpha = 0.25, scale = T, norm = F, single = T) {
    row_names = rownames(X)
    num_genes = nrow(X)
    gene_sets = lapply(gene_sets, function(genes) {which(row_names %in% genes)})

    # Ranks for genes
    R = matrixStats::colRanks(X, preserveShape = T, ties.method = 'average')

    # Calculate enrichment score (es) for each sample (column)
    es = apply(R, 2, function(R_col) {
        gene_ranks = order(R_col, decreasing = TRUE)

        # Calc es for each gene set
        es_sample = sapply(gene_sets, function(gene_set_idx) {
            # pos: match (within the gene set)
            # neg: non-match (outside the gene set)
            indicator_pos = gene_ranks %in% gene_set_idx
            indicator_neg = !indicator_pos

            rank_alpha  = (R_col[gene_ranks] * indicator_pos) ^ alpha

            step_cdf_pos = cumsum(rank_alpha)    / sum(rank_alpha)
            step_cdf_neg = cumsum(indicator_neg) / sum(indicator_neg)

            step_cdf_diff = step_cdf_pos - step_cdf_neg

            # Normalize by gene number
            if (scale) step_cdf_diff = step_cdf_diff / num_genes

            # Use ssGSEA or not
            if (single) {
                sum(step_cdf_diff)
            } else {
                step_cdf_diff[which.max(abs(step_cdf_diff))]
            }
        })
        unlist(es_sample)
    })

    if (length(gene_sets) == 1) es = matrix(es, nrow = 1)

    # Normalize by absolute diff between max and min
    if (norm) es = es / diff(range(es))

    # Prepare output
    rownames(es) = names(gene_sets)
    colnames(es) = colnames(X)
    return(es)
}

### Find DE genes for clusters, adopted from Zemke 2023 Nature
### need to be tested
run.edgeR <- function(counts, sample, group, prefix, pairwise = FALSE){
    library("edgeR")

    ### counts: a gene by group matrix. Can obtained through AggregateExpression() 
    ### sample: a dataframe recording sample and group information. 
    ### group: name of group information in sample (group %in% colnames(sample))
    ### prefix: saving file as...
    cat("now creating DGE class for edgeR...\n")
    group.list <- factor(sample[, group])
    dge = DGEList(counts = counts, group = group.list, samples = sample)
    keep = filterByExpr(dge)
    dge = dge[keep, ,keep.lib.sizes=FALSE]
    dge = calcNormFactors(dge, 'TMM')
    design = model.matrix(~0 + group, data = dge$samples)
    colnames(design) <- make.names(gsub(group, "", colnames(design)))
    dge = estimateDisp(dge, trend.method = 'locfit', design = design)
    fit = glmFit(dge, design = design)
    
    cat("now performing glmLRT for individual clusters...\n")
    coef <- unique(group.list)
    for (f in coef){
        ### Perform an ANOVA like test for all groups
        cat(paste0("cluster: ", f, "\n"))
        one_v_rest = c(rep(-1/(length(coef)-1), times = length(coef)))
        one_v_rest[which(coef == f)] = 1  
        vrest = glmLRT(fit, contrast = one_v_rest)
        write.table(topTags(vrest, n = "Inf")$table, file = paste0(prefix, "_edgeR_", f, "_vs_all_groups.xls"),
                    sep = "\t", quote = F)
        cat(paste0("cluster: ", f, " finish.\n"))
        if(pairwise){
            cat(paste0("Now calculate pairwise test for cluster: ", f, "\n"))
            for (s in coef){
            ### Perform pairwise test for each group
                if (s == f){
                        next
                }
                one_v_one = c(rep(0, times = length(coef)))
                one_v_one[which(coef == f)] = 1  
                one_v_one[which(coef == s)] = -1 
                vtest = glmLRT(fit, contrast = one_v_one)
                write.table(topTags(vtest, n = "Inf")$table, file = paste0(prefix, "_edgeR_", f, "_vs_", s, ".xls"),
                            sep = "\t", quote = F)
            }
        }
    }
}

### add a continous option for edgeR
### also enable spline fit
edgeR_ovr_pairwise <- function (counts, sample, group, prefix, pairwise = FALSE,
        mode = c("group", "continuous"),
        spline_df = NULL,
        shape_flat_tol = log(2)) {  # min fold-change across depth to be "non-flat"
  stopifnot(is.matrix(counts) || inherits(counts, "dgCMatrix"))
  stopifnot(group %in% colnames(sample))
  mode <- match.arg(mode)

  message("Checking sample order matches counts columns...")
  if (!all(colnames(counts) %in% rownames(sample)))
      stop("All count columns must exist as rownames in 'sample'.")
  sample <- sample[colnames(counts), , drop = FALSE]

  suppressPackageStartupMessages(library(edgeR))

  message("Creating DGEList...")
  dge <- DGEList(counts = counts, samples = sample)

  if (mode == "group") {
      grp <- factor(sample[[group]])
      if (nlevels(grp) < 2L) stop("Need at least two groups for DE.")
      keep <- filterByExpr(dge, group = grp)
  } else {
      covar <- sample[[group]]
      if (!is.numeric(covar)) stop("For mode = 'continuous', 'group' must be numeric.")
      keep <- filterByExpr(dge)
  }

  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge, method = "TMM")

  # ---- design ----
  if (mode == "group") {
      design <- model.matrix(~0 + grp)
      colnames(design) <- levels(grp)
  } else if (is.null(spline_df)) {
      design <- model.matrix(~covar)
      colnames(design) <- c("Intercept", group)
  } else {
      stopifnot(is.numeric(spline_df), spline_df >= 2)
      suppressPackageStartupMessages(library(splines))
      ns_mat <- ns(covar, df = spline_df)
      design <- model.matrix(~ns_mat)
      colnames(design) <- c("Intercept",
                            paste0(group, "_ns", seq_len(ncol(ns_mat))))
  }

  dge <- estimateDisp(dge, design = design, robust = TRUE)
  fit <- glmQLFit(dge, design = design, robust = TRUE)

  # ---- shape classifier (used only in spline mode) ----
  classify_spline_shape <- function(curve, flat_tol) {
      if (diff(range(curve)) < flat_tol) return("flat")
      d <- diff(curve)
      # Ignore derivative changes smaller than 10% of the peak slope:
      # these are noise-level wiggles, not real extrema.
      slope_tol <- 0.10 * max(abs(d))
      d_signed  <- ifelse(abs(d) < slope_tol, 0, sign(d))
      d_signed  <- d_signed[d_signed != 0]
      if (length(d_signed) == 0) return("flat")
      sc <- sum(abs(diff(d_signed)) == 2)   # only count true +/- flips
      if (sc == 0) return(if (d_signed[1] > 0) "monotone_up" else "monotone_down")
      if (sc == 1) {
          endpoint_mean <- mean(curve[c(1, length(curve))])
          up   <- max(curve) - endpoint_mean
          down <- endpoint_mean - min(curve)
          return(if (up > down) "peaked" else "U")
      }
      "multi_modal"
  }

  run_contrast <- function(w, outname) {
      qlf <- glmQLFTest(fit, contrast = w)
      tab <- topTags(qlf, n = Inf)$table
      write.table(tab, file = paste0(prefix, "_edgeR_", outname, ".tsv"),
                  sep = "\t", quote = FALSE, row.names = TRUE)
  }

  # ---- tests ----
  if (mode == "group") {
      levs <- levels(grp); k <- length(levs)
      message("Running one-vs-rest contrasts...")
      for (f in levs) {
          w <- rep(-1/(k - 1), k); names(w) <- levs; w[f] <- 1
          run_contrast(w, paste0(f, "_vs_all"))
      }
      if (pairwise) {
          message("Running pairwise contrasts...")
          for (i in seq_len(k)) for (j in seq_len(k)) {
              if (i >= j) next
              f <- levs[i]; s <- levs[j]
              w <- setNames(rep(0, k), levs); w[f] <- 1; w[s] <- -1
              run_contrast(w, paste0(f, "_vs_", s))
          }
      }
  } else if (is.null(spline_df)) {
      message("Running continuous variable test (linear)...")
      qlf <- glmQLFTest(fit, coef = group)
      tab <- topTags(qlf, n = Inf)$table
      write.table(tab,
                  file = paste0(prefix, "_edgeR_", group, "_continuous.tsv"),
                  sep = "\t", quote = FALSE, row.names = TRUE)
  } else {
      message(sprintf("Running continuous variable test (natural spline, df=%d)...",
                      spline_df))
      spline_cols <- 2:ncol(design)
      qlf <- glmQLFTest(fit, coef = spline_cols)
      tab <- topTags(qlf, n = Inf)$table

      message("Classifying per-gene spline shapes...")
      grid     <- seq(min(covar), max(covar), length.out = 200)
      ns_grid  <- predict(ns_mat, grid)                         # same knots as the fit
      coef_sp  <- fit$coefficients[, spline_cols, drop = FALSE] # genes x spline_df
      curves   <- coef_sp %*% t(ns_grid)                         # genes x grid
      shapes   <- apply(curves, 1, classify_spline_shape,
                        flat_tol = shape_flat_tol)

      # peak / trough location in depth units (useful for interpretation)
      peak_depth <- grid[apply(curves, 1, which.max)]
      trough_depth <- grid[apply(curves, 1, which.min)]
      depth_range <- apply(curves, 1, function(v) diff(range(v)))

      tab$shape        <- shapes[rownames(tab)]
      tab$peak_depth   <- peak_depth[rownames(tab)]
      tab$trough_depth <- trough_depth[rownames(tab)]
      tab$depth_range  <- depth_range[rownames(tab)]

      write.table(tab,
                  file = paste0(prefix, "_edgeR_", group, "_spline", spline_df, ".tsv"),
                  sep = "\t", quote = FALSE, row.names = TRUE)
  }

  message("Done.")
}


### Calculate tissue / cell type specificity of marker genes
### adopted from hicat (https://rdrr.io/github/AllenInstitute/scrattch.hicat/src/R/util.R#sym-calc_tau)
calc_tau <- function(m) {  
  row_maxes <- apply(m, 1, max)
  
  m <- m / row_maxes
  tau <- Matrix::rowSums(1 - m) / (ncol(m) - 1)
  tau[is.na(tau)] <- 0
  return(tau)
}
