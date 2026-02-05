suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(FNN))
suppressPackageStartupMessages(library(RSpectra))

run_centroid_prediction <- function(ref_emb, query_emb, rmeta, coord_cols = c("x", "y"), k = 50) {  
  message("2. Finding ", k, " nearest neighbors...")
  # KNN Search
  knn_res <- FNN::get.knnx(data = ref_emb, query = query_emb, k = k)
  idx_mat <- knn_res$nn.index
  
  message("3. Calculating centroids...")
  uniform_weight <- 1 / k
  n_query <- nrow(query_emb)
  n_ref <- nrow(ref_emb)

  sparse_i <- rep(1:n_query, k)
  sparse_j <- as.vector(idx_mat)
  sparse_x <- rep(uniform_weight, length(sparse_j))
  
  TransformMat <- sparseMatrix(i = sparse_i, j = sparse_j, x = sparse_x, dims = c(n_query, n_ref))
  ref_coords_mat <- as.matrix(rmeta[, coord_cols])
  
  pred_coords <- TransformMat %*% ref_coords_mat
  
  # --- Output ---
  results <- as.data.frame(as.matrix(pred_coords))
  colnames(results) <- paste0("pred_", coord_cols)
  results$cell_id <- rownames(query_emb)
  
  return(results)
}

run_merfish_classification <- function(ref_emb, query_emb, rmeta, ref_label_col, k = 50, sd = 1, prob_threshold = 0.4) {
  ref_labels <- rmeta[[ref_label_col]]
  names(ref_labels) <- rownames(rmeta)
  
  # --- Step 3: Neighborhood-based Classifier ---
  message("Running neighborhood-based classifier...")
  
  # 1. Identify nearest s (k=50) cells in the joint embedding space 
  # dist(c, a_i) calculation 
  knn_res <- FNN::get.knnx(data = ref_emb, query = query_emb, k = k)
  
  dist_matrix <- knn_res$nn.dist   # The Euclidean distances
  index_matrix <- knn_res$nn.index # The indices of the neighbors in ref_emb
  
  # 2. Calculate Weighted Distance D_{c,i}
  # "D_{c,i} = 1 - dist(c, a_i) / dist(c, a_s)" 
  # dist(c, a_s) is the distance to the 50th neighbor (last column of dist_matrix)
  dist_cs <- dist_matrix[, k]
  
  # Create a matrix for D_ci (broadcasting the division by the k-th distance)
  D_ci <- 1 - (dist_matrix / dist_cs)
  
  # 3. Convert to Weighted Similarity S_{c,i}
  # "S_{c,i} = 1 - e^(-D_{c,i} / (2/sd)^2)" where sd=1 
  denominator <- (2 / sd)^2
  S_ci <- 1 - exp(-D_ci / denominator)
  
  # 4. Normalize Similarity W_{c,i}
  # "W_{c,i} = S_{c,i} / Sum(S_{c,j})" 
  row_sums <- rowSums(S_ci)
  W_ci <- S_ci / row_sums
  
  # --- Step 4: Prediction ---
  # "We then compute label predictions for query cells as P^l = WL" 
  
  unique_clusters <- unique(ref_labels)
  prob_matrix <- matrix(0, nrow = nrow(query_emb), ncol = length(unique_clusters))
  colnames(prob_matrix) <- unique_clusters
  rownames(prob_matrix) <- rownames(query_emb)
  
  # Map neighbor indices to actual labels
  # Reshape vector of labels into matrix matching the KNN structure
  neighbor_labels <- matrix(ref_labels[index_matrix], ncol = k)
  
  # Calculate probabilities (Weighted vote)
  for (cluster in unique_clusters) {
    # Binary matrix: 1 if neighbor is in cluster, 0 otherwise [cite: 245]
    is_cluster <- (neighbor_labels == cluster) * 1
    
    # Sum weights for this cluster
    prob_matrix[, cluster] <- rowSums(W_ci * is_cluster)
  }
  
  # --- Step 5: Assignment and Filtering ---
  # "Each cell... was assigned to the cluster label that had the maximum probability." [cite: 250]
  max_probs <- apply(prob_matrix, 1, max)
  predicted_labels <- colnames(prob_matrix)[apply(prob_matrix, 1, which.max)]
  
  # "Cells with the maximum assignment probability less than 0.4 was removed" 
  final_labels <- predicted_labels
  final_labels[max_probs < prob_threshold] <- NA
  
  # Return results
  results <- data.frame(
    cell_id = rownames(query_emb),
    predicted_cluster = final_labels,
    max_probability = max_probs,
    stringsAsFactors = FALSE
  )
  
  return(results)
}

getNN <- function(ref_emb, query_emb, k = 50, sd = 1, weight_scheme = c("fang", "uniform")) {
  weight_scheme <- match.arg(weight_scheme)
  
  # Validation: Ensure dimensions match
  if (ncol(ref_emb) != ncol(query_emb)) {
    stop("Error: Reference and Query embeddings must have the same number of dimensions (columns).")
  }
  
  # --- Step 1: Find Neighbors ---
  message("1. Finding ", k, " nearest neighbors...")
  knn_res <- FNN::get.knnx(data = ref_emb, query = query_emb, k = k)
  
  idx_mat <- knn_res$nn.index # Integer indices
  dist_mat <- knn_res$nn.dist # Euclidean distances
  
  # --- Step 2: Calculate Weights ---
  message("2. Calculating ", weight_scheme, " weights...")
  
  if (weight_scheme == "fang") {
    # Fang et al. (2022) Exponential Decay Formula
    dist_k <- dist_mat[, k] # Distance to k-th neighbor
    
    # Avoid division by zero
    dist_k[dist_k == 0] <- 1e-6 
    
    D_ci <- 1 - (dist_mat / dist_k)
    denom <- (2 / sd)^2
    S_ci <- 1 - exp(-D_ci / denom)
    
    # Normalize rows to sum to 1
    row_sums <- rowSums(S_ci)
    row_sums[row_sums == 0] <- 1
    weights_mat <- S_ci / row_sums
    
  } else {
    # Uniform Weights (1/k)
    weights_mat <- matrix(1/k, nrow = nrow(dist_mat), ncol = ncol(dist_mat))
  }
  
  # --- Step 3: Format Output ---
  # Map integer indices back to Reference Cell IDs
  if (!is.null(rownames(ref_emb))) {
    ref_names <- rownames(ref_emb)
    neighbor_names <- matrix(ref_names[idx_mat], nrow = nrow(idx_mat), ncol = ncol(idx_mat))
  } else {
    warning("Reference embeddings do not have rownames. Returning integer indices.")
    neighbor_names <- idx_mat
  }
  
  # Assign Query Cell IDs to rows
  if (!is.null(rownames(query_emb))) {
    rownames(neighbor_names) <- rownames(query_emb)
    rownames(weights_mat) <- rownames(query_emb)
    rownames(dist_mat) <- rownames(query_emb)
  }
  
  return(list(
    neighbors = neighbor_names,
    weights = weights_mat,
    distances = dist_mat
  ))
}

calculate_label_consistency <- function(neighbor_res, ref_labels, query_labels) {
  
  # 1. inputs setup
  neighbors_mat <- neighbor_res$neighbors
  weights_mat <- neighbor_res$weights
  query_ids <- rownames(neighbors_mat)
  
  # Ensure query_labels are in the same order as the result matrix
  if (!all(query_ids %in% names(query_labels))) {
    stop("Some query cells in the result are missing from query_labels.")
  }
  current_query_labels <- query_labels[query_ids]
  
  # 2. Map Neighbor IDs to Neighbor Labels
  # We create a matrix of labels with the same shape as neighbors_mat
  # This relies on ref_labels being a named vector
  neighbor_labels_mat <- matrix(ref_labels[neighbors_mat], 
                                nrow = nrow(neighbors_mat), 
                                ncol = ncol(neighbors_mat))
  
  # 3. Compare Neighbor Labels to Query Labels
  # Create a logical matrix: TRUE if neighbor matches query, FALSE otherwise
  is_match_mat <- neighbor_labels_mat == current_query_labels
  
  # Handle NAs (if some neighbors didn't have labels)
  is_match_mat[is.na(is_match_mat)] <- FALSE
  
  # 4. Calculate Scores
  # Simple Purity: Fraction of neighbors that match
  raw_purity <- rowMeans(is_match_mat)
  
  # Weighted Purity: Weighted sum of matches (using Fang weights)
  # rowSums(Weight * 1 if match, 0 if not)
  weighted_purity <- rowSums(weights_mat * is_match_mat)
  
  # 5. Compile Results
  results <- data.frame(
    cell_id = query_ids,
    cell_type = current_query_labels,
    purity_score = raw_purity,
    weighted_score = weighted_purity,
    stringsAsFactors = FALSE
  )
  
  return(results)
}

#' Propagate Labels with Distance Weighting
#'
#' @param obj Seurat object.
#' @param label_col String. Metadata column with known labels.
#' @param image_col String. Metadata column for Image ID.
#' @param coord_cols Character vector. Coordinates.
#' @param k Integer. Neighbors.
#' @param layer_order Character vector. Biological order for depth calculation.
#' @param weighted Boolean. If TRUE, weights neighbors by distance (Gaussian decay).
#' @param sigma_scale Numeric. Controls weight decay. Higher = smoother/wider influence. 
#'                    Default 1.0 means sigma is the mean distance to the k-th neighbor.
#'
#' @return Seurat object with propagated labels, weighted scores, and depth.
propagate_layer_labels_weighted <- function(obj, label_col = "layer_label", 
                                            image_col = "image_id",
                                            coord_cols = c("x", "y"), k = 20,
                                            layer_order = NULL,
                                            weighted = TRUE, sigma_scale = 1.0) {
  
  meta <- obj@meta.data
  
  # 1. Coordinate Validation
  if (!all(coord_cols %in% colnames(meta))) {
    coords <- GetTissueCoordinates(obj)
    colnames(coords)[1:2] <- c("x", "y")
    meta <- cbind(meta, coords[rownames(meta), ])
  }
  
  # Setup result vectors
  all_cells <- rownames(meta)
  res_layer <- as.character(meta[[label_col]])
  names(res_layer) <- all_cells
  
  res_score <- rep(NA, length(all_cells))
  names(res_score) <- all_cells
  res_score[!is.na(res_layer) & res_layer != "Unknown"] <- 1.0
  
  res_depth <- rep(NA, length(all_cells))
  names(res_depth) <- all_cells
  
  # Anchor depths
  if (!is.null(layer_order)) {
    anchor_depths <- match(res_layer, layer_order)
    res_depth[!is.na(res_layer)] <- anchor_depths[!is.na(res_layer)]
  }
  
  # 2. Iterate Images
  unique_images <- unique(meta[[image_col]])
  message("Processing ", length(unique_images), " images...")
  
  for (img in unique_images) {
    img_meta <- meta[meta[[image_col]] == img, ]
    
    is_unknown <- is.na(img_meta[[label_col]]) | img_meta[[label_col]] == "Unknown"
    ref_cells <- img_meta[!is_unknown, ]
    query_cells <- img_meta[is_unknown, ]
    
    if (nrow(query_cells) == 0 || nrow(ref_cells) == 0) next
    
    # 3. Spatial KNN
    current_k <- min(k, nrow(ref_cells))
    knn_res <- FNN::get.knnx(data = as.matrix(ref_cells[, coord_cols]), 
                             query = as.matrix(query_cells[, coord_cols]), 
                             k = current_k)
    
    idx_mat <- knn_res$nn.index
    dist_mat <- knn_res$nn.dist
    ref_labels <- as.character(ref_cells[[label_col]])
    
    # 4. Calculate Weights (Gaussian Kernel)
    if (weighted) {
      # Sigma: Adaptive kernel width based on local density
      # We set sigma = median distance to k-th neighbor * scale
      # This ensures the kernel scales if one FOV is denser than another
      local_sigma <- median(dist_mat[, current_k]) * sigma_scale
      
      # Gaussian Weight: exp(-d^2 / (2*sigma^2))
      # Adding small epsilon to sigma to avoid division by zero
      weights_mat <- exp(- (dist_mat^2) / (2 * (local_sigma + 1e-6)^2))
    } else {
      weights_mat <- matrix(1, nrow = nrow(dist_mat), ncol = ncol(dist_mat))
    }
    
    # Normalize weights so rows sum to 1 (for interpretation as probability)
    row_sums <- rowSums(weights_mat)
    weights_mat <- weights_mat / row_sums
    
    # 5. Apply Logic (Weighted Vote)
    img_results <- lapply(1:nrow(query_cells), function(i) {
      indices <- idx_mat[i, ]
      current_weights <- weights_mat[i, ]
      current_labels <- ref_labels[indices]
      
      # A. Weighted Vote
      # Sum weights for each unique label
      unique_labs <- unique(current_labels)
      weighted_counts <- sapply(unique_labs, function(lbl) {
        sum(current_weights[current_labels == lbl])
      })
      
      # Winner is the one with highest Total Weight
      winner <- names(which.max(weighted_counts))
      
      # B. Prediction Score (Weighted Confidence)
      # This effectively represents the probability P(Label | Neighbors)
      score <- max(weighted_counts) # Since we normalized weights to 1, this is 0.0-1.0
      
      # C. Continuous Depth (Weighted Mean)
      depth_val <- NA
      if (!is.null(layer_order)) {
        numeric_labels <- match(current_labels, layer_order)
        # We need to filter out NAs but keep weights aligned
        valid_mask <- !is.na(numeric_labels)
        
        if (any(valid_mask)) {
          vals <- numeric_labels[valid_mask]
          wts <- current_weights[valid_mask]
          # Weighted Average
          depth_val <- sum(vals * wts) / sum(wts)
        }
      }
      
      return(c(winner, score, depth_val))
    })
    
    # Unpack results
    img_results_mat <- do.call(rbind, img_results)
    
    query_ids <- rownames(query_cells)
    res_layer[query_ids] <- img_results_mat[, 1]
    res_score[query_ids] <- as.numeric(img_results_mat[, 2])
    res_depth[query_ids] <- as.numeric(img_results_mat[, 3])
  }
  
  # 6. Store
  obj$propagated_layer <- res_layer
  obj$prediction_score <- res_score
  obj$continuous_depth <- res_depth
  
  return(obj)
}

### recap of ENVI analysis but with CCA embedding
run_diffusion_maps_R <- function(
  data_mat,
  n_components = 10,
  knn = 30,
  alpha = 0
) {
  # data_mat: cells × dimensions (PCA / integrated embedding)
  
  N <- nrow(data_mat)
  stopifnot(knn < N)

  ## -----------------------------
  ## 1. kNN graph
  ## -----------------------------
  knn_res <- get.knn(data_mat, k = knn)
  idx <- knn_res$nn.index
  dist <- knn_res$nn.dist

  ## -----------------------------
  ## 2. Adaptive bandwidth (self-tuning)
  ## -----------------------------
  adaptive_k <- floor(knn / 3)
  sigma <- dist[, adaptive_k]

  ## -----------------------------
  ## 3. Construct anisotropic kernel
  ## -----------------------------
  i <- rep(seq_len(N), knn)
  j <- as.vector(idx)
  d <- as.vector(dist)

  d_scaled <- d / sigma[i]
  W <- sparseMatrix(
    i = i,
    j = j,
    x = exp(-d_scaled),
    dims = c(N, N)
  )

  # Symmetrize
  K <- W + t(W)

  ## -----------------------------
  ## 4. Alpha normalization
  ## -----------------------------
  D <- rowSums(K)

  if (alpha > 0) {
    D_alpha <- D^(-alpha)
    D_alpha[!is.finite(D_alpha)] <- 0
    Dmat <- Diagonal(x = D_alpha)
    K <- Dmat %*% K %*% Dmat
    D <- rowSums(K)
  }

  ## -----------------------------
  ## 5. Markov normalization
  ## -----------------------------
  D_inv <- 1 / D
  D_inv[!is.finite(D_inv)] <- 0
  Tmat <- Diagonal(x = D_inv) %*% K

  ## -----------------------------
  ## 6. Eigendecomposition
  ## -----------------------------
  eig <- eigs(
    Tmat,
    k = n_components,
    which = "LR",
    tol = 1e-4,
    maxit = 1000
  )

  values <- Re(eig$values)
  vectors <- Re(eig$vectors)

  # sort eigenvalues descending
  ord <- order(values, decreasing = TRUE)
  values <- values[ord]
  vectors <- vectors[, ord]

  # normalize eigenvectors
  vectors <- apply(vectors, 2, function(v) v / sqrt(sum(v^2)))

  ## Drop DC0
  list(
    eigenvalues = values[-1],
    diffusion_components = vectors[, -1]
  )
}
