require(dplyr)
require(tidyr)
require(ggplot2)
require(readr)
require(parallel)
require(Seurat)
require(WGCNA)
require(pbmcapply)



# Negative Binomial Distribution Model 
MLE_NB <- function(obs, ref, count_obs_size,top_cell_type_for_fine_tune = NULL,dispersion_min = 1e-5,cell_id=NULL,top_correlated_ref_cell_names=NULL) {
    
    # calculate the mean and dispersion of each gene
    means <- rowMeans(ref)
    vars <- apply(ref, 1, var)
    
    # handle the dispersion calculation, ensure it is positive  
    dispersions <- pmax((vars - means) / means^2, dispersion_min)  # set the minimum dispersion
    dispersions[is.na(dispersions)] <- dispersion_min
    
    # calculate the log-likelihood of the negative binomial distribution
    G <- length(obs)
    K <- ncol(ref)
    # 1) replicate into G×K matrices
    counts_mat <- matrix(obs, nrow = G, ncol = K)
    mu_mat     <- ref * count_obs_size
    size_mat   <- matrix(1 / dispersions,      nrow = G, ncol = K, byrow = TRUE)
    # 2) compute one big matrix of log-probs
    logp_mat <- dnbinom(counts_mat, size = size_mat, mu = mu_mat, log = TRUE)
    # handle the NaN value
    logp_mat[is.na(logp_mat)] <- -Inf
    # handle the Inf value
    logp_mat[is.infinite(logp_mat)] <- 0
    # 3) sum over genes (rows) to get one log-likelihood per column
    log_likelihoods <- colSums(logp_mat,na.rm = T)
    # if all the log-likelihoods are -Inf, return the default value
    if(all(is.infinite(log_likelihoods))|all(is.na(log_likelihoods))) {
      return(data.frame(
        best_prob = 0,
        cell_id = cell_id,
        ref_cell_type = top_correlated_ref_cell_names[1],
        top_cor_celltype = top_correlated_ref_cell_names[1],
        ci_95_matches = 0,
        log_likelihood = -Inf,
        best_cell_types = paste0(top_correlated_ref_cell_names,collapse = ';'),
        stringsAsFactors = FALSE
      ))
    }
    
    # calculate the probability
    max_loglik <- max(log_likelihoods)
    likelihood_ratios <- exp(log_likelihoods - max_loglik)
    probs <- likelihood_ratios / sum(likelihood_ratios)
    
    # find the best match
    max_idx <- which.max(log_likelihoods)
    best_prob <- probs[max_idx]
    
    # calculate the confidence interval  
    cum_probs <- cumsum(probs)
    lower_idx <- which(cum_probs >= 0.025)[1]
    upper_idx <- which(cum_probs >= 0.975)[1]
    ci_95 <- upper_idx - lower_idx + 1
    
    # extract the most fit cell types
    candidate_finetune_cell_number = max(ci_95, top_cell_type_for_fine_tune)
    best_cell_types <- colnames(ref)[order(probs, decreasing = TRUE)][1:candidate_finetune_cell_number]
    
    # return a single row data frame
    return(data.frame(
      best_prob = best_prob,
      cell_id = cell_id,
      ref_cell_type = colnames(ref)[max_idx],
      top_cor_celltype = top_correlated_ref_cell_names[1],
      ci_95_matches = ci_95,
      log_likelihood = log_likelihoods[max_idx],
      best_cell_types = paste0(best_cell_types, collapse = ';'),
      stringsAsFactors = FALSE
    ))
}

# dynamic marker gene selection
Select_Dynamic_Marker_Genes <- function(obs, ref, top_candidates, fine_tune_marker_select_num = 10) {
    # 1. extract the subset of data
    ref_subset <- ref[, top_candidates]
    ref_subset <- ref_subset[!grepl("MT-|^MT|^RN|^LINC|^MIR|^AC|^ENSG|^RPL|^RPS|\\.|-IT|-AS",rownames(ref_subset)),]
    
    # 2. calculate the statistics of each gene
    means <- rowMeans(ref_subset)
    sds <- rowSds(as.matrix(ref_subset))
    cv <- sds / means
    
    # 3. calculate the expression difference of each gene in different cell types
    # use combn to generate all the cell type pairs
    type_pairs <- combn(1:ncol(ref_subset), 2)
    
    # calculate all the differences at once
    diff_matrix <- abs(ref_subset[, type_pairs[1,]] - ref_subset[, type_pairs[2,]])
    
    # 4. calculate the comprehensive score of each gene
    max_diff <- apply(diff_matrix, 1, max)
    score <- max_diff * cv
    
    # 5. select the candidate genes
    candidate_genes <- names(score[score > quantile(score, 0.95,na.rm=T)])
    
    # 6. ensure the candidate genes exist in obs
    candidate_genes <- candidate_genes[candidate_genes %in% names(obs)]
    
    # 7. if the candidate genes are too many, select the highest score ones
    if(length(candidate_genes) > fine_tune_marker_select_num) {
      candidate_genes <- names(sort(score[candidate_genes], decreasing = TRUE)[1:fine_tune_marker_select_num])
    }
    
    return(candidate_genes)
}


# get_cell_type_by_MLE_without_mixture
# cell_id: the cell id to be annotated
# cell_expression_counts_full: the expression counts of the cell
# ref_expr: the expression counts of the reference, it should be normalized by total sum of the reference cell. We will normalize it in the function using single cell size factor.
# top_correlated_ref_cell_names: the top correlated reference cell names, performed in last step of MLE_cor.
# common_vars: the common variables, performed in last step of MLE_cor.
# iter_num_max: the maximum number of iterations, default is 3.
# fine_tune_marker_select_num: the number of fine-tune marker genes, default is 10.
# top_cell_type_for_fine_tune: the number of top cell types to be fine-tuned, default is 5.

Get_Cell_Type_by_MLE <- function(
    cell_id, 
    cell_expression_counts_full, 
    ref_expr, 
    count_obs,
    top_correlated_ref_cell_names, 
    common_vars,
    iter_num_max = 3,
    fine_tune_marker_select_num = 10,
    top_cell_type_for_fine_tune = 5,
    probability_threshold = 0.99,
    ci_95_matches_threshold = 5,
    dispersion_min = 1e-5
) {
  # 
  count_obs = sum(cell_expression_counts_full)

  ref_expr = ref_expr[, top_correlated_ref_cell_names]

  # initial MLE result
  iter_num = 0 
  iterative_result <- MLE_NB(
    obs = cell_expression_counts_full[common_vars],
    ref = ref_expr[common_vars, top_correlated_ref_cell_names],
    count_obs_size = count_obs,
    top_cell_type_for_fine_tune = top_cell_type_for_fine_tune,
    dispersion_min = dispersion_min,
    cell_id = cell_id,
    top_correlated_ref_cell_names = top_correlated_ref_cell_names
  ) %>% dplyr::mutate(iter_num = iter_num)
  
  # if the confidence interval is too wide, perform the iteration optimization
  
  while(iterative_result$best_prob < probability_threshold & (iterative_result$ci_95_matches > ci_95_matches_threshold) & iter_num < iter_num_max) {
    # select the marker genes
    top_candidates_iter = stringr::str_split(iterative_result$best_cell_types,';')[[1]] %>% 
      unlist() %>% 
      as.character()
    
    markers <- Select_Dynamic_Marker_Genes(
      obs = cell_expression_counts_full, 
      ref = ref_expr[,top_candidates_iter], 
      top_candidates = top_candidates_iter,
      fine_tune_marker_select_num = fine_tune_marker_select_num
    )
    
    iter_num = iter_num + 1 
    if(sum(markers %in% names(cell_expression_counts_full)) > 0) {
      markers <- markers[markers %in% names(cell_expression_counts_full)]
      # use the marker genes to recalculate
      if(length(markers) > 0) {
        iterative_result <- MLE_NB(
          obs = cell_expression_counts_full[markers],
          ref = ref_expr[markers,top_candidates_iter],
          count_obs_size = count_obs,
          top_cell_type_for_fine_tune = top_cell_type_for_fine_tune,
          dispersion_min = dispersion_min,
          cell_id = cell_id,
          top_correlated_ref_cell_names = top_correlated_ref_cell_names
        ) %>% dplyr::mutate(iter_num = iter_num)
      }
    } else {
      return(iterative_result)
    }
  }
  return(iterative_result)
}


# function to perform pre-processing

Preprocess_for_MLE <- function(input_seurat_obj,reference_data_matrix,nfeatures_for_variable_in_obj=3000,nfeatures_for_variable_in_ref=40,chunk_size=5000,ncores=1){

    # input_seurat_obj <- seurat_obj
    # reference_data_matrix <- ref_mtx
    # nfeatures_for_variable_in_obj <- 3000
    # nfeatures_for_variable_in_ref <- 40
    # chunk_size <- 5000
    # ncores <- 3

  if(sum(grepl("RNA",input_seurat_obj@assays %>% names())) == 0){
    stop("The input seurat object should have an RNA assay")
  }

  cat(as.vector(date()),'preprocess start \n')
  seurat_obj <- input_seurat_obj
  ref_mtx <- reference_data_matrix

  GetAssayData(seurat_obj, slot = 'count',assay = 'RNA') -> count_data
  # count_data_norm <- sweep(count_data, 2, colSums(count_data), FUN = "/")
  count_data_norm <- t(t(count_data)/colSums(count_data))

  # combine variable features to make a coherent highly variable gene set
  # compute variable features in the reference dataset 
  rowVars(ref_mtx) -> vars_per_gene
  # select top 5% variable gene in the reference set 
  variable_in_ref <- rownames(ref_mtx)[vars_per_gene>quantile(vars_per_gene,seq(0,1,0.01))[95]] 
  # remove the ribosomal and mitochondrial and LINC etc genes 
  variable_in_ref <- variable_in_ref[!grepl("^RNU|^RF|^SNOR|^AC|^RPL|^RPS|^MT-|^SNOR|^MIR|^MT|^LINC|-AS|^RN\\d|^RNA\\d|-IT",variable_in_ref)]
  # quickly calculate top-N DEG for each group
  t(ref_mtx) %>% scale() %>% t() -> scaled_ref_mtx
  scaled_ref_mtx <- scaled_ref_mtx[!grepl('^RNU|^RF|^SNOR|^AC|^RPL|^RPS|^MT-|^SNOR|^MIR|^MT|^LINC|-AS|^RN\\d|^RNA\\d|-IT',rownames(scaled_ref_mtx)),]
  top_per_column <- lapply(colnames(scaled_ref_mtx), function(id,N=nfeatures_for_variable_in_ref) {
    # get the top N values of each column
    x=scaled_ref_mtx[,id]
    y = order(x, decreasing = TRUE)
    top_indices <- y[1:N]
    top_names = rownames(scaled_ref_mtx)[top_indices]
    top_Z = x[top_indices]
    # return the top N values
    return(data.frame(gene=top_names,Z=top_Z,ref_id=id))
  }) %>% bind_rows()

  DefaultAssay(seurat_obj) <- 'RNA'
  seurat_obj <- seurat_obj %>% NormalizeData() %>% FindVariableFeatures(nfeatures = nfeatures_for_variable_in_obj) 
  VariableFeatures(seurat_obj) -> variable_in_obj 
  # combine the variable features into one 
  variable_features_combine <- variable_in_obj[variable_in_obj %in% c(top_per_column$gene,variable_in_ref)] %>% unique()
  # clean-up those features not in query dataset 
  common_vars <- variable_features_combine
  common_vars <- intersect(rownames(ref_mtx), common_vars)
  common_vars <- intersect(rownames(seurat_obj[['RNA']]), common_vars)
  # clean up useless ones 
  common_vars <- common_vars[!grepl('^RNU|^RF|^SNOR|^AC|^RPL|^RPS|^MT-|^SNOR|^MIR|^MT|^LINC|-AS|^RN\\d|^RNA\\d|-IT|^\\d',common_vars)]

  # use correlation to find the most similar cells
  if(ncol(count_data_norm) > 5000){
    # subset the thing to multiple chunks 
    n_chunks <- ceiling(ncol(count_data_norm) / chunk_size)
    pbmcapply::pbmclapply(1:n_chunks, function(i){
      start_idx <- (i - 1) * chunk_size + 1
      end_idx <- min(i * chunk_size, ncol(count_data_norm))
      chunk_data <- count_data_norm[common_vars, start_idx:end_idx]
      WGCNA::cor(chunk_data[common_vars,] %>% log1p, ref_mtx[common_vars,] %>% log1p, method = 'spearman')
    }, mc.cores = ncores) -> cor_list
    do.call(rbind, cor_list) -> cor_mtx
  }else{
    WGCNA::cor(count_data_norm[common_vars,] %>% log1p, ref_mtx[common_vars,] %>% log1p, method = 'spearman') -> cor_mtx
  }

  cat(as.vector(date()),'preprocess end \n')
  return(list(common_vars=common_vars,cor_mtx=cor_mtx,count_data=count_data,ref_mtx=ref_mtx))
}

# predict the cell type by MLE
# input_seurat_obj: the seurat object
# reference_data_matrix: the reference data matrix, it needs to be normalized.
# If it is a count matrix, normalize it by the sequencing sum size. 
# If it is normalized by sum, then no further normalization is performed. 
# If it is a log1p matrix, then the matrix should be converted to a count matrix and then normalized by the sequencing sum size.
# preprocess_result: the preprocess result, if not provided, the function will perform the pre-processing
# ncores: the number of cores to use
# top_correlation_percentile_for_MLE: the percentile of the correlation to use for the MLE, between 0 and 1, default is 0.95
# return: a data frame with the cell type prediction result

Predict_MLE <- function(input_seurat_obj,reference_data_matrix,preprocess_result=NULL,top_correlation_percentile_for_MLE = 0.95,ncores=3,MLE_hit_probability_threshold = 0.99,MLE_hit_ci_95_matches_threshold = 5,compute_cells_limit=NULL,dispersion_min = 1e-5){

  if(ncores > 1){
    require(pbmcapply)
  }

  if(is.null(preprocess_result)){
    # perform the pre-processing
    pre_process_result <- Preprocess_for_MLE(input_seurat_obj,reference_data_matrix) 
  }
  
  common_vars <- preprocess_result$common_vars
  cor_mtx <- preprocess_result$cor_mtx
  count_data <- preprocess_result$count_data
  ref_mtx <- preprocess_result$ref_mtx

  cat(as.vector(date()),'start \n')

  if(!is.null(compute_cells_limit)){
    count_data <- count_data[,1:compute_cells_limit]
  }

  cat (length(colnames(count_data)), "cells to be annotated \n")

  if(ncores > 1){
    results <- pbmcapply::pbmclapply(colnames(count_data),
        function(cell_id, common_vars, cor_mtx, count_data, ref_mtx, 
           top_correlation_percentile_for_MLE, MLE_hit_probability_threshold, 
           MLE_hit_ci_95_matches_threshold, dispersion_min) {
        cell_expression_counts_full <- count_data[,cell_id] %>% as.vector()
        names(cell_expression_counts_full) <- rownames(count_data)
        cor_vector <- cor_mtx[cell_id,]
        test_use_ids <- min(10,max(2,sum(cor_vector>quantile(cor_vector,top_correlation_percentile_for_MLE,na.rm = T))))
        top_correlated_ref_cell_names <- colnames(cor_mtx)[order(cor_mtx[cell_id,], decreasing = TRUE)][1:test_use_ids]
        # compute the MLE result
        df1 <- Get_Cell_Type_by_MLE(cell_id = cell_id,cell_expression_counts_full = cell_expression_counts_full,ref_expr = ref_mtx,count_obs = sum(cell_expression_counts_full),top_correlated_ref_cell_names = top_correlated_ref_cell_names,common_vars = common_vars,iter_num_max = 3,   fine_tune_marker_select_num = 10,top_cell_type_for_fine_tune = 5,probability_threshold = MLE_hit_probability_threshold,ci_95_matches_threshold = MLE_hit_ci_95_matches_threshold,dispersion_min = dispersion_min) %>% dplyr::mutate(compare_set_size=test_use_ids)
        return(df1)
      },    common_vars = common_vars,
            cor_mtx = cor_mtx,
            count_data = count_data,
            ref_mtx = ref_mtx,
            top_correlation_percentile_for_MLE = top_correlation_percentile_for_MLE,
            MLE_hit_probability_threshold = MLE_hit_probability_threshold,
            MLE_hit_ci_95_matches_threshold = MLE_hit_ci_95_matches_threshold,
            dispersion_min = dispersion_min,
            mc.cores = ncores 
        ) %>% bind_rows() 
  } else {
        results <- lapply(colnames(count_data), function(cell_id, common_vars, cor_mtx, count_data, ref_mtx, 
           top_correlation_percentile_for_MLE, MLE_hit_probability_threshold, 
           MLE_hit_ci_95_matches_threshold, dispersion_min) {
        cell_expression_counts_full <- count_data[,cell_id] %>% as.vector()
        names(cell_expression_counts_full) <- rownames(count_data)
        cor_vector <- cor_mtx[cell_id,]
        test_use_ids <- min(10,max(2,sum(cor_vector>quantile(cor_vector,top_correlation_percentile_for_MLE,na.rm = T))))
        top_correlated_ref_cell_names <- colnames(cor_mtx)[order(cor_mtx[cell_id,], decreasing = TRUE)][1:test_use_ids]
        # compute the MLE result
        df1 <- Get_Cell_Type_by_MLE(cell_id = cell_id,cell_expression_counts_full = cell_expression_counts_full,ref_expr = ref_mtx,count_obs = sum(cell_expression_counts_full),top_correlated_ref_cell_names = top_correlated_ref_cell_names,common_vars = common_vars,iter_num_max = 3,   fine_tune_marker_select_num = 10,top_cell_type_for_fine_tune = 5,probability_threshold = MLE_hit_probability_threshold,ci_95_matches_threshold = MLE_hit_ci_95_matches_threshold,dispersion_min = dispersion_min) %>% dplyr::mutate(compare_set_size=test_use_ids)
        return(df1)
      },    common_vars = common_vars,
            cor_mtx = cor_mtx,
            count_data = count_data,
            ref_mtx = ref_mtx,
            top_correlation_percentile_for_MLE = top_correlation_percentile_for_MLE,
            MLE_hit_probability_threshold = MLE_hit_probability_threshold,
            MLE_hit_ci_95_matches_threshold = MLE_hit_ci_95_matches_threshold,
            dispersion_min = dispersion_min
        ) %>% bind_rows() 
  }

  cat(as.vector(date()),'end \n')

  return(list(MLE_results=results,common_vars=common_vars,cor_mtx=cor_mtx,count_data=count_data,ref_mtx=ref_mtx))
}

## Tree-based distance measure functions 


# Compute_Transcriptional_Span_by_MLE_95CI
# Compute the transcriptional span of the reference cells
# ref_mtx: reference expression matrix. Note this should be the similar matrix as the one for estimating MLE result. 
# MLE_result_list: MLE result list. This should be the result from the same matrix as ref_mtx.
# clustering_method: the clustering method for the reference cells. Default is 'average'. Suggest not change.
# correlation_method_for_ref_mtx: the correlation method for the reference matrix. Default is 'spearman' Suggest not change.
# ncore: the number of cores
# return: a data frame with the cell id and the transcriptional span


Compute_Transcriptional_Span_by_MLE_95CI <- function(ref_mtx,MLE_result_list,clustering_method='average',correlation_method_for_ref_mtx='spearman',ncore=3){

  require(WGCNA)
  require(dplyr)
  require(ape)
  require(stringr)

  if(ncore > 1){
    require(pbmcapply)
  }
  
  # 1. compute the correlation matrix of the reference cells
  ref_cor <- WGCNA::cor(ref_mtx, method = correlation_method_for_ref_mtx)  # Spearman or pearson
  ref_dist <- as.dist(1 - ref_cor)

  # 2. hierarchical clustering
  ref_hc <- hclust(ref_dist, method = clustering_method)  # 你也可以用"ward.D2"等

  # 3. Compute the cophenetic distance matrix
  cophenetic_mat <- cophenetic(ref_hc)
  cophenetic_mat <- as.matrix(cophenetic_mat)

  # 4. calculate the distance of MLE 95%CI hits on tree
  if(ncore > 1){
    sum_dist_MLE_95CI <- pbmcapply::pbmclapply(1:nrow(MLE_result_list$MLE_results), function(x){
      candidate_mle_95ci_names_for_cell <- MLE_result_list$MLE_results$best_cell_types[x] %>% stringr::str_split(pattern=';',simplify = T) %>% unlist()
      mle_95ci_names_for_cell <- candidate_mle_95ci_names_for_cell[1:min(length(candidate_mle_95ci_names_for_cell),MLE_result_list$MLE_results$ci_95_matches[x]+1)]
      sum_dist_MLE_95CI <- Compute_Individual_Cell_Transcriptional_Span_on_Tree(ref_hc, mle_95ci_names_for_cell,cophenetic_mat )
      sum_dist_MLE_95CI_df <-data.frame(cell_id=MLE_result_list$MLE_results$cell_id[x],sum_dist_MLE_95CI=sum_dist_MLE_95CI) 
      return(sum_dist_MLE_95CI_df)
    },mc.cores = ncore) %>% dplyr::bind_rows() 
  }else{
    sum_dist_MLE_95CI <- lapply(1:nrow(MLE_result_list$MLE_results), function(x){
      mle_95ci_names_for_cell <- MLE_result_list$MLE_results$best_cell_types[x] %>% stringr::str_split(pattern=';',simplify = T) %>% unlist()
      sum_dist_MLE_95CI <- Compute_Individual_Cell_Transcriptional_Span_on_Tree(ref_hc, mle_95ci_names_for_cell,cophenetic_mat )
      sum_dist_MLE_95CI_df <-data.frame(cell_id=MLE_result_list$MLE_results$cell_id[x],sum_dist_MLE_95CI=sum_dist_MLE_95CI) 
      return(sum_dist_MLE_95CI_df) 
    }) %>% dplyr::bind_rows() 
  }
  return(sum_dist_MLE_95CI)
}


# Compute_Individual_Cell_Transcriptional_Span_on_Tree
# Compute the sum of the upper triangular part of the cophenetic distance matrix for a given set of nodes
# Example: sum_dist_topN <- get_tree_distance_sum(ref_hc, node_names=top_correlated_ref_cell_names,cophenetic_mat) 
# tree: a hclust object
# node_names: the names of the nodes to compute the distance
# cophenetic_mat: the cophenetic matrix built from the hclust object

Compute_Individual_Cell_Transcriptional_Span_on_Tree <- function(tree,node_names,cophenetic_mat) {
  # tree: a hclust object
  # node_names: the names of the nodes to compute the distance
  # cophenetic_mat: the cophenetic matrix built from the hclust object
  # 1. get the indices of these nodes in the tree
  idx <- match(node_names, tree$label)
  # 2. Extract the submatrix of the cophenetic matrix for these nodes
  sub_dist <- cophenetic_mat[idx, idx]
  # 3. Sum the upper triangular part of the submatrix
  result <- sum(sub_dist[upper.tri(sub_dist)],na.rm = TRUE)
  # 4. Return the result
  return(result)
}