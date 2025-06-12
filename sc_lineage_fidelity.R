
sc_lineage_fidelity <- function (scdata, ref_data=NULL,numbers_plot = 3,ncores=10,top_n_candidate_for_fidelity=10,species='human') 
{
  require(WGCNA)
  require(dplyr)
  require(reshape2)
  require(parallel)
  gettissue <- function (x, Num = 3) 
  {
    top_markers <- order(x, decreasing = T)[1:Num]
    return(top_markers)
  }
  
  if(is.null(ref_data) & species == 'human'){
    ref_data = readRDS('/gpfs/output/ECS_Research/data/scRNA/HPA_scRNA_ref/ref_expr_human.HPA.rds')
  }
  if(is.null(ref_data) & species == 'mouse'){
    ref_data = readRDS('/gpfs/output/ECS_Research/data/scRNA/HPA_scRNA_ref/scMCA_mouse_ref_expr.rds')
  }
  ref.expr = ref_data 
  common_genes <- rownames(ref.expr)[rownames(ref.expr) %in% rownames(scdata)]
  ref.expr <- ref.expr[common_genes, ]
  tst.expr <- scdata[common_genes,]  
  
  tst.expr[is.na(tst.expr)] <- 0
  tst.expr <- as.matrix(t(t(tst.expr)/colSums(tst.expr)) * 
                          1e+05)
  tst.expr <- log(tst.expr + 1)
  #parallel code 
    ref.mat <- log(ref.expr + 1)
    # Get the number of columns in tst.expr
    n_tst <- ncol(tst.expr)
    # Split the column indices of tst.expr into ncores groups
    split_indices <- split(1:n_tst, cut(1:n_tst, breaks = ncores, labels = FALSE))
    # Perform parallel correlation computations
    cors_list <- parallel::mclapply(split_indices, function(idx) {
      # Subset tst.expr for the current chunk
      sub_tst <- tst.expr[, idx, drop = FALSE]
      # Compute the correlation matrix between all columns of ref.mat
      # and the subset of tst.expr columns; add use="pairwise.complete.obs"
      # if needed to handle any missing computations.
      WGCNA::cor(ref.mat, sub_tst, use = "pairwise.complete.obs")
    }, mc.cores = ncores)
  # Combine the correlation matrices (column-wise) to form the full result
    cors <- do.call(cbind, cors_list)
  #
  cors_index <- apply(cors, 2, gettissue, numbers_plot)
  cors_index <- sort(unique(as.integer(cors_index)))
  scblast.result <- apply(cors, 2, function(x) rownames(cors)[which.max(x)])
  cors_in = cors[cors_index, ]
  colnames(cors_in) <- colnames(scdata)
  cors_out = reshape2::melt(cors_in)[c(2, 1, 3)]
  colnames(cors_out) <- c("Cell", "Cell type", "Score")
  cors_out <- as.data.frame(cors_out %>% group_by(Cell) %>% 
                              top_n(n = numbers_plot, wt = Score))
  lineage_fidelity <- cors_out %>% dplyr::group_by(Cell) %>% dplyr::top_n(n=top_n_candidate_for_fidelity,wt=Score) %>% dplyr::summarise(lineage_fidelity=sd(Score)/mean(Score),lineage_fidelity_var=var(Score))
  result <- list()
  result[["cors_matrix"]] <- cors
  result[["top_cors"]] <- numbers_plot
  result[["scMCA"]] <- scblast.result
  result[["scMCA_probility"]] <- cors_out
  result[["lineage_fidelity"]] <- lineage_fidelity
  return(result)
}