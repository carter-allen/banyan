#' Fit Banyan community detection model for spatially resolved single cell data
#'
#' This function allows you to fit the Bayesian multi-layer SBM for integration of spatial and gene expression data for identifying cell sub-populations
#' @param seurat_obj A Seurat object with PCA reduction and spatial coordinates. If provided, the exp and coords arguments are ignored
#' @param exp Gene expression dimension reduction feature matrix (e.g., PCs) with rows as cells and columns as features
#' @param coords A matrix or data frame with rows as cells and 2 columns for coordinates. Rows should be ordered the same as in exp.
#' @param K The number of sub-populations to infer
#' @param n_pcs The number of principal components to use from the Seurat object
#' @param R A length 2 vector of integers for the number of neighbors to use. R[1] corresponds to the number of neighbors in gene expression network and R[2] for spatial.
#' @param z_init Logical for whether or not to initialize the community allocation vector using the Louvain algorithm applied to the gene expression layer.
#' @param a0 Dirichlet prior parameter (shared across all communities)
#' @param b10 Beta prior number of connections
#' @param b20 Beta prior number on non-edges
#' @param n_iter The number of total MCMC iterations to run.
#' @param burn The number of MCMC iterations to discard as burn-in. A total of n_iter - burn iterations will be saved.
#' @param verbose Whether or not to print cluster allocations at each iteration
#' @param s Louvain resolution parameter
#'
#' @keywords SBM MLSBM Gibbs Bayesian networks spatial gene expression
#' @importFrom mlsbm fit_mlsbm
#' @importFrom scran buildKNNGraph
#' @export
#' @return A list of MCMC samples, including the MAP estimate of cluster indicators (z)
#' @examples
#' 
fit_banyan <- function(seurat_obj = NULL,
                       exp = NULL,
                       coords = NULL,
                       K,
                       n_pcs = 16,
                       R = NULL,
                       z_init = TRUE,
                       a0 = 2,
                       b10 = 1,
                       b20 = 1,
                       n_iter = 1000,
                       burn = 100, 
                       verbose = TRUE,
                       s = 1.2)
{
  # case 1: user provides a Seurat object
  if(!is.null(seurat_obj))
  {
    # check PCA reductions exist in the Seurat object
    if(!is.null(seurat_obj@reductions$pca))
    {
      exp <- seurat_obj@reductions$pca@cell.embeddings[,1:2]
    }
    else
    {
      return("Error: No PCA reductions found in the supplied Seurat object. Try passing features through using the exp argument, or add PCA reductions to seurat_obj$reductions$pca.")
    }
    # check coordinates exist in the Seurat object
    if(length(seurat_obj@images) != 1)
    {
      return("Error: Please provide Seurat object with exactly 1 image slot to access spatial coordinates from.")
    }
    else
    {
      coords_x <- seurat_obj@images[[1]]@coordinates$col
      coords_y <- seurat_obj@images[[1]]@coordinates$row
      coords <- data.frame(x = coords_x,
                           y = coords_y)
    }
  }
  
  # case 2: user provided exp and coords directly
  # check dimensions of everything
  if(nrow(exp) != nrow(coords))
  {
    return("Error: Number of rows (cells) should match between exp and coords.")
  }
  else
  {
    n <- nrow(exp)
  }
  
  if(nrow(exp) < ncol(exp))
  {
    return("Error: Number of features (columns) should be less than number of cells (rows) of exp.")
  }
  
  if(ncol(coords) != 2)
  {
    return("Error: Coordinates should be an N x 2 matrix")
  }
  else
  {
    colnames(coords) <- c("x","y")
  }
  
  # set number of neighbors
  if(is.null(R))
  {
    r1 <- r2 <- round(sqrt(n))
  }
  else
  {
    # user provided some invalid but non-null # of neighbors
    if(length(R) != 2)
    {
      return("Error: Please set R to either NULL or a length 2 vector of integer number of neighbors for each layer.")
    }
    else
    {
      # number of neighbors for gene expression
      r1 <- round(R[1])
      # number of neighbors for spatial
      r2 <- round(R[2])
    }
  }
  
  # build networks
  # gene expression
  exp <- as.matrix(exp)
  G1 <- scran::buildKNNGraph(exp,r1,transposed = TRUE)
  A1 <- igraph::as_adjacency_matrix(G1,sparse = FALSE)
  # spatial
  coords <- as.matrix(coords)
  G2 <- scran::buildKNNGraph(coords,r2,transposed = TRUE)
  A2 <- igraph::as_adjacency_matrix(G2,sparse = FALSE)
  
  AL <- list(A1,A2)
  
  fit <- mlsbm::fit_mlsbm(A = AL,
                          K = K,
                          z_init = z_init,
                          a0 = a0,
                          b10 = b10,
                          b20 = b20, 
                          n_iter = n_iter,
                          burn = burn,
                          verbose = verbose,
                          r = s)
  return(fit)
}