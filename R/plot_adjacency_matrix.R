#' Plot cell-cell adjacency matrix sorted by labels from fit_banyan()
#'
#' This function allows you to visualize the inferred adjacency matrix after running fit_banyan()
#' @param fit A list returned by fit_banyan(). Must have slots named "A" and "z".
#' @param layer Either "expression" (default) or "spatial" for which adjacency network to use
#'
#' @keywords SBM MLSBM Gibbs Bayesian networks spatial gene expression
#' @import ggplot2
#' @importFrom ggplotify as.ggplot
#' @importFrom pheatmap pheatmap
#' @export
#' @return A ggplot object
#' 
plot_adjacency_matrix <- function(fit,layer = "expression")
{
  if(layer == "expression") A_mat <- fit$A[[1]]
  else if(layer == "spatial") A_mat <- fit$A[[2]]
  else
  {
    message("Error: Invalid layer argument. Should be either expression or spatial. Using expression.")
    A_mat <- fit$A[[1]]
  }
  z_map <- fit$z
  g_ret <- ggplotify::as.ggplot(pheatmap::pheatmap(A_mat[order(z_map),order(z_map)],
                                                   cluster_rows = F,
                                                   cluster_cols = F,
                                                   legend = F,
                                                   color = c("white","black"))) 
  return(g_ret)
}