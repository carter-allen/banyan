#' Plot tissue labels from fit_banyan()
#'
#' This function allows you to visualize the inferred cell sub-populations after running fit_banyan()
#' @param fit A list returned by fit_banyan(). Must have slots named "coords" and "z".
#'
#' @keywords SBM MLSBM Gibbs Bayesian networks spatial gene expression
#' @import ggplot2
#' @importFrom rlang .data
#' @export
#' @return A ggplot object
#' 
plot_labels <- function(fit)
{
  coords = fit$coords
  z_map = fit$z
  coords$label = as.factor(z_map)
  g = ggplot(data = coords, aes(x = .data$x, y = .data$y, color = .data$label)) + 
    geom_point() + 
    theme_classic() + 
    xlab(NULL) + 
    ylab(NULL)
  return(g)
}