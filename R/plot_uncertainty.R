#' Plot uncertainty of tissue labels from fit_banyan()
#'
#' This function allows you to visualize the uncertainty levels inferred cell sub-populations after running fit_banyan() and get_scores()
#' @param fit A list returned by fit_banyan(). Must have slots named "coords" and "U_scores".
#'
#' @keywords SBM MLSBM Gibbs Bayesian networks spatial gene expression
#' @import ggplot2
#' @importFrom rlang .data
#' @export
#' @return A ggplot object
#' 
plot_uncertainty <- function(fit)
{
  coords = fit$coords
  coords = as.data.frame(coords)
  colnames(coords) = c("x","y")
  Uncertainty = fit$U_scores
  coords$Uncertainty = Uncertainty
  g = ggplot(data = coords, aes(x = .data$x, y = .data$y, color = .data$Uncertainty)) + 
    geom_point() + 
    theme_classic() + 
    xlab(NULL) + 
    ylab(NULL) + 
    scale_color_viridis_c(option = "A")
  return(g)
}