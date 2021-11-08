#' Plot uncertainty of tissue labels from fit_banyan()
#'
#' This function allows you to visualize the uncertainty levels inferred cell sub-populations after running fit_banyan() and get_scores()
#' @param fit A list returned by fit_banyan(). Must have slots named "coords" and "U_scores".
#'
#' @keywords SBM MLSBM Gibbs Bayesian networks spatial gene expression
#' @import ggplot2
#' @export
#' @return A ggplot object
#' @examples
#' 
plot_uncertainty <- function(fit)
{
  coords = fit$coords
  Uncertainty = fit$U_scores
  coords$Uncertainty = Uncertainty
  g = ggplot(data = coords, aes(x = x, y = y, color = Uncertainty)) + 
    geom_point() + 
    theme_classic() + 
    xlab(NULL) + 
    ylab(NULL) + 
    scale_color_viridis_c(option = "A")
  return(g)
}