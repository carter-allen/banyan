#' Plot continuous phenotype scores from fit_banyan()
#'
#' This function allows you to visualize continuous propensity of each cell spot towards each sub-population (e.g., continuous phenoptypes) after running fit_banyan() and get_scores()
#' @param fit A list returned by fit_banyan(). Must have slots named "coords" and "C_scores".
#' @param k Which cell sub-population to compute propensity towards
#'
#' @keywords SBM MLSBM Gibbs Bayesian networks spatial gene expression
#' @import ggplot2
#' @export
#' @return A ggplot object
#' @examples
#' 
plot_propensity <- function(fit, k = 1)
{
  coords = fit$coords
  Propensity = fit$C_scores[,k]
  coords$Propensity = Propensity
  g = ggplot(data = coords, aes(x = x, y = y, color = Propensity)) + 
    geom_point() + 
    theme_classic() + 
    xlab(NULL) + 
    ylab(NULL) + 
    scale_color_viridis_c(option = "A")
  return(g)
}