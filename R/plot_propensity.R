#' Plot continuous phenotype scores from fit_banyan()
#'
#' This function allows you to visualize continuous propensity of each cell spot towards each sub-population (e.g., continuous phenoptypes) after running fit_banyan() and get_scores()
#' @param fit A list returned by fit_banyan(). Must have slots named "coords" and "C_scores".
#' @param k Which cell sub-population to compute propensity towards
#'
#' @keywords SBM MLSBM Gibbs Bayesian networks spatial gene expression
#' @import ggplot2
#' @importFrom rlang .data
#' @export
#' @return A ggplot object
#' 
plot_propensity <- function(fit, k = 1)
{
  coords = fit$coords
  coords = as.data.frame(coords)
  colnames(coords) = c("x","y")
  Propensity = fit$C_scores[,k]
  coords$Propensity = Propensity
  g = ggplot(data = coords, aes(x = .data$x, y = .data$y, color = .data$Propensity)) + 
    geom_point() + 
    theme_classic() + 
    xlab(NULL) + 
    ylab(NULL) + 
    scale_color_viridis_c(option = "A")
  return(g)
}