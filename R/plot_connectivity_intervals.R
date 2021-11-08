#' Plot community structure of cell sub-populations as posterior credible intervals
#'
#' This function allows you to visualize the community structure of cell sub-populations in posterior credible interval format via the connectivity parameters of the BANYAN model
#' @param fit A list returned by fit_banyan().
#'
#' @keywords SBM MLSBM Gibbs Bayesian networks spatial gene expression
#' @import ggplot2 dplyr
#' @importFrom tidyr pivot_longer
#' @importFrom tidyselect everything
#' @export
#' @return A ggplot object
#' @examples
#' 
plot_connectivity_intervals <- function(fit)
{
  K = fit$K
  Kc = choose(K,2) + K
  THETA = fit$PM
  n_sim = dim(THETA)[1]
  
  thetas_df <- matrix(0,nrow = n_sim,ncol = Kc)
  k_count = 1
  t_names <- NULL
  for(k1 in 1:K)
  {
    for(k2 in k1:K)
    {
      thetas_df[,k_count] = THETA[,k1,k2]
      t_names <- c(t_names,paste0("theta_",k1,k2))
      k_count = k_count + 1
    }
  }
  colnames(thetas_df) <- t_names
  thetas_df <- as.data.frame(thetas_df)
  
  thetas_df_long <- thetas_df %>% 
    tidyr::pivot_longer(cols = tidyselect::everything(),
                        names_to = "theta",
                        values_to = "value")
  
  g <- thetas_df_long %>%
    mutate(pair = substr(theta,7,8),
           Type = ifelse(substr(theta,7,7) == substr(theta,8,8),
                         "Within Community",
                         "Between Community")) %>%
    group_by(pair,Type) %>%
    summarize(Connectivity = median(value),
              LB = quantile(value,probs = 0.025),
              UB = quantile(value,probs = 0.975)) %>%
    ggplot(.,aes(x = stats::reorder(pair,-Connectivity),y = Connectivity,color = Type)) + 
    geom_point() + 
    geom_errorbar(aes(ymin = LB, ymax = UB)) + 
    theme_classic() + 
    theme(axis.text.x = element_text(family = "serif",size = 12),
          axis.text.y = element_text(family = "serif",size = 12),
          text = element_text(family = "serif"),
          axis.title.x = element_text(family = "serif", face = "bold"),
          axis.title.y = element_text(family = "serif",angle = 90, face = "bold")) + 
    xlab("Cell Sub-Population Pair") + 
    ylab("Connectivity") 
  return(g)
}