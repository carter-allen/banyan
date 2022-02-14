#' Plot community structure of cell sub-populations as posterior credible intervals
#'
#' This function allows you to visualize the community structure of cell sub-populations in posterior credible interval format via the connectivity parameters of the BANYAN model
#' @param fit A list returned by fit_banyan().
#'
#' @keywords SBM MLSBM Gibbs Bayesian networks spatial gene expression
#' @import ggplot2 dplyr
#' @import patchwork
#' @importFrom rlang .data
#' @importFrom tidyr pivot_longer separate
#' @importFrom tidyselect everything
#' @importFrom stats median quantile reorder
#' @export
#' @return A ggplot object
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
      t_names <- c(t_names,paste0("theta_",k1,"-",k2))
      k_count = k_count + 1
    }
  }
  colnames(thetas_df) <- t_names
  thetas_df <- as.data.frame(thetas_df)
  
  thetas_df_long <- thetas_df %>% 
    tidyr::pivot_longer(cols = tidyselect::everything(),
                        names_to = "theta",
                        values_to = "value")
  
  g_df = thetas_df_long %>%
    separate(col = .data$theta, 
             sep = "_",
             into = c("param","comb"), 
             remove = FALSE) %>%
    separate(col = .data$comb,
             sep = "-",
             into = c("x_val","y_val"),
             remove = FALSE) %>%
    mutate(x_val = as.numeric(.data$x_val),
           y_val = as.numeric(.data$y_val)) %>%
    mutate(pair = paste0(.data$x_val,.data$y_val),
           Type = ifelse(.data$x_val == .data$y_val,
                         "Within Community",
                         "Between Community")) %>%
    group_by(.data$pair,.data$Type) %>%
    summarize(Connectivity = median(.data$value),
              LB = quantile(.data$value,probs = 0.025),
              UB = quantile(.data$value,probs = 0.975))
  
  g_df_within = g_df %>%
    filter(.data$Type == "Within Community")
  
  g_df_between <- g_df %>%
    filter(.data$Type == "Between Community")
  
  g1 = ggplot(g_df_within,aes(x = stats::reorder(.data$pair,-.data$Connectivity),
                      y = .data$Connectivity)) + 
    geom_point() + 
    geom_errorbar(aes(ymin = .data$LB, ymax = .data$UB)) + 
    theme_classic() + 
    # scale_color_viridis_c(option = "A") + 
    theme(axis.text.x = element_text(family = "serif",size = 12),
          axis.text.y = element_text(family = "serif",size = 12),
          text = element_text(family = "serif"),
          axis.title.x = element_text(family = "serif", face = "bold"),
          axis.title.y = element_text(family = "serif",angle = 90, face = "bold")) + 
    xlab("Cell Sub-Population Pair") + 
    ylab("Connectivity") +
    ggtitle("Within Community Connectivity")
  
  g2 = ggplot(g_df_between,aes(x = stats::reorder(.data$pair,-.data$Connectivity),
                              y = .data$Connectivity)) + 
    geom_point() + 
    geom_errorbar(aes(ymin = .data$LB, ymax = .data$UB)) + 
    theme_classic() + 
    # scale_color_viridis_c(option = "A") + 
    theme(axis.text.x = element_text(family = "serif",size = 12),
          axis.text.y = element_text(family = "serif",size = 12),
          text = element_text(family = "serif"),
          axis.title.x = element_text(family = "serif", face = "bold"),
          axis.title.y = element_text(family = "serif",angle = 90, face = "bold")) + 
    xlab("Cell Sub-Population Pair") + 
    ylab("Connectivity") +
    ggtitle("Between Community Connectivity")
  
  return(g1 + g2)
}