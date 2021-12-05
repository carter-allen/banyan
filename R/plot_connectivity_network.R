#' Plot community structure parameters as a K x K network
#'
#' This function allows you to visualize the inferred community structure as a community-community connectivity network
#' @param fit A list returned by fit_banyan()
#'
#' @keywords SBM MLSBM Gibbs Bayesian networks spatial gene expression
#' @import ggplot2 ggraph dplyr
#' @importFrom tidyr pivot_longer
#' @importFrom tidyselect everything
#' @importFrom rlang .data
#' @importFrom stats median
#' @importFrom igraph graph_from_data_frame
#' @export
#' @return A ggplot object
#' 
plot_connectivity_network <- function(fit)
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
  
  g_df = thetas_df_long %>%
    mutate(x_val = substr(.data$theta,7,7),
           y_val = substr(.data$theta,8,8),
           Type = ifelse(substr(.data$theta,7,7) == substr(.data$theta,8,8),
                         "Within Community",
                         "Between Community")) %>%
    group_by(.data$x_val,.data$y_val,.data$Type) %>%
    summarize(Connectivity = median(.data$value))
  
  g_df_within <- g_df %>%
    ungroup() %>%
    filter(.data$Type == "Within Community") %>%
    select(.data$x_val, .data$Connectivity) %>%
    rename(Community = .data$x_val,
           Within = .data$Connectivity) %>%
    as.data.frame()
  
  z_map = fit$z
  size = as.vector(table(z_map))
  g_df_size = data.frame(
    Community = as.character(1:K),
    Size = size
  )
  g_df_within = inner_join(g_df_within,g_df_size, by = "Community")
  
  
  g_df_between <- g_df %>%
    ungroup() %>%
    filter(.data$Type == "Between Community") %>%
    select(.data$x_val, .data$y_val, .data$Connectivity) %>%
    rename(from = .data$x_val, to = .data$y_val, Between = .data$Connectivity)
  
  G <- graph_from_data_frame(g_df_between,
                             directed = FALSE,
                             vertices = g_df_within)
  G_layout <- create_layout(G, layout = "igraph", algorithm = "nicely")
  
  g_ret <- ggraph(G_layout) +
    geom_edge_link(aes(width = .data$Between, alpha = .data$Between)) + 
    geom_node_point(aes(fill = .data$Within, size = .data$Size), shape = 21) +
    theme_graph() + 
    scale_size(range = c(12,22)) +
    geom_node_label(aes(label = .data$name), repel=FALSE) + 
    scale_fill_viridis_c(option = "A") 
  
  return(g_ret)
}