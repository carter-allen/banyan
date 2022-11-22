#' Get ACC parameters
#'
#' Conduct ACC based on user-defined tissue architecture labels
#'
#' @param z User-defined cell spot labels 
#' @param AL List of adjacency matrices
#' @param b10 beta prior parameter 1
#' @param b20 beta prior parameter 2
#'
#' @return an adjacency matrix
#' @export
#' @importFrom utils setTxtProgressBar txtProgressBar
get_ACC <- function(z,AL,b10,b20)
{
  K = length(unique(z))
  L = length(AL)
  n = length(z)
  
  beta_mat <- matrix(0, nrow = 2, ncol = choose(K,2) + K)
  
  pb <- txtProgressBar(min = 0, max = choose(K,2) + K, style = 3)
  k_count = 1
  beta_mat_names = NULL
  for(a in 1:K)
  {
    for(b in a:K)
    {
      sum_A_ij = 0
      sum_one_minus_A_ij = 0
      
      for(i_sample in 1:(n-1))
      {
        for(j_sample in (i_sample+1):n)
        {
          if(z[i_sample] == a && z[j_sample] == b)
          {
            for(l in 1:L)
            {
              Al = AL[[l]]
              sum_A_ij = sum_A_ij + Al[i_sample,j_sample]
              sum_one_minus_A_ij = sum_one_minus_A_ij + (1-Al[i_sample,j_sample])
            }
          }
        }
      }
      
      beta1_star = b10 + sum_A_ij 
      beta2_star = b20 + sum_one_minus_A_ij
      
      beta_mat[,k_count] = c(beta1_star,beta2_star)
      beta_mat_names = c(beta_mat_names,paste0(a,"-",b))
      
      setTxtProgressBar(pb, k_count)
      k_count = k_count + 1
    }
  }
  
  colnames(beta_mat) <- beta_mat_names
  return(beta_mat)
}