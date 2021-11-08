#' Calculate continuous uncertainty scores
#'
#' This function allows you to augment the discrete cell type assignments with continuous propensity and uncertainty scores
#' @param fit A list returned by fit_banyan()
#'
#' @keywords SBM MLSBM Gibbs Bayesian networks spatial gene expression
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
#' @return A list with populated entries C_scores (N x K matrix for cell type propensities) and U_scores (N x 1 vector of uncertainty scores)
#' 
get_scores <- function(fit)
{
  n = length(fit$z)
  K = fit$K
  AL = fit$A
  L = length(AL)
  z_map = fit$z
  pi_post = fit$pi
  P_post = fit$P
  
  C_scores <- matrix(0,nrow = n,ncol = K)
  U_scores <- rep(0,n)
  
  print("Calculating continuous phenotype and uncertainty scores")
  pb <- txtProgressBar(min = 0, max = n, style = 3)
  for(i in 1:n)
  {
    setTxtProgressBar(pb, i)
    pi_star <- rep(0,K)
    for(k in 1:K)
    {
      pi_star[k] = log(pi_post[k])
      for(j in 1:n)
      {
        if(i != j)
        {
          for(l in 1:L)
          {
            Al = AL[[l]]
            pi_star[k] = ((pi_star[k]) + 
                            log(P_post[k,z_map[j]]^Al[i,j]) + 
                            log((1 - P_post[k,z_map[j]])^(1-Al[i,j])))
          }
        }
      }
    }
    C_scores[i,] <- exp(pi_star)
    U_scores[i] <- 1 - exp(pi_star[z_map[i]])
  }
  close(pb)
  
  fit$C_scores <- C_scores
  fit$U_scores <- U_scores
  return(fit)
}