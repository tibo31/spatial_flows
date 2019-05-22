NE_impact <- function (AW, AW_Wd = NULL, AW_Wo = NULL, 
                           all_dest, all_origin) {
  
  # initialization 
  res_ne1 <- 0
  res_ne2 <- 0
  res_ne3 <- 0
  res_ne4 <- 0 
  
  unique_dest <- unique(all_dest)
  n_d <- length(unique_dest)
  unique_origin <- unique(all_origin)
  n_o <- length(unique_origin) 
  
  for (k in 1:n_d) {
    origin_k <- unique_origin[k]
    ind_row_A_ij_kl <- (all_origin != origin_k) & (all_dest != origin_k)
    # We first compute the 1st part of the equation
    ind_col_A_ij_kl <- all_origin == origin_k  
    res_ne1 <- res_ne1 + sum(AW[ind_row_A_ij_kl, ind_col_A_ij_kl])
    
    # Then, we compute the 2nd part of the equation
    ind_col_A_ij_lk <- all_dest == origin_k 
    res_ne2 <- res_ne2 + sum(AW[ind_row_A_ij_kl, ind_col_A_ij_lk])
    
    # Then, we compute the 3rd part of the equation
    if (!is.null(AW_Wo))
      res_ne3 <- res_ne3 + sum(AW_Wo[ind_row_A_ij_kl, ind_col_A_ij_kl]) 
    
    if (!is.null(AW_Wd))
      res_ne4 <- res_ne4 + sum(AW_Wd[ind_row_A_ij_kl, ind_col_A_ij_lk]) 
  }
  
  res <- list(res_ne1, res_ne2, res_ne3, res_ne4)
  return(res)
}