intra_impact <- function (AW, AW_Wd = NULL, AW_Wo = NULL,
                           all_dest, all_origin) {
  
  # initialization 
  res_i1 <- 0
  res_i2 <- 0
  res_i3 <- 0
  res_i4 <- 0
  
  unique_dest <- unique(all_dest)
  n_d <- length(unique_dest)
  unique_origin <- unique(all_origin)
  n_o <- length(unique_origin) 
  
  unique_obs <- intersect(unique_dest, unique_origin)
  n <- length(unique_obs)
    
  for (k in 1:n) {
    origin_k <- unique_origin[k]
    ind_row_A_kk_kl <- (all_origin == origin_k) & (all_dest == origin_k)
    # We first compute the 1st part of the equation
    ind_col_A_kk_kl <- all_origin == origin_k  
    res_i1 <- res_i1 + sum(AW[ind_row_A_kk_kl, ind_col_A_kk_kl])
    
    # Then, we compute the 2nd part of the equation
    ind_col_A_kk_lk <- all_dest == origin_k 
    res_i2 <- res_i2 + sum(AW[ind_row_A_kk_kl, ind_col_A_kk_lk])
    
    # Then, we compute the 3rd part of the equation
    if (!is.null(AW_Wo))
      res_i3 <- res_i3 + sum(AW_Wo[ind_row_A_kk_kl, ind_col_A_kk_kl]) 
 
    if (!is.null(AW_Wd))
      res_i4 <- res_i4 + sum(AW_Wd[ind_row_A_kk_kl, ind_col_A_kk_lk])    
  }
  
  res <- list(res_i1, res_i2, res_i3, res_i3)
  return(res)
}