OE_impact <- function (AW, AW_Wd = NULL, AW_Wo = NULL, 
                           all_dest, all_origin) {
  
  # initialization 
  res_o1 <- 0
  res_o2 <- 0
  res_o3 <- 0
  res_o4 <- 0
  
  unique_dest <- unique(all_dest)
  n_d <- length(unique_dest)
  unique_origin <- unique(all_origin)
  n_o <- length(unique_origin) 
  
  for (k in 1:n_o) {
    origin_k <- unique_origin[k]
    ind_k <- all_origin == origin_k
    ind_row_A_kj_kl <- ind_k & (all_dest != origin_k)
    # We first compute the 1st part of the equation
    ind_col_A_kj_kl <- ind_k 
    res_o1 <- res_o1 + sum(AW[ind_row_A_kj_kl, ind_col_A_kj_kl])
    
    # Then, we compute the 2nd part of the equation
    ind_col_A_kj_lk <- all_dest == origin_k 
    res_o2 <- res_o2 + sum(AW[ind_row_A_kj_kl, ind_col_A_kj_lk])
    
    # Then, we compute the 3rd part of the equation
    if (!is.null(AW_Wo))
      res_o3 <- res_o3 + sum(AW_Wo[ind_row_A_kj_kl, ind_col_A_kj_kl]) 
    
    # Then, we compute the 4th part of the equation
    if (!is.null(AW_Wd))
      res_o4 <- res_o4 + sum(AW_Wd[ind_row_A_kj_kl, ind_col_A_kj_lk]) 
  }
  
  res <- list(res_o1, res_o2, res_o3, res_o4)
  return(res)
}