DE_impact <- function (AW, AW_Wd = NULL, AW_Wo = NULL,
                           all_dest, all_origin) {
  
  # initialization 
  res_d1 <- 0
  res_d2 <- 0
  res_d3 <- 0
  res_d4 <- 0
  
  unique_dest <- unique(all_dest)
  n_d <- length(unique_dest)
  unique_origin <- unique(all_origin)
  n_o <- length(unique_origin) 

  for (k in 1:n_d) {
    dest_k <- unique_dest[k]
    ind_k <- all_dest == dest_k
    ind_row_A_ik_kl <- (all_origin != dest_k) & ind_k
    # We first compute the 1st part of the equation
    ind_col_A_ik_kl <- all_origin == dest_k  
    res_d1 <- res_d1 + sum(AW[ind_row_A_ik_kl, ind_col_A_ik_kl])
    
    # Then, we compute the 2nd part of the equation
    ind_col_A_ik_lk <- ind_k 
    res_d2 <- res_d2 + sum(AW[ind_row_A_ik_kl, ind_col_A_ik_lk])
    
    # Then, we compute the 3rd part of the equation
    if (!is.null(AW_Wo))
      res_d3 <- res_d3 + sum(AW_Wo[ind_row_A_ik_kl, ind_col_A_ik_kl]) 
  
    if (!is.null(AW_Wd))
      res_d4 <- res_d4 + sum(AW_Wd[ind_row_A_ik_kl, ind_col_A_ik_lk]) 
    
  }
  
  res <- list(res_d1, res_d2, res_d3, res_d4)
  return(res)
}