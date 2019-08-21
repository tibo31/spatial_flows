gravity_model_geo <- function(Y, G, x_d, x_o, 
                              ind_d = NULL, ind_o = NULL,
                              lagged = F, centered = F, 
                              DW = NULL, OW = NULL) {
  
  # verification  
  if (lagged & (is.null(DW) | is.null(OW))){
    stop("DW and OW must be given for lagging x variable")
  }

  # initialization
  x_matrix_d <- as(x_d, "matrix")
  x_matrix_o <- as(x_o, "matrix")
  n_d <- nrow(x_matrix_d)
  n_o <- nrow(x_matrix_o)
  # another verification 
  stopifnot(n_d == nrow(Y), n_o ==ncol(Y))
  stopifnot(n_d == nrow(G), n_o ==ncol(G))  
  if (lagged) {
    stopifnot(n_d == nrow(DW), n_o == ncol(OW)) 
  }
  N <- n_d * n_o
  p_d <- ncol(x_matrix_d)
  p_o <- ncol(x_matrix_o)
  iota_n_d <- rep(1, n_d)
  iota_n_o <- rep(1, n_o)
  
  # we compute constants for the computation of the effects  
  sum_G <- sum(G)
  sum_Y <- sum(Y)
  # check concerning the coices of variables
  all_variables <- ifelse(is.null(ind_d) & is.null(ind_o), T, F) 

  # construct the lagged variables 
  if (lagged) {
    x_matrix_d_lagged <- DW %*% x_matrix_d  
    x_matrix_o_lagged <- OW %*% x_matrix_o
  } 
  
  # names of the variable at the origin
  names_x_o <- colnames(x_matrix_o)
  names_x_d <- colnames(x_matrix_d)
  if (is.null(names_x_o)) {
    names_x_o <- paste0("x_o_", seq_len(ncol(x_o)))
  }
  if (is.null(names_x_d)) {
    names_x_d <- paste0("x_d_", seq_len(ncol(x_d)))
  }
  
  # Select variables
  if (is.null(ind_d)) {
    ind_d <- 1:p_d
  }
  
  if (is.null(ind_o)) {
    ind_o <- 1:p_o
  }
  
  if (lagged) {
    # create lagged variables
    x_matrix_d <- cbind(x_matrix_d[, ind_d],
                        x_matrix_d_lagged[, ind_d])
    x_matrix_o <- cbind(x_matrix_o[, ind_o],
                        x_matrix_o_lagged[, ind_o])
    names_x_d <- c(names_x_d[ind_d], 
                        paste0("lagged_", names_x_d[ind_d]))
    names_x_o <- c(names_x_o[ind_o], 
                        paste0("lagged_", names_x_o[ind_o]))    
    } else {
      x_matrix_d <- as.matrix(x_matrix_d[, ind_d]) # as.matrix if there is only 1 variable
      x_matrix_o <- as.matrix(x_matrix_o[, ind_o])
      names_x_d <- names_x_d[ind_d]
      names_x_o <- names_x_o[ind_o]  
    }
    
    # dimension 
    R_d <- length(names_x_d)
    R_o <- length(names_x_o)
    
    # names of the result variable
    names_x_od <- c(names_x_d, names_x_o) 
    
    if (centered) {
      mean_x_matrix_d <- apply(x_matrix_d, 2, mean)
      mean_x_matrix_o <- apply(x_matrix_o, 2, mean)
      x_matrix_d <- x_matrix_d - matrix(rep(mean_x_matrix_d, 
                                            each = n_d), n_d, ncol(x_matrix_d))
      x_matrix_o <- x_matrix_o - matrix(rep(mean_x_matrix_o, 
                                            each = n_o), n_o, ncol(x_matrix_o))
      # compute some blocks for the estimates
      mean_d <- rep(0, R_d)
      mean_o <- rep(0, R_o)
      outer_x_d_x_o <- matrix(0, R_d, R_o) 
    } else {
      mean_d <- n_o * colSums(x_matrix_d)
      mean_o <- n_d * colSums(x_matrix_o)
      outer_x_d_x_o <- matrix(colSums(kronecker(x_matrix_d, x_matrix_o)), 
                              R_d, R_o, byrow = T)
    }
    
    # we compute block by block 
    crossprod_x_d <- crossprod(x_matrix_d)
    crossprod_x_o <- crossprod(x_matrix_o)
    G_x_d <- crossprod(x_matrix_d, G)  %*% iota_n_o
    G_x_o <- crossprod(x_matrix_o, crossprod(G, iota_n_d))
    
    # We code the formula 83.3 in Thomas and LeSage (2014)
    zprime_z <- rbind(
      c(N, mean_d, mean_o, sum_G),
      cbind(mean_d, n_o * crossprod_x_d, outer_x_d_x_o, G_x_d),
      cbind(mean_o, t(outer_x_d_x_o), n_d * crossprod_x_o, G_x_o),
      cbind(sum_G, t(G_x_d), t(G_x_o), sum(diag(crossprod(G))))
    )
    
    # We code the formula 83.4 in Thomas and LeSage (2014) :
    z_prime_y <- rbind(sum_Y,
                       crossprod(x_matrix_d, Y %*% iota_n_o),
                       crossprod(x_matrix_o, t(Y) %*% iota_n_d),
                       sum(diag(t(G) %*% Y)))
    
    
  
  # The estimates of the coefficients is given by formula 83.5 in Thomas and LeSage (2014) :
  hat_matrix <- qr.solve(zprime_z)
  hat_delta <- hat_matrix %*% z_prime_y
  n_K <- length(hat_delta)
  
  # computation of the significance
  
  x_matrix_d <- x_matrix_d * matrix(hat_delta[2:(R_d + 1)], 
                                    n_d, R_d, byrow = T)
  x_matrix_o <- x_matrix_o * matrix(hat_delta[(R_d + 2):(n_K - 1)], 
                                n_o, R_o, byrow = T) 
  # prediction
  residus <- Y - (hat_delta[n_K] * G + hat_delta[1])
  rowSums_o <- rowSums(x_matrix_o)
  for (k in 1:n_d) {
    residus[k, ] <- residus[k, ] - (sum(x_matrix_d[k, ]) + rowSums_o)
  }
  s2 <- sum(residus ^ 2) / (N - n_K)
  hat_matrix <- s2 * hat_matrix
  t_statistic <- hat_delta / sqrt(diag(hat_matrix)) 
  
  # correction of the intercept in case of centering
  if (centered) {
    hat_delta[1] <- hat_delta[1] - sum(
      hat_delta[-c(1, length(hat_delta))] *
        c(mean_x_matrix_d, mean_x_matrix_o)) 
    t_statistic[1] <- NA
  }
  
  rownames(hat_delta) <- c("(Intercept)", names_x_od, "Distance")
  res <- cbind(hat_delta, t_statistic)
  dimnames(res)[[2]] <- c("Estimate", "t value")
  return(res)
}
