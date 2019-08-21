gravity_model <- function(x, Y, G, ind_d = NULL, ind_o = NULL,
                          lagged = F, centered = F, w = NULL, 
                          intra_x = F) {
  
  # verification
  if ((nrow(Y) != ncol(Y)) | (nrow(G) != ncol(G)) | 
      (nrow(Y) != ncol(G))) {
    stop("Y and G should be square matrices with same dimension")
  }

  if (lagged) {
    if (is.null(w))
      stop("w must be given for lagging x variable")
    if ((ncol(w) != nrow(w)) | (ncol(Y) != nrow(w)))
      stop("the dimension of w is incorrect")
  }

  # initialization
  x_matrix <- as(x, "matrix")
  n <- nrow(x)
  # another verification 
  stopifnot(n == nrow(Y))
  N <- n * n
  p <- ncol(x)
  iota_n <- rep(1, n)

  # we compute constants for the computation of the effects  
  sum_G <- sum(G)
  sum_Y <- sum(Y)
  # check concerning the coices of variables
  all_variables <- ifelse(is.null(ind_d) & is.null(ind_o), T, F) 

  # names of the variable 
  names_x <- colnames(x)
  if (is.null(names_x)) {
    names_x <- paste0("x_", seq_len(ncol(x)))
  }
  
  # separate intraregional case and normal case 
  if (!intra_x) {  
    ###############################################
    # normal case
    if (is.null(ind_d)) {
      ind_d <- 1:p
    }

    if (is.null(ind_o)) {
      ind_o <- 1:p
    }
  
    if (lagged) {
      # construct the lagged variables 
      x_matrix_lagged <- w %*% x_matrix

      # create lagged variables
      x_matrix_d <- cbind(x_matrix[, ind_d],
                          x_matrix_lagged[, ind_d])
      x_matrix_o <- cbind(x_matrix[, ind_o],
                          x_matrix_lagged[, ind_o])
      names_x_d <- paste0("Dest_", c(names_x[ind_d], 
                          paste0("lagged_", names_x[ind_d])))
      names_x_o <- paste0("Origin_", c(names_x[ind_o], 
                          paste0("lagged_", names_x[ind_o])))    
    } else {
      x_matrix_d <- as.matrix(x_matrix[, ind_d]) # as.matrix if there is only 1 variable
      x_matrix_o <- as.matrix(x_matrix[, ind_o])
      names_x_d <- paste0("Dest_", names_x[ind_d])
      names_x_o <- paste0("Origin_", names_x[ind_o])   
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
                                            each = n), n, ncol(x_matrix_d))
      x_matrix_o <- x_matrix_o - matrix(rep(mean_x_matrix_o, 
                                            each = n), n, ncol(x_matrix_o))
      # compute some blocks for the estimates
      mean_d <- rep(0, R_d)
      mean_o <- rep(0, R_o)
      outer_x_d_x_o <- matrix(0, R_d, R_o) 
    } else {
      mean_d <- n * colSums(x_matrix_d)
      mean_o <- n * colSums(x_matrix_o)
      outer_x_d_x_o <- matrix(colSums(kronecker(x_matrix_d, x_matrix_o)), 
                              R_d, R_o, byrow = T)
    }
    
    # we compute block by block 
    crossprod_x_d <- crossprod(x_matrix_d)
    crossprod_x_o <- crossprod(x_matrix_o)
    G_x_d <- crossprod(x_matrix_d, G %*% iota_n)
    G_x_o <- crossprod(x_matrix_o, G %*% iota_n)
    
    # We code the formula 83.3 in Thomas and LeSage (2014)
    zprime_z <- rbind(
      c(N, mean_d, mean_o, sum_G),
      cbind(mean_d, n * crossprod_x_d, outer_x_d_x_o, G_x_d),
      cbind(mean_o, t(outer_x_d_x_o), n * crossprod_x_o, G_x_o),
      cbind(sum_G, t(G_x_d), t(G_x_o), sum(diag(crossprod(G))))
    )
    
    # We code the formula 83.4 in Thomas and LeSage (2014) :
    z_prime_y <- rbind(sum_Y,
                       crossprod(x_matrix_d, Y %*% iota_n),
                       crossprod(x_matrix_o, t(Y) %*% iota_n),
                       sum(diag(G %*% Y)))
    
    
  } else {
    ###############################################
    # intraregional case
    if (!all.equal(ind_d, ind_o)) {
      warning("In the case of an intraregional model, all variables must be both at origin and destination")
    }
    ind_od <- intersect(ind_d, ind_o)
    if (is.null(ind_od)) {
      ind_od <- 1:p
    }
    names_x <- names_x[ind_od]
    if (lagged) {
      # construct the lagged variables 
      x_matrix_lagged <- w %*% x_matrix
      x_matrix <- cbind(x_matrix[, ind_od],
                          x_matrix_lagged[, ind_od])
      names_x_od <- c(names_x, paste0("lagged_", names_x))
      names_x_d <- paste0("Dest_", names_x_od)
      names_x_o <- paste0("Origin_", names_x_od) 
      names_x_i <- paste0("Intra_", names_x_od) 
    } else {
      x_matrix <- as.matrix(x_matrix[, ind_od])
      names_x_d <- paste0("Dest_", names_x) 
      names_x_o <- paste0("Origin_", names_x) 
      names_x_i <- paste0("Intra_", names_x) 
    }
    
    # dimension 
    R <- length(names_x_d)
    
    # names of the result variable
    names_x_od <- c(names_x_d, names_x_o, names_x_i) 

    # we fix the block :
    sums_x <- colSums(x_matrix)
    mean_d <- mean_o <- (n - 1) * sums_x
    zeros_d <- zeros_o <- numeric(R)
    mean_i <- sums_x
    zeros_i_i <- matrix(0, R, R)
    zeros_i <- numeric(R)
    mean_o_d <- matrix(0, R, R)
    crossprod_x <- crossprod(x_matrix)
    outer_x <- matrix(colSums(kronecker(x_matrix, x_matrix)), 
                            R, R, byrow = T)
    G_x <- crossprod(x_matrix, G %*% iota_n)
    sum_G_x <- rep(sum(G_x), R)
    
    zprime_z <- rbind(
      c(N, n, mean_d, mean_o, mean_i, sum_G),
      c(n, n, zeros_d, zeros_o, mean_i, 0),
      cbind(mean_d, zeros_d, (n - 1) * crossprod_x, 
            outer_x - crossprod_x, zeros_i_i, G_x),
      cbind(mean_o, zeros_o, outer_x - crossprod_x, (n - 1) * crossprod_x, zeros_i_i, G_x), 
      cbind(mean_i, mean_i, zeros_i_i, zeros_i_i, crossprod_x, zeros_i),
      c(sum_G, 0, t(G_x), t(G_x), zeros_i, sum(diag(crossprod(G))))
    )
    
    # We code the formula 83.4 in Thomas and LeSage (2014) :
    z_prime_y <- rbind(sum_Y,
                       sum(diag(Y)),
                       crossprod(x_matrix, Y %*% iota_n) - crossprod(x_matrix, diag(Y)),
                       crossprod(x_matrix, t(Y) %*% iota_n) - crossprod(x_matrix, diag(Y)),
                       crossprod(x_matrix, diag(Y)),
                       sum(diag(G %*% Y)))
  }
  
  # The estimates of the coefficients is given by formula 83.5 in Thomas and LeSage (2014) :
  hat_matrix <- qr.solve(zprime_z)
  hat_delta <- hat_matrix %*% z_prime_y
  n_K <- length(hat_delta)
  
  # computation of the significance
  if (!intra_x) {
    x_matrix_d <- x_matrix_d * matrix(hat_delta[2:(R_d + 1)], 
                                n, R_d, byrow = T)
    x_matrix_o <- x_matrix_o * matrix(hat_delta[(R_d + 2):(n_K - 1)], 
                                n, R_o, byrow = T) 
    # prediction
    residus <- Y - (hat_delta[n_K] * G + hat_delta[1])
    rowSums_o <- rowSums(x_matrix_o)
    for (k in 1:n) {
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
  } else {
    x_matrix_d <- x_matrix * matrix(hat_delta[3:(2 + R)], 
                                    n, R, byrow = T)
    rowSums_o <- x_matrix %*% hat_delta[(R + 3):(2 * R + 2)] 
 
    # prediction
    residus <- Y - (hat_delta[n_K] * G + hat_delta[1] + diag(n) * hat_delta[2])
    for (k in 1:n) {
      residus[k, -k] <- (residus[k, ] - (sum(x_matrix_d[k, ]) + rowSums_o))[-k]
    }
    diag(residus) <- diag(residus) -  x_matrix %*% hat_delta[(2 * R + 3):(n_K - 1)]
    s2 <- sum(residus ^ 2) / (N - n_K)
    hat_matrix <- s2 * hat_matrix
    t_statistic <- hat_delta / sqrt(diag(hat_matrix)) 
    rownames(hat_delta) <- c("(Intercept)", "c_i", names_x_od, "Distance")
    res <- cbind(hat_delta, t_statistic)
    dimnames(res)[[2]] <- c("Estimate", "t value")
  }
  
  return(res)
}