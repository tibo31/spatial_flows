gravity_model <- function(x, Y, G, ind_d = NULL, ind_o = NULL,
                          lagged = F, w = NULL) {
  
  # verification
  if (lagged & is.null(w))
    stop("w must be gived for lagging x variable")
  
  # initialization
  n <- nrow(x)
  N <- n * n
  p <- ncol(x)
  iota_n <- rep(1, n)
  # we center the data: 
  x_matrix <- as(x, "matrix")
  x_matrix_centered <- x_matrix - matrix(rep(apply(x_matrix, 2, mean), 
                                             each = n), n, p)
  if (lagged) {
    x_matrix_lagged <- w %*% x_matrix
    x_matrix_centered_lagged <- x_matrix_lagged - matrix(rep(apply(x_matrix_lagged, 2, mean), 
                                                      each = n), n, p)
  } 
  
  # names of the variable 
  names_x <- colnames(x)
  
  
  if (is.null(ind_d)) {
    ind_d <- 1:p
    names_x_d <- paste0(names_x, "_d") 
  } else {
    names_x_d <- paste0(names_x[ind_d], "_d") 
  }
  
  if (is.null(ind_o)) {
    ind_o <- 1:p
    names_x_o <- paste0(names_x, "_o") 
  } else {
    names_x_o <- paste0(names_x[ind_o], "_o") 
  }
  
  if (lagged) {
    names_x_d <- c(names_x_d, paste0("lagged_", names_x_d))
    names_x_o <- c(names_x_o, paste0("lagged_", names_x_o))    
  } 
  
  # names of the result variable
  names_x_od <- c(names_x_d, names_x_o) 
  
  R_d <- length(names_x_d)
  R_o <- length(names_x_o)
  
  # R <- R_d + R_o
  zero_k_d <- rep(0, R_d)
  zero_k_d_o <- matrix(0, R_d, R_o)
  zero_k_o <- rep(0, R_o)
  zero_k_o_d <- matrix(0, R_o, R_d)
  # zero_k_prime_zero_k <- crossprod(t(zero_k))
  
  if (lagged) {
    x_matrix_centered_d <- cbind(x_matrix_centered[, ind_d],
                               x_matrix_centered_lagged[, ind_d])
    x_matrix_centered_o <- cbind(x_matrix_centered[, ind_o],
                               x_matrix_centered_lagged[, ind_o])
  } else {
    x_matrix_centered_d <- x_matrix_centered[, ind_d]
    x_matrix_centered_o <- x_matrix_centered[, ind_o]
  }
  
  # We code the formula 83.3 in Thomas and LeSage (2014) 
  zprime_z <- rbind(c(N, zero_k_d, zero_k_o, iota_n %*% G %*% iota_n),
                    cbind(zero_k_d, 
                          n * crossprod(x_matrix_centered_d), 
                          zero_k_d_o, 
                          crossprod(x_matrix_centered_d, G %*% iota_n)),
                    cbind(zero_k_o, 
                          zero_k_o_d, 
                          n * crossprod(x_matrix_centered_o), 
                          crossprod(x_matrix_centered_o, G %*% iota_n)),
                    cbind(iota_n %*% G %*% iota_n, 
                          iota_n %*% crossprod(G, x_matrix_centered_d), 
                          iota_n %*% crossprod(G, x_matrix_centered_o), 
                          sum(diag(crossprod(G))))
  )
  
  # We code the formula 83.4 in Thomas and LeSage (2014) :
  z_prime_y <- rbind(iota_n %*% Y %*% iota_n,
                     crossprod(x_matrix_centered_d, Y %*% iota_n),
                     crossprod(x_matrix_centered_o, t(Y) %*% iota_n),
                     sum(diag(G %*% Y)))
  
  # The estimates of the coefficients is given by formula 83.5 in Thomas and LeSage (2014) :
  hat_delta <- solve(1 / n^2 * zprime_z, 1 / n^2 * z_prime_y)
  rownames(hat_delta) <- c("(Intercept)", names_x_od, "g")
  return(hat_delta)
}