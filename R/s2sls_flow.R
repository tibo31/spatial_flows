s2sls_flow <- function(x_d, x_o, y, g, 
                       instru_x_d = !logical(ncol(x_d)),
                       instru_x_o = !logical(ncol(x_o)),
                       W_d = NULL, W_o = NULL, W_w = NULL,
                       constant = TRUE) {

  # we add the distance variable to x
  x_d_matrix <- as(x_d, "matrix")
  x_o_matrix <- as(x_o, "matrix")
  x_matrix <- cbind(x_d_matrix, x_o_matrix, g)
  n <- nrow(x_matrix)
  # centered data
  # x_matrix <- x_matrix - matrix(rep(apply(x_matrix, 2, mean), 
  #                                  each = n), n, ncol(x_matrix))  
  
  # we add the constant if not included
  if (constant) {
    if(! any(apply(x_matrix, 2, function(x) all(x == 1)))) {
      x_matrix <- cbind(1, x_matrix)
      colnames(x_matrix)[1] <- "(intercept)"
    }
  }
  
  nvars <- ncol(x_matrix)  # number of x + distance
  names_x <- colnames(x_matrix)
  
  # initialisation
  if (!is.null(W_d)) {
    W_X_d <- W_d %*% cbind(x_d_matrix[, instru_x_d], g)
    W2_X_d <- W_d %*% W_X_d
    W_Y_d <- W_d %*% y
    names_rho_d <- "rho_d"
  } else {
    W_X_d <- NULL
    W2_X_d <- NULL
    W_Y_d <- NULL
    names_rho_d <- NULL
  }  
    
  if (!is.null(W_o)) {
    W_X_o <- W_o %*% cbind(x_o_matrix[, instru_x_o], g)
    W2_X_o <- W_o %*% W_X_o
    W_Y_o <- W_o %*% y
    names_rho_o <- "rho_o"
  } else {
    W_X_o <- NULL
    W2_X_o <- NULL
    W_Y_o <- NULL
    names_rho_o <- NULL
  }  
  
  if (!is.null(W_w)) {
    W_X_w <- W_w %*%  g
    W2_X_w <- W_w %*% W_X_w
    W_Y_w <- W_w %*% y
    names_rho_w <- "rho_w"
  } else {
    W_X_w <- NULL
    W2_X_w <- NULL
    W_Y_w <- NULL
    names_rho_w <- NULL
  }  
  
  W_X <- cbind(W_X_d, W_X_o, W_X_w)
  W2_X <- cbind(W2_X_d, W2_X_o, W2_X_w)
  H_n <- cbind(x_matrix, W_X, W2_X)
  W_Y <- cbind(W_Y_d, W_Y_o, W_Y_w)
  names_rho <- c(names_rho_d, names_rho_o, names_rho_w)
  
  
  P_H <- H_n %*% chol2inv(qr(H_n)$qr) %*% t(H_n)
  Z_tilde <- cbind(x_matrix, P_H %*% W_Y)

  cste <- chol2inv(qr(Z_tilde)$qr) %*% t(Z_tilde)

  res_k <- cste %*% y
  res_beta <- res_k[1:nvars]
  RHO <- res_k[(nvars + 1):length(res_k)]  
  hat_Uk <-  y - x_matrix %*% res_beta - W_Y %*% RHO

  SIGMA <- t(hat_Uk) %*% hat_Uk / n
  
  # compute the sd
  Z <- cbind(x_matrix, W_Y)
  H <- as.matrix(H_n)
  zhat <- H %*% (MASS::ginv(crossprod(H)) %*% crossprod(H, as.matrix(Z)))
  beta_rho <- c(res_beta, RHO)
  yhat <- Z %*% beta_rho
  res <- y - yhat
  resvar <- sum(res ^ 2) / length(y)
  varcoef <- resvar * MASS::ginv(crossprod(zhat))
  std <- sqrt(diag(varcoef))
  st <- beta_rho/std
  
  names(res_beta) <- names_x
  sd_beta <- st[1:ncol(x_matrix)]
  names(sd_beta) <- names_x
  sd_rho <- st[(ncol(x_matrix) + 1):length(st)]
  names(sd_rho) <- names_rho  
  names(RHO) <- names_rho
  
  return(list(res_beta = res_beta,
              RHO = RHO,
              SIGMA = SIGMA,
              sd_beta = sd_beta,
              sd_rho = sd_rho)
  )
}



