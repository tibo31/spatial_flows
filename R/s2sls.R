s2sls_flow <- function(x, y, g, W) {
  
  # initialization
  n <- nrow(x)
  
  # we add the distance variable to x
  x_matrix <- as(x, "matrix")
  x_matrix <- cbind(x_matrix, g)
  
  # centered data
  x_matrix <- x_matrix - matrix(rep(apply(x_matrix, 2, mean), 
                                    each = n), n, ncol(x_matrix))  
  
  # we add the constant if not included
  if(! any(apply(x_matrix, 2, function(x) all(x == 1)))) {
    x_matrix <- cbind(1, x_matrix)
    colnames(x_matrix)[1] <- "(intercept)"
  }
  
  nvars <- ncol(x_matrix)  # number of x + distance
  names_x <- colnames(x_matrix)
  
  # initialisation
  W_X <- W %*% x_matrix[, -1]
  W2_X <- W %*% W_X
  H_n <- cbind(x_matrix, W_X, W2_X)
  W_Y <- W %*% y
  
  P_H <- H_n %*% chol2inv(qr(H_n)$qr) %*% t(H_n)
  Z_tilde <- cbind(x_matrix, P_H %*% W %*% y)
  
  cste <- chol2inv(qr(Z_tilde)$qr) %*% t(Z_tilde)

  res_k <- cste %*% y
  res_beta <- res_k[1:nvars]
  RHO <- res_k[nvars + 1]  
  hat_Uk <-  y - x_matrix %*% res_beta - W_Y %*% RHO

  SIGMA <- t(hat_Uk) %*% hat_Uk / n
  
  return(list(res_beta = res_beta,
              RHO = RHO,
              SIGMA = SIGMA))
}
