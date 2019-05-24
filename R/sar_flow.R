sar_flow <- function(x, Y, G, w, ind_d = NULL, ind_o = NULL, model = "model_9",
                     lagged = F) {
  
  # x, a data_frame or a matrix with explanatory variable
  # Y, the matrix of flows
  # G, the matrix of distance
  # w, spatial weight matrix observed on the n sample 
  # model, the model choosen
  
  # initialization
  n <- nrow(x)
  p <- ncol(x)
  iota_n <- rep(1, n)
  
  # names of the variable 
  names_x <- colnames(x)
  
  if (is.null(ind_d)) {
    ind_d <- 1:p
    names_x_d <- paste0(names_x, "_d") 
  } else {
    names_x_d <- paste0(names_x[ind_d], "_d") 
  }
  R_d <- length(names_x_d) # number of dest variable
  
  if (is.null(ind_o)) {
    ind_o <- 1:p
    names_x_o <- paste0(names_x, "_o") 
  } else {
    names_x_o <- paste0(names_x[ind_o], "_o") 
  }
  R_o <- length(names_x_o) # number of origin variable
  
  # we center the data: 
  x_matrix <- as(x, "matrix")
  x_matrix_centered <- x_matrix - matrix(rep(apply(x_matrix, 2, mean), 
                                             each = n), n, p)

  if (lagged) {
    x_matrix_lagged <- w %*% x_matrix
    x_matrix_centered_lagged <- x_matrix_lagged - matrix(rep(apply(x_matrix_lagged, 2, mean), 
                                                             each = n), n, p)
    x_matrix_centered_d <- cbind(x_matrix_centered[, ind_d],
                                 x_matrix_centered_lagged[, ind_d])
    x_matrix_centered_o <- cbind(x_matrix_centered[, ind_o],
                                 x_matrix_centered_lagged[, ind_o])
    # rename with the lagged variable
    names_x_d <- c(names_x_d, paste0("lagged_", names_x_d))
    names_x_o <- c(names_x_o, paste0("lagged_", names_x_o))    
  } else {
    x_matrix_centered_d <- x_matrix_centered[, ind_d]
    x_matrix_centered_o <- x_matrix_centered[, ind_o]
  }
  
  R_d_with_lagged <- length(names_x_d)
  R_o_with_lagged <- length(names_x_o)
  
  # names of the result variable
  names_x_od <- c(names_x_d, names_x_o) 
  
  if (!lagged) {
    nvars <- R_d + R_o + 2 # + 2 for intercept and distance
  } else {
    nvars <- 2 * (R_d + R_o) + 2
  }
  # we center the matrix of distance
  G_dot <- G - mean(G)
  
  # number of rho 
  names_rho <- switch(substr(model, 7, 7),
                      "9" = c("rho_d", "rho_o", "rho_w"),
                      "8" = c("rho_d", "rho_o", "rho_w"),
                      "7" = c("rho_d", "rho_o"),
                      "6" = "rho_odw",
                      "5" = "rho_od",
                      "4" = "rho_w",
                      "3" = "rho_o",
                      "2" = "rho_d",
                      "1" = NULL)
  
  nb_rho <- length(names_rho)
  
  if (nb_rho == 0) {
    return(gravity_model(x, Y, G, lagged = lagged, w = w))
  } 
  
  # initialization of rho
  pvec <- runif(nb_rho)
  rho <- 0.7 * pvec / sum(pvec) 
  
  if (model == "model_8") {
    rho[3] <- - rho[1] * rho[2]
  }
  
  # initialization of bayesian parameters
  # number of replicates
  ndraw <- 5500
  # to store the results
  bsave <- matrix(0, ndraw, nvars)
  psave <- matrix(0, ndraw, nb_rho)
  ssave <- rep(0, ndraw)
  acc_rate <- matrix(0, ndraw, nb_rho)
  # Parameters related to the prior distribution
  # V <- matrix(1, n, n)
  sige <- 1
  
  cc1 <- 0.2  # tuning parameter, equation (5.28) in Lesage book
  acc1 <- 0
  
  if (nb_rho == 2 | nb_rho == 3) {
    cc2 <- cc1
    acc2 <- 0
  }
  
  if (nb_rho == 3 & model == "model_9") {
    cc3 <- cc1
    acc3 <- 0
  }
  
  rmin <- -1
  rmax <- 1
  # Definition of some matrices:
  WY <- w %*% Y
  YW <- Y %*% t(w)
  WYW <- w %*% Y %*% t(w)
  # we store some computational matrices
  zpzt <- matrix(0, nvars, nvars)
  zpzty1 <- zpzty2 <- zpzty3 <- zpzty4 <- numeric(nvars)
  
  for (i in 1:n) {
    
    z <- cbind(iota_n, 
               cbind(x_matrix_centered_d, 
                     matrix(x_matrix_centered_o[i, ], n, R_o_with_lagged, byrow = T),
                     G_dot[i, ])
               )
    
    fy <- switch(substr(model, 7, 7),   # (8.15) in LeSage book                                
                 "9" = cbind(Y[, i], WY[, i], YW[, i], WYW[, i]),
                 "8" = cbind(Y[, i], WY[, i], YW[, i], WYW[, i]),
                 "7" = cbind(Y[, i], WY[, i], YW[, i]),
                 "6" = cbind(Y[, i], (WY[, i] + YW[, i] + WYW[, i])/3),
                 "5" = cbind(Y[, i], 0.5 * (WY[, i] + YW[, i])),
                 "4" = cbind(Y[, i], WYW[, i]),
                 "3" = cbind(Y[, i], YW[, i]),
                 "2" = cbind(Y[, i], WY[, i]),
                 "1" = NULL)
    
    zt <- z 
    zpzt <- zpzt + crossprod(z)
    zpzty1 <- zpzty1 + crossprod(zt, fy[, 1])
    zpzty2 <- zpzty2 + crossprod(zt, fy[, 2]) 
    if (ncol(fy) > 2) 
      zpzty3 <- zpzty3 + crossprod(zt, fy[, 3])
    if (ncol(fy) > 3)
      zpzty4 <- zpzty4 + crossprod(zt, fy[, 4])     
  }
  
  zpzti <- solve(zpzt)
  
  # Bayesian algorithm
  pb <- txtProgressBar(min = 0, max = ndraw, initial = 0, 
                       char = "=", width = (getOption("width")), 
                       style = 1)
  counter <- 0
  
  for (iter in 1:ndraw) {
    tau <- c(1, -rho)
    
    bdraw1 <- as.numeric(mvtnorm::rmvnorm(1, sigma = sige * zpzti, 
                                          method = "chol")) + zpzti %*% zpzty1
    bdraw2 <- as.numeric(mvtnorm::rmvnorm(1, sigma = sige * zpzti, 
                                          method = "chol")) + zpzti %*% zpzty2
    bdraw <- cbind(bdraw1, bdraw2)
    
    if (nb_rho == 2) {
      bdraw3 <- as.numeric(mvtnorm::rmvnorm(1, sigma = sige * zpzti, 
                                            method = "chol")) + zpzti %*% zpzty3 
      bdraw <- cbind(bdraw, bdraw3)
    } else {
      if (nb_rho == 3) {
        bdraw3 <- as.numeric(mvtnorm::rmvnorm(1, sigma = sige * zpzti, 
                                              method = "chol")) + zpzti %*% zpzty3 
        bdraw4 <- as.numeric(mvtnorm::rmvnorm(1, sigma = sige * zpzti, 
                                              method = "chol")) + zpzti %*% zpzty4 
        bdraw <- cbind(bdraw, bdraw3, bdraw4)
      }
    }
    
    
    # update for beta ends here
      beta_draw <- bdraw %*% tau
    
    # update for sige starts here
    for (ij in 1:(nb_rho + 1)) {
      beta_draw_temp <- switch(as.character(ij),
                               "1" = bdraw1,
                               "2" = bdraw2,
                               "3" = bdraw3,
                               "4" = bdraw4)
      alpha <- beta_draw_temp[1]
      bd <- beta_draw_temp[2:(2 + R_d_with_lagged - 1)]
      bo <- beta_draw_temp[(2 + R_d_with_lagged):(2 + R_d_with_lagged + R_o_with_lagged - 1)]
      gamma_coeff <- beta_draw_temp[length(beta_draw_temp)]
      
      xdb <- kronecker(iota_n, as(x_matrix_centered_d, "matrix") %*% bd) 
      xob <- kronecker(as(x_matrix_centered_o, "matrix") %*% bo, iota_n)
     
      y_hat <- alpha + xob + xdb + as.vector(G_dot) * gamma_coeff
      
      if (ij == 1)
        E1 <- matrix(as.vector(Y) - y_hat, n, n)
      
      if (ij == 2)
        E2 <- matrix(as.vector(WY) - y_hat, n, n)
      
      if (ij == 3)
        E3 <- matrix(as.vector(YW) - y_hat, n, n)
      
      if (ij == 4)
        E4 <- matrix(as.vector(WYW) - y_hat, n, n)
    }
    
    Q <- matrix(0, nb_rho + 1, nb_rho + 1)
    Q[1, 1] <- sum(E1^2)
    Q[1, 2] <- Q[2, 1] <- sum(E1 * E2)  
    Q[2, 2] <- sum(E2^2)  
    
    if (nb_rho == 2 | nb_rho == 3) {
      Q[1, 3] <- Q[3, 1] <- sum(E1 * E3)
      Q[2, 3] <- Q[3, 2] <- sum(E2 * E3)
      Q[3, 3] <- sum(E3^2)
    }
    
    if (nb_rho == 3) {
      Q[1, 4] <- Q[4, 1] <- sum(E1 * E4)
      Q[2, 4] <- Q[4, 2] <- sum(E2 * E4)
      Q[3, 4] <- Q[4, 3] <- sum(E3 * E4)
      Q[4, 4] <- sum(E4^2)
    }
    
    epe <- t(tau) %*% Q %*% tau # Sum of Square residuals (p. 222, in Lesage book) 
    
    nu <- 0 
    d0 <- 0
    nu1 <- n * n + 2 * nu
    d1 <- 2 * d0 + epe
    chi <- rgamma(1, nu1 * 0.5) * 2
    sige <- as.numeric(d1/chi) 
    # update for sige ends here
    
    # update for rho1, rho2, rho3 starts here
    # update rho1 using metropolis hastings
    rho_x <- c_sarf(rho, sige, Q, traces, n, nvars)
    accept <- T
    rho1_c <- rho[1] + cc1*rnorm(1)
    
    while (accept) {
      if ((rho1_c > rmin) & (rho1_c < rmax))
        accept <- F
      else
        rho1_c <- rho[1] + cc1*rnorm(1)
    }
    
    if (nb_rho == 1)
      rho_temp <- rho1_c
    
    if (nb_rho == 2)
      rho_temp <- c(rho1_c, rho[2])
    
    if (nb_rho == 3)
      rho_temp <- c(rho1_c, rho[2], rho[3])
    
    if (model == "model_8")
      rho_temp[3] <- - rho_temp[1] * rho_temp[2]
    
    rho_y <- c_sarf(rho_temp, sige, Q, traces, n, nvars)
    ru <- runif(1)
    
    if ((rho_y - rho_x) > exp(1)) {
      p <- 1
    } else {
      ratio <- exp(rho_y - rho_x)
      p <- min(1, ratio)
    }
    
    if (ru < p) {
      rho[1] <- rho1_c
      acc1 <- acc1 + 1
    }
    
    acc_rate[iter, 1] <- acc1/iter
    
    if (acc_rate[iter, 1] < 0.4) {
      cc1 <- cc1 / 1.1
    } else {
      if (acc_rate[iter, 1] > 0.6)
        cc1 <- cc1 * 1.1
    }
    
    if (model == "model_8")
      rho[3] <- - rho[1] * rho[2]
    
    # update rho2 using metropolis hastings
    if (nb_rho == 2 | nb_rho == 3) {
      rho_x <- c_sarf(rho, sige, Q, traces, n, nvars)
      accept <- T
      rho2_c <- rho[2] + cc2*rnorm(1)
      
      while (accept) {
        if ((rho2_c > rmin) & (rho2_c < rmax))
          accept <- F
        else
          rho2_c <- rho[2] + cc2*rnorm(1)
      }
      
      if (nb_rho == 2)
        rho_temp <- c(rho[1], rho2_c)
      
      if (nb_rho == 3)
        rho_temp <- c(rho[1], rho2_c, rho[3])
      
      if (model == "model_8")
        rho_temp[3] <- - rho_temp[1] * rho_temp[2]
      
      rho_y <- c_sarf(rho_temp, sige, Q, traces, n, nvars)
      ru <- runif(1)
      
      if ((rho_y - rho_x) > exp(1)) {
        p <- 1
      } else {
        ratio <- exp(rho_y - rho_x)
        p <- min(1, ratio)
      }
      
      if (ru < p) {
        rho[2] <- rho2_c
        acc2 <- acc2 + 1
      }
      
      acc_rate[iter, 2] <- acc2/iter
      
      if (acc_rate[iter, 2] < 0.4) {
        cc2 <- cc2 / 1.1
      } else {
        if (acc_rate[iter, 2] > 0.6)
          cc2 <- cc2 * 1.1
      }
    }
    
    if (model == "model_8")
      rho[3] <- - rho[1] * rho[2]
    
    # update rho3 using metropolis hastings
    if (nb_rho == 3 & model != "model_8") {
      rho_x <- c_sarf(rho, sige, Q, traces, n, nvars)
      accept <- T
      rho3_c <- rho[3] + cc3*rnorm(1)
      
      while (accept) {
        if ((rho3_c > rmin) & (rho3_c < rmax))
          accept <- F
        else
          rho3_c <- rho[3] + cc3*rnorm(1)
      }
      
      rho_temp <- c(rho[1], rho[2], rho3_c)
      rho_y <- c_sarf(rho_temp, sige, Q, traces, n, nvars)
      ru <- runif(1)
      
      if ((rho_y - rho_x) > exp(1)) {
        p <- 1
      } else {
        ratio <- exp(rho_y - rho_x)
        p <- min(1, ratio)
      }
      
      if (ru < p) {
        rho[3] <- rho3_c
        acc3 <- acc3 + 1
      }
      
      acc_rate[iter, 3] <- acc3/iter
      
      if (acc_rate[iter, 3] < 0.4) {
        cc3 <- cc3 / 1.1
      } else {
        if (acc_rate[iter, 3] > 0.6)
          cc3 <- cc3 * 1.1
      }
      
    }
    
    # save the results
    bsave[iter, ] <- beta_draw
    ssave[iter] <- sige
    psave[iter, ] <- rho
    
    counter <- counter + 1
    setTxtProgressBar(pb, counter)
  }
  
  res_sar <- data.frame(mean = c(apply(psave, 2, mean), apply(bsave, 2, mean)),
                        lower_05 = c(apply(psave, 2, quantile, 0.05), 
                                     apply(bsave, 2, quantile, 0.05)),
                        lower_95 = c(apply(psave, 2, quantile, 0.95), 
                                     apply(bsave, 2, quantile, 0.95)),
                        row.names = c(names_rho, "(intercept)", 
                                      names_x_od, "g"))
  res_sar$t_stat <- res_sar$mean/c(apply(psave, 2, sd), apply(bsave, 2, sd)) 
  
  return(res_sar)
}