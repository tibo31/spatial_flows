sar_flow_2 <- function(x, y, g, W_d = NULL, W_o = NULL, W_w = NULL, 
                       ind_d = NULL, ind_o = NULL, model = "", centered = T) {
  
  # x, a data_frame or a matrix with explanatory variable of full size
  # Y, the vector of flows
  # G, the vector of distances
  # w, spatial weight matrix observed on the N sample 
  # model, the model choosen
  
  # initialization
  n <- nrow(x)
  
  # we add the distance variable to x
  x_matrix <- as(x, "matrix")
  x_matrix <- cbind(x_matrix, g)
  
  # centered data
  if (centered) {
    x_matrix <- x_matrix - matrix(rep(apply(x_matrix, 2, mean), 
                                    each = n), n, ncol(x_matrix))  
  }
  
  # we add the constant if not included
  if(! any(apply(x_matrix, 2, function(x) all(x == 1)))) {
    x_matrix <- cbind(1, x_matrix)
    colnames(x_matrix)[1] <- "(intercept)"
  }
  
  nvars <- ncol(x_matrix)  # number of x + distance
  names_x <- colnames(x_matrix)
  
  # verification
  stopifnot(n == length(y), n == length(g))
  
  # determine the good model 
  if (is.null(W_d) && is.null(W_o) && is.null(W_w))
    names_rho <- NULL
  
  if (!is.null(W_d) && is.null(W_o) && is.null(W_w)) {
    names_rho <- "rho_d"
    y2 <- W_d %*% y
    zpzty2 <- crossprod(x_matrix, y2) 
    W_2 <- W_d
  }
  
  if (is.null(W_d) && !is.null(W_o) && is.null(W_w)) {
    names_rho <- "rho_o" 
    y2 <- W_o %*% y
    zpzty2 <- crossprod(x_matrix, y2)
    W_2 <- W_o
  }
  
  if (is.null(W_d) && is.null(W_o) && !is.null(W_w)) {
    names_rho <- "rho_w"  
    y2 <- W_w %*% y 
    zpzty2 <- crossprod(x_matrix, y2)  
    W_2 <- W_w
  }
  
  if (!is.null(W_d) && !is.null(W_o) && is.null(W_w)) {
    if (model != "model_5") {
      names_rho <- c("rho_d", "rho_o")
      y2 <- W_d %*% y
      zpzty2 <- crossprod(x_matrix, y2)
      W_2 <- W_d
      y3 <- W_o %*% y
      zpzty3 <- crossprod(x_matrix, y3)
      W_3 <- W_o
    } else {
      names_rho <- c("rho_od")
      W_2 <- 0.5 * (W_d + W_o)
      y2 <- W_2 %*% y
      zpzty2 <- crossprod(x_matrix, y2)
    }
  }
  
  if (!is.null(W_d) && !is.null(W_o) && !is.null(W_w)) {
    if (model == "model_6") {
      names_rho <- c("rho_odw")
      W_2 <- 1 / 3 * (W_d + W_o + W_w)
      y2 <- W_2 %*% y
      zpzty2 <- crossprod(x_matrix, y2)
    } else {
      names_rho <- c("rho_d", "rho_o", "rho_w")
      y2 <- W_d %*% y
      zpzty2 <- crossprod(x_matrix, y2) 
      W_2 <- W_d
      y3 <- W_o %*% y
      zpzty3 <- crossprod(x_matrix, y3)
      W_3 <- W_o
      y4 <- W_w %*% y 
      zpzty4 <- crossprod(x_matrix, y4) 
      W_4 <- W_w
    }
  }
  
  # number of rho 
  nb_rho <- length(names_rho)
  
  if (nb_rho == 0) {
    return(lm(y ~ x + g))
  }
  
  # initialization of rho
  pvec <- runif(nb_rho)
  rho <- 0.7 * pvec / sum(pvec) 
  
  if (model == "model_8")
    rho[3] <- - rho[1] * rho[2]
  
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
  
  if (nb_rho == 3 && model != "model_8") {
    cc3 <- cc1
    acc3 <- 0
  }
  
  rmin <- -1
  rmax <- 1
  
  # x-prime*w
  x_prime_x <- crossprod(x_matrix)
  zpzti <- solve(x_prime_x)
  zpzty1 <- crossprod(x_matrix, y)
  
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
    
    beta_draw <- as(bdraw, "matrix") %*% tau
    # update for beta ends here
    
    # update for sige starts here
    y_hat <- x_matrix %*% bdraw[, 1] 
    E1 <- y - y_hat
    y_hat <- x_matrix %*% bdraw[, 2]     
    E2 <- y2 - y_hat
    
    Q <- matrix(0, nb_rho + 1, nb_rho + 1)
    Q[1, 1] <- sum(E1^2)
    Q[1, 2] <- Q[2, 1] <- sum(E1 * E2)  
    Q[2, 2] <- sum(E2^2)  
    
    if (nb_rho == 2 | nb_rho == 3) {
      y_hat <- x_matrix %*% bdraw[, 3] 
      E3 <- y3 - y_hat
      Q[1, 3] <- Q[3, 1] <- sum(E1 * E3)
      Q[2, 3] <- Q[3, 2] <- sum(E2 * E3)
      Q[3, 3] <- sum(E3^2)
    }
    
    if (nb_rho == 3) {
      y_hat <- x_matrix %*% bdraw[, 4] 
      E4 <- y4 - y_hat
      Q[1, 4] <- Q[4, 1] <- sum(E1 * E4)
      Q[2, 4] <- Q[4, 2] <- sum(E2 * E4)
      Q[3, 4] <- Q[4, 3] <- sum(E3 * E4)
      Q[4, 4] <- sum(E4^2)
    }
    
    epe <- t(tau) %*% Q %*% tau # Sum of Square residuals (p. 222, in Lesage book) 
    
    nu <- 0 
    d0 <- 0
    nu1 <- n + 2 * nu
    d1 <- 2 * d0 + epe
    chi <- rgamma(1, nu1 * 0.5) * 2
    sige <- as.numeric(d1/chi) 
    # update for sige ends here
    
    # update for rho1, rho2, rho3 starts here
    # update rho1 using metropolis hastings
    logdet <- as.numeric(determinant(
      switch(as.character(nb_rho), 
             "1" = Diagonal(n) -  rho[1] * W_2, 
             "2" = Diagonal(n) -  rho[1] * W_2 - rho[2] * W_3,
             "3" = Diagonal(n) -  rho[1] * W_2 - rho[2] * W_3 - rho[3] * W_4), 
      logarithm = TRUE)$modulus)
    
    rho_x <- logdet - epe/(2*sige)
    
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
    
    logdet_y <- as.numeric(determinant(
      switch(as.character(nb_rho), 
             "1" = Diagonal(n) -  rho_temp[1] * W_2, 
             "2" = Diagonal(n) -  rho_temp[1] * W_2 - rho_temp[2] * W_3,
             "3" = Diagonal(n) -  rho_temp[1] * W_2 - rho_temp[2] * W_3 - rho_temp[3] * W_4), 
      logarithm = TRUE)$modulus)
    
    tau_y <- c(1, -rho_temp)
    epe_y <- t(tau_y) %*% Q %*% tau_y
    rho_y <- logdet_y - epe_y/(2*sige)
    
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
      logdet_2 <- as.numeric(determinant(
        switch(as.character(nb_rho), 
               "2" = Diagonal(n) -  rho[1] * W_2 - rho[2] * W_3,
               "3" = Diagonal(n) -  rho[1] * W_2 - rho[2] * W_3 - rho[3] * W_4), 
        logarithm = TRUE)$modulus) 
      
      
      
      
      tau_x <- c(1, -rho)
      epe_x <- t(tau_x) %*% Q %*% tau_x
      rho_x <- logdet_2 - epe_x/(2 * sige)
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
      
      logdet_y <- as.numeric(determinant(
        switch(as.character(nb_rho), 
               "2" = Diagonal(n) -  rho_temp[1] * W_2 - rho_temp[2] * W_3,
               "3" = Diagonal(n) -  rho_temp[1] * W_2 - rho_temp[2] * W_3 - rho_temp[3] * W_4), 
        logarithm = TRUE)$modulus) 
      
      tau_y <- c(1, -rho_temp)
      epe_y <- t(tau_y) %*% Q %*% tau_y
      rho_y <- logdet_y - epe_y/(2*sige)
      
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
    if (nb_rho == 3 && model != "model_8") {
      logdet_2 <- as.numeric(determinant(
        Diagonal(n) -  rho[1] * W_2 - rho[2] * W_3 - rho[3] * W_4, 
        logarithm = TRUE)$modulus) 
      
      tau_x <- c(1, -rho)
      epe_x <- t(tau_x) %*% Q %*% tau_x
      rho_x <- logdet_2 - epe_x/(2 * sige)
      accept <- T
      rho3_c <- rho[3] + cc3*rnorm(1)
      
      while (accept) {
        if ((rho3_c > rmin) & (rho3_c < rmax))
          accept <- F
        else
          rho3_c <- rho[3] + cc3*rnorm(1)
      }
      
      rho_temp <- c(rho[1], rho[2], rho3_c)
      logdet_y <- as.numeric(determinant(
        Diagonal(n) -  rho_temp[1] * W_2 - rho_temp[2] * W_3 - rho_temp[3] * W_4, 
        logarithm = TRUE)$modulus) 
      
      tau_y <- c(1, -rho_temp)
      epe_y <- t(tau_y) %*% Q %*% tau_y
      rho_y <- logdet_y - epe_y/(2*sige)
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
    bsave[iter, ] <- as.numeric(beta_draw)
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
                        row.names = c(names_rho, names_x))
  res_sar$t_stat <- res_sar$mean/c(apply(psave, 2, sd), apply(bsave, 2, sd)) 
  res_sar
}