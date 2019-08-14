DGP_flow_sdm <- function(z, delta, rho = NULL,
                         W_d = NULL, W_o = NULL, W_w = NULL,
                         seed = NULL, sigma = 1, message = F,
                         model = "") {
  
  # verification
  if (model == "model_2")
    stopifnot(!is.null(W_d), length(rho) == 1)
  
  if (model == "model_3")
    stopifnot(!is.null(W_o), length(rho) == 1)
              
  if (model == "model_4")
    stopifnot(!is.null(W_w), length(rho) == 1)
  
  if (model == "model_5")
    stopifnot(!is.null(W_o) & !is.null(W_d), length(rho) == 1) 
  
  if (model == "model_6") {
    stopifnot(!is.null(W_o) & !is.null(W_d) & !is.null(W_w), length(rho) == 1) 
  }
    
  if (model == "model_7")
    stopifnot(!is.null(W_o) & !is.null(W_d), length(rho) == 2) 
  
  if (model == "model_8") {
    stopifnot(!is.null(W_o) & !is.null(W_d) & !is.null(W_w), length(rho) == 2)   
    rho[3] <- - rho[1] * rho[2] 
  }
  
  if (model == "model_9")
    stopifnot(!is.null(W_o) & !is.null(W_d) & !is.null(W_w), length(rho) == 3)  
  
  # initialisation
  N <- nrow(z)
  z_delta <- z %*% delta 

  # Model 9
  a_w <- switch(model,
              "model_9" = rho[1] * W_d + rho[2] * W_o + rho[3] * W_w,
              "model_8" = rho[1] * W_d + rho[2] * W_o + rho[3] * rho[2] * W_w,
              "model_7" = rho[1] * W_d + rho[2] * W_o, 
              "model_6" = rho[1] * (W_d + W_o + W_w)/3, 
              "model_5" = rho[1] * (W_d + W_o)/2, 
              "model_4" = rho[1] * W_w, 
              "model_3" = rho[1] * W_o, 
              "model_2" = rho[1] * W_d, 
              "model_1" = 0)
   
  # inversion
  a_w <- solve(diag(N) - a_w)

  # compute the signal
  signal <- a_w %*% z_delta 
  sd_signal <- sd(signal)
    
  # simulation of the noise
  repeat_simu <- function (x) {
    set.seed(x)
    eps <- rnorm(N, sd = sigma) 
    noise <- a_w %*% eps
    sd_noise <- sd(noise)
    return(list(y_hat = signal + noise, sd = sd_noise))
  }
  
  if (!is.null(seed)) {
    res_1 <- lapply(seed, repeat_simu)
    res_2 <- sapply(res_1, function(x) x$y_hat)
    sd_noise <- mean(sapply(res_1, function(x) x$sd))
  } else {
    res_1 <- repeat_simu(seed)
    res_2 <- res_1$y_hat
    sd_noise <- res_1$sd
  }
  
  if (message) {
    cat("sd(noise)/sd(signal):\n")
    cat(sd_noise/sd_signal, "\n")
  }
  
  return(res_2)
}