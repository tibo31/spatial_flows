DGP_flow_sdm <- function(z, delta, rho, W_d, W_o, W_w,
                         seed = NULL, sigma = 1, message = F) {
  
  # initialisation
  set.seed(seed)
  N <- nrow(z)
  eps <- rnorm(N, sd = sigma) 
  z_delta <- z %*% delta 
  # gravity 
  delta_gravity <- delta
  delta_gravity[4:5] <- 0 
  z_delta_gravity <- z %*% delta_gravity 
  
  # Model 9 
  irW_9 <- solve(diag(N) - (rho[1] * W_d + rho[2] * W_o + rho[3] * W_w))
  signal_9 <- irW_9 %*% z_delta
  bruit_9 <- irW_9 %*% eps
  y_9 <- signal_9 + bruit_9
  
  irW_8 <- solve(diag(N) - rho[1] * W_d - 
                      rho[2] * W_o + 
                      rho[3] * rho[2] * W_w)
  signal_8 <- irW_8 %*% z_delta
  bruit_8 <- irW_8 %*% eps
  y_8 <- signal_8 + bruit_8
  
  irW_7 <- solve(diag(N) - rho[1] * W_d - 
                      rho[2] * W_o)
  signal_7 <- irW_7 %*% z_delta
  bruit_7 <- irW_7 %*% eps
  y_7 <- signal_7 + bruit_7
  
  irW_6 <- solve(diag(N) - rho[1] * 
                      (W_d + W_o + W_w)/3)
  signal_6 <- irW_6 %*% z_delta
  bruit_6 <- irW_6 %*% eps
  y_6 <- signal_6 + bruit_6
  
  irW_5 <- solve(diag(N) - rho[1] * 
                      (W_d + W_o)/2)
  signal_5 <- irW_5 %*% z_delta
  bruit_5 <- irW_5 %*% eps
  y_5 <- signal_5 + bruit_5
  
  irW_4 <- solve(diag(N) - rho[3] * W_w)
  signal_4 <- irW_4 %*% z_delta
  bruit_4 <- irW_4 %*% eps
  y_4 <- signal_4 + bruit_4
  
  irW_3 <- solve(diag(N) - rho[2] * W_o)
  signal_3 <- irW_3 %*% z_delta
  bruit_3 <- irW_3 %*% eps
  y_3 <- signal_3 + bruit_3
  
  irW_2 <- solve(diag(N) - rho[1] * W_d)
  signal_2 <- irW_2 %*% z_delta
  bruit_2 <- irW_2 %*% eps
  y_2 <- signal_2 + bruit_2
  
  signal_1 <- z_delta
  bruit_1 <- eps
  y_1 <- signal_1 + bruit_1
  
  signal_0 <- z_delta_gravity
  bruit_0 <- eps
  y_0 <- signal_0 + bruit_0
  
  if (message) {
    cat("signal/noise:\n")
    cat("y9: ", mean(signal_9/abs(bruit_9)), "\n")
    cat("y8: ", mean(signal_8/abs(bruit_8)), "\n")
    cat("y7: ", mean(signal_7/abs(bruit_7)), "\n")
    cat("y6: ", mean(signal_6/abs(bruit_6)), "\n")
    cat("y5: ", mean(signal_5/abs(bruit_5)), "\n")
    cat("y4: ", mean(signal_4/abs(bruit_4)), "\n")
    cat("y3: ", mean(signal_3/abs(bruit_3)), "\n")
    cat("y2: ", mean(signal_2/abs(bruit_2)), "\n")
    cat("y1: ", mean(signal_1/abs(bruit_1)), "\n")
    cat("y0: ", mean(signal_0/abs(bruit_0)), "\n")
  }
  
  return(cbind(y_9, y_8, y_7, y_6, y_5, y_4, y_3, y_2, y_1, y_0))
}