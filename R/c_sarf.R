c_sarf <- function(rho, sige, Q, traces, n, nvars) {
  # nvars not necessary
  
  nb_rho <- length(rho) 
  
  # check if we estimate the model (8)
  if (nb_rho == 3) {
    if (rho[3] == - rho[1] * rho[2])
      nb_rho <- 2
  }
  
  logdet <- switch(as.character(nb_rho),
                   "3" = fodet1(rho, traces, n), 
                   "2" = lndetmc(rho[1], traces, n) + lndetmc(rho[2], traces, n),
                   "1" = lndetmc(rho, traces, n))
  
  tau <- c(1, -rho)
  epe <- t(tau) %*% Q %*% tau
  res <- logdet - epe/(2*sige)
  
  return(res)
}