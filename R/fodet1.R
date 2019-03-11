fodet1 <- function(parms, traces, n) {
  
  # initialization
  nb_rho <- length(parms)
  # condition 
  stopifnot(nb_rho == 3, all(parms != 0))
  
  # initialization
  titer <- 800
  rez <- numeric(titer)
  
  # precision machine
  eps <- 10^(-50)
  
  # initialisation
  miter <- length(traces)
  tw <- c(traces, n)
  op <- rep(1, 3)
  
  p <- parms
  vvv <- p
  
  for (miteri in 2:miter) {
    vprod <- rep(1, 3^miteri)
    vis <- rep(0, 3^miteri)
    
    for (k in 1:miteri) {
      vvv_temp <- rep(parms, each = 3^(miteri - k))
      vprod <- vprod * vvv_temp
      vis <- vis + as.numeric(vvv_temp == p[1])
    }
    
    vks <- rev(vis)
    vjs <- miteri - vis - vks 
    
    va <- vis + vks
    vb <- vjs + vks
    
    vaa <- va + (va == 0) * (miter + 1)
    vbb <- vb + (vb == 0) * (miter + 1) 
    
    ts <- tw[vaa] * tw[vbb]
    
    tparts <- vprod * ts
    
    tracei <- sum(tparts)
    rez[miteri] <- tracei/miteri
  }
  
  scalarparm <- sum(p)
  ti <- tracei
  
  for (miteri in (miter + 1):titer) {
    ti <- scalarparm * ti
    rez[miteri] <- ti/miteri
    if (abs(rez[miteri] - rez[miteri - 1]) < eps)
      break
  }
  
  return(-sum(rez))
}