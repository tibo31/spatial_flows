powerWeights <- function (W, rho, order = 250, X, tol = .Machine$double.eps^(3/5)) {
  timings <- list()
  #   .ptime_start <- proc.time()
  n <- dim(W)[1]
  dX <- dim(X)
  if (dX[1] == n) 
    side <- "R"
  else if (dX[2] == n) 
    side <- "L"
  else stop("W and X non-conformant")
  aW <- rho * W
  if (side == "R") 
    last <- aW %*% X
  else last <- X %*% aW
  acc <- X + last
  conv <- FALSE
  iter <- 1
  series <- numeric(order)
  while (iter < order) {
    #  cat(iter)
    if (side == "R") {
      last <- aW %*% last
      acc <- acc + last
    }
    else {
      last <- last %*% aW
      acc <- acc + last
    }
    series[iter] <- mean(last)
    if (series[iter] < tol) {
      conv <- TRUE
      break
    }
    iter <- iter + 1
  }
  if (!conv) 
    warning("not converged within order iterations")
  # timings[["make_power_sum"]] <- proc.time() - .ptime_start
  attr(acc, "internal") <- list(series = series, order = order, 
                                tol = tol, iter = iter, conv = conv)
  attr(acc, "timings") <- do.call("rbind", timings)[, c(1, 
                                                        3)]
  acc
}