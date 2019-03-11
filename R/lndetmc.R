lndetmc <- function(parms, traces, n) {
  
  # condition
  stopifnot(length(parms) == 1)
  
  # initialization
  m <- length(traces)
  
  res <- traces * parms^(1:m)/(1:m) 
  
  return(-n * sum(res))
}