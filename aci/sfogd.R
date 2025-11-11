sfogd <- function(betas, alpha, gamma = 0.01, alphaInit = alpha) {
  T <- length(betas)
  alphaSeq    <- numeric(T)
  errSeqAdapt <- integer(T)

  a  <- alphaInit
  g2 <- 0
  
  for (t in seq_len(T)) {
    b <- betas[t]
    alphaSeq[t] <- a
    
    # (1) Compute Err_t
    miss <- (a > b)
    errSeqAdapt[t] <- as.integer(miss)
    
    # (2) Compute Gradient := -(α - 1{a > b})
    grad <- if (miss) (1 - alpha) else (-alpha)
    
    # (3) Compute Scale-Free Step Size
    g2 <- g2 + grad^2
    step <- gamma / sqrt(g2 + 1e-12)
    
    # (4) Update Alpha
    a <- a - step * grad
    if (a < 1e-6) a <- 1e-6
    if (a > 1 - 1e-6) a <- 1 - 1e-6
  }
  
  list(
    alphaSeq    = alphaSeq,
    errSeqAdapt = errSeqAdapt
  )
}
