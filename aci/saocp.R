vecZeroMin <- Vectorize(function(x) min(x, 0))
pinball <- function(u, alpha) alpha*u - vecZeroMin(u)

L_time <- function(i, g = 8){
  if (any(i <= 0L)) stop("i must be positive")
  n <- 0L; x <- i
  while (x %% 2L == 0L) { x <- x %/% 2L; n <- n + 1L }
  g * (2^n)
}

saocp <- function(betas, alpha, gamma = 1/ sqrt(3),
                  g = 8, alphaInit = alpha) {
  
  T <- length(betas)
  
  alphaSeq     <- numeric(T)  
  errSeqAdapt  <- integer(T)  
  
  experts <- list()
  theta_prev <- alphaInit
  
  for (t in 1:T) {
    beta_t <- betas[t]
    
    experts[[length(experts) + 1L]] <- list(
      theta = theta_prev, start = t, life = L_time(t, g),
      w = 0.0, g2 = 0.0, sum_g = 0.0, sum_wg = 0.0
    )
    
    experts <- Filter(function(e) e$start > (t - e$life) && e$start <= t, experts)
    m <- length(experts); idx <- seq_len(m)
    
    starts <- vapply(experts, function(e) e$start, 0L)
    pi_raw <- (starts^(-2)) / (1 + floor(log(starts, 2)))
    pi_raw[!is.finite(pi_raw)] <- 0
    if (sum(pi_raw) == 0) pi_raw[] <- 1
    pi <- pi_raw / sum(pi_raw)
    
    wts  <- vapply(experts, function(e) e$w, 0.0)
    phat <- pi * pmax(0, wts)
    p <- if (sum(phat) > 0) phat / sum(phat) else pi
    
    thetas_now <- vapply(experts, function(e) e$theta, 0.0)
    theta_t <- if (m == 0L) 0 else as.numeric(sum(p * thetas_now))
    if (t == 1L) theta_t <- 0
    
    alphaSeq[t]     <- theta_t
    errSeqAdapt[t]  <- as.integer(theta_t > beta_t)
    
    L_star <- pinball(beta_t - theta_t, alpha)
    
    for (i in idx) {
      e <- experts[[i]]
      
      L_i   <- pinball(beta_t - e$theta, alpha)
      diffL <- L_star - L_i
      
      p_i_t <- e$w 
      g_i <- if (p_i_t > 0) diffL else max(diffL, 0)
      
      e$sum_g  <- e$sum_g  + g_i
      e$sum_wg <- e$sum_wg + p_i_t * g_i
      
      age <- t - e$start + 1
      e$w <- (1 / age) * e$sum_g * (1 + e$sum_wg)
      
      miss <- (e$theta > beta_t)
      grad_theta <- if (miss) (1 - alpha) else (-alpha)
      e$g2    <- e$g2 + grad_theta^2
      step    <- gamma / sqrt(e$g2)
      e$theta <- e$theta - step * grad_theta
      
      experts[[i]] <- e
    }
    
    theta_prev <- theta_t
  }
  
  list(
    alphaSeq     = alphaSeq,
    errSeqAdapt  = errSeqAdapt
  )
}
