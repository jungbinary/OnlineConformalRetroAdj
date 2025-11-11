# =============================================================================
# conformalRetroAdj(): Online Conformal Inference with Retrospective Adjustment
# =============================================================================

conformalRetroAdj <- function(
    X, y,
    alpha    = 0.1,
    t_init   = 250,
    lambda   = NULL,
    sigma    = NULL,
    kernel   = "rbf",  # "rbf" or "ntk"
    d        = NULL,   # data window (sliding KRR); defaults to t_init inside
    downdate = TRUE,   # enable sliding-window downdate
    methods  = "dtaci",# "dtaci","agaci","aci","sfogd","saocp"
    aci_gam  = 0.005,
    gammas   = c(0.001, 0.002, 0.004, 0.008, 0.016, 0.032, 0.064, 0.128)
) {
  
  ## 0. Safeguard
  n <- length(y)
  p <- dim(X)[2]
  stopifnot(t_init < n)
  
  ## 1. Data
  X_tot <- X
  y_tot <- y
  
  win_size <- if (is.null(d)) {
    t_init
  } else if (d < t_init) {
    d
  } else {
    d
  }
  
  ## 2. KRR Initialization
  if (!is.null(d) && d < t_init) {
    model <- krr_init(
      X_tot[(t_init - d + 1):t_init, , drop = FALSE],
      y_tot[(t_init - d + 1):t_init],
      lambda = lambda,
      kernel = kernel
    )
  } else {
    model <- krr_init(
      X_tot[1:t_init, , drop = FALSE],
      y_tot[1:t_init],
      lambda = lambda,
      kernel = kernel
    )
  }
  
  ## 3. Online Loop
  L_list <- vector("list", n - t_init)
  U_list <- vector("list", n - t_init)
  betas  <- numeric(n - t_init)
  
  for (tt in seq_len(n - t_init)) {
    idx   <- t_init + tt
    x_new <- X_tot[idx, , drop = FALSE]
    y_new <- y_tot[idx]
    
    # Jackknife+ residuals -> beta_t
    resid     <- compute_resid_jk(model, x_new)
    betas[tt] <- compute_beta_jk(y_new, resid$L, resid$U)
    L_list[[tt]] <- resid$L
    U_list[[tt]] <- resid$U
    
    if (downdate && (length(model$y) >= win_size)) {
      model <- krr_downdate(model)
    }
    model <- krr_update(model, x_new, y_new)
  }
  
  ## 4. DtACI
  if (methods == "dtaci") {
    aci_out <- dtaci(betas = betas, alpha = alpha, gammas = gammas, sigma = 1/(n - t_init), etaAdapt = TRUE, etaLookback = (n-t_init))
    alphas  <- pmin(1, pmax(0, aci_out[[1]]))
  } else if (methods == "agaci") {
    aci_out <- agaci(betas = betas, alpha = alpha, gammas = gammas, alphaInit = alpha)
    alphas  <- pmin(1, pmax(0, aci_out[[1]]))
  } else if (methods == "aci") {
    aci_out <- dtaci(betas = betas, alpha = alpha, gammas = aci_gam)
    alphas  <- pmin(1, pmax(0, aci_out[[1]]))
  } else if (methods == "sfogd") {
    aci_out <- sfogd(betas, alpha, gamma = 0.005, alphaInit = alpha)
    alphas  <- pmin(1, pmax(0, aci_out[[1]]))
  } else if (methods == "saocp") {
    aci_out <- saocp(betas, alpha, gamma = 0.005, g = 8)
    alphas  <- pmin(1, pmax(0, aci_out[[1]]))
  }
  
  ## 5. Interval Construction
  L_vec <- vapply(seq_along(L_list), function(i) {
    n_i <- length(L_list[[i]])
    quantile(L_list[[i]],
             probs = min(1, floor((n_i + 1) * alphas[i]) / n_i),
             type  = 1)
  }, numeric(1))
  
  U_vec <- vapply(seq_along(U_list), function(i) {
    n_i <- length(U_list[[i]])
    quantile(U_list[[i]],
             probs = min(1, ceiling((n_i + 1) * (1 - alphas[i])) / n_i),
             type  = 1)
  }, numeric(1))
  
  ## 6. Evaluation
  y_test   <- y_tot[(t_init + 1):n]
  coverage <- mean((L_vec <= y_test) & (y_test <= U_vec))
  len_mean <- mean(U_vec - L_vec)
  err_t    <- as.numeric((L_vec > y_test) | (y_test > U_vec))
  
  ## 7. Return
  invisible(list(
    L        = L_vec,
    U        = U_vec,
    alpha_t  = alphas,
    beta_t   = betas,
    err_t    = err_t,
    coverage = coverage,
    mean_len = len_mean,
    model    = model,
    R        = resid$R
  ))
}


