# =======================================================================
# forwardMOA(): Online Conformal Inference with KRR
# =======================================================================

forwardKRR <- function(
    X, y,
    alpha   = 0.1,
    t_init  = 250,
    w       = 250,  # calibration window (residuals for OCI)
    d       = NULL, # data window (sliding KRR)
    lambda  = NULL,
    sigma   = NULL,
    kernel  = "rbf", # kernel: "rbf" or "ntk"
    gammas  = c(0.001, 0.002, 0.004, 0.008, 0.016, 0.032, 0.064, 0.128),
    aci_gam = 0.005,
    methods = "dtaci" # "dtaci", "agaci", "aci", "sfogd", "saocp"
){
  
  n <- nrow(X); stopifnot(t_init < n, w > 0, length(y) == n)
  if (is.null(d)) d <- t_init
  stopifnot(d >= t_init, d <= n)
  
  X_tot <- X; y_tot <- y
  
  ## KRR Init
  model <- krr_init(X_tot[1:t_init, , drop = FALSE],
                    y_tot[1:t_init],
                    lambda = lambda, kernel = kernel)
  
  ## Warm Up
  {
    x_seed      <- X_tot[t_init, , drop = FALSE]
    y_seed      <- y_tot[t_init]
    k_seed      <- compute_kernel_mat(x_seed, model$X, sigma = model$sigma, kernel = model$kernel) 
    yhat_seed   <- drop(k_seed %*% model$Q %*% model$y)
    seed_resid  <- abs(y_seed - yhat_seed)    
  }
  
  ## Index Sets
  pred_start <- t_init + 1
  idx_pred   <- pred_start:n
  T_pred     <- length(idx_pred)
  
  beta_start <- t_init + 1
  idx_beta   <- beta_start:n
  T_beta     <- length(idx_beta); stopifnot(T_beta > 0)
  
  y_hat <- rep(NA_real_, n)
  resid <- rep(NA_real_, n)
  
  for (t in idx_pred) {
    x_new      <- X_tot[t, , drop = FALSE]
    y_new      <- y_tot[t]
    k_x        <- compute_kernel_mat(x_new, model$X, sigma = model$sigma, kernel = model$kernel)
    y_hat[t]   <- drop(k_x %*% model$Q %*% model$y)
    resid[t]   <- abs(y_hat[t] - y_new)
    
    if (t > d) model <- krr_downdate(model)
    model <- krr_update(model, x_new, y_new)
  }
  
  betas <- numeric(T_beta)
  for (j in seq_along(idx_beta)) {
    t <- idx_beta[j]
    win_start <- max(pred_start, t - w)
    win_idx   <- win_start:(t - 1)
    window    <- resid[win_idx]
    
    if (length(window) == 0L) {
      window <- seed_resid
    } else if (win_start == pred_start) {
      window <- c(seed_resid, window)
    }
    
    if (length(window) > w) window <- tail(window, w)
    
    betas[j] <- pmin(1, pmax(0, 1 - mean(resid[t] < window, na.rm = TRUE)))
  }
  
  if (methods == "dtaci") {
    aci_out <- dtaci(betas = betas, alpha = alpha, gammas = gammas, sigma = 1/(n - t_init), etaAdapt = TRUE, etaLookback = (n-t_init))
    alphas  <- pmin(1, pmax(0, aci_out[[1]]))
  } else if (methods == "agaci") {
    aci_out <- agaci(betas = betas, alpha = alpha, gammas = gammas, alphaInit = alpha)
    alphas  <- pmin(1, pmax(0, aci_out[[1]]))
  } else if (methods == "aci") {
    aci_out <- dtaci(betas = betas, alpha = alpha, gammas = aci_gam, etaAdapt = TRUE)
    alphas  <- pmin(1, pmax(0, aci_out[[1]]))
  } else if (methods == "sfogd") {
    aci_out <- sfogd(betas, alpha, gamma = 0.005, alphaInit = alpha)
    alphas  <- pmin(1, pmax(0, aci_out[[1]]))
  } else if (methods == "saocp") {
    aci_out <- saocp(betas, alpha, gamma = 0.005, g = 8)
    alphas  <- pmin(1, pmax(0, aci_out[[1]]))
  }

    L_vec <- U_vec <- numeric(T_beta)
  for (j in seq_along(idx_beta)) {
    t <- idx_beta[j]
    win_start <- max(pred_start, t - w)
    win_idx   <- win_start:(t - 1)
    window    <- resid[win_idx]
    
    if (length(window) == 0L) {
      window <- seed_resid
    } else if (win_start == pred_start) {
      window <- c(seed_resid, window)
    }
    if (length(window) > w) window <- tail(window, w)
    
    w_eff <- length(window)
    q <- stats::quantile(window,
                         probs = min(1, ceiling((w_eff + 1) * (1 - alphas[j])) / w_eff),
                         type = 1, na.rm = TRUE)
    L_vec[j] <- y_hat[t] - q
    U_vec[j] <- y_hat[t] + q
  }
  
  y_test   <- y_tot[idx_beta]
  err_t    <- as.numeric((L_vec > y_test) | (y_test > U_vec))
  coverage <- mean((L_vec <= y_test) & (y_test <= U_vec))
  mean_len <- mean(U_vec - L_vec)
  
  invisible(list(
    L = L_vec, U = U_vec,
    alpha_t = alphas, beta_t = betas,
    coverage = coverage, mean_len = mean_len,
    model = model, y_hat = y_hat, resid = resid,
    d = d, err_t = err_t
  ))
}
