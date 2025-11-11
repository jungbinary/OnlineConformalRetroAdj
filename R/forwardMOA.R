# =======================================================================
# forwardMOA(): Online Conformal Inference with MOA
# =======================================================================

forwardMOA <- function(
    X, y,
    alpha   = 0.1,
    t_init  = 250,
    w       = 250,  # calibration window (residuals for OCI)
    gammas  = c(0.001, 0.002, 0.004, 0.008, 0.016, 0.032, 0.064, 0.128),
    moa_model  = "FIMTDD", # base learner: FIMT-DD
    methods    = "dtaci",  # "dtaci", "agaci", "aci", "sfogd", "saocp"
    moa_control = list()   # additional MOA parameters
){
  
  stopifnot(is.matrix(X) || is.data.frame(X))
  X <- as.data.frame(X)
  n <- nrow(X); p <- ncol(X)
  stopifnot(t_init < n, w > 0, length(y) == n)
  
  if (is.null(colnames(X))) colnames(X) <- paste0("x", seq_len(p))
  xcols <- colnames(X)
  
  dat_init <- data.frame(y = y[1:t_init], X[1:t_init, , drop = FALSE])
  ds_init  <- RMOA::datastream_dataframe(dat_init)
  
  base_mod <- RMOA::MOA_regressor(model = moa_model, control = moa_control)
  fit <- RMOA::trainMOA(
    model     = base_mod,
    formula   = y ~ .,
    data      = ds_init,
    chunksize = nrow(dat_init),
    reset     = TRUE,
    trace     = FALSE
  ) 
  
  newx_pred <- data.frame(y = NA_real_, X[1, , drop = FALSE])
  one_row   <- data.frame(y = NA_real_, X[1, , drop = FALSE])
  
  newx_pred[1, xcols] <- X[t_init, , drop = FALSE]
  yhat_seed  <- as.numeric(predict(fit, newdata = newx_pred, type = "response"))
  seed_resid <- abs(y[t_init] - yhat_seed)

  pred_start <- t_init + 1
  idx_pred   <- pred_start:n
  T_pred     <- length(idx_pred)
  
  beta_start <- t_init + 1
  idx_beta   <- beta_start:n
  T_beta     <- length(idx_beta)
  stopifnot(T_beta > 0)
  
  y_hat <- rep(NA_real_, n)
  resid <- rep(NA_real_, n)
  
  for (t in idx_pred) {
    newx_pred[1, xcols] <- X[t, , drop = FALSE]
    y_hat[t]  <- as.numeric(predict(fit, newdata = newx_pred, type = "response"))
    resid[t]  <- abs(y_hat[t] - y[t])
    
    one_row[1, "y"]     <- y[t]
    one_row[1, xcols]   <- X[t, , drop = FALSE]
    ds_one <- RMOA::datastream_dataframe(one_row) 
    
    fit <- RMOA::trainMOA(
      model     = fit$model, 
      formula   = y ~ .,
      data      = ds_one,
      chunksize = 1L,
      reset     = FALSE,
      trace     = FALSE
    )
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
    
    betas[j] <- 1 - mean(resid[t] < window, na.rm = TRUE) 
  }
  
  if (methods == "dtaci") {
    aci_out <- dtaci(betas = betas, alpha = alpha, gammas = gammas, sigma = 1/(n - t_init), etaAdapt = TRUE, etaLookback = (n-t_init))
    alphas  <- pmin(1, pmax(0, aci_out[[1]]))
  } else if (methods == "agaci") {
    aci_out <- agaci(betas = betas, alpha = alpha, gammas = gammas, alphaInit = alpha)
    alphas  <- pmin(1, pmax(0, aci_out[[1]]))
  } else if (methods == "aci") {
    aci_out <- dtaci(betas = betas, alpha = alpha, gammas = 0.005, etaAdapt = TRUE)
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
    q <- stats::quantile(
      window,
      probs = min(1, ceiling((w_eff + 1) * (1 - alphas[j])) / w_eff),
      type = 1, na.rm = TRUE
    )
    L_vec[j] <- y_hat[t] - q
    U_vec[j] <- y_hat[t] + q
  }
  
  y_test   <- y[idx_beta]
  err_t    <- as.numeric((L_vec > y_test) | (y_test > U_vec))
  coverage <- mean((L_vec <= y_test) & (y_test <= U_vec))
  mean_len <- mean(U_vec - L_vec)
  
  invisible(list(
    L = L_vec, U = U_vec,
    alpha_t = alphas, beta_t = betas, err_t = err_t,
    coverage = coverage, mean_len = mean_len,
    model = fit, y_hat = y_hat, resid = resid
  ))
}
