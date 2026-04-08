# =============================================================================
# conformalRetroAdjNaiveRefit():
# Separate naive retrospective-adjustment implementation for runtime comparison.
# =============================================================================

if (!exists("compute_kernel_mat", mode = "function") ||
    !exists("krr_predict", mode = "function")) {
  stop("Source R/baseKRR.R before sourcing conformalRetroAdj_naive.R")
}

krr_fit_fixed <- function(X_fit, y_fit, lambda, sigma, kernel = "rbf") {
  stopifnot(!is.null(lambda), !is.null(sigma))
  stopifnot(is.matrix(X_fit) || is.data.frame(X_fit))

  X_fit <- as.matrix(X_fit)
  n_fit <- nrow(X_fit)
  stopifnot(length(y_fit) == n_fit, n_fit >= 1L)

  k_mat <- compute_kernel_mat(X_fit, X_fit, sigma = sigma, kernel = kernel)
  q_mat <- solve(k_mat + lambda * diag(n_fit))

  list(
    X = X_fit,
    y = y_fit,
    Q = q_mat,
    lambda = lambda,
    sigma = sigma,
    kernel = kernel,
    cv_mse = NA_real_
  )
}

compute_resid_jk_naive <- function(model, x_new) {
  x_train <- model$X
  y_train <- model$y
  n_train <- nrow(x_train)

  signed_resid <- numeric(n_train)
  f_new_loo <- numeric(n_train)

  for (i in seq_len(n_train)) {
    keep_idx <- seq_len(n_train) != i
    loo_model <- krr_fit_fixed(
      X_fit = x_train[keep_idx, , drop = FALSE],
      y_fit = y_train[keep_idx],
      lambda = model$lambda,
      sigma = model$sigma,
      kernel = model$kernel
    )

    y_hat_i <- krr_predict(loo_model, x_train[i, , drop = FALSE])
    signed_resid[i] <- y_train[i] - y_hat_i
    f_new_loo[i] <- krr_predict(loo_model, x_new)
  }

  list(
    R = abs(signed_resid),
    L = f_new_loo - abs(signed_resid),
    U = f_new_loo + abs(signed_resid)
  )
}

conformalRetroAdjNaiveRefit <- function(
    X, y,
    alpha = 0.1,
    t_init = 250,
    lambda = NULL,
    sigma = NULL,
    kernel = "rbf",
    d = NULL,
    downdate = TRUE,
    methods = "dtaci",
    aci_gam = 0.005,
    gammas = c(0.001, 0.002, 0.004, 0.008, 0.016, 0.032, 0.064, 0.128)
) {
  n <- length(y)
  stopifnot(t_init < n)

  X_tot <- X
  y_tot <- y

  win_size <- if (is.null(d)) {
    t_init
  } else if (d < t_init) {
    d
  } else {
    d
  }

  if (!is.null(d) && d < t_init) {
    model <- krr_init(
      X_tot[(t_init - d + 1L):t_init, , drop = FALSE],
      y_tot[(t_init - d + 1L):t_init],
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

  l_list <- vector("list", n - t_init)
  u_list <- vector("list", n - t_init)
  betas <- numeric(n - t_init)

  for (tt in seq_len(n - t_init)) {
    idx <- t_init + tt
    x_new <- X_tot[idx, , drop = FALSE]
    y_new <- y_tot[idx]

    resid <- compute_resid_jk_naive(model, x_new)
    betas[tt] <- compute_beta_jk(y_new, resid$L, resid$U)
    l_list[[tt]] <- resid$L
    u_list[[tt]] <- resid$U

    if (downdate && (length(model$y) >= win_size)) {
      model <- krr_downdate(model)
    }
    model <- krr_update(model, x_new, y_new)
  }

  if (methods == "dtaci") {
    aci_out <- dtaci(
      betas = betas,
      alpha = alpha,
      gammas = gammas,
      sigma = 1 / (n - t_init),
      etaAdapt = TRUE,
      etaLookback = (n - t_init)
    )
    alphas <- pmin(1, pmax(0, aci_out[[1]]))
  } else if (methods == "agaci") {
    aci_out <- agaci(betas = betas, alpha = alpha, gammas = gammas, alphaInit = alpha)
    alphas <- pmin(1, pmax(0, aci_out[[1]]))
  } else if (methods == "aci") {
    aci_out <- dtaci(betas = betas, alpha = alpha, gammas = aci_gam)
    alphas <- pmin(1, pmax(0, aci_out[[1]]))
  } else if (methods == "sfogd") {
    aci_out <- sfogd(betas, alpha, gamma = 0.005, alphaInit = alpha)
    alphas <- pmin(1, pmax(0, aci_out[[1]]))
  } else if (methods == "saocp") {
    aci_out <- saocp(betas, alpha, gamma = 0.005, g = 8)
    alphas <- pmin(1, pmax(0, aci_out[[1]]))
  } else {
    stop("Unknown ACI method: ", methods)
  }

  l_vec <- vapply(seq_along(l_list), function(i) {
    n_i <- length(l_list[[i]])
    quantile(
      l_list[[i]],
      probs = min(1, floor((n_i + 1L) * alphas[i]) / n_i),
      type = 1
    )
  }, numeric(1))

  u_vec <- vapply(seq_along(u_list), function(i) {
    n_i <- length(u_list[[i]])
    quantile(
      u_list[[i]],
      probs = min(1, ceiling((n_i + 1L) * (1 - alphas[i])) / n_i),
      type = 1
    )
  }, numeric(1))

  y_test <- y_tot[(t_init + 1L):n]
  coverage <- mean((l_vec <= y_test) & (y_test <= u_vec))
  len_mean <- mean(u_vec - l_vec)
  err_t <- as.numeric((l_vec > y_test) | (y_test > u_vec))

  invisible(list(
    L = l_vec,
    U = u_vec,
    alpha_t = alphas,
    beta_t = betas,
    err_t = err_t,
    coverage = coverage,
    mean_len = len_mean,
    model = model,
    R = resid$R
  ))
}
