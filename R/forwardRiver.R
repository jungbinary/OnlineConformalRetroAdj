# =======================================================================
# forwardRiver(): Online Conformal Inference with River
# =======================================================================

.forwardRiver_locate_module <- function(py_module_path = NULL) {
  candidates <- character(0)
  
  if (!is.null(py_module_path)) {
    if (dir.exists(py_module_path)) {
      candidates <- c(candidates, file.path(py_module_path, "river_helpers.py"))
    } else {
      candidates <- c(candidates, py_module_path)
    }
  } else {
    wd <- getwd()
    candidates <- c(
      file.path(wd, "python", "river_helpers.py"),
      file.path(wd, "..", "python", "river_helpers.py"),
      file.path(wd, "..", "..", "python", "river_helpers.py")
    )
    recursive_hits <- list.files(
      wd,
      pattern = "^river_helpers\\.py$",
      recursive = TRUE,
      full.names = TRUE
    )
    candidates <- c(candidates, recursive_hits)
  }
  
  candidates <- unique(candidates[file.exists(candidates)])
  if (!length(candidates)) {
    stop(
      "Could not locate 'python/river_helpers.py'. ",
      "Set `py_module_path` to the helper file or its containing directory.",
      call. = FALSE
    )
  }
  
  normalizePath(candidates[[1]], winslash = "/", mustWork = TRUE)
}

.forwardRiver_import_helper <- function(py_module_path = NULL) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' is required for `forwardRiver()`.", call. = FALSE)
  }
  
  module_file <- .forwardRiver_locate_module(py_module_path)
  module_dir  <- dirname(module_file)
  module_name <- tools::file_path_sans_ext(basename(module_file))
  
  tryCatch(
    {
      reticulate::import("importlib", delay_load = FALSE)$invalidate_caches()
      reticulate::import_from_path(
        module = module_name,
        path = module_dir,
        delay_load = FALSE,
        convert = FALSE
      )
    },
    error = function(e) {
      stop(
        "Failed to import Python helper module '", module_name,
        "' from '", module_dir, "': ", conditionMessage(e),
        call. = FALSE
      )
    }
  )
}

.forwardRiver_py_to_r <- function(x) {
  reticulate::py_to_r(x)
}

.forwardRiver_row_to_list <- function(df_row, xcols) {
  row_list <- as.list(df_row[1, xcols, drop = FALSE])
  names(row_list) <- xcols
  
  lapply(row_list, function(value) {
    value <- value[[1]]
    if (is.factor(value)) {
      as.character(value)
    } else if (inherits(value, "integer64")) {
      as.numeric(value)
    } else {
      value
    }
  })
}

forwardRiver <- function(
    X, y,
    alpha   = 0.1,
    t_init  = 250,
    w       = 250,
    gammas  = c(0.001, 0.002, 0.004, 0.008, 0.016, 0.032, 0.064, 0.128),
    river_model   = "ARFRegressor",
    methods       = "dtaci",
    river_control = list(),
    py_module_path = NULL
){
  
  stopifnot(is.matrix(X) || is.data.frame(X))
  X <- as.data.frame(X)
  n <- nrow(X); p <- ncol(X)
  stopifnot(t_init < n, w > 0, length(y) == n)
  
  if (is.null(colnames(X))) colnames(X) <- paste0("x", seq_len(p))
  xcols <- colnames(X)
  
  if (!is.list(river_control)) {
    stop("`river_control` must be a list.", call. = FALSE)
  }
  if (length(river_control) && is.null(names(river_control))) {
    stop("`river_control` must be a named list.", call. = FALSE)
  }
  
  river_helper <- .forwardRiver_import_helper(py_module_path = py_module_path)
  river_control_arg <- if (length(river_control)) river_control else NULL
  
  dat_init <- data.frame(y = y[1:t_init], X[1:t_init, , drop = FALSE])
  x_init   <- lapply(
    seq_len(nrow(dat_init)),
    function(i) .forwardRiver_row_to_list(dat_init[i, , drop = FALSE], xcols)
  )
  
  base_mod <- .forwardRiver_py_to_r(
    river_helper$create_model(
      model_name = river_model,
      control = river_control_arg
    )
  )
  fit <- .forwardRiver_py_to_r(
    river_helper$fit_rows(
      model = base_mod,
      rows = x_init,
      y = as.numeric(dat_init$y)
    )
  )
  
  newx_pred <- data.frame(y = NA_real_, X[1, , drop = FALSE])
  one_row   <- data.frame(y = NA_real_, X[1, , drop = FALSE])
  
  newx_pred[1, xcols] <- X[t_init, , drop = FALSE]
  yhat_seed  <- as.numeric(.forwardRiver_py_to_r(
    river_helper$predict_one(fit, .forwardRiver_row_to_list(newx_pred, xcols))
  ))
  if (length(yhat_seed) != 1L || is.na(yhat_seed)) {
    stop(
      "River helper returned a non-scalar or NA prediction during warm-up.",
      call. = FALSE
    )
  }
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
    yhat_t <- as.numeric(.forwardRiver_py_to_r(
      river_helper$predict_one(fit, .forwardRiver_row_to_list(newx_pred, xcols))
    ))
    if (length(yhat_t) != 1L || is.na(yhat_t)) {
      stop(
        "River helper returned a non-scalar or NA prediction at t = ", t, ".",
        call. = FALSE
      )
    }
    y_hat[t]  <- yhat_t
    resid[t]  <- abs(y_hat[t] - y[t])
    
    one_row[1, "y"]     <- y[t]
    one_row[1, xcols]   <- X[t, , drop = FALSE]
    
    fit <- .forwardRiver_py_to_r(
      river_helper$learn_one(
        model = fit,
        row = .forwardRiver_row_to_list(one_row, xcols),
        y = y[t]
      )
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
