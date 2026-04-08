suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(scales)
})

# Benchmark defaults ----------------------------------------------------

runtime_t_init <- 250L
runtime_alpha <- 0.1
runtime_kernel <- "rbf"
runtime_aci_method <- "dtaci"
runtime_scale_windows <- c(100L, 250L, 500L, 750L, 1000L)

runtime_end_to_end_window <- 100L
runtime_end_to_end_reps <- 5L
runtime_end_to_end_timeout_sec <- 3600

runtime_scale_block_updates <- 50L
runtime_scale_warmup_updates <- 1L
runtime_scale_reps <- 5L
runtime_scale_timeout_sec <- 300
runtime_scale_data_window_krr <- runtime_t_init

default_runtime_workers <- max(
  1L,
  min(6L, parallel::detectCores(logical = FALSE))
)

# Output paths ----------------------------------------------------------

runtime_output_dir <- file.path("experiments", "results", "appendix", "runtime_benchmark")
runtime_notes_dir <- file.path(runtime_output_dir, "notes")
runtime_results_dir <- file.path(runtime_output_dir, "results")
runtime_tables_dir <- file.path(runtime_output_dir, "tables")
runtime_figures_dir <- file.path(runtime_output_dir, "figures")

end_to_end_raw_path <- file.path(runtime_results_dir, "runtime_end_to_end_raw.csv")
update_only_raw_path <- file.path(runtime_results_dir, "runtime_update_only_raw.csv")
end_to_end_summary_path <- file.path(runtime_tables_dir, "runtime_end_to_end_summary.csv")
update_only_summary_path <- file.path(runtime_tables_dir, "runtime_update_only_summary.csv")
timeout_summary_path <- file.path(runtime_tables_dir, "runtime_timeout_summary.csv")
main_table_path <- file.path(runtime_tables_dir, "runtime_main_table.csv")
main_table_tex_path <- file.path(runtime_tables_dir, "runtime_main_table.tex")
scale_table_tex_path <- file.path(runtime_tables_dir, "runtime_scalability_table.tex")
protocol_note_path <- file.path(runtime_notes_dir, "runtime_protocol.txt")
protocol_tex_path <- file.path(runtime_notes_dir, "runtime_protocol.tex")
scalability_plot_pdf_path <- file.path(runtime_figures_dir, "runtime_scalability.pdf")
scalability_plot_png_path <- file.path(runtime_figures_dir, "runtime_scalability.png")
session_info_path <- file.path(runtime_notes_dir, "session_info.txt")

# Method metadata -------------------------------------------------------

runtime_method_specs <- list(
  retro_eff = list(
    label = "RetroAdj",
    family = "retro",
    residual_fun = compute_resid_jk
  ),
  retro_naive = list(
    label = "RetroAdj-naive",
    family = "retro",
    residual_fun = compute_resid_jk_naive
  ),
  fw_krr = list(
    label = "FW-KRR",
    family = "fw_krr"
  ),
  fw_amf = list(
    label = "FW-AMF",
    family = "fw_river",
    river_model = "AMFRegressor",
    river_control = list(seed = 1L)
  ),
  fw_fimtdd = list(
    label = "FW-FIMT-DD",
    family = "fw_moa",
    moa_model = "FIMTDD"
  ),
  fw_amrules = list(
    label = "FW-AMRules",
    family = "fw_moa",
    moa_model = "AMRulesRegressor"
  )
)
runtime_method_keys <- names(runtime_method_specs)
runtime_method_order <- vapply(runtime_method_specs, `[[`, character(1), "label")

runtime_method_palette <- c(
  "RetroAdj" = "#EE6A70",
  "RetroAdj-naive" = "#4C78A8",
  "FW-KRR" = "#B8A9C9",
  "FW-AMF" = "#9DBCE3",
  "FW-FIMT-DD" = "#90C8AC",
  "FW-AMRules" = "#E3C57B"
)

runtime_base_theme <- theme_bw(base_size = 15) +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.text = element_text(size = 13),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 14),
    plot.margin = margin(t = 8, r = 8, b = 8, l = 8)
  )

dir.create(runtime_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(runtime_notes_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(runtime_results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(runtime_tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(runtime_figures_dir, recursive = TRUE, showWarnings = FALSE)

# Setup helpers ---------------------------------------------------------

configure_single_thread_math <- function() {
  Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1"
  )
  compiler::enableJIT(3L)
  invisible(NULL)
}

source_runtime_dependencies <- function() {
  suppressPackageStartupMessages(library(RMOA))
  if (!nzchar(Sys.getenv("RETICULATE_PYTHON"))) {
    py_candidates <- unique(c(Sys.which("python")))
    py_candidates <- py_candidates[nzchar(py_candidates) & file.exists(py_candidates)]
    if (length(py_candidates)) {
      Sys.setenv(
        RETICULATE_PYTHON = normalizePath(py_candidates[[1]], winslash = "/", mustWork = TRUE)
      )
    }
  }
  source("R/baseKRR.R", encoding = "UTF-8")
  source("aci/dtaci.R", encoding = "UTF-8")
  source("aci/agaci.R", encoding = "UTF-8")
  source("aci/sfogd.R", encoding = "UTF-8")
  source("aci/saocp.R", encoding = "UTF-8")
  source("R/conformalRetroAdj.R", encoding = "UTF-8")
  source("R/forwardKRR.R", encoding = "UTF-8")
  source("R/forwardMOA.R", encoding = "UTF-8")
  source("R/forwardRiver.R", encoding = "UTF-8")
  source("experiments/appendix/runtime_benchmark/conformalRetroAdj_naive.R", encoding = "UTF-8")
  invisible(NULL)
}

resolve_existing_path <- function(candidates) {
  for (candidate in candidates) {
    if (file.exists(candidate)) {
      return(candidate)
    }
  }
  stop("Could not find any of: ", paste(candidates, collapse = ", "))
}

load_elec2_runtime_data <- function(csv_path = "data/real/elec2.csv") {
  dat <- read.csv(csv_path)
  dat <- na.omit(dat)
  y_all <- dat[["nswprice"]]

  lag_k <- 10L
  x_embed <- embed(y_all, lag_k + 1L)
  g <- rbind(
    x_embed[1:250, , drop = FALSE],
    x_embed[(nrow(x_embed) - 2750 + 1L):nrow(x_embed), , drop = FALSE]
  )

  list(
    dataset = "Elec2",
    dataset_slug = "elec2",
    X = g[, -1, drop = FALSE],
    y = g[, 1],
    n_obs = nrow(g)
  )
}

load_aig_runtime_data <- function(
    csv_path = resolve_existing_path(c("data/real/AIG.csv", file.path(getwd(), "data", "real", "AIG.csv")))
) {
  dat <- read.csv(csv_path)
  idx <- (9012L - 3000L):(9012L + 3000L)
  y_all <- dat[idx, "Close"]

  lag_k <- 10L
  x_embed <- embed(y_all, lag_k + 1L)

  list(
    dataset = "AIG",
    dataset_slug = "aig",
    X = x_embed[, -1, drop = FALSE],
    y = x_embed[, 1],
    n_obs = nrow(x_embed)
  )
}

load_runtime_dataset <- function(dataset_slug) {
  switch(
    tolower(dataset_slug),
    elec2 = load_elec2_runtime_data(),
    aig = load_aig_runtime_data(),
    stop("Unknown dataset: ", dataset_slug)
  )
}

.runtime_krr_hyperparams <- new.env(parent = emptyenv())

get_fixed_krr_hyperparams <- function(dataset_slug) {
  key <- tolower(dataset_slug)
  if (exists(key, envir = .runtime_krr_hyperparams, inherits = FALSE)) {
    return(get(key, envir = .runtime_krr_hyperparams, inherits = FALSE))
  }

  dat <- load_runtime_dataset(key)
  base_model <- krr_init(
    dat$X[1:runtime_t_init, , drop = FALSE],
    dat$y[1:runtime_t_init],
    kernel = runtime_kernel
  )

  params <- list(lambda = base_model$lambda, sigma = base_model$sigma)
  assign(key, params, envir = .runtime_krr_hyperparams)
  params
}

method_spec <- function(method_key) {
  spec <- runtime_method_specs[[method_key]]
  if (is.null(spec)) {
    stop("Unknown method_key: ", method_key)
  }
  spec
}

# Timing utilities ------------------------------------------------------

time_with_timeout <- function(expr_fun, timeout_sec = Inf) {
  gc(verbose = FALSE)
  start_elapsed <- Sys.time()
  status <- "completed"
  message_text <- NA_character_
  value <- NULL

  on.exit({
    setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
  }, add = TRUE)

  if (is.finite(timeout_sec)) {
    setTimeLimit(elapsed = timeout_sec, transient = TRUE)
  }

  value <- tryCatch(
    expr_fun(),
    error = function(e) {
      cond_msg <- conditionMessage(e)
      if (grepl("elapsed time limit", cond_msg, fixed = TRUE)) {
        status <<- "timeout"
      } else {
        status <<- "error"
      }
      message_text <<- cond_msg
      NULL
    }
  )

  elapsed <- as.numeric(difftime(Sys.time(), start_elapsed, units = "secs"))
  list(status = status, elapsed = elapsed, value = value, message = message_text)
}

time_update_loop_with_timeout <- function(expr_fun, timeout_sec = Inf) {
  gc(verbose = FALSE)
  start_elapsed <- Sys.time()
  status <- "completed"
  message_text <- NA_character_

  on.exit({
    setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
  }, add = TRUE)

  if (is.finite(timeout_sec)) {
    setTimeLimit(elapsed = timeout_sec, transient = TRUE)
  }

  tryCatch(
    expr_fun(),
    error = function(e) {
      cond_msg <- conditionMessage(e)
      if (grepl("elapsed time limit", cond_msg, fixed = TRUE)) {
        status <<- "timeout"
      } else {
        status <<- "error"
      }
      message_text <<- cond_msg
      NULL
    }
  )

  elapsed <- as.numeric(difftime(Sys.time(), start_elapsed, units = "secs"))
  list(status = status, elapsed = elapsed, message = message_text)
}

empty_runtime_raw_df <- function() {
  data.frame(
    benchmark = character(),
    dataset = character(),
    dataset_slug = character(),
    method = character(),
    method_key = character(),
    window = integer(),
    repetition = integer(),
    status = character(),
    timeout_sec = numeric(),
    total_runtime_sec = numeric(),
    runtime_per_step_ms = numeric(),
    timed_updates = integer(),
    updates_per_block = integer(),
    timing_repeats = integer(),
    message = character(),
    stringsAsFactors = FALSE
  )
}

make_runtime_row <- function(
    benchmark,
    dat,
    method_key,
    window,
    repetition,
    timed_updates,
    status,
    timeout_sec,
    elapsed,
    updates_per_block = NA_integer_,
    timing_repeats = NA_integer_,
    message = NA_character_
) {
  spec <- method_spec(method_key)
  data.frame(
    benchmark = benchmark,
    dataset = dat$dataset,
    dataset_slug = dat$dataset_slug,
    method = spec$label,
    method_key = method_key,
    window = as.integer(window),
    repetition = as.integer(repetition),
    status = status,
    timeout_sec = timeout_sec,
    total_runtime_sec = elapsed,
    runtime_per_step_ms = if (!is.na(elapsed) && !is.na(timed_updates) && timed_updates > 0L) {
      (elapsed * 1000) / timed_updates
    } else {
      NA_real_
    },
    timed_updates = as.integer(timed_updates),
    updates_per_block = as.integer(updates_per_block),
    timing_repeats = as.integer(timing_repeats),
    message = ifelse(is.na(message), "", message),
    stringsAsFactors = FALSE
  )
}

# Update-only benchmark -------------------------------------------------

retro_interval_quantiles <- function(L_vals, U_vals, alpha_fixed = runtime_alpha) {
  n_i <- length(L_vals)
  lower_q <- stats::quantile(
    L_vals,
    probs = min(1, floor((n_i + 1L) * alpha_fixed) / n_i),
    type = 1
  )
  upper_q <- stats::quantile(
    U_vals,
    probs = min(1, ceiling((n_i + 1L) * (1 - alpha_fixed)) / n_i),
    type = 1
  )
  c(lower = lower_q, upper = upper_q)
}

forward_interval_quantile <- function(window_resid, alpha_fixed = runtime_alpha) {
  w_eff <- length(window_resid)
  stats::quantile(
    window_resid,
    probs = min(1, ceiling((w_eff + 1L) * (1 - alpha_fixed)) / w_eff),
    type = 1,
    na.rm = TRUE
  )
}

prepare_retro_state <- function(dat, window, benchmark_steps, hyperparams) {
  timed_end <- dat$n_obs
  timed_start <- timed_end - benchmark_steps + 1L
  train_end <- timed_start - 1L
  train_start <- train_end - window + 1L
  stopifnot(train_start >= 1L)

  model <- krr_init(
    dat$X[train_start:train_end, , drop = FALSE],
    dat$y[train_start:train_end],
    lambda = hyperparams$lambda,
    sigma = hyperparams$sigma,
    kernel = runtime_kernel
  )

  list(
    model = model,
    timed_indices = seq.int(timed_start, timed_end),
    active_window = as.integer(window)
  )
}

run_retro_update_block <- function(dat, method_key, window, benchmark_steps, hyperparams, timeout_sec) {
  spec <- method_spec(method_key)
  state <- prepare_retro_state(dat, window, benchmark_steps, hyperparams)
  model <- state$model
  timed_indices <- state$timed_indices

  timing <- time_update_loop_with_timeout(
    expr_fun = function() {
      for (idx in timed_indices) {
        x_new <- dat$X[idx, , drop = FALSE]
        y_new <- dat$y[idx]
        resid <- spec$residual_fun(model, x_new)
        interval <- retro_interval_quantiles(resid$L, resid$U, alpha_fixed = runtime_alpha)
        beta <- compute_beta_jk(y_new, resid$L, resid$U)

        if (length(model$y) >= window) {
          model <- krr_downdate(model)
        }
        model <- krr_update(model, x_new, y_new)
      }
    },
    timeout_sec = timeout_sec
  )

  if (exists("resid", inherits = FALSE)) {
    rm(resid)
  }
  if (exists("interval", inherits = FALSE)) {
    rm(interval)
  }
  if (exists("beta", inherits = FALSE)) {
    rm(beta)
  }
  if (exists("model", inherits = FALSE)) {
    rm(model)
  }
  gc(verbose = FALSE)
  timing
}

prepare_fw_krr_state <- function(dat, window, benchmark_steps, hyperparams) {
  timed_end <- dat$n_obs
  timed_start <- timed_end - benchmark_steps + 1L
  warmup_end <- timed_start - 1L
  warmup_start <- warmup_end - window + 1L
  init_end <- warmup_start - 1L
  init_start <- init_end - runtime_scale_data_window_krr + 1L
  stopifnot(init_start >= 1L)

  model <- krr_init(
    dat$X[init_start:init_end, , drop = FALSE],
    dat$y[init_start:init_end],
    lambda = hyperparams$lambda,
    sigma = hyperparams$sigma,
    kernel = runtime_kernel
  )

  residual_buffer <- numeric(0)
  for (idx in seq.int(warmup_start, warmup_end)) {
    x_new <- dat$X[idx, , drop = FALSE]
    y_new <- dat$y[idx]
    y_hat <- krr_predict(model, x_new)
    resid_t <- abs(y_hat - y_new)
    residual_buffer <- c(residual_buffer, resid_t)
    if (length(residual_buffer) > window) {
      residual_buffer <- tail(residual_buffer, window)
    }

    if (length(model$y) >= runtime_scale_data_window_krr) {
      model <- krr_downdate(model)
    }
    model <- krr_update(model, x_new, y_new)
  }

  list(
    model = model,
    residual_buffer = residual_buffer,
    timed_indices = seq.int(timed_start, timed_end)
  )
}

run_fw_krr_update_block <- function(dat, window, benchmark_steps, hyperparams, timeout_sec) {
  state <- prepare_fw_krr_state(dat, window, benchmark_steps, hyperparams)
  model <- state$model
  residual_buffer <- state$residual_buffer
  timed_indices <- state$timed_indices

  timing <- time_update_loop_with_timeout(
    expr_fun = function() {
      for (idx in timed_indices) {
        x_new <- dat$X[idx, , drop = FALSE]
        y_new <- dat$y[idx]
        q <- forward_interval_quantile(residual_buffer, alpha_fixed = runtime_alpha)
        y_hat <- krr_predict(model, x_new)
        lower <- y_hat - q
        upper <- y_hat + q
        resid_t <- abs(y_hat - y_new)
        beta <- pmin(1, pmax(0, 1 - mean(resid_t < residual_buffer, na.rm = TRUE)))

        residual_buffer <- c(residual_buffer, resid_t)
        if (length(residual_buffer) > window) {
          residual_buffer <- tail(residual_buffer, window)
        }

        if (length(model$y) >= runtime_scale_data_window_krr) {
          model <- krr_downdate(model)
        }
        model <- krr_update(model, x_new, y_new)
      }
    },
    timeout_sec = timeout_sec
  )

  if (exists("q", inherits = FALSE)) {
    rm(q)
  }
  if (exists("lower", inherits = FALSE)) {
    rm(lower)
  }
  if (exists("upper", inherits = FALSE)) {
    rm(upper)
  }
  if (exists("resid_t", inherits = FALSE)) {
    rm(resid_t)
  }
  if (exists("beta", inherits = FALSE)) {
    rm(beta)
  }
  if (exists("model", inherits = FALSE)) {
    rm(model)
  }
  if (exists("residual_buffer", inherits = FALSE)) {
    rm(residual_buffer)
  }
  gc(verbose = FALSE)
  timing
}

prepare_fw_moa_state <- function(dat, window, benchmark_steps, moa_model) {
  X_df <- as.data.frame(dat$X)
  p <- ncol(X_df)
  if (is.null(colnames(X_df))) {
    colnames(X_df) <- paste0("x", seq_len(p))
  }
  xcols <- colnames(X_df)

  timed_end <- dat$n_obs
  timed_start <- timed_end - benchmark_steps + 1L
  warmup_end <- timed_start - 1L
  warmup_start <- warmup_end - window + 1L
  init_end <- warmup_start - 1L
  init_start <- init_end - runtime_t_init + 1L
  stopifnot(init_start >= 1L)

  dat_init <- data.frame(y = dat$y[init_start:init_end], X_df[init_start:init_end, , drop = FALSE])
  ds_init <- RMOA::datastream_dataframe(dat_init)
  base_mod <- RMOA::MOA_regressor(model = moa_model, control = list())
  fit <- RMOA::trainMOA(
    model = base_mod,
    formula = y ~ .,
    data = ds_init,
    chunksize = nrow(dat_init),
    reset = TRUE,
    trace = FALSE
  )

  newx_pred <- data.frame(y = NA_real_, X_df[1, , drop = FALSE])
  one_row <- data.frame(y = NA_real_, X_df[1, , drop = FALSE])
  residual_buffer <- numeric(0)

  for (idx in seq.int(warmup_start, warmup_end)) {
    newx_pred[1, xcols] <- X_df[idx, , drop = FALSE]
    y_hat <- as.numeric(predict(fit, newdata = newx_pred, type = "response"))
    resid_t <- abs(y_hat - dat$y[idx])
    residual_buffer <- c(residual_buffer, resid_t)
    if (length(residual_buffer) > window) {
      residual_buffer <- tail(residual_buffer, window)
    }

    one_row[1, "y"] <- dat$y[idx]
    one_row[1, xcols] <- X_df[idx, , drop = FALSE]
    ds_one <- RMOA::datastream_dataframe(one_row)
    fit <- RMOA::trainMOA(
      model = fit$model,
      formula = y ~ .,
      data = ds_one,
      chunksize = 1L,
      reset = FALSE,
      trace = FALSE
    )
  }

  list(
    fit = fit,
    X_df = X_df,
    xcols = xcols,
    newx_pred = newx_pred,
    one_row = one_row,
    residual_buffer = residual_buffer,
    timed_indices = seq.int(timed_start, timed_end)
  )
}

run_fw_moa_update_block <- function(dat, window, benchmark_steps, moa_model, timeout_sec) {
  state <- prepare_fw_moa_state(dat, window, benchmark_steps, moa_model)
  fit <- state$fit
  X_df <- state$X_df
  xcols <- state$xcols
  newx_pred <- state$newx_pred
  one_row <- state$one_row
  residual_buffer <- state$residual_buffer
  timed_indices <- state$timed_indices

  timing <- time_update_loop_with_timeout(
    expr_fun = function() {
      for (idx in timed_indices) {
        newx_pred[1, xcols] <- X_df[idx, , drop = FALSE]
        q <- forward_interval_quantile(residual_buffer, alpha_fixed = runtime_alpha)
        y_hat <- as.numeric(predict(fit, newdata = newx_pred, type = "response"))
        lower <- y_hat - q
        upper <- y_hat + q
        resid_t <- abs(y_hat - dat$y[idx])
        beta <- pmin(1, pmax(0, 1 - mean(resid_t < residual_buffer, na.rm = TRUE)))

        residual_buffer <- c(residual_buffer, resid_t)
        if (length(residual_buffer) > window) {
          residual_buffer <- tail(residual_buffer, window)
        }

        one_row[1, "y"] <- dat$y[idx]
        one_row[1, xcols] <- X_df[idx, , drop = FALSE]
        ds_one <- RMOA::datastream_dataframe(one_row)
        fit <- RMOA::trainMOA(
          model = fit$model,
          formula = y ~ .,
          data = ds_one,
          chunksize = 1L,
          reset = FALSE,
          trace = FALSE
        )
      }
    },
    timeout_sec = timeout_sec
  )

  if (exists("q", inherits = FALSE)) {
    rm(q)
  }
  if (exists("lower", inherits = FALSE)) {
    rm(lower)
  }
  if (exists("upper", inherits = FALSE)) {
    rm(upper)
  }
  if (exists("resid_t", inherits = FALSE)) {
    rm(resid_t)
  }
  if (exists("beta", inherits = FALSE)) {
    rm(beta)
  }
  if (exists("fit", inherits = FALSE)) {
    rm(fit)
  }
  if (exists("residual_buffer", inherits = FALSE)) {
    rm(residual_buffer)
  }
  gc(verbose = FALSE)
  timing
}

prepare_fw_river_state <- function(dat, window, benchmark_steps, river_model, river_control) {
  X_df <- as.data.frame(dat$X)
  p <- ncol(X_df)
  if (is.null(colnames(X_df))) {
    colnames(X_df) <- paste0("x", seq_len(p))
  }
  xcols <- colnames(X_df)

  timed_end <- dat$n_obs
  timed_start <- timed_end - benchmark_steps + 1L
  warmup_end <- timed_start - 1L
  warmup_start <- warmup_end - window + 1L
  init_end <- warmup_start - 1L
  init_start <- init_end - runtime_t_init + 1L
  stopifnot(init_start >= 1L)

  river_helper <- .forwardRiver_import_helper(py_module_path = "python")
  river_control_arg <- if (length(river_control)) river_control else NULL

  dat_init <- data.frame(y = dat$y[init_start:init_end], X_df[init_start:init_end, , drop = FALSE])
  x_init <- lapply(
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

  newx_pred <- data.frame(y = NA_real_, X_df[1, , drop = FALSE])
  one_row <- data.frame(y = NA_real_, X_df[1, , drop = FALSE])
  residual_buffer <- numeric(0)

  for (idx in seq.int(warmup_start, warmup_end)) {
    newx_pred[1, xcols] <- X_df[idx, , drop = FALSE]
    y_hat <- as.numeric(.forwardRiver_py_to_r(
      river_helper$predict_one(fit, .forwardRiver_row_to_list(newx_pred, xcols))
    ))
    resid_t <- abs(y_hat - dat$y[idx])
    residual_buffer <- c(residual_buffer, resid_t)
    if (length(residual_buffer) > window) {
      residual_buffer <- tail(residual_buffer, window)
    }

    one_row[1, "y"] <- dat$y[idx]
    one_row[1, xcols] <- X_df[idx, , drop = FALSE]
    fit <- .forwardRiver_py_to_r(
      river_helper$learn_one(
        model = fit,
        row = .forwardRiver_row_to_list(one_row, xcols),
        y = dat$y[idx]
      )
    )
  }

  list(
    fit = fit,
    river_helper = river_helper,
    X_df = X_df,
    xcols = xcols,
    newx_pred = newx_pred,
    one_row = one_row,
    residual_buffer = residual_buffer,
    timed_indices = seq.int(timed_start, timed_end)
  )
}

run_fw_river_update_block <- function(dat, window, benchmark_steps, river_model, river_control, timeout_sec) {
  state <- prepare_fw_river_state(dat, window, benchmark_steps, river_model, river_control)
  fit <- state$fit
  river_helper <- state$river_helper
  X_df <- state$X_df
  xcols <- state$xcols
  newx_pred <- state$newx_pred
  one_row <- state$one_row
  residual_buffer <- state$residual_buffer
  timed_indices <- state$timed_indices

  timing <- time_update_loop_with_timeout(
    expr_fun = function() {
      for (idx in timed_indices) {
        newx_pred[1, xcols] <- X_df[idx, , drop = FALSE]
        q <- forward_interval_quantile(residual_buffer, alpha_fixed = runtime_alpha)
        y_hat <- as.numeric(.forwardRiver_py_to_r(
          river_helper$predict_one(fit, .forwardRiver_row_to_list(newx_pred, xcols))
        ))
        lower <- y_hat - q
        upper <- y_hat + q
        resid_t <- abs(y_hat - dat$y[idx])
        beta <- pmin(1, pmax(0, 1 - mean(resid_t < residual_buffer, na.rm = TRUE)))

        residual_buffer <- c(residual_buffer, resid_t)
        if (length(residual_buffer) > window) {
          residual_buffer <- tail(residual_buffer, window)
        }

        one_row[1, "y"] <- dat$y[idx]
        one_row[1, xcols] <- X_df[idx, , drop = FALSE]
        fit <- .forwardRiver_py_to_r(
          river_helper$learn_one(
            model = fit,
            row = .forwardRiver_row_to_list(one_row, xcols),
            y = dat$y[idx]
          )
        )
      }
    },
    timeout_sec = timeout_sec
  )

  if (exists("q", inherits = FALSE)) {
    rm(q)
  }
  if (exists("lower", inherits = FALSE)) {
    rm(lower)
  }
  if (exists("upper", inherits = FALSE)) {
    rm(upper)
  }
  if (exists("resid_t", inherits = FALSE)) {
    rm(resid_t)
  }
  if (exists("beta", inherits = FALSE)) {
    rm(beta)
  }
  if (exists("fit", inherits = FALSE)) {
    rm(fit)
  }
  if (exists("residual_buffer", inherits = FALSE)) {
    rm(residual_buffer)
  }
  gc(verbose = FALSE)
  timing
}

run_update_only_once <- function(dataset_slug, method_key, window, repetition) {
  dat <- load_runtime_dataset(dataset_slug)
  spec <- method_spec(method_key)
  hyperparams <- if (spec$family %in% c("retro", "fw_krr")) {
    get_fixed_krr_hyperparams(dataset_slug)
  } else {
    NULL
  }

  warmup_result <- switch(
    spec$family,
    retro = run_retro_update_block(
      dat,
      method_key,
      window,
      runtime_scale_warmup_updates,
      hyperparams,
      timeout_sec = runtime_scale_timeout_sec
    ),
    fw_krr = run_fw_krr_update_block(
      dat,
      window,
      runtime_scale_warmup_updates,
      hyperparams,
      timeout_sec = runtime_scale_timeout_sec
    ),
    fw_moa = run_fw_moa_update_block(
      dat,
      window,
      runtime_scale_warmup_updates,
      spec$moa_model,
      timeout_sec = runtime_scale_timeout_sec
    ),
    fw_river = run_fw_river_update_block(
      dat,
      window,
      runtime_scale_warmup_updates,
      spec$river_model,
      spec$river_control,
      timeout_sec = runtime_scale_timeout_sec
    )
  )

  if (!identical(warmup_result$status, "completed")) {
    return(make_runtime_row(
      benchmark = "update_only",
      dat = dat,
      method_key = method_key,
      window = window,
      repetition = repetition,
      timed_updates = NA_integer_,
      status = warmup_result$status,
      timeout_sec = runtime_scale_timeout_sec,
      elapsed = warmup_result$elapsed,
      updates_per_block = runtime_scale_block_updates,
      timing_repeats = 1L,
      message = paste0("warm-up failed: ", warmup_result$message)
    ))
  }

  result <- switch(
    spec$family,
    retro = run_retro_update_block(
      dat,
      method_key,
      window,
      runtime_scale_block_updates,
      hyperparams,
      timeout_sec = runtime_scale_timeout_sec
    ),
    fw_krr = run_fw_krr_update_block(
      dat,
      window,
      runtime_scale_block_updates,
      hyperparams,
      timeout_sec = runtime_scale_timeout_sec
    ),
    fw_moa = run_fw_moa_update_block(
      dat,
      window,
      runtime_scale_block_updates,
      spec$moa_model,
      timeout_sec = runtime_scale_timeout_sec
    ),
    fw_river = run_fw_river_update_block(
      dat,
      window,
      runtime_scale_block_updates,
      spec$river_model,
      spec$river_control,
      timeout_sec = runtime_scale_timeout_sec
    )
  )

  if (identical(result$status, "completed")) {
    timed_updates <- as.integer(runtime_scale_block_updates)
    elapsed <- result$elapsed
    timing_repeats <- 1L
  } else {
    timed_updates <- NA_integer_
    elapsed <- result$elapsed
    timing_repeats <- 1L
  }

  make_runtime_row(
    benchmark = "update_only",
    dat = dat,
    method_key = method_key,
    window = window,
    repetition = repetition,
    timed_updates = timed_updates,
    status = result$status,
    timeout_sec = runtime_scale_timeout_sec,
    elapsed = elapsed,
    updates_per_block = runtime_scale_block_updates,
    timing_repeats = timing_repeats,
    message = result$message
  )
}

run_end_to_end_once <- function(dataset_slug, method_key, repetition) {
  dat <- load_runtime_dataset(dataset_slug)
  total_updates <- dat$n_obs - runtime_t_init
  set.seed(1000L + repetition)

  result <- time_with_timeout(
    expr_fun = function() {
      switch(
        method_key,
        retro_eff = conformalRetroAdj(
          X = dat$X,
          y = dat$y,
          alpha = runtime_alpha,
          t_init = runtime_t_init,
          d = runtime_end_to_end_window,
          kernel = runtime_kernel,
          methods = runtime_aci_method
        ),
        retro_naive = conformalRetroAdjNaiveRefit(
          X = dat$X,
          y = dat$y,
          alpha = runtime_alpha,
          t_init = runtime_t_init,
          d = runtime_end_to_end_window,
          kernel = runtime_kernel,
          methods = runtime_aci_method
        ),
        fw_krr = forwardKRR(
          X = dat$X,
          y = dat$y,
          alpha = runtime_alpha,
          t_init = runtime_t_init,
          w = runtime_end_to_end_window,
          d = runtime_scale_data_window_krr,
          kernel = runtime_kernel,
          methods = runtime_aci_method
        ),
        fw_amf = forwardRiver(
          X = dat$X,
          y = dat$y,
          alpha = runtime_alpha,
          t_init = runtime_t_init,
          w = runtime_end_to_end_window,
          river_model = "AMFRegressor",
          methods = runtime_aci_method,
          river_control = list(seed = 1L),
          py_module_path = "python"
        ),
        fw_fimtdd = forwardMOA(
          X = dat$X,
          y = dat$y,
          alpha = runtime_alpha,
          t_init = runtime_t_init,
          w = runtime_end_to_end_window,
          moa_model = "FIMTDD",
          methods = runtime_aci_method
        ),
        fw_amrules = forwardMOA(
          X = dat$X,
          y = dat$y,
          alpha = runtime_alpha,
          t_init = runtime_t_init,
          w = runtime_end_to_end_window,
          moa_model = "AMRulesRegressor",
          methods = runtime_aci_method
        )
      )
    },
    timeout_sec = runtime_end_to_end_timeout_sec
  )

  make_runtime_row(
    benchmark = "end_to_end",
    dat = dat,
    method_key = method_key,
    window = runtime_end_to_end_window,
    repetition = repetition,
    timed_updates = total_updates,
    status = result$status,
    timeout_sec = runtime_end_to_end_timeout_sec,
    elapsed = result$elapsed,
    updates_per_block = total_updates,
    timing_repeats = 1L,
    message = result$message
  )
}

# Output helpers --------------------------------------------------------

read_runtime_csv_or_empty <- function(path) {
  if (!file.exists(path)) {
    return(empty_runtime_raw_df())
  }

  out <- suppressMessages(readr::read_csv(path, show_col_types = FALSE))
  as.data.frame(out, stringsAsFactors = FALSE)
}

sd_or_na <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) <= 1L) {
    return(NA_real_)
  }
  stats::sd(x)
}

mean_or_na <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) {
    return(NA_real_)
  }
  mean(x)
}

format_mean_sd_ascii <- function(mean_val, sd_val, digits = 2L) {
  if (is.na(mean_val)) {
    return("timeout")
  }
  if (is.na(sd_val)) {
    return(sprintf(paste0("%.", digits, "f +/- NA"), mean_val))
  }
  sprintf(
    paste0("%.", digits, "f +/- %.", digits, "f"),
    mean_val,
    sd_val
  )
}

format_mean_sd_tex <- function(mean_val, sd_val, digits = 2L) {
  if (is.na(mean_val)) {
    return("timeout")
  }
  if (is.na(sd_val)) {
    return(sprintf(paste0("%.", digits, "f $\\\\pm$ NA"), mean_val))
  }
  sprintf(
    paste0("%.", digits, "f $\\\\pm$ %.", digits, "f"),
    mean_val,
    sd_val
  )
}

summarise_runtime_raw <- function(raw_df) {
  if (nrow(raw_df) == 0L) {
    return(data.frame())
  }

  raw_df %>%
    group_by(
      benchmark,
      dataset,
      dataset_slug,
      method,
      method_key,
      window
    ) %>%
    summarise(
      reps_requested = n(),
      completed_reps = sum(status == "completed"),
      timeout_reps = sum(status == "timeout"),
      error_reps = sum(status == "error"),
      mean_total_runtime_sec = mean_or_na(total_runtime_sec[status == "completed"]),
      sd_total_runtime_sec = sd_or_na(total_runtime_sec[status == "completed"]),
      mean_runtime_per_step_ms = mean_or_na(runtime_per_step_ms[status == "completed"]),
      sd_runtime_per_step_ms = sd_or_na(runtime_per_step_ms[status == "completed"]),
      ms_per_update_display = format_mean_sd_ascii(
        mean_or_na(runtime_per_step_ms[status == "completed"]),
        sd_or_na(runtime_per_step_ms[status == "completed"])
      ),
      ms_per_update_display_tex = format_mean_sd_tex(
        mean_or_na(runtime_per_step_ms[status == "completed"]),
        sd_or_na(runtime_per_step_ms[status == "completed"])
      ),
      .groups = "drop"
    ) %>%
    mutate(
      benchmark = factor(
        benchmark,
        levels = c("end_to_end", "update_only"),
        labels = c("End-to-end", "Update-only")
      ),
      dataset = factor(dataset, levels = c("Elec2", "AIG")),
      method = factor(method, levels = runtime_method_order)
    ) %>%
    arrange(benchmark, dataset, method, window)
}

build_runtime_main_table <- function(end_summary_df) {
  if (nrow(end_summary_df) == 0L) {
    return(data.frame())
  }

  elec2_df <- end_summary_df %>%
    filter(dataset == "Elec2") %>%
    transmute(
      method = as.character(method),
      elec2_ms_per_update = ms_per_update_display
    )

  aig_df <- end_summary_df %>%
    filter(dataset == "AIG") %>%
    transmute(
      method = as.character(method),
      aig_ms_per_update = ms_per_update_display
    )

  out <- data.frame(
    method = runtime_method_order,
    stringsAsFactors = FALSE
  ) %>%
    left_join(elec2_df, by = "method") %>%
    left_join(aig_df, by = "method") %>%
    mutate(
      benchmark = sprintf(
        "End-to-end full stream, w = %d, mean +/- sd ms/update",
        runtime_end_to_end_window
      )
    ) %>%
    select(benchmark, method, elec2_ms_per_update, aig_ms_per_update)

  out
}

build_runtime_timeout_summary <- function(update_summary_df) {
  if (nrow(update_summary_df) == 0L) {
    return(data.frame())
  }

  update_summary_df %>%
    filter(timeout_reps > 0L | completed_reps < reps_requested) %>%
    transmute(
      dataset = as.character(dataset),
      method = as.character(method),
      window = window,
      reps_requested = reps_requested,
      completed_reps = completed_reps,
      timeout_reps = timeout_reps,
      error_reps = error_reps,
      reported_ms_per_update = ms_per_update_display
    ) %>%
    arrange(dataset, match(method, runtime_method_order), window)
}

build_scalability_plot <- function(update_summary_df) {
  plot_df <- update_summary_df %>%
    filter(completed_reps > 0L, !is.na(mean_runtime_per_step_ms)) %>%
    mutate(
      dataset = factor(as.character(dataset), levels = c("Elec2", "AIG")),
      method = factor(as.character(method), levels = runtime_method_order)
    )

  ggplot(
    plot_df,
    aes(
      x = window,
      y = mean_runtime_per_step_ms,
      color = method,
      group = method
    )
  ) +
    geom_line(linewidth = 1.1, lineend = "round") +
    geom_point(size = 2.5) +
    geom_errorbar(
      aes(
        ymin = pmax(mean_runtime_per_step_ms - sd_runtime_per_step_ms, 1e-6),
        ymax = mean_runtime_per_step_ms + sd_runtime_per_step_ms
      ),
      width = 18,
      linewidth = 0.45
    ) +
    scale_color_manual(values = runtime_method_palette) +
    scale_x_continuous(breaks = runtime_scale_windows) +
    scale_y_log10(labels = label_number(accuracy = 0.1, big.mark = ",")) +
    facet_wrap(~dataset, nrow = 1) +
    labs(
      x = "Window size (w)",
      y = "Runtime per update (ms)"
    ) +
    runtime_base_theme +
    guides(
      color = guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(linewidth = 2.8)
      )
    ) +
    theme(panel.spacing.x = grid::unit(14, "pt"))
}

save_scalability_plot <- function(plot_obj) {
  ggsave(
    scalability_plot_pdf_path,
    plot_obj,
    width = 12,
    height = 4.4,
    units = "in",
    bg = "white"
  )
  ggsave(
    scalability_plot_png_path,
    plot_obj,
    width = 12,
    height = 4.4,
    units = "in",
    dpi = 400,
    bg = "white"
  )
  invisible(list(pdf = scalability_plot_pdf_path, png = scalability_plot_png_path))
}

write_runtime_main_table_tex <- function(main_table_df, path = main_table_tex_path) {
  if (nrow(main_table_df) == 0L) {
    return(invisible(NULL))
  }

  lines <- c(
    "\\begin{table}[t]",
    "\\centering",
    sprintf(
      "\\caption{End-to-end wall-clock runtime in mean $\\\\pm$ SD milliseconds per update for the full stream benchmark at $w=%d$.}",
      runtime_end_to_end_window
    ),
    "\\label{tab:runtime-main}",
    "\\begin{tabular}{lll}",
    "\\toprule",
    "Method & Elec2 & AIG \\\\",
    "\\midrule"
  )

  row_lines <- apply(main_table_df, 1, function(row) {
    paste(
      row[c("method", "elec2_ms_per_update", "aig_ms_per_update")],
      collapse = " & "
    )
  })
  row_lines <- paste0(row_lines, " \\\\")

  lines <- c(
    lines,
    row_lines,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}"
  )
  writeLines(lines, con = path, useBytes = TRUE)
  invisible(path)
}

write_runtime_scalability_table_tex <- function(update_summary_df, path = scale_table_tex_path) {
  if (nrow(update_summary_df) == 0L) {
    return(invisible(NULL))
  }

  latex_df <- update_summary_df %>%
    transmute(
      Dataset = as.character(dataset),
      Method = as.character(method),
      Window = window,
      `Mean $\\pm$ SD (ms/update)` = ms_per_update_display_tex,
      `Completed/Requested` = sprintf("%d/%d", completed_reps, reps_requested),
      `Timeout Reps` = timeout_reps
    )

  lines <- c(
    "\\begin{table}[t]",
    "\\centering",
    paste0(
      "\\caption{Update-only scalability benchmark. Hyperparameter search and one-time setup are excluded. Each repetition performs one untimed warm-up update followed by a timed block of ",
      runtime_scale_block_updates,
      " updates; values report mean $\\\\pm$ SD milliseconds per update over completed repetitions, and timeout counts are shown explicitly.}"
    ),
    "\\label{tab:runtime-scalability}",
    "\\begin{tabular}{lllrrr}",
    "\\toprule",
    "Dataset & Method & Window & Mean $\\pm$ SD (ms/update) & Completed/Requested & Timeout Reps \\\\",
    "\\midrule"
  )

  row_lines <- apply(latex_df, 1, function(row) {
    paste(row, collapse = " & ")
  })
  row_lines <- paste0(row_lines, " \\\\")

  lines <- c(
    lines,
    row_lines,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}"
  )
  writeLines(lines, con = path, useBytes = TRUE)
  invisible(path)
}

write_runtime_protocol_notes <- function(
    path_txt = protocol_note_path,
    path_tex = protocol_tex_path
) {
  txt_lines <- c(
    "Runtime benchmarks use wall-clock elapsed time measured from Sys.time() differences and are reported in milliseconds per update.",
    sprintf(
      "End-to-end benchmark: full stream execution at w = %d with %d repetitions. The timed section includes each method's native initialization, online updates, adaptive alpha updates, interval construction, and final evaluation.",
      runtime_end_to_end_window,
      runtime_end_to_end_reps
    ),
    paste0(
      "Update-only scalability benchmark: timed from a pre-initialized state on a fixed contiguous block of updates. Hyperparameter search and one-time setup are excluded. "
    ),
    paste0(
      "For RetroAdj and RetroAdj-naive, the timed section includes Jackknife+ residual generation, beta computation, interval quantiles, sliding-window downdate, and model update; the two methods differ only in the residual-generation routine."
    ),
    paste0(
      "For FW-KRR, the timed section includes prediction, conformal quantile computation from the residual buffer, beta computation, residual-buffer maintenance, and KRR update; initial KRR fitting and residual warm-up are excluded."
    ),
    paste0(
      "For FW-AMF, the timed section includes River prediction, conformal quantile computation from the residual buffer, beta computation, residual-buffer maintenance, and the single River online update; initial River model fitting and residual warm-up are excluded."
    ),
    paste0(
      "For FW-FIMT-DD and FW-AMRules, the timed section includes prediction, conformal quantile computation from the residual buffer, beta computation, residual-buffer maintenance, and the single MOA online update; initial MOA training and residual warm-up are excluded."
    ),
    sprintf(
      "Each update-only repetition performs one untimed warm-up update followed by a fixed %d-update timed block, and results are averaged over %d repetitions. Timeout results are reported explicitly when the timed update loop exceeds %.0f seconds.",
      runtime_scale_block_updates,
      runtime_scale_reps,
      runtime_scale_timeout_sec
    ),
    "All workers use single-thread BLAS settings (OMP, OpenBLAS, MKL threads set to 1), and the R JIT level is fixed at 3."
  )

  tex_lines <- c(
    "\\paragraph{Runtime protocol.}",
    sprintf(
      "We measured wall-clock elapsed time from \\texttt{Sys.time()} differences and report milliseconds per update. In the end-to-end benchmark, each method was run on the full stream at $w=%d$ for %d repetitions, and the timed section includes native initialization, online updates, adaptive $\\alpha_t$ updates, interval construction, and final evaluation.",
      runtime_end_to_end_window,
      runtime_end_to_end_reps
    )
  )

  tex_lines <- c(
    tex_lines,
    "In the update-only scalability benchmark, hyperparameter search and one-time setup were excluded. Each method was initialized once on a fixed contiguous pre-benchmark segment, after which only the benchmarked online updates were timed.",
    "For RetroAdj and RetroAdj-naive, the timed section includes Jackknife+ residual generation, $\\beta_t$ computation, interval quantiles, sliding-window downdate, and model update; the two methods differ only in how Jackknife+ residuals are computed.",
    "For FW-KRR, the timed section includes prediction, residual-window quantile computation, $\\beta_t$ computation, residual-buffer maintenance, and KRR update, excluding initial KRR fitting and residual warm-up.",
    "For FW-AMF, the timed section includes River prediction, residual-window quantile computation, $\\beta_t$ computation, residual-buffer maintenance, and the single River online update, excluding initial River model fitting and residual warm-up.",
    "For FW-FIMT-DD and FW-AMRules, the timed section includes prediction, residual-window quantile computation, $\\beta_t$ computation, residual-buffer maintenance, and the single MOA online update, excluding initial MOA training and residual warm-up.",
    sprintf(
      "Each update-only repetition performs one untimed warm-up update followed by a fixed %d-update timed block, and results are averaged over %d repetitions. Repetitions whose timed update loop exceeds %.0f seconds are recorded as timeouts.",
      runtime_scale_block_updates,
      runtime_scale_reps,
      runtime_scale_timeout_sec
    ),
    "All workers used single-thread BLAS settings (OMP, OpenBLAS, and MKL threads set to 1), and the R JIT level was fixed at 3."
  )

  writeLines(txt_lines, con = path_txt, useBytes = TRUE)
  writeLines(tex_lines, con = path_tex, useBytes = TRUE)
  invisible(list(txt = path_txt, tex = path_tex))
}

save_runtime_outputs <- function(
    end_to_end_raw,
    update_only_raw
) {
  end_summary <- summarise_runtime_raw(end_to_end_raw)
  update_summary <- summarise_runtime_raw(update_only_raw)
  main_table <- build_runtime_main_table(end_summary)
  timeout_summary <- build_runtime_timeout_summary(update_summary)

  readr::write_csv(end_to_end_raw, end_to_end_raw_path)
  readr::write_csv(update_only_raw, update_only_raw_path)
  readr::write_csv(end_summary, end_to_end_summary_path)
  readr::write_csv(update_summary, update_only_summary_path)
  readr::write_csv(main_table, main_table_path)
  readr::write_csv(timeout_summary, timeout_summary_path)
  write_runtime_main_table_tex(main_table, main_table_tex_path)
  write_runtime_scalability_table_tex(update_summary, scale_table_tex_path)
  write_runtime_protocol_notes()
  writeLines(capture.output(sessionInfo()), con = session_info_path, useBytes = TRUE)

  if (nrow(update_summary) > 0L) {
    plot_obj <- build_scalability_plot(update_summary)
    save_scalability_plot(plot_obj)
  }

  invisible(list(
    end_summary = end_summary,
    update_summary = update_summary,
    main_table = main_table,
    timeout_summary = timeout_summary
  ))
}
