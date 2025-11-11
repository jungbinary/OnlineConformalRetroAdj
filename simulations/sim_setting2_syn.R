g_bump <- function(X, center) {
  r <- sqrt(rowSums((sweep(X, 2, center))^2))
  s <- pmax(1 - r, 0)
  (s^6) * (35 * r^2 + 18 * r + 3)
}

gen_wendland <- function(n = 1000, t_init = 250,
                         a1 = 1.0, a2 = -1.0,
                         b1 = 0.0, b2 = 0.4,
                         c1 = c(0.25, 0.25, 0.25),
                         c2 = c(0.3, 0.4, 0.5),
                         sigma = 0.5, seed = NULL) {
  stopifnot(t_init < n)
  if (!is.null(seed)) set.seed(seed)
  
  X <- matrix(runif(n * 3), ncol = 3)
  
  idx1 <- 1:t_init
  idx2 <- (t_init + 1):n
  
  f <- numeric(n)
  f[idx1] <- a1 * g_bump(X[idx1, , drop = FALSE], c1) + b1
  f[idx2] <- a2 * g_bump(X[idx2, , drop = FALSE], c2) + b2
  
  y <- f + rnorm(n, sd = sigma)
  
  list(X = X, y = y,
       f = f,
       params = list(t_init = t_init, a1 = a1, a2 = a2,
                     b1 = b1, b2 = b2, c1 = c1, c2 = c2, sigma = sigma))
}


one_experiment_sim2 <- function(alpha, t_init = 250, d, kernel, methods, seed, n = 1000){

  run_safe <- function(expr){
    t0 <- Sys.time()
    out <- tryCatch(eval.parent(substitute(expr)), error = function(e) e)
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    list(out = out, time = elapsed)
  }
  
  dat <- gen_wendland(n = n, t_init = t_init, seed = seed)
  X <- dat$X ; y <- dat$y
  
  r_retro <- run_safe({
    conformalRetroAdj(X, y, alpha = alpha, t_init = t_init, d = d,
                      kernel = kernel, methods = methods)
  })
  r_k  <- run_safe({
    forwardKRR(X, y, alpha = alpha, t_init = t_init, d = d,
                 kernel = kernel, methods = methods)
  })
  r_am <- run_safe({
    forwardMOA(X, y, alpha = alpha, t_init = t_init, w = d,
                 moa_model = "AMRulesRegressor", methods = methods)
  })
  r_tr <- run_safe({
    forwardMOA(X, y, alpha = alpha, t_init = t_init, w = d,
                 moa_model = "FIMTDD", methods = methods)
  })
  
  get_metric <- function(res, name){
    if (inherits(res$out, "error")) NA_real_ else as.numeric(res$out[[name]])
  }
  
  # coverage / length / time vectors -------------------------------------------
  cov  <- c(
    RetroAdj = get_metric(r_retro, "coverage"),
    KRR      = get_metric(r_k ,  "coverage"),
    AMRules  = get_metric(r_am,  "coverage"),
    FIMTDD   = get_metric(r_tr,  "coverage")
  )
  len  <- c(
    RetroAdj = get_metric(r_retro, "mean_len"),
    KRR      = get_metric(r_k ,  "mean_len"),
    AMRules  = get_metric(r_am,  "mean_len"),
    FIMTDD   = get_metric(r_tr,  "mean_len")
  )
  time <- c(
    RetroAdj = r_retro$time,
    KRR      = r_k$time,
    AMRules  = r_am$time,
    FIMTDD   = r_tr$time
  )
  
  list(
    cov = cov, len = len, time = time
  )
}

library(future.apply)
run_sim_parallel_sim2 <- function(alpha = 0.1, t_init = 250, d = 250,
                                  R = 50, seed = 2025, kernel, methods,
                                  n = 1000,
                                  workers = max(1, parallel::detectCores() - 1)) {
  set.seed(seed)
  plan(multisession, workers = workers)
  on.exit(plan(sequential), add = TRUE)
  
  out <- future_lapply(
    seq_len(R),
    function(i) {
      one_experiment_sim2(alpha = alpha, t_init = t_init, d = d,
                     kernel = kernel, methods = methods,
                     seed = seed + i, n = n)
    },
    future.seed = TRUE
  )
  
  cov_mat  <- do.call(rbind, lapply(out, `[[`, "cov"))
  len_mat  <- do.call(rbind, lapply(out, `[[`, "len"))
  time_mat <- do.call(rbind, lapply(out, `[[`, "time"))
  
  colnames(cov_mat)  <- names(out[[1]]$cov)
  colnames(len_mat)  <- names(out[[1]]$len)
  colnames(time_mat) <- names(out[[1]]$time)
  rownames(cov_mat)  <- rownames(len_mat) <- rownames(time_mat) <-
    paste0("iter_", seq_len(R))
  
  list(coverage = cov_mat, len = len_mat, time = time_mat)
}

# run per-ACI variant ----------------------------------------------------------
result_aci_2   <- run_sim_parallel_sim2(alpha=0.1, t_init=250, d=250, R = 50, seed = 2025,
                                 kernel="rbf", methods= "aci")
result_sfogd_2 <- run_sim_parallel_sim2(alpha=0.1, t_init=250, d=250, R = 50, seed = 2025,
                                 kernel="rbf", methods= "sfogd")
result_saocp_2 <- run_sim_parallel_sim2(alpha=0.1, t_init=250, d=250, R = 50, seed = 2025,
                                 kernel="rbf", methods= "saocp")
result_dtaci_2 <- run_sim_parallel_sim2(alpha=0.1, t_init=250, d=250, R = 50, seed = 2025,
                                 kernel="rbf", methods= "dtaci")
result_agaci_2 <- run_sim_parallel_sim2(alpha=0.1, t_init=250, d=250, R = 50, seed = 2025,
                                 kernel="rbf", methods= "agaci")

library(ggplot2)

# 1) Collect results (sim2) ----------------------------------------------------
res_sim2 <- list(
  ACI   = result_aci_2,
  AgACI = result_agaci_2,
  DtACI = result_dtaci_2,
  SFOGD = result_sfogd_2,
  SAOCP = result_saocp_2
)

res <- res_sim2

# 2) Summary helper: mean + 95% CI --------------------------------------------
summarise_mat <- function(mat, method, metric) {
  m   <- colMeans(mat, na.rm = TRUE)
  sdv <- apply(mat, 2, sd, na.rm = TRUE)
  n   <- nrow(mat)
  sem <- sdv / sqrt(n)
  data.frame(method = method, model = names(m), metric = metric,
             mean = unname(m), lo = unname(m - 1.96*sem), hi = unname(m + 1.96*sem))
}

# 3) Build tidy data frame -----------------------------------------------------
df <- do.call(
  rbind,
  unlist(lapply(names(res), function(meth) {
    r <- res[[meth]]
    list(
      summarise_mat(r$coverage, meth, "Coverage"),
      summarise_mat(r$len,      meth, "Width"),
      summarise_mat(log10(r$time), meth, "log(Time)")
    )
  }), recursive = FALSE)
)

df$model  <- factor(
  df$model,
  levels = c("FIMTDD","AMRules","KRR","RetroAdj"),
  labels = c("FIMTDD","AMRules","KRR","RetroAdj")
)
df$method <- factor(df$method, levels = c("ACI","AgACI","DtACI","SFOGD","SAOCP"))

# 4) Palette & plotting helpers -----------------------------------------------
cols <- c(
  "FIMTDD"  = "#90C8AC",
  "AMRules" = "#E3C57B",
  "KRR"     = "#B8A9C9",
  "RetroAdj"= "#EE6A70"
)

base_theme <- theme_bw(base_size = 15) +
  theme(plot.title=element_text(face="bold", size=14),
        plot.subtitle=element_text(size=15),
        axis.title.x  = element_text(size = 15),
        axis.title.y  = element_text(size = 15),
        legend.position="bottom")

plot_metric <- function(df, metric, alpha_target = 0.1) {
  g <- ggplot(df[df$metric == metric, ],
              aes(x = method, y = mean, fill = model)) +
    geom_col(position = position_dodge(width = .7), width = .6) +
    geom_errorbar(aes(ymin = lo, ymax = hi),
                  position = position_dodge(width = .7), width = .2) +
    scale_fill_manual(
      name   = "Methods",
      values = cols,
      labels = c("FW-FIMTDD", "FW-AMRules", "FW-KRR", "RetroAdj")
    ) +
    labs(
      x = "ACI Methods",
      y = metric
    ) +
    base_theme +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE))  # single horizontal row
  
  if (tolower(metric) == "coverage") {
    g <- g + geom_hline(yintercept = 1 - alpha_target, linetype = "dashed")
  }
  g
}

# 5) Make plots ----------------------------------------------------------------
p_cov  <- plot_metric(df, "Coverage", alpha_target = 0.1) +
  coord_cartesian(ylim = c(0.8, 1))
p_len  <- plot_metric(df, "Width") +
  coord_cartesian(ylim = c(1.5, 3))

library(patchwork)

# 12x4 -------------------------------------------------------------------------
p_2 <- (p_cov | p_len)
p_2
