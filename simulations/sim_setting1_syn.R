# Data Generation ======================================================
gen_beta_shift_only <- function(n = 500,
                                p = 10,
                                t_init = 500,
                                sigma = 0.5,
                                beta0 = NULL,
                                beta1 = NULL,
                                seed = 2025) {
  stopifnot(t_init < n, p >= 1)
  set.seed(seed)
  
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("x", seq_len(p))
  
  if (is.null(beta0)) {
    beta0 <- c(1.0, 0.8, 0.0, 0.0, 0.5, 0.0, 0.3, 0.0, 0.0, 0.2)[seq_len(p)]
  }
  if (is.null(beta1)) {
    beta1 <- c(0.0, -1.2, 0.7, 0.4, 0.0, 0.0, 0.9, 0.0, -0.6, 0.0)[seq_len(p)]
  }
  
  f_true <- numeric(n)
  f_true[1:t_init]       <- X[1:t_init, ] %*% beta0
  f_true[(t_init+1):n]   <- X[(t_init+1):n, ] %*% beta1
  
  eps <- rnorm(n, sd = sigma)
  y <- as.numeric(f_true + eps)
  
  list(
    X = X, y = y,
    beta0 = as.numeric(beta0),
    beta1 = as.numeric(beta1),
    f_true = as.numeric(f_true),
    t_init = t_init,
    info = list(n = n, p = p, sigma = sigma, seed = seed)
  )
}


one_experiment <- function(alpha, t_init = 250, d, kernel, methods, seed, n = 1000){

  run_safe <- function(expr){
    t0 <- Sys.time()
    out <- tryCatch(eval.parent(substitute(expr)), error = function(e) e)
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    list(out = out, time = elapsed)
  }
  
  dat <- gen_beta_shift_only(n = n, t_init = t_init, seed = seed)
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

run_sim_parallel <- function(alpha = 0.1, t_init = 250, d = 250,
                             R = 10, seed = 2025, kernel, methods,
                             n = 1000,
                             workers = max(1, parallel::detectCores() - 1)) {
  set.seed(seed)
  plan(multisession, workers = workers)
  on.exit(plan(sequential), add = TRUE)
  
  out <- future_lapply(
    seq_len(R),
    function(i) {
      one_experiment(alpha = alpha, t_init = t_init, d = d,
                     kernel = kernel, methods = methods,
                     seed = seed + i, n = n)
    },
    future.seed = TRUE
  )
  
  # bind results into matrices -------------------------------------------------
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
result_aci_1   <- run_sim_parallel(alpha=0.1, t_init=250, d=250, R = 50, seed = 2025,
                                 kernel="rbf", methods= "aci")
result_sfogd_1 <- run_sim_parallel(alpha=0.1, t_init=250, d=250, R = 50, seed = 2025,
                                 kernel="rbf", methods= "sfogd")
result_saocp_1 <- run_sim_parallel(alpha=0.1, t_init=250, d=250, R = 50, seed = 2025,
                                 kernel="rbf", methods= "saocp")
result_dtaci_1 <- run_sim_parallel(alpha=0.1, t_init=250, d=250, R = 50, seed = 2025,
                                 kernel="rbf", methods= "dtaci")
result_agaci_1 <- run_sim_parallel(alpha=0.1, t_init=250, d=250, R = 50, seed = 2025,
                                 kernel="rbf", methods= "agaci")

library(ggplot2)

# 1) Collect results -----------------------------------------------------------
res_sim1 <- list(
  ACI   = result_aci_1,
  AgACI = result_agaci_1,
  DtACI = result_dtaci_1,
  SFOGD = result_sfogd_1,
  SAOCP = result_saocp_1
)

res <- res_sim1

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
  coord_cartesian(ylim = c(1, 7))

p_1 <- ( p_cov | p_len )
p_1

# 12x4 -------------------------------------------------------------------------
row_title <- function(txt, size = 16, bottom_margin = 4){
  cowplot::ggdraw() +
    cowplot::draw_label(
      txt, x = 0, y = 1, hjust = 0, vjust = 1,
      fontface = "bold", size = size
    ) +
    theme(plot.margin = margin(t = 0, r = 0, b = bottom_margin, l = 0))
}


# Run After sim2syn.R
title1 <- row_title("Setting 1 (Linear)",     size = 16, bottom_margin = 6)
title2 <- row_title("Setting 2 (Non-Linear)", size = 16, bottom_margin = 6)


final_plot <- (title1 / p_1 / title2 / p_2) +
  plot_layout(
    ncol    = 1,
    heights = c(0.09, 1, 0.09, 1),  
    guides  = "collect"
  ) &
  theme(
    legend.position  = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 15),
    plot.margin      = margin(t = 8, r = 8, b = 8, l = 8) 
  ) &
  guides(
    color = guide_legend(nrow = 1, byrow = TRUE,
                         override.aes = list(size = 2.2))
  )

final_plot
