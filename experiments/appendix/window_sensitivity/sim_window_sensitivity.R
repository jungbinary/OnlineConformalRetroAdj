suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(parallel)
  library(patchwork)
})

source("R/baseKRR.R", encoding = "UTF-8")
source("aci/dtaci.R", encoding = "UTF-8")
source("aci/agaci.R", encoding = "UTF-8")
source("aci/sfogd.R", encoding = "UTF-8")
source("aci/saocp.R", encoding = "UTF-8")
source("R/conformalRetroAdj.R", encoding = "UTF-8")
source("R/forwardKRR.R", encoding = "UTF-8")

output_dir <- file.path("experiments", "results", "appendix", "window_sensitivity")
figure_stem <- "window_sensitivity"
rep_csv_path <- file.path(output_dir, "window_sensitivity_rep_results.csv")
summary_csv_path <- file.path(output_dir, "window_sensitivity_summary.csv")
summary_tex_path <- file.path(output_dir, "window_sensitivity_table.tex")
protocol_txt_path <- file.path(output_dir, "window_sensitivity_protocol.txt")

window_grid <- c(100L, 250L, 500L, 750L, 1000L)
t_init_n <- 250L
alpha0 <- 0.1
kernel_name <- "rbf"
aci_method <- "dtaci"
n_total <- 1000L
n_rep <- 50L
seed_base <- 2025L
n_workers <- max(1L, min(4L, parallel::detectCores(logical = FALSE)))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

method_levels <- c("RetroAdj", "FW-KRR")
method_palette <- c(
  "RetroAdj" = "#EE6A70",
  "FW-KRR" = "#B8A9C9"
)

base_theme <- theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    legend.position = "bottom"
  )

configure_single_thread_math <- function() {
  Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1"
  )
  invisible(NULL)
}

gen_beta_shift_only <- function(
    n = 500L,
    p = 10L,
    t_init = 500L,
    sigma = 0.5,
    beta0 = NULL,
    beta1 = NULL,
    seed = 2025L
) {
  stopifnot(t_init < n, p >= 1L)
  set.seed(seed)

  x <- matrix(rnorm(n * p), n, p)
  colnames(x) <- paste0("x", seq_len(p))

  if (is.null(beta0)) {
    beta0 <- c(1.0, 0.8, 0.0, 0.0, 0.5, 0.0, 0.3, 0.0, 0.0, 0.2)[seq_len(p)]
  }
  if (is.null(beta1)) {
    beta1 <- c(0.0, -1.2, 0.7, 0.4, 0.0, 0.0, 0.9, 0.0, -0.6, 0.0)[seq_len(p)]
  }

  f_true <- numeric(n)
  f_true[1:t_init] <- x[1:t_init, , drop = FALSE] %*% beta0
  f_true[(t_init + 1L):n] <- x[(t_init + 1L):n, , drop = FALSE] %*% beta1

  y <- as.numeric(f_true + rnorm(n, sd = sigma))
  list(X = x, y = y)
}

mean_or_zero_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1L) {
    return(0)
  }
  stats::sd(x)
}

run_one_rep <- function(rep_id) {
  seed_value <- seed_base + rep_id
  dat <- gen_beta_shift_only(
    n = n_total,
    t_init = t_init_n,
    seed = seed_value
  )
  x <- dat$X
  y <- dat$y

  base_model <- krr_init(
    x[1:t_init_n, , drop = FALSE],
    y[1:t_init_n],
    kernel = kernel_name
  )

  bind_rows(lapply(window_grid, function(window_size) {
    retro_rst <- conformalRetroAdj(
      X = x,
      y = y,
      alpha = alpha0,
      t_init = t_init_n,
      lambda = base_model$lambda,
      sigma = base_model$sigma,
      kernel = kernel_name,
      d = window_size,
      methods = aci_method
    )

    fw_rst <- forwardKRR(
      X = x,
      y = y,
      alpha = alpha0,
      t_init = t_init_n,
      lambda = base_model$lambda,
      sigma = base_model$sigma,
      kernel = kernel_name,
      w = window_size,
      methods = aci_method
    )

    bind_rows(
      data.frame(
        rep = rep_id,
        seed = seed_value,
        window = window_size,
        method = "RetroAdj",
        coverage = retro_rst$coverage,
        rep_mean_width = retro_rst$mean_len,
        stringsAsFactors = FALSE
      ),
      data.frame(
        rep = rep_id,
        seed = seed_value,
        window = window_size,
        method = "FW-KRR",
        coverage = fw_rst$coverage,
        rep_mean_width = fw_rst$mean_len,
        stringsAsFactors = FALSE
      )
    )
  }))
}

summarise_reps <- function(rep_df) {
  rep_df %>%
    mutate(method = factor(method, levels = method_levels)) %>%
    group_by(method, window) %>%
    summarise(
      n_rep = n(),
      mean_coverage = mean(coverage),
      sd_coverage = mean_or_zero_sd(coverage),
      mean_width = mean(rep_mean_width),
      sd_width = mean_or_zero_sd(rep_mean_width),
      coverage_gap = abs(mean_coverage - (1 - alpha0)),
      coverage_pm = sprintf("%.4f +/- %.4f", mean_coverage, sd_coverage),
      width_pm = sprintf("%.4f +/- %.4f", mean_width, sd_width),
      .groups = "drop"
    ) %>%
    arrange(method, window)
}

plot_metric <- function(summary_df, metric = c("coverage", "width")) {
  metric <- match.arg(metric)

  if (metric == "coverage") {
    y_col <- "mean_coverage"
    sd_col <- "sd_coverage"
    y_lab <- "Coverage"
    p <- ggplot(summary_df, aes(x = window, y = .data[[y_col]], color = method)) +
      geom_hline(
        yintercept = 1 - alpha0,
        color = "black",
        linewidth = 0.6
      )
  } else {
    y_col <- "mean_width"
    sd_col <- "sd_width"
    y_lab <- "Width"
    p <- ggplot(summary_df, aes(x = window, y = .data[[y_col]], color = method))
  }

  plot_df <- summary_df %>%
    mutate(
      ymin = .data[[y_col]] - .data[[sd_col]],
      ymax = .data[[y_col]] + .data[[sd_col]]
    )

  p +
    geom_vline(
      xintercept = 250,
      color = "gray40",
      linetype = "dashed",
      linewidth = 0.6
    ) +
    geom_ribbon(
      data = plot_df,
      aes(
        x = window,
        ymin = ymin,
        ymax = ymax,
        fill = method,
        group = method
      ),
      alpha = 0.18,
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    geom_line(
      data = plot_df,
      aes(
        x = window,
        y = .data[[y_col]],
        color = method,
        group = method
      ),
      linewidth = 1.2,
      lineend = "round",
      inherit.aes = FALSE
    ) +
    geom_point(
      data = plot_df,
      aes(x = window, y = .data[[y_col]], color = method),
      size = 2.8,
      inherit.aes = FALSE
    ) +
    scale_color_manual(values = method_palette) +
    scale_fill_manual(values = method_palette) +
    scale_x_continuous(breaks = window_grid) +
    labs(x = "Window Size (w)", color = "Method") +
    ylab(y_lab) +
    base_theme +
    theme(legend.position = "bottom")
}

save_summary_tex <- function(summary_df, output_path) {
  latex_df <- summary_df %>%
    transmute(
      Method = as.character(method),
      Window = window,
      `Mean Coverage $\\pm$ SD` = sprintf("%.4f $\\\\pm$ %.4f", mean_coverage, sd_coverage),
      `Mean Width $\\pm$ SD` = sprintf("%.4f $\\\\pm$ %.4f", mean_width, sd_width)
    )

  lines <- c(
    "\\begin{table}[t]",
    "\\centering",
    paste0(
      "\\caption{Synthetic Setting 1 window-size sensitivity over 50 replications. ",
      "Each row reports the mean $\\\\pm$ SD across replications.}"
    ),
    "\\label{tab:syn-setting1-window-tradeoff}",
    "\\begin{tabular}{llll}",
    "\\toprule",
    "Method & Window & Mean Coverage $\\pm$ SD & Mean Width $\\pm$ SD \\\\",
    "\\midrule"
  )

  row_lines <- apply(latex_df, 1, function(row) paste(row, collapse = " & "))
  row_lines <- paste0(row_lines, " \\\\")

  lines <- c(
    lines,
    row_lines,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}"
  )

  writeLines(lines, con = output_path, useBytes = TRUE)
  invisible(output_path)
}

configure_single_thread_math()
rep_ids <- seq_len(n_rep)

if (n_workers <= 1L) {
  rep_results <- bind_rows(lapply(rep_ids, run_one_rep))
} else {
  cl <- parallel::makePSOCKcluster(n_workers)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterCall(cl, function() {
    suppressPackageStartupMessages(library(dplyr))
    source("R/baseKRR.R", encoding = "UTF-8")
    source("aci/dtaci.R", encoding = "UTF-8")
    source("aci/agaci.R", encoding = "UTF-8")
    source("aci/sfogd.R", encoding = "UTF-8")
    source("aci/saocp.R", encoding = "UTF-8")
    source("R/conformalRetroAdj.R", encoding = "UTF-8")
    source("R/forwardKRR.R", encoding = "UTF-8")
    Sys.setenv(
      OMP_NUM_THREADS = "1",
      OPENBLAS_NUM_THREADS = "1",
      MKL_NUM_THREADS = "1"
    )
    invisible(NULL)
  })
  parallel::clusterExport(
    cl,
    varlist = c(
      "aci_method",
      "alpha0",
      "gen_beta_shift_only",
      "kernel_name",
      "n_total",
      "run_one_rep",
      "seed_base",
      "t_init_n",
      "window_grid"
    ),
    envir = environment()
  )
  rep_results <- bind_rows(parallel::parLapplyLB(cl, rep_ids, run_one_rep))
}

summary_df <- summarise_reps(rep_results)

write.csv(rep_results, rep_csv_path, row.names = FALSE)
write.csv(summary_df, summary_csv_path, row.names = FALSE)
save_summary_tex(summary_df, summary_tex_path)

protocol_lines <- c(
  "Synthetic Setting 1 window-size sensitivity protocol",
  sprintf("Target miscoverage level alpha = %.2f.", alpha0),
  sprintf("ACI update method = %s.", aci_method),
  sprintf("Window grid = {%s}.", paste(window_grid, collapse = ", ")),
  sprintf("Number of replications = %d.", n_rep),
  sprintf("Replication seeds = %d, %d, ..., %d.", seed_base + 1L, seed_base + 2L, seed_base + n_rep),
  sprintf("Data generation follows Setting 1 with n = %d and t_init = %d.", n_total, t_init_n),
  "For each replication, the KRR hyperparameters lambda and sigma are tuned once on the initial 250 observations and then reused across both methods and all window sizes.",
  "For RetroAdj, the active-data window is set to d = w. For FW-KRR, the residual calibration window is set to w while the KRR update window follows the default implementation used in the paper.",
  "Reported standard deviations are empirical SDs across the 50 synthetic replications."
)
writeLines(protocol_lines, con = protocol_txt_path, useBytes = TRUE)

p_cov <- plot_metric(summary_df, metric = "coverage")
p_width <- plot_metric(summary_df, metric = "width")
final_plot <- p_cov + p_width + patchwork::plot_layout(ncol = 2, guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 13)
  ) &
  guides(color = guide_legend(nrow = 1, byrow = TRUE))

pdf_path <- file.path(output_dir, paste0(figure_stem, ".pdf"))
png_path <- file.path(output_dir, paste0(figure_stem, ".png"))

ggsave(
  filename = pdf_path,
  plot = final_plot,
  width = 12,
  height = 4.5,
  units = "in",
  bg = "white"
)
ggsave(
  filename = png_path,
  plot = final_plot,
  width = 12,
  height = 4.5,
  units = "in",
  dpi = 400,
  bg = "white"
)

print(summary_df)
cat("Saved repetition-level results to:", rep_csv_path, "\n")
cat("Saved summary to:", summary_csv_path, "\n")
cat("Saved LaTeX table to:", summary_tex_path, "\n")
cat("Saved protocol to:", protocol_txt_path, "\n")
cat("Saved figures to:", pdf_path, "and", png_path, "\n")
