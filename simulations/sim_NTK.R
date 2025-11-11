suppressPackageStartupMessages({
  library(dynaTree)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(colorspace)
  library(scales)
  library(zoo)
  library(patchwork)
})

load_elec2 <- function(csv_path = "./simulations/elec2.csv") {
  if (exists("elec2", where = as.environment("package:dynaTree"), inherits = FALSE)) {
    data("elec2", package = "dynaTree")
    as.data.frame(elec2)
  } else if (file.exists(csv_path)) {
    read.csv(csv_path)
  } else stop("elec2 dataset not found.")
}

dat <- load_elec2()
dat <- na.omit(dat)
target <- "nswprice"
y_all <- dat[[target]]

lag_k <- 10L
X_embed <- embed(y_all, lag_k + 1)

G <- rbind(
  X_embed[1:250, , drop = FALSE],
  X_embed[(nrow(X_embed) - 2750 + 1):nrow(X_embed), , drop = FALSE]
)

y <- G[, 1]
X <- G[, -1, drop = FALSE]

t_init_n <- 250
alpha0 <- 0.1

rst_rbf <- conformalRetroAdj(X, y, t_init = t_init_n, alpha = alpha0, kernel = "rbf", methods = "dtaci")
rst_ntk <- conformalRetroAdj(X, y, t_init = t_init_n, alpha = alpha0, kernel = "ntk", methods = "dtaci")

# Two-way comparison: RetroAdj (RBF) vs RetroAdj (NTK)
method_labels <- c(rbf = "RetroAdj (RBF)", ntk = "RetroAdj (NTK)")
pal2 <- c("RetroAdj (RBF)" = "#EE6A70", "RetroAdj (NTK)" = "#4C78A8")
method_levels2 <- c("RetroAdj (RBF)", "RetroAdj (NTK)")

label_and_shift_2 <- function(df_wide, t_init, w_keep, value_name){
  df_wide %>%
    pivot_longer(c("rbf","ntk"), names_to = "Method_key", values_to = value_name) %>%
    mutate(Method = recode(Method_key, !!!method_labels),
           Method = factor(Method, levels = method_levels2),
           x_plot = t + t_init) %>%
    filter(t >= w_keep)
}

build_long_2 <- function(sim_rbf, sim_ntk, t_init, w_cov, w_len){
  df_cov <- tibble(
    t = seq_along(sim_rbf$err_t),
    rbf = local_coverage(sim_rbf$err_t, w_cov),
    ntk = local_coverage(sim_ntk$err_t, w_cov)
  )
  df_cov_long <- label_and_shift_2(df_cov, t_init, w_cov, "LocalCov")
  
  len_rbf <- sim_rbf$U - sim_rbf$L
  len_ntk <- sim_ntk$U - sim_ntk$L
  df_len <- tibble(
    t = seq_along(len_rbf),
    rbf = local_mean(len_rbf, w_len),
    ntk = local_mean(len_ntk, w_len)
  )
  df_len_long <- label_and_shift_2(df_len, t_init, w_len, "LocalLen")
  
  list(df_cov_long = df_cov_long, df_len_long = df_len_long)
}

plot_cov <- function(df_cov_long){
  ggplot(df_cov_long, aes(x = x_plot, y = LocalCov, color = Method)) +
    geom_hline(yintercept = 0.9, color = "black", linetype = "solid", linewidth = 0.6) +
    geom_line(aes(alpha = ifelse(Method == "RetroAdj (NTK)", 1, 0.85),
                  linewidth = ifelse(Method == "RetroAdj (NTK)", 1.2, 1.0)),
              lineend = "round") +
    scale_color_manual(values = pal2) +
    scale_alpha_identity() + scale_linewidth_identity() +
    labs(x = latex2exp::TeX("$t$"), y = "Local Coverage", color = "Method") +
    coord_cartesian(ylim = c(0.8, 1.0)) + base_theme
}

plot_len <- function(df_len_long){
  ggplot(df_len_long, aes(x = x_plot, y = LocalLen, color = Method)) +
    geom_line(aes(alpha = ifelse(Method == "RetroAdj (NTK)", 1, 0.85),
                  linewidth = ifelse(Method == "RetroAdj (NTK)", 1.2, 1.0)),
              lineend = "round") +
    scale_color_manual(values = pal2) +
    scale_alpha_identity() + scale_linewidth_identity() +
    labs(x = latex2exp::TeX("$t$"), y = "Local Width", color = "Method") +
    coord_cartesian(ylim = c(0, 0.10)) + base_theme
}

t_init <- 250; w_cov <- 250; w_len <- 250
lg_two <- build_long_2(rst_rbf, rst_ntk, t_init, w_cov, w_len)
df_cov_long_two <- lg_two$df_cov_long
df_len_long_two <- lg_two$df_len_long

p_cov_two <- plot_cov(df_cov_long_two)
p_len_two <- plot_len(df_len_long_two)

plot_elec2_ntk <- (p_cov_two | p_len_two)

final_plot_real_elec2 <- (p_cov_two | p_len_two) +
  patchwork::plot_layout(ncol = 2, guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 15),
    plot.margin = margin(t = 8, r = 8, b = 8, l = 8)
  ) &
  guides(color = guide_legend(nrow = 1, byrow = TRUE,
                              override.aes = list(linewidth = 3)))

final_plot_real_elec2

