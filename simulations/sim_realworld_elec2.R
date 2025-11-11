suppressPackageStartupMessages({
  library(dynaTree)
  library(ggplot2)
  library(dplyr)
  library(zoo)
  library(tidyr)
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

rst_1 <- conformalRetroAdj(X, y, t_init = t_init_n, alpha = alpha0, kernel = "rbf", methods = "dtaci")
rst_2 <- forwardKRR(X, y, t_init = t_init_n, alpha = alpha0, kernel = "rbf", methods = "dtaci")
rst_3 <- forwardMOA(X, y, t_init = t_init_n, w = 250, alpha = alpha0, methods = "dtaci")
rst_4 <- forwardMOA(X, y, t_init = t_init_n, w = 250, alpha = alpha0, moa_model = "AMRulesRegressor", methods = "dtaci")
rst_add <- conformalRetroAdj(X, y, t_init = t_init_n, alpha = alpha0, kernel = "ntk", methods = "dtaci")

safe_len_mean <- function(rst) {
  if (!is.null(rst$U) && !is.null(rst$L)) mean(rst$U - rst$L, na.rm = TRUE)
  else if (!is.null(rst$mean_len)) rst$mean_len
  else NA_real_
}

print_results <- function(..., digits = 4) {
  lst <- list(...)
  nms <- names(lst)
  if (is.null(nms) || any(nms == "")) nms <- paste0("rst_", seq_along(lst))
  covs <- sapply(lst, function(x) if (!is.null(x$coverage)) x$coverage else NA_real_)
  lens <- sapply(lst, safe_len_mean)
  df <- data.frame(
    Model = nms,
    Coverage = round(covs, digits),
    Mean_Length = round(lens, digits),
    check.names = FALSE
  )
  print(df, row.names = FALSE)
  invisible(df)
}

print_results(
  "RetroAdj" = rst_1,
  "FW-KRR" = rst_2,
  "FW-FIMTDD" = rst_3,
  "FW-AMRules" = rst_4
)

local_coverage <- function(err_vec, w) 1 - rollmean(err_vec, k = w, fill = NA, align = "right")
local_mean     <- function(x, w)        rollmean(x, k = w, fill = NA, align = "right")

label_and_shift <- function(df_wide, t_init, w_keep, value_name){
  df_wide %>%
    pivot_longer(c("retro","okrr","otree","orule"),
                 names_to = "Method_key", values_to = value_name) %>%
    mutate(
      Method = recode(Method_key,
                      retro = "RetroAdj", okrr = "FW-KRR",
                      otree = "FW-FIMTDD",   orule = "FW-AMRules"),
      Method = factor(Method, levels = c("FW-FIMTDD","FW-AMRules","FW-KRR","RetroAdj")),
      x_plot = t + t_init
    ) %>%
    filter(t >= w_keep)
}

pal <- c("FW-FIMTDD"="#90C8AC","FW-AMRules"="#E3C57B","FW-KRR"="#B8A9C9","RetroAdj"="#EE6A70")
base_theme <- theme_bw(base_size = 15) +
  theme(plot.title=element_text(face="bold", size=14),
        plot.subtitle=element_text(size=15),
        axis.title.x  = element_text(size = 15),
        axis.title.y  = element_text(size = 15),
        legend.position="bottom")

build_long <- function(sim_retro, sim_okrr, sim_otree, sim_orule, t_init, w_cov, w_len){
  # coverage
  df_cov <- tibble(
    t     = seq_along(sim_retro$err_t),
    retro = local_coverage(sim_retro$err_t, w_cov),
    okrr  = local_coverage(sim_okrr$err_t,  w_cov),
    otree = local_coverage(sim_otree$err_t, w_cov),
    orule = local_coverage(sim_orule$err_t, w_cov)
  )
  df_cov_long <- label_and_shift(df_cov, t_init, w_cov, "LocalCov")
  # length
  len_retro <- sim_retro$U - sim_retro$L
  len_okrr  <- sim_okrr$U  - sim_okrr$L
  len_otree <- sim_otree$U - sim_otree$L
  len_orule <- sim_orule$U - sim_orule$L
  df_len <- tibble(
    t     = seq_along(len_retro),
    retro = local_mean(len_retro, w_len),
    okrr  = local_mean(len_okrr,  w_len),
    otree = local_mean(len_otree, w_len),
    orule = local_mean(len_orule, w_len)
  )
  df_len_long <- label_and_shift(df_len, t_init, w_len, "LocalLen")
  list(df_cov_long=df_cov_long, df_len_long=df_len_long)
}

plot_cov <- function(df_cov_long){
  ggplot(df_cov_long, aes(x=x_plot, y=LocalCov, color=Method)) +
    geom_hline(yintercept=0.9, color="black", linetype="solid", size=0.6) +
    geom_line(aes(alpha=ifelse(Method=="RetroAdj",1,0.75),
                  size =ifelse(Method=="RetroAdj",1.2,1)), lineend="round") +
    scale_color_manual(values=pal) + scale_alpha_identity() + scale_size_identity() +
    labs(x=TeX("$t$"), y="Local Coverage", color="Method") +
    coord_cartesian(ylim=c(0.7,1.0)) + base_theme
}

plot_len <- function(df_len_long){
  ggplot(df_len_long, aes(x=x_plot, y=LocalLen, color=Method)) +
    geom_line(aes(alpha=ifelse(Method=="RetroAdj",1.2,0.75),
                  size =ifelse(Method=="RetroAdj",1.2,1)), lineend="round") +
    scale_color_manual(values=pal) + scale_alpha_identity() + scale_size_identity() +
    labs(x=TeX("$t$"), y="Local Width", color="Method") +
    coord_cartesian(ylim=c(0, 0.15)) + base_theme
}

t_init <- 250; w_cov <- 250; w_len <- 250; method = "dtaci"
lg_real <- build_long(rst_1, rst_2, rst_3, rst_4, t_init, w_cov, w_len)
df_cov_long_real <- lg_real$df_cov_long; df_len_long_real <- lg_real$df_len_long

p_cov_real <- plot_cov(df_cov_long_real)
p_len_real <- plot_len(df_len_long_real)

row_real   <- (p_cov_real | p_len_real)

final_plot_real <- row_real +
  plot_layout(ncol = 2, heights = c(0.09, 1), guides = "collect") &
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        plot.margin = margin(t=8, r=8, b=8, l=8),
        legend.text = element_text(size = 15)) &
  guides(color = guide_legend(nrow = 1, byrow = TRUE,
                              override.aes = list(size = 3)))

final_plot_real  # 12 x 8
