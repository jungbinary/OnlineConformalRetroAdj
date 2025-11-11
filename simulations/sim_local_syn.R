# ---- Packages ----
library(dplyr)
library(tidyr)
library(ggplot2)
library(zoo)
library(latex2exp)
library(patchwork)
library(gridExtra) 

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

dat <- gen_beta_shift_only(n = 1000, t_init = 250, seed = 2020)
X <- dat$X ; y <- dat$y

sim1_retro <- conformalRetroAdj(X = X, y = y, alpha = 0.1, t_init = 250, d = 250, kernel = "rbf", methods = "dtaci")
sim1_okrr <- forwardKRR(X = X, y = y, alpha = 0.1, t_init = 250, w = 250, kernel = "rbf", methods = "dtaci")
sim1_otree <- forwardMOA(X= X, y = y, alpha = 0.1, t_init = 250, w = 250, moa_model = "FIMTDD", methods = "dtaci")
sim1_orule <- forwardMOA(X = X, y = y, alpha = 0.1, t_init = 250, w = 250, moa_model = "AMRulesRegressor", methods = "dtaci")

# ---- Params ----
t_init <- 250
w_cov  <- 250
w_len  <- 250

local_coverage <- function(err_vec, w) 1 - rollmean(err_vec, k = w, fill = NA, align = "right")
local_mean     <- function(x, w)        rollmean(x, k = w, fill = NA, align = "right")

label_and_shift <- function(df_wide, t_init, w_keep, value_name){
  df_wide %>%
    pivot_longer(c("retro","okrr","otree","orule"),
                 names_to = "Method_key", values_to = value_name) %>%
    mutate(
      Method = recode(Method_key,
                      retro  = "RetroAdj",
                      okrr  = "FW-KRR",
                      otree = "FW-FIMTDD",
                      orule = "FW-AMRules"),
      Method = factor(Method, levels = c("FW-FIMTDD","FW-AMRules","FW-KRR", "RetroAdj")),
      x_plot = t + t_init
    ) %>%
    filter(t >= w_keep)
}

pal <- c(
  "FW-FIMTDD"         = "#90C8AC",
  "FW-AMRules"        = "#E3C57B",
  "FW-KRR"            = "#B8A9C9",
  "RetroAdj" = "#EE6A70"
)

base_theme <- theme_bw(base_size = 15) +
  theme(plot.title=element_text(face="bold", size=14),
        plot.subtitle=element_text(size=15),
        axis.title.x  = element_text(size = 15),
        axis.title.y  = element_text(size = 15),
        legend.position="bottom")

# ---- Coverage data ----
df_cov <- tibble(
  t     = seq_along(sim1_retro$err_t),
  retro  = local_coverage(sim1_retro$err_t, w_cov),
  okrr  = local_coverage(sim1_okrr$err_t,  w_cov),
  otree = local_coverage(sim1_otree$err_t, w_cov),
  orule = local_coverage(sim1_orule$err_t, w_cov)
)

df_cov_long <- label_and_shift(df_cov, t_init, w_cov, "LocalCov")

# ---- Length data ----
len_retro <- sim1_retro$U - sim1_retro$L
len_okrr  <- sim1_okrr$U - sim1_okrr$L
len_otree <- sim1_otree$U - sim1_otree$L
len_orule <- sim1_orule$U - sim1_orule$L

df_len <- tibble(
  t     = seq_along(len_retro),
  retro  = local_mean(len_retro,  w_len),
  okrr  = local_mean(len_okrr,  w_len),
  otree = local_mean(len_otree, w_len),
  orule = local_mean(len_orule, w_len)
)

df_len_long <- label_and_shift(df_len, t_init, w_len, "LocalLen")

# ---- Plots ----
p_cov_1 <-
  ggplot(df_cov_long, aes(x = x_plot, y = LocalCov, color = Method)) +
  geom_hline(yintercept = 0.9, color = "black", linetype = "solid", size = 0.6) +
  geom_line(aes(
    alpha = ifelse(Method == "RetroAdj", 1, 0.85),
    size  = ifelse(Method == "RetroAdj", 1.2, 1)
  ),
  lineend = "round"
  ) +
  scale_color_manual(values = pal) +
  scale_alpha_identity() +
  scale_size_identity() +
  labs(x = TeX("$t$"), y = "Local Coverage", color = "Method") +
  coord_cartesian(ylim = c(0.8, 1.0)) +
  base_theme

p_len_1 <-
  ggplot(df_len_long, aes(x = x_plot, y = LocalLen, color = Method)) +
  geom_line(aes(
    alpha = ifelse(Method == "RetroAdj", 1.2, 0.85),
    size  = ifelse(Method == "RetroAdj", 1.2, 1)
  ),
  lineend = "round"
  ) +
  scale_color_manual(values = pal) +
  scale_alpha_identity() +
  scale_size_identity() +
  labs(x = TeX("$t$"), y = "Local Width", color = "Method")+
  base_theme

(p_cov_1 | p_len_1 )


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

dat <- gen_wendland(n = 1000, seed = 10)
X <- dat$X ; y <- dat$y

sim2_retro <- conformalRetroAdj(X = X, y = y, alpha = 0.1, t_init = 250, d = 250, kernel = "rbf", methods = "dtaci")
sim2_okrr <- forwardKRR(X = X, y = y, alpha = 0.1, t_init = 250, w = 250, kernel = "rbf", methods = "dtaci")
sim2_otree <- forwardMOA(X= X, y = y, alpha = 0.1, t_init = 250, w = 250, moa_model = "FIMTDD", methods = "dtaci")
sim2_orule <- forwardMOA(X = X, y = y, alpha = 0.1, t_init = 250, w = 250, moa_model = "AMRulesRegressor", methods = "dtaci")


# ---- Params ----
t_init <- 250
w_cov  <- 250
w_len  <- 250

# ---- Coverage data ----
df_cov <- tibble(
  t     = seq_along(sim2_retro$err_t),
  retro  = local_coverage(sim2_retro$err_t, w_cov),
  okrr  = local_coverage(sim2_okrr$err_t,  w_cov),
  otree = local_coverage(sim2_otree$err_t, w_cov),
  orule = local_coverage(sim2_orule$err_t, w_cov)
)

df_cov_long <- label_and_shift(df_cov, t_init, w_cov, "LocalCov")

# ---- Length data ----
len_retro <- sim2_retro$U - sim2_retro$L
len_okrr  <- sim2_okrr$U - sim2_okrr$L
len_otree <- sim2_otree$U - sim2_otree$L
len_orule <- sim2_orule$U - sim2_orule$L

df_len <- tibble(
  t     = seq_along(len_retro),
  retro  = local_mean(len_retro,  w_len),
  okrr  = local_mean(len_okrr,  w_len),
  otree = local_mean(len_otree, w_len),
  orule = local_mean(len_orule, w_len)
)

df_len_long <- label_and_shift(df_len, t_init, w_len, "LocalLen")

# ---- Plots ----
p_cov_2 <-
  ggplot(df_cov_long, aes(x = x_plot, y = LocalCov, color = Method)) +
  geom_hline(yintercept = 0.9, color = "black", linetype = "solid", size = 0.6) +
  geom_line(aes(
    alpha = ifelse(Method == "RetroAdj", 1, 0.85),
    size  = ifelse(Method == "RetroAdj", 1.2, 1)
  ),
  lineend = "round"
  ) +
  scale_color_manual(values = pal) +
  scale_alpha_identity() +
  scale_size_identity() +
  labs(x = TeX("$t$"), y = "Local Coverage", color = "Method") +
  coord_cartesian(ylim = c(0.8, 1.0)) +
  base_theme

p_len_2 <-
  ggplot(df_len_long, aes(x = x_plot, y = LocalLen, color = Method)) +
  geom_line(aes(
    alpha = ifelse(Method == "RetroAdj", 1.2, 0.85),
    size  = ifelse(Method == "RetroAdj", 1.2, 1)
  ),
  lineend = "round"
  ) +
  scale_color_manual(values = pal) +
  scale_alpha_identity() +
  scale_size_identity() +
  labs(x = TeX("$t$"), y = "Local Width", color = "Method")+
  base_theme

(p_cov_2 | p_len_2 )


################################################################################

# 15 X 6

plot_1 <- (p_cov_1 | p_len_1) +
  plot_layout(guides = "collect") + 
  plot_annotation(
    title = "Setting 1 (Linear)",
    theme = theme(
      plot.title = element_text(face = "bold", size = 15))
  ) &
  theme(legend.position = "bottom")


plot_2 <- (p_cov_2 | p_len_2) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Setting 2 (Non-Linear)",
    theme = theme(
      plot.title = element_text(face = "bold", size = 15))
  ) &
  theme(legend.position = "bottom")

plot_1


###############################################################################

library(patchwork)
library(cowplot) 

row_title <- function(txt, size = 16, bottom_margin = 4){
  cowplot::ggdraw() +
    cowplot::draw_label(
      txt, x = 0, y = 1, hjust = 0, vjust = 1,
      fontface = "bold", size = size
    ) +
    theme(plot.margin = margin(t = 0, r = 0, b = bottom_margin, l = 0))
}

title1 <- row_title("Setting 1 (Linear)",     size = 16, bottom_margin = 6)
title2 <- row_title("Setting 2 (Non-Linear)", size = 16, bottom_margin = 6)

row1 <- (p_cov_1 | p_len_1)
row2 <- (p_cov_2 | p_len_2)

final_plot <- (title1 / row1 / title2 / row2) +
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
