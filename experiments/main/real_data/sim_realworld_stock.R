suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(zoo)
  library(tidyr)
  library(patchwork)
  library(ggforce)
  library(latex2exp)
  library(scales)
})

source_if_missing <- function(fun_name, path) {
  if (!exists(fun_name, mode = "function")) {
    source(path, encoding = "UTF-8")
  }
}

source_if_missing("dtaci", "aci/dtaci.R")
source_if_missing("agaci", "aci/agaci.R")
source_if_missing("sfogd", "aci/sfogd.R")
source_if_missing("saocp", "aci/saocp.R")
source_if_missing("krr_init", "R/baseKRR.R")
source_if_missing("conformalRetroAdj", "R/conformalRetroAdj.R")
source_if_missing("forwardKRR", "R/forwardKRR.R")
source_if_missing("forwardMOA", "R/forwardMOA.R")
source_if_missing("forwardRiver", "R/forwardRiver.R")

if (!nzchar(Sys.getenv("RETICULATE_PYTHON"))) {
  py_candidates <- unique(c(
    Sys.which("python"),
    file.path(Sys.getenv("USERPROFILE"), "anaconda3", "python.exe")
  ))
  py_candidates <- py_candidates[nzchar(py_candidates) & file.exists(py_candidates)]
  if (length(py_candidates)) {
    Sys.setenv(RETICULATE_PYTHON = normalizePath(py_candidates[[1]], winslash = "/", mustWork = TRUE))
  }
}

py_module_dir <- file.path(getwd(), "python")
results_root <- file.path("experiments", "results")
figure_dir <- file.path(results_root, "FIGURE")
csv_dir <- file.path(results_root, "CSV")
base_theme <- theme_bw(base_size = 15) +
  theme(plot.title=element_text(face="bold", size=14),
        plot.subtitle=element_text(size=15),
        axis.title.x  = element_text(size = 15),
        axis.title.y  = element_text(size = 15),
        legend.position="bottom")

resolve_existing_path <- function(candidates) {
  for (candidate in candidates) {
    if (file.exists(candidate)) {
      return(candidate)
    }
  }
  stop("Could not find any of: ", paste(candidates, collapse = ", "))
}

csv_path <- resolve_existing_path(c(
  "data/real/AIG.csv",
  file.path(getwd(), "data", "real", "AIG.csv")
))

dat <- read.csv(csv_path)
head(dat)
y_all <- dat[(9012-3000):(9012+3000),'Close']

idx <- (9012 - 3000):(9012 + 3000)
dates_all   <- as.Date(dat[idx, "Date"])
lag_k <- 10L
dates_embed <- dates_all[(lag_k + 250 + 1):length(dates_all)]

df_plot <- data.frame(
  Date  = as.Date(dat[idx, "Date"]),  
  Close = dat[idx, "Close"]
)

crisis_date <- as.Date("2008-09-15")

ggplot(df_plot, aes(x = Date, y = Close)) +
  geom_line(color = "orange", linewidth = 1.1) +
  geom_vline(xintercept = as.numeric(crisis_date),
             linetype = "dashed", color = "gray40", linewidth = 0.8) +
  annotate("text",
           x = crisis_date + 60, 
           y = max(df_plot$Close) * 0.95,
           label = "2008-09-15", color = "gray40", size = 3.8, hjust = 0) +
  labs(
    x = "Date",
    y = "Closing Price (USD)"
  ) +
  scale_x_date(
    date_labels = "%Y-%m",      
    date_breaks = "5 year"   
  ) +
  base_theme

start_date <- as.Date(dat[(9012 - 3000 + 250 + 250):(9012+3000), "Date"])
tail(start_date)

X_embed <- embed(y_all, lag_k + 1)

G <- X_embed
y <- G[, 1]
X <- G[, -1, drop = FALSE]

t_init_n <- 250
alpha0 <- 0.1

rst_1 <- conformalRetroAdj(X, y, t_init = t_init_n, alpha = alpha0, kernel = "rbf", methods = "dtaci")
rst_2 <- forwardKRR(X, y, t_init = t_init_n, alpha = alpha0, kernel = "rbf", methods = "dtaci")
rst_3 <- forwardMOA(X, y, t_init = t_init_n, w = 250, alpha = alpha0, moa_model = "FIMTDD", methods = "dtaci")
rst_4 <- forwardMOA(X, y, t_init = t_init_n, w = 250, alpha = alpha0, moa_model = "AMRulesRegressor", methods = "dtaci")
rst_5 <- forwardRiver(
  X, y,
  t_init = t_init_n,
  w = 250,
  alpha = alpha0,
  river_model = "AMFRegressor",
  methods = "dtaci",
  river_control = list(seed = 1L),
  py_module_path = py_module_dir
)

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

result_summary <- print_results(
  "RetroAdj" = rst_1,
  "FW-KRR" = rst_2,
  "FW-FIMTDD" = rst_3,
  "FW-AMRules" = rst_4,
  "FW-AMF" = rst_5
)

local_coverage <- function(err_vec, w) 1 - rollmean(err_vec, k = w, fill = NA, align = "right")
local_mean     <- function(x, w)        rollmean(x, k = w, fill = NA, align = "right")

label_and_shift <- function(df_wide, t_init, w_keep, value_name){
  df_wide %>%
    pivot_longer(c("retro","okrr","otree","orule","oamf"),
                 names_to = "Method_key", values_to = value_name) %>%
    mutate(
      Method = recode(Method_key,
                      retro = "RetroAdj", okrr = "FW-KRR",
                      otree = "FW-FIMTDD",   orule = "FW-AMRules", oamf = "FW-AMF"),
      Method = factor(Method, levels = c("FW-FIMTDD","FW-AMRules","FW-AMF","FW-KRR","RetroAdj")),
      x_plot = t + t_init
    ) %>%
    filter(t >= w_keep)
}

pal <- c("FW-FIMTDD"="#90C8AC","FW-AMRules"="#E3C57B","FW-AMF"="#9DBCE3","FW-KRR"="#B8A9C9","RetroAdj"="#EE6A70")
legend_order <- c("RetroAdj", "FW-KRR", "FW-AMF", "FW-AMRules", "FW-FIMTDD")

build_long <- function(sim_retro, sim_okrr, sim_otree, sim_orule, sim_oamf, t_init, w_cov, w_len, date){
  
  df_cov <- tibble(
    t     = date[251:5741],
    retro = local_coverage(sim_retro$err_t, w_cov)[251:5741],
    okrr  = local_coverage(sim_okrr$err_t,  w_cov)[251:5741],
    otree = local_coverage(sim_otree$err_t, w_cov)[251:5741],
    orule = local_coverage(sim_orule$err_t, w_cov)[251:5741],
    oamf  = local_coverage(sim_oamf$err_t,  w_cov)[251:5741]
  )
  df_cov_long <- label_and_shift(df_cov, t_init, w_cov, "LocalCov")
  
  len_retro <- sim_retro$U - sim_retro$L
  len_okrr  <- sim_okrr$U  - sim_okrr$L
  len_otree <- sim_otree$U - sim_otree$L
  len_orule <- sim_orule$U - sim_orule$L
  len_oamf  <- sim_oamf$U - sim_oamf$L
  df_len <- tibble(
    t     = date[251:5741],
    retro = local_mean(len_retro, w_len)[251:5741],
    okrr  = local_mean(len_okrr,  w_len)[251:5741],
    otree = local_mean(len_otree, w_len)[251:5741],
    orule = local_mean(len_orule, w_len)[251:5741],
    oamf  = local_mean(len_oamf,  w_len)[251:5741]
  )
  df_len_long <- label_and_shift(df_len, t_init, w_len, "LocalLen")
  list(df_cov_long=df_cov_long, df_len_long=df_len_long)
}

plot_cov <- function(df_cov_long){
  ggplot(df_cov_long, aes(x=x_plot, y=LocalCov, color=Method)) +
    geom_hline(yintercept=0.9, color="black", linetype="solid", size=0.6) +
    geom_line(data = subset(df_cov_long, Method != "RetroAdj"),
              alpha = 0.75, size = 1, lineend="round") +
    geom_line(data = subset(df_cov_long, Method == "RetroAdj"),
              alpha = 1, size = 1.2, lineend="round") +
    scale_color_manual(values=pal, breaks = legend_order) +
    labs(x=TeX("$t$"), y="Local Coverage", color="Method") +
    scale_x_date(
      date_labels = "%Y-%m",    
      date_breaks = "5 year"   
    ) +
    coord_cartesian(ylim=c(0.7,1.0)) + base_theme
}

plot_len <- function(df_len_long){
  if (!inherits(df_len_long$x_plot, "Date")) {
    df_len_long$x_plot <- as.Date(df_len_long$x_plot)
  }
  
  ggplot(df_len_long, aes(x = x_plot, y = LocalLen, color = Method, group = Method)) +
    geom_line(data = subset(df_len_long, Method != "RetroAdj"),
              alpha = 0.75, size = 1, lineend = "round") +
    geom_line(data = subset(df_len_long, Method == "RetroAdj"),
              alpha = 1.0, size = 1.2, lineend = "round") +
    scale_color_manual(values = pal, breaks = legend_order) +
    labs(x = latex2exp::TeX("$t$"), y = "Local Width", color = "Method") +
    scale_x_date(date_labels = "%Y-%m", date_breaks = "5 years") +
    coord_cartesian(ylim = c(0, 120)) +
    base_theme
}

t_init <- 250; w_cov <- 250; w_len <- 250; method = "dtaci"
lg_real <- build_long(rst_1, rst_2, rst_3, rst_4, rst_5, t_init, w_cov, w_len, dates_embed)
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

dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(csv_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(
  result_summary,
  file = file.path(csv_dir, "aig_summary_amf.csv"),
  row.names = FALSE
)
write.csv(
  df_cov_long_real,
  file = file.path(csv_dir, "aig_local_coverage_amf.csv"),
  row.names = FALSE
)
write.csv(
  df_len_long_real,
  file = file.path(csv_dir, "aig_local_width_amf.csv"),
  row.names = FALSE
)
ggsave(
  filename = file.path(figure_dir, "aig_plot_amf.png"),
  plot = final_plot_real,
  width = 12,
  height = 4,
  dpi = 300,
  bg = "white"
)
ggsave(
  filename = file.path(figure_dir, "aig_plot_amf.pdf"),
  plot = final_plot_real,
  width = 12,
  height = 4,
  bg = "white"
)

final_plot_real

