suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(zoo)
  library(tidyr)
  library(patchwork)
  library(latex2exp)
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

resolve_existing_path <- function(candidates) {
  for (candidate in candidates) {
    if (file.exists(candidate)) {
      return(candidate)
    }
  }
  stop("Could not find any of: ", paste(candidates, collapse = ", "))
}

data_path <- resolve_existing_path(c(
  "data/real/communities_data.csv",
  file.path(getwd(), "data", "real", "communities_data.csv")
))
demographic_cols <- c(3, 4, 5, 6)  
set.seed(2026)

crime_data <- read.csv(data_path, header = FALSE)
n_total <- nrow(crime_data)
col_y <- ncol(crime_data) ; col_y 

selected_col <- demographic_cols[2]

crime_sorted <- crime_data[order(-crime_data[[selected_col]]), ]

crime_mixed <- rbind(
  crime_sorted[1:250, ],
  crime_sorted[seq(n_total, 251, by = -1), ]
)

y_all <- crime_mixed[[col_y]]
X_all <- as.matrix(crime_mixed[, -c(col_y, demographic_cols), drop = FALSE])

crime_1 <- conformalRetroAdj(X_all, y_all, alpha = 0.1, methods = 'dtaci')
crime_2 <- forwardKRR(X_all, y_all, alpha = 0.1, t_init = 250, kernel = "rbf", methods = "dtaci")
crime_3 <- forwardMOA(X_all, y_all, alpha = 0.1, t_init = 250, w = 250, moa_model = "FIMTDD", methods = "dtaci")
crime_4 <- forwardMOA(X_all, y_all, alpha = 0.1, t_init = 250, w = 250, moa_model = "AMRulesRegressor", methods = "dtaci")
crime_5 <- forwardRiver(
  X_all, y_all,
  alpha = 0.1,
  t_init = 250,
  w = 250,
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
  "RetroAdj" = crime_1,
  "FW-KRR" = crime_2,
  "FW-FIMTDD" = crime_3,
  "FW-AMRules" = crime_4,
  "FW-AMF" = crime_5
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
base_theme <- theme_bw(base_size = 15) +
  theme(plot.title=element_text(face="bold", size=14),
        plot.subtitle=element_text(size=15),
        axis.title.x  = element_text(size = 15),
        axis.title.y  = element_text(size = 15),
        legend.position="bottom")

build_long <- function(sim_retro, sim_okrr, sim_otree, sim_orule, sim_oamf, t_init, w_cov, w_len){
  
  df_cov <- tibble(
    t     = seq_along(sim_retro$err_t),
    retro = local_coverage(sim_retro$err_t, w_cov),
    okrr  = local_coverage(sim_okrr$err_t,  w_cov),
    otree = local_coverage(sim_otree$err_t, w_cov),
    orule = local_coverage(sim_orule$err_t, w_cov),
    oamf  = local_coverage(sim_oamf$err_t,  w_cov)
  )
  df_cov_long <- label_and_shift(df_cov, t_init, w_cov, "LocalCov")
  
  len_retro <- sim_retro$U - sim_retro$L
  len_okrr  <- sim_okrr$U  - sim_okrr$L
  len_otree <- sim_otree$U - sim_otree$L
  len_orule <- sim_orule$U - sim_orule$L
  len_oamf  <- sim_oamf$U - sim_oamf$L
  df_len <- tibble(
    t     = seq_along(len_retro),
    retro = local_mean(len_retro, w_len),
    okrr  = local_mean(len_okrr,  w_len),
    otree = local_mean(len_otree, w_len),
    orule = local_mean(len_orule, w_len),
    oamf  = local_mean(len_oamf,  w_len)
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
    coord_cartesian(ylim=c(0.8,1.0)) + base_theme
}

plot_len <- function(df_len_long){
  ggplot(df_len_long, aes(x=x_plot, y=LocalLen, color=Method)) +
    geom_line(data = subset(df_len_long, Method != "RetroAdj"),
              alpha = 0.75, size = 1, lineend="round") +
    geom_line(data = subset(df_len_long, Method == "RetroAdj"),
              alpha = 1, size = 1.2, lineend="round") +
    scale_color_manual(values=pal, breaks = legend_order) +
    labs(x=TeX("$t$"), y="Local Width", color="Method") +
    coord_cartesian(ylim=c(0, 2.0)) + base_theme
}

t_init <- 250; w_cov <- 250; w_len <- 250; method = "dtaci"
lg_crime <- build_long(crime_1, crime_2, crime_3, crime_4, crime_5, t_init, w_cov, w_len)
df_cov_long_crime <- lg_crime$df_cov_long; df_len_long_crime <- lg_crime$df_len_long

p_cov_crime <- plot_cov(df_cov_long_crime)
p_len_crime <- plot_len(df_len_long_crime)


row_crime   <- (p_cov_crime | p_len_crime)

final_plot_crime <- row_crime +
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
  file = file.path(csv_dir, "crime_summary_amf.csv"),
  row.names = FALSE
)
write.csv(
  df_cov_long_crime,
  file = file.path(csv_dir, "crime_local_coverage_amf.csv"),
  row.names = FALSE
)
write.csv(
  df_len_long_crime,
  file = file.path(csv_dir, "crime_local_width_amf.csv"),
  row.names = FALSE
)
ggsave(
  filename = file.path(figure_dir, "crime_plot_amf.png"),
  plot = final_plot_crime,
  width = 12,
  height = 4,
  dpi = 300,
  bg = "white"
)
ggsave(
  filename = file.path(figure_dir, "crime_plot_amf.pdf"),
  plot = final_plot_crime,
  width = 12,
  height = 4,
  bg = "white"
)

final_plot_crime 

