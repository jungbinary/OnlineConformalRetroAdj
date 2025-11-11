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

## ------------------------------------------------------------
## Test for Robustness
## ------------------------------------------------------------
## Settings
gammaCandidates <- c(0.001, 0.002, 0.004, 0.008, 0.016, 0.032, 0.064, 0.128, 0.256, 0.512, 1.024, 2.048)
totalRep <- length(gammaCandidates)

alpha_0 = 0.1
t_init_n = 250
dat <- gen_beta_shift_only(n = 1000, t_init = t_init_n, seed = 2026) ; X <- dat$X ; y <- dat$y

RetroAdjCov <- vector(length = totalRep)
RetroAdjLen <- vector(length = totalRep)
FWCov <- vector(length = totalRep)
FWLen <- vector(length = totalRep)

for (i in 1:totalRep){
  RetroAdjRst <- conformalRetroAdj(X, y, alpha = alpha_0, d = 250, t_init = t_init_n, methods = "aci", aci_gam = gammaCandidates[i])
  RetroAdjCov[i] <- RetroAdjRst$coverage
  RetroAdjLen[i] <- RetroAdjRst$mean_len
  FWRst <- original_krr(X, y, alpha = alpha_0, w = 250, t_init = t_init_n, methods = "aci", aci_gam = gammaCandidates[i])
  FWCov[i] <- FWRst$coverage
  FWLen[i] <- FWRst$mean_len
}

RetroAdjRst.DtACI <- conformalRetroAdj(X, y, alpha = alpha_0, d = 250, t_init = t_init_n, methods = "dtaci")
FWRst.DtACI <- original_krr(X, y, alpha = alpha_0, t_init_n, d = 250, methods = "dtaci")


## ------------------------------------------------------------
## Theme
## ------------------------------------------------------------
base_theme <- theme_bw(base_size = 15) +
  theme(
    plot.title   = element_text(face = "bold", size = 14),
    plot.subtitle= element_text(size = 15),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    legend.position = "bottom"
  )

## ------------------------------------------------------------
## Data frames (assumes the vectors/objects already exist)
## ------------------------------------------------------------
df <- data.frame(
  gamma = gammaCandidates,
  RetroAdj_ACI   = RetroAdjCov,
  FWKRR_ACI      = FWCov,
  RetroAdj_Length= RetroAdjLen,
  FWKRR_Length   = FWLen
)

df_cov <- df %>%
  select(gamma, RetroAdj_ACI, FWKRR_ACI) %>%
  pivot_longer(-gamma, names_to = "Method", values_to = "Coverage")

df_len <- df %>%
  transmute(
    gamma,
    RetroAdj_ACI = RetroAdj_Length,
    FWKRR_ACI    = FWKRR_Length
  ) %>%
  pivot_longer(-gamma, names_to = "Method", values_to = "Width")

## ------------------------------------------------------------
## Legend levels, labels, styles
## ------------------------------------------------------------
cols <- c(
  "RetroAdj_ACI"               = "#F8766D",
  "FWKRR_ACI"                  = "#B8A9C9",
  "RetroAdj (DtACI baseline)"  = darken("#F8766D", amount = 0.5),
  "FWKRR (DtACI baseline)"     = darken("#B8A9C9", amount = 0.5)
)

ltys <- c(
  "RetroAdj_ACI"               = "solid",
  "FWKRR_ACI"                  = "solid",
  "RetroAdj (DtACI baseline)"  = "twodash",
  "FWKRR (DtACI baseline)"     = "twodash"
)

levels_all <- c(
  "RetroAdj_ACI",
  "FWKRR_ACI",
  "RetroAdj (DtACI baseline)",
  "FWKRR (DtACI baseline)"
)

labels_all <- c(
  "RetroAdj_ACI"               = "RetroAdj (ACI)",
  "FWKRR_ACI"                  = "FW-KRR (ACI)",
  "RetroAdj (DtACI baseline)"  = "RetroAdj (DtACI)",
  "FWKRR (DtACI baseline)"     = "FW-KRR (DtACI)"
)

## ------------------------------------------------------------
## Horizontal lines for baselines
## ------------------------------------------------------------
df_hcov <- data.frame(
  y     = c(RetroAdjRst.DtACI$coverage, FWRst.DtACI$coverage),
  label = c("RetroAdj (DtACI baseline)", "FWKRR (DtACI baseline)")
)

df_hlen <- data.frame(
  y     = c(RetroAdjRst.DtACI$mean_len, FWRst.DtACI$mean_len),
  label = c("RetroAdj (DtACI baseline)", "FWKRR (DtACI baseline)")
)

## ------------------------------------------------------------
## Coverage plot
## ------------------------------------------------------------
p_cov <- ggplot() +
  geom_line(
    data = df_cov,
    aes(x = gamma, y = Coverage, color = Method, linetype = Method),
    linewidth = 1.5, alpha = 1.2
  ) +
  geom_point(
    data = df_cov,
    aes(x = gamma, y = Coverage, color = Method),
    size = 3
  ) +
  scale_x_log10(
    breaks = c(0.001, 0.004, 0.016, 0.064, 0.256, 1.024),
    labels = number_format(accuracy = 0.001, trim = TRUE)
  ) +
  ## Target coverage line (no legend)
  geom_hline(
    yintercept = 1 - alpha_0,
    color      = "black",
    linewidth  = 0.5,
    show.legend = FALSE
  ) +
  ## DtACI baselines (with legend)
  geom_hline(
    data = df_hcov,
    aes(yintercept = y, color = label, linetype = label),
    linewidth = 1.2, alpha = 0.6, show.legend = TRUE
  ) +
  scale_color_manual(
    name   = "Method",
    values = cols,
    breaks = levels_all,
    labels = labels_all[levels_all]
  ) +
  scale_linetype_manual(
    name   = "Method",
    values = ltys,
    breaks = levels_all,
    labels = labels_all[levels_all]
  ) +
  labs(x = expression(gamma), y = "Coverage") +
  base_theme +
  coord_cartesian(ylim = c(0.8, 1))  

## ------------------------------------------------------------
## Width plot
## ------------------------------------------------------------
p_len <- ggplot() +
  geom_line(
    data = df_len,
    aes(x = gamma, y = Width, color = Method, linetype = Method),
    linewidth = 1.5, alpha = 1.2
  ) +
  geom_point(
    data = df_len,
    aes(x = gamma, y = Width, color = Method),
    size = 3
  ) +
  scale_x_log10(
    breaks = c(0.001, 0.004, 0.016, 0.064, 0.256, 1.024),
    labels = number_format(accuracy = 0.001, trim = TRUE)
  ) +
  geom_hline(
    data = df_hlen,
    aes(yintercept = y, color = label, linetype = label),
    linewidth = 1.2, alpha = 0.6, show.legend = TRUE
  ) +
  scale_color_manual(
    name   = "Method",
    values = cols,
    breaks = levels_all,
    labels = labels_all[levels_all]
  ) +
  scale_linetype_manual(
    name   = "Method",
    values = ltys,
    breaks = levels_all,
    labels = labels_all[levels_all]
  ) +
  labs(x = expression(gamma), y = "Width") +
  base_theme


(p_cov | p_len) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        legend.text = element_text(size = 15))

