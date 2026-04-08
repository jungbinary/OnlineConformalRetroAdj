# Install and attach the main R dependencies used in this repository,
# then source the core method files under `R/` and `aci/`.

pkgs <- c(
  "RMOA", "rJava", "reticulate", "sgd",
  "dplyr", "ggplot2", "zoo", "patchwork", "readr", "scales",
  "latex2exp", "cowplot", "tidyr", "future.apply", "dynaTree",
  "colorspace", "gridExtra", "ggforce"
)

load_or_install <- function(pkg, repos = "https://cloud.r-project.org") {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = repos)
  }
  suppressPackageStartupMessages(
    library(pkg, character.only = TRUE)
  )
  invisible(pkg)
}

source_file <- function(path) {
  if (!file.exists(path)) {
    stop("Missing source file: ", path, call. = FALSE)
  }
  source(path, echo = FALSE)
  invisible(path)
}

core_files <- c(
  file.path("R", "baseKRR.R"),
  file.path("aci", "dtaci.R"),
  file.path("aci", "agaci.R"),
  file.path("aci", "sfogd.R"),
  file.path("aci", "saocp.R"),
  file.path("R", "conformalRetroAdj.R"),
  file.path("R", "forwardKRR.R"),
  file.path("R", "forwardMOA.R"),
  file.path("R", "forwardRiver.R")
)

invisible(lapply(pkgs, load_or_install))
invisible(lapply(core_files, source_file))

if (!file.exists(file.path("python", "river_helpers.py"))) {
  warning(
    "python/river_helpers.py was not found. River-based experiments will not run.",
    call. = FALSE
  )
}

message("Dependencies loaded and core files sourced from R/ and aci/.")
