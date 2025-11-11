# Call up the dependencies all at once.

pkgs <- c("RMOA", "sgd", "ggplot2", "dplyr")

load_or_install <- function(pkg){
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

lapply(pkgs, load_or_install)
lapply(list.files("FUNCTIONS", pattern = "\\.R$", full.names = TRUE),
       source, echo = TRUE)

lapply(list.files("ACI", pattern = "\\.R$", full.names = TRUE),
       source, echo = TRUE)
