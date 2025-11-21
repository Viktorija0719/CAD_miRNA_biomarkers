# -----------------------------
# CRAN packages
# -----------------------------
install.packages(c(
  "tidyr",
  "readr",
  "caret",
  "ggplot2",
  "ggrepel",
  "scales",
  "dplyr",
  "stringr",
  "optparse",
  "tibble",
  "purrr",
  "broom",
  "splines"
), repos = "https://cran.r-project.org")

# -----------------------------
# Bioconductor packages
# -----------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cran.r-project.org")

BiocManager::install(c(
  "edgeR",
  "limma",
  "DESeq2",
  "EnhancedVolcano"  # for volcano plots
), ask = FALSE, update = TRUE)

# -----------------------------
# Optional: validate installation
# -----------------------------
BiocManager::valid()