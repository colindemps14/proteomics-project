# Run once in the RStudio Console from the project folder.
cran_packages <- c("readxl", "ggplot2", "knitr", "BiocManager", "weightedGCM")
missing <- setdiff(cran_packages, rownames(installed.packages()))
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
bioc_packages <- c("graph", "RBGL")
missing <- setdiff(bioc_packages, rownames(installed.packages()))
if (length(missing)) BiocManager::install(missing, ask = FALSE, update = FALSE)
if (!requireNamespace("pcalg", quietly = TRUE)) {
  install.packages("pcalg", repos = "https://cloud.r-project.org")
}
