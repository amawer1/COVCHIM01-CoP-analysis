# Installs the R packages needed by the scripts in "R scripts/".
# Only packages that are not already installed are installed.

options(repos = c(CRAN = "https://cloud.r-project.org"))

cran_packages <- c(
  "tidyverse", "patchwork", "glmnet", "logistf", "pROC", "purrr",
  "ggrepel", "writexl", "svglite"
)
bioc_packages <- c("ComplexHeatmap")

missing <- function(pkgs) pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]

to_install <- missing(cran_packages)
if (length(to_install) > 0) install.packages(to_install)

to_install <- missing(bioc_packages)
if (length(to_install) > 0) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(to_install, update = FALSE, ask = FALSE)
}

still_missing <- missing(c(cran_packages, bioc_packages))
if (length(still_missing) > 0) {
  stop("Failed to install: ", paste(still_missing, collapse = ", "))
}
message("All required R packages are installed.")
