# COVCHIM01-CoP-analysis

Analysis of baseline and post inoculation immunology data and its association with virology data in the COVCHIM01 SARS-CoV-2 Controlled Human Infection Model.

Associated manuscript:  'Correlates of protection against SARS-CoV-2 shedding in a controlled human infection study in immune-experienced participants'

Scripts:

i) `Baseline immune markers heatmap and nasal serum correlation plot.r`: An R script for generating pairwise correlation between baseline immune markers using hierarchical clustered heatmaps across specific comparable groups of interest.

ii) `Repeated EN and firth multivariable models.r`: An R script for performing repeated Elastic Net feature selection analysis of important immune markers as well as Firth's penalised regression of quantifiable viral shedding based on the best selected markers and additional descriptive summary plots.

iii) `Single marker logistic regression models.r` An R script for generating additional Volcano and forest plots based on estimated odds ratios following immune marker-specific logistic regression (profile likelihood CIs)

iV) `Ridge regression.r` An R script for peforming ridge regression the anti-spike binding antibody immune markers


Requirements: 

Tested on Windows 11 version 24H2 and and R version 4.4.1. The required packages are listed at the top of each script


Setup:

Install R (version 4.4 or later) and make sure `Rscript` is on your `PATH`. Then, from the repository root, run the command for your operating system:

| OS | Command | What it does |
|---|---|---|
| macOS | `make install_mac` | Installs `nlopt` via [Homebrew](https://brew.sh), then the R packages |
| Linux (Debian/Ubuntu) | `make install_linux` | Installs the system libraries needed to build the R packages via `apt-get` (requires `sudo`), then the R packages |
| Windows | `make install_windows` | Installs the R packages (CRAN provides pre-built binaries, so no system libraries are needed) |

All three run `install_packages.R`, which installs any missing packages from CRAN and Bioconductor. If you don't have `make`, you can run `Rscript install_packages.R` directly; on macOS and Linux, install the system libraries listed in the `Makefile` first.

`make` is not installed on Windows by default. It is included with [Rtools](https://cran.r-project.org/bin/windows/Rtools/) (add Rtools' `usr/bin` folder to your `PATH`), or can be installed with `winget install GnuWin32.Make` or `choco install make`.

The scripts use paths relative to the repository root, so run them from there, e.g. `Rscript "R scripts/Ridge regression.r"`, or open `COVCHIM01-CoP-analysis.Rproj` in RStudio and source them.


Data:

Data are available from the same repository, use Combined_nasal_serum_mastersheet_final.csv in the data folder.


Expected output:

The scripts reproduce figures 2,3 and 4 available witihn the manuscript.