# Installation

## Installation

Follow the step-by-step guidelines to install MeTime

### Prerequisites

- R version: `>= 3.5.0`

- OS tested: `Linux/macOS/Windows`

- System dependencies (Debian libraries specified):
  `libhdf5-dev, libcurl4-gnutls-dev, libssl-dev, libxml2-dev, libpng-dev, libxt-dev, zlib1g-dev, libbz2-dev, liblzma-dev, libglpk40, libgit2-dev`

- Follow the installation steps below

### 1. Install dependencies manually

``` r

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools", repos = "https://cloud.r-project.org")
}

# Use both CRAN + Bioconductor repositories during dependency resolution
options(repos = BiocManager::repositories())
```

### 2. Install package from GitHub

``` r

# Install MeTime and dependencies
# (set dependencies = TRUE to install Imports/Suggests)
devtools::install_github("compneurobio/MeTime", dependencies = TRUE)
```

### 3. Alternative installation

If you prefer automatic dependency solving across CRAN and Bioconductor,
`pak` is often the easiest route:

``` r

if (!requireNamespace("pak", quietly = TRUE)) {
  install.packages("pak", repos = "https://cloud.r-project.org")
}
pak::pak("compneurobio/MeTime")
```

### 4. Verify installation

``` r

library(MeTime)
packageVersion("MeTime")
```

Expected output:

``` text
Loading required package: tidyverse
── Attaching core tidyverse packages ──────────────────────────────────────── tidyverse 2.0.0 ──
✔ dplyr     1.1.4     ✔ readr     2.1.6
✔ forcats   1.0.1     ✔ stringr   1.6.0
✔ ggplot2   4.0.1     ✔ tibble    3.3.1
✔ lubridate 1.9.4     ✔ tidyr     1.3.2
✔ purrr     1.1.0     
── Conflicts ────────────────────────────────────────────────────────── tidyverse_conflicts() ──
✖ dplyr::filter() masks stats::filter()
✖ dplyr::lag()    masks stats::lag()
ℹ Use the conflicted package to force all conflicts to become errors

[1] '1.0.0'
```

### Troubleshooting

- **Error:** package build fails due to missing system libraries.  
  **Fix:** install OS-level dependencies listed in prerequisites, then
  reinstall.

### Next

- Continue to
  [`quickstart.md`](https://compneurobio.github.io/MeTime/articles/quickstart.md).
