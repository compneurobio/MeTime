# Installation

Follow the step-by-step guidelines to install MeTime 

## Prerequisites

- R version: `>= 4.3.0`
- OS tested: `Linux/macOS/Windows`
- System dependencies (if any): `<list>`


## Install dependencies manually (optional)

```r
install.packages(c("tidyverse", "knitr"))
BiocManager::install()
```
## Install from GitHub

```r
devtools::install_github("compneurobio/MeTime")
```

## Verify installation

```r
library(MeTime)
packageVersion("MeTime")
```

Expected output:

```text
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

## Troubleshooting

- **Error:** `<error text>`  
  **Fix:** `<resolution>`

## Next

- Continue to [`quickstart.md`](quickstart.md).
