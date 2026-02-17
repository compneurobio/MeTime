# Installation

> One-line summary: how to install MeTime and verify it works.

## Prerequisites

- R version: `>= <VERSION>`
- OS tested: `<Linux/macOS/Windows>`
- System dependencies (if any): `<list>`

## Install from GitHub

```r
# install.packages("pak")
pak::pak("<org>/MeTime")
```

## Install dependencies manually (optional)

```r
install.packages(c("tidyverse", "knitr"))
```

## Verify installation

```r
library(MeTime)
packageVersion("MeTime")
```

Expected output:

```text
[1] '<VERSION>'
```

## Troubleshooting

- **Error:** `<error text>`  
  **Fix:** `<resolution>`

## Next

- Continue to [`quickstart.md`](quickstart.md).
