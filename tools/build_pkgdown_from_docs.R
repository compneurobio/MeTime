#!/usr/bin/env Rscript

# Build pkgdown site from docs/ without keeping vignettes/ in the repository.

source(file.path("tools", "sync_docs_to_vignettes.R"))

if (!requireNamespace("pkgdown", quietly = TRUE)) {
  stop("pkgdown must be installed to build the site.")
}

root <- normalizePath(".", mustWork = TRUE)
vignette_dir <- file.path(root, "vignettes")
created_vignettes <- FALSE

if (!dir.exists(vignette_dir)) {
  dir.create(vignette_dir, recursive = TRUE)
  created_vignettes <- TRUE
}

on.exit({
  if (created_vignettes && dir.exists(vignette_dir)) {
    unlink(vignette_dir, recursive = TRUE, force = TRUE)
    message("Removed temporary vignettes/ directory.")
  }
}, add = TRUE)

sync_docs_to_vignettes(root = root, vin_root = vignette_dir)
pkgdown::build_site()
