#!/usr/bin/env Rscript

# Sync markdown docs into vignettes/*.Rmd so pkgdown can render them as articles.

root <- normalizePath('.', mustWork = TRUE)
docs_root <- file.path(root, 'docs')
vin_root <- file.path(root, 'vignettes')

if (!dir.exists(docs_root)) {
  stop('docs/ folder not found. Nothing to sync.')
}

if (!dir.exists(vin_root)) {
  dir.create(vin_root, recursive = TRUE)
}

md_files <- list.files(
  docs_root,
  pattern = '\\.(md|MD)$',
  recursive = TRUE,
  full.names = TRUE
)

# Exclude docs hub from article generation.
exclude <- c('README.md')
md_files <- md_files[!basename(md_files) %in% exclude]

if (length(md_files) == 0) {
  message('No markdown files found to sync.')
  quit(status = 0)
}

slugify <- function(x) {
  x <- tolower(x)
  x <- gsub('[^a-z0-9]+', '-', x)
  x <- gsub('(^-|-$)', '', x)
  x
}

for (f in md_files) {
  rel <- sub(paste0('^', docs_root, '/?'), '', f)
  stem <- tools::file_path_sans_ext(rel)
  out_name <- paste0(slugify(gsub('/', '-', stem)), '.Rmd')
  out <- file.path(vin_root, out_name)

  lines <- readLines(f, warn = FALSE)
  title <- basename(stem)

  if (length(lines) > 0 && grepl('^#\\s+', lines[1])) {
    title <- sub('^#\\s+', '', lines[1])
  }

  header <- c(
    '---',
    paste0('title: "', title, '"'),
    'output: rmarkdown::html_vignette',
    'vignette: >',
    paste0('  %\\VignetteIndexEntry{', title, '}'),
    '  %\\VignetteEngine{knitr::rmarkdown}',
    '  %\\VignetteEncoding{UTF-8}',
    '---',
    ''
  )

  writeLines(c(header, lines), out)
  message('Synced: ', rel, ' -> vignettes/', out_name)
}

message('Done.')
