# Connect `docs/` content to pkgdown (website setup)

This guide shows how to publish a pkgdown website while keeping `docs/` as your editing workspace.

## Option A (recommended): `docs/` as source + auto-sync to `vignettes/`

`pkgdown` builds article pages from package vignettes/articles (`vignettes/*.Rmd`), not directly from arbitrary markdown files.

Use this workflow:

1. Author docs in `docs/**/*.md`.
2. Convert selected docs pages into `vignettes/*.Rmd` with the script below.
3. Build pkgdown site.

### Sync script

Use `tools/sync_docs_to_vignettes.R` (added in this repo):

```bash
Rscript tools/sync_docs_to_vignettes.R
```

This script copies markdown pages into `vignettes/` as article `.Rmd` files with pkgdown-compatible YAML headers.

## Option B: Keep long-form docs in GitHub only

If you do not want duplicated article files, keep docs in `docs/` and add direct navbar links to GitHub pages from `_pkgdown.yml`.

---

## 1) Minimal package setup

Install tools:

```r
install.packages(c("pkgdown", "rmarkdown", "knitr"))
```

## 2) Configure `_pkgdown.yml`

The current config now includes a docs-oriented navbar and article structure. Edit titles and paths as docs mature.

## 3) Build website locally

```bash
Rscript tools/sync_docs_to_vignettes.R
Rscript -e "pkgdown::build_site(preview = FALSE)"
```

Generated static site appears under `docs/` by default (pkgdown output directory).

> Note: this repo now tracks `docs/` source files. Prefer deploying built site artifacts through GitHub Actions (`gh-pages`) instead of committing generated site output to the main branch.

## 4) Publish with GitHub Pages (recommended)

Use GitHub Actions and deploy from `gh-pages` branch:

```r
usethis::use_github_action("pkgdown")
```

Then in GitHub repo settings:

- Pages source: `Deploy from a branch`
- Branch: `gh-pages`

## 5) Suggested release workflow

1. Edit docs templates under `docs/`.
2. Run sync script.
3. Build/test site locally.
4. Commit source docs + vignette updates.
5. Push to trigger pkgdown deployment.

