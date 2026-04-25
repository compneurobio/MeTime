# Contributing to MeTime

Thanks for your interest in improving MeTime.

## Ways to contribute

- Report bugs and confusing documentation.
- Propose enhancements and new analyses.
- Improve examples, vignettes, and tutorials.
- Submit pull requests for fixes and features.

## Development workflow

1.  Fork the repository.
2.  Create a focused branch (`feat/...` or `fix/...`).
3.  Make your change with matching documentation updates.
4.  Run local checks.
5.  Open a pull request with context and a reproducible example.

## Suggested local setup

``` r
install.packages(c("devtools", "roxygen2", "testthat"))
```

If available in your environment:

``` r
devtools::document()
devtools::check()
```

## Pull request checklist

Scope is focused and explained clearly.

User-facing docs are updated (`README`, `docs/`, or `vignettes/`).

Any API change is reflected in examples.

`NEWS.md` updated for user-visible changes.

Backward compatibility considered.

## Reporting issues

Please include:

- Minimal reproducible example.
- Full error output.
- [`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html).
- Input shape details.

## Code style notes

- Keep functions modular and composable.
- Prefer explicit argument names in examples.
- Preserve naming conventions (`add_*`, `mod_*`, `calc_*`, `meta_*`,
  `get_*`, `write_*`).

## Code of conduct

By participating, you agree to follow the project Code of Conduct in
`CODE_OF_CONDUCT.md`.
