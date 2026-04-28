# MeTime Docs Hub

Use this folder as the single source of truth for user-facing documentation pages.

## Sections

- `getting-started/`: first 10-minute onboarding flow ([README](getting-started/README.md)).
- `user-guides/`: task-oriented usage guidance ([README](user-guides/README.md)).
- `reference/`: architecture and API mental model ([README](reference/README.md)).
- `case-studies/`: narrative, proof-oriented examples ([README](case-studies/README.md)).
- `developer/`: contribution and extension workflows ([README](developer/README.md)).

## Writing guidelines

1. Duplicate the right template page.
2. Replace placeholders.
3. Keep sections in the same order for consistency.
4. Link each page to related docs where useful.

## Website publishing

The pkgdown website is generated through CI from this repository's documentation and vignettes.

## Installation fallback

If local installation fails, use the container path in [getting-started/docker-fallback.md](getting-started/docker-fallback.md).
