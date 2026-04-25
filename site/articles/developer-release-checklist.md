# Release Checklist

## Release Checklist

> One-line summary: pre-release QA for package and docs website.

### Package

`DESCRIPTION` version bumped

`NEWS.md` updated

`R CMD check` passes

### Documentation

README quickstart still valid

docs templates updated where needed

pkgdown site builds without broken links

### Publish

Tag release

Publish GitHub release notes

Deploy pkgdown site

Publish Docker image to GHCR (if tagging a release)

Verify package appears under GitHub Packages
